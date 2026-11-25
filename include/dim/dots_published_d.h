/* This software is developed by caoweiquan322 OR DynamicFatty.
 * All rights reserved.
 *
 * Author: caoweiquan322
 */

/**
  * @file
  * @brief DotsSimplifier.h defines the DotsSimplifier class.
  * @author caoweiquan322
  * @date 2015/12/14
  * @version 1.0
  */
#ifndef DOTSSIMPLIFIERD_H
#define DOTSSIMPLIFIERD_H

#include "qtcompat.h"
#include "DotsException.h"

/**
 * @brief The DotsSimplifier class implements the trajectory simplification algorithm DOTS.
 *
 * DOTS is logogram of "Directed acyclic graph based Online Trajectory Simplification algorithm".
 * Please refer to my paper for details. The DOTS algorithm has several advantages:
 *
 * First, it works in an online manner. So it would adapts to streaming data well.
 * Second, it's a near-optimal algorithm. As a result, the accuracy is compared to those batch mode algorithms
 * like TS/MRPA.
 * Third, it solves not only the min-# problem but also (partially) the min-e problem.
 * Forth, the time cost is relatively low. The time complexity is O(N/M) for each input point. Note that N/M
 * represents the simplification or compression rate.
 */
class DotsSimplifierD : public QObject
{
    Q_OBJECT
public:
    // Regular member methods.
    /**
     * @brief DotsSimplifier is the default constructor.
     * @param parent is the QT parent object.
     */
    explicit DotsSimplifierD(QObject *parent = 0, DotsSimplifierD *cascadeRoot = NULL);

    /**
     * @brief setParameters specifies DOTS settings for trajectory simplification.
     * @param lssdTh is the LSSD threshold for DAG searching.
     * @param k is the factor for upper bound. LSSD that exceeds lssdTh*k would be ignored by DAG searching.
     * @param maxVkSize is the maximum size of each layer of DAG tree.
     */
    void setParameters(double lssdTh, double k = 2.0, int maxVkSize = 1e6);

    /**
     * @brief resetInternalData resets all internal data structures for DOTS algorithm.
     */
    void resetInternalData();

    /**
     * @brief feedData feeds a 2D spatio temporary point to DOTS.
     * @param x is the x position.
     * @param y is the y position.
     * @param t is the timestamp.
     */
    static void batchDotsByIndexD(const QVector< QVector<double> > &coord, const QVector<double> &t,
                           QVector<int> &simplifiedIndexOut, double lssdThreshold);

    inline void feedData(const QVector<double>& coord, double t)
    {
        dim_ = coord.count();
        if (finished)
            DotsException("Feeding data is NOT allowed after the simplifier finished. "\
                          "Suggest calling resetInternalData() first.").raise();
        if (!isCascadeRoot)
            DotsException("We can only feed index to non-root simplifier. "\
                          "Try feedIndex() instead.").raise();

        // Store data.
        
        
        ptt.append(t);
        ptIndex.append(ptIndex.count());
        // Update internal data.
        if (ptIndex.count() == 1)
        {
            if (qAbs(t) > 3600*24*365)
            {
                DotsException("The first timestamp seems so big that DOTS would potentially fail "\
                              "due to numerical errors. Would you please consider normalizing the "\
                              "input data properly first?").raise();
            }
            if (qAbs(coord[0]) > 2000*1000 || (coord.count()>1 && qAbs(coord[1]) > 2000*1000))
            {
                DotsException("The first position seems so big that DOTS would potentially fail "\
                              "due to numerical errors. Would you please consider normalizing the "\
                              "input data properly first?").raise();
            }

            ptc.append(coord);
            tSum.append(t);
            t2Sum.append(t*t);
            for (int k = 0; k < dim_; ++k) {
                axSum.append(QVector<double>());
                ax2Sum.append(QVector<double>());
                axTSum.append(QVector<double>());
                axSum[k].append(coord[k]);
                ax2Sum[k].append(coord[k]*coord[k]);
                axTSum[k].append(coord[k]*t);
            }

            // Setup the initial vK set {0}.
            vK.append(0);
            terminated.append(false);
            numTerminated = 0;

            // Setup the following-path for viterbi decoding.
            QVector<int> path;
            path.append(0);
            pathK.append(path);

            // Set input/output queue.
            inputCount = 1;
            outputCount = 0;
            simplifiedIndex.append(0);
        }
        else
        {
            ptc.append(coord);
            int count = tSum.count()-1;
            tSum.append(tSum[count]+t);
            t2Sum.append(t2Sum[count]+t*t);
            for (int k = 0; k < dim_; ++k) {
                axSum[k].append(axSum[k][count]+coord[k]);
                ax2Sum[k].append(ax2Sum[k][count]+coord[k]*coord[k]);
                axTSum[k].append(axTSum[k][count]+coord[k]*t);
            }
        }
        // Initialize issed&parents.
        issed.append(0);
        parents.append(-1);
    }

    inline void feedIndex(int index)
    {
        if (finished)
            DotsException("Feeding data is NOT allowed after the simplifier finished. "\
                          "Suggest calling resetInternalData() first.").raise();
        if (isCascadeRoot)
            DotsException("We can only feed data to the cascade root simplifier. "\
                          "Try feedData() instead.").raise();

        ptIndex.append(index);
        // Update internal data.
        if (ptIndex.count() == 1)
        {
            // Setup the initial vK set {0}.
            vK.append(0);
            terminated.append(false);
            numTerminated = 0;

            // Setup the following-path for viterbi decoding.
            QVector<int> path;
            path.append(0);
            pathK.append(path);

            // Set input/output queue.
            inputCount = 1;
            outputCount = 0;
            simplifiedIndex.append(ptIndex.at(0));
        }
        // Initialize issed&parents.
        issed.append(0);
        parents.append(-1);
    }

    /**
     * @brief readOutputData checks if the simplifier outputs any data after the recent feeds. Output data will be
     * stored in corresponding parameters if returned true.
     * @param x the x value to output.
     * @param y the y value to output.
     * @param t the timestamp to output.
     * @return true if there's output data, false otherwise.
     */
    inline bool readOutputData(double &x, double &y, double &t)
    {
        int index = -1;
        if (!readOutputIndex(index))
            return false;

        QVector<double> coord = pC->at(index);
        t = pT->at(index);
        return true;
    }

    /**
     * @brief readOutputIndex checks if the simplifier outputs any data after the recent feeds. Output data will be
     * stored in corresponding parameters if returned true.
     * @param index the selected index to output.
     * @return true if there's output data, false otherwise.
     */
    inline bool readOutputIndex(int &index)
    {
        // Run DAG search to produce potentially more output data.
        if (!finished && outputCount >= simplifiedIndex.count())
            directedAcyclicGraphSearch();

        // Retrieve one data.
        if (outputCount < simplifiedIndex.count())
        {
            index = ptIndex.at(simplifiedIndex.at(outputCount));
            ++outputCount;
            return true;
        }
        // No output yet.
        return false;
    }

    /**
     * @brief getSimplifiedIndex retrieves index of the i-th point of simplified trajectory.
     * @param i the point number of simplified trajectory to retrieve.
     * @return index of the i-th point.
     */
    inline int getSimplifiedIndex(int i)
    {
        if (i<0 || i>=simplifiedIndex.count())
            DotsException(QString("Index %1 is out of range [0, %2)").arg(i).arg(simplifiedIndex.count())).raise();

        return ptIndex.at(simplifiedIndex.at(i));
    }

    /**
     * @brief finish sets the finish flag for DOTS algorithm. No more data could be feeded after calling this method.
     */
    void finish();

    /**
     * @brief getLssdThreshold retrieves the LSSD threshold used in this simplifier.
     * @return the LSSD threshold used in this simplifier.
     */
    double getLssdThreshold();

    /**
     * @brief getAverageSED gets the average SED error. Note that this method must be called after finish().
     * @return the average SED error.
     */
    double getAverageSED();

    /**
     * @brief getMaxLSSD retrieves the maximum LSSD generated during the simplification. Note that this method must
     * be called after finish().
     * @return the maximum LSSD.
     */
    double getMaxLSSD();

    /**
     * @brief batchDots provides a batched simplification utility by invoking the online DOTS simplifier.
     * @param x
     * @param y
     * @param t
     * @param ox
     * @param oy
     * @param ot
     * @param lssdThreshold
     */
    static void batchDots(const QVector<double> &x, const QVector<double> &y, const QVector<double> &t,
                          QVector<double> &ox, QVector<double> &oy, QVector<double> &ot,
                          double lssdThreshold);

    static void batchDotsByIndex(const QVector<double> &x, const QVector<double> &y, const QVector<double> &t,
                                 QVector<int> &simplifiedIndex, double lssdThreshold);

    /**
     * @brief batchDotsCascade provides a batched simplification utility by invoking the online DOTS simplifier in
     * cascade mode. The cascade mode is not as accurate as normal mode but is much faster. Refer to this method to
     * learn how to use cascaded DOTS online.
     * @param x
     * @param y
     * @param t
     * @param ox
     * @param oy
     * @param ot
     * @param lssdThreshold
     */
    static void batchDotsCascade(const QVector<double> &x, const QVector<double> &y, const QVector<double> &t,
                                 QVector<double> &ox, QVector<double> &oy, QVector<double> &ot,
                                 double lssdThreshold);

    static void batchDotsCascadeByIndex(const QVector<double> &x, const QVector<double> &y, const QVector<double> &t,
                                        QVector<int> &simplifiedIndex, double lssdThreshold);

protected:
    /**
     * @brief directedAcyclicGraphSearch does a DAG search among the feeded spatio-temporal 2D data. The output queue
     * would be updated when appropriate.
     */
    inline void directedAcyclicGraphSearch()
    {
        int numPoints = ptIndex.count();
        if (finished)
        {
            // Construct the DAG completely.
            while (true)
            {
                if (inputCount>=numPoints)
                    break;

                // Move v* points to vL.
                for (int i=inputCount; i<numPoints; ++i)
                {
                    // Update vK if parent was not assigned yet.
                    if (parents.at(i) < 0)
                    {
                        for (int j=0;j<vK.count(); ++j)
                        {
                            int jIndex = vK.at(j);
                            if (!terminated.at(j))
                            {
                                double distance = getLSSD(jIndex,i);
                                if (distance < lssdTh)
                                {
                                    vL.append(i);
                                    issed[i] = issed[jIndex]+distance;
                                    parents[i] = jIndex;
                                    break;
                                }
                                else if (distance > lssdUpperBound)
                                {
                                    terminated[j] = true;
                                    ++numTerminated;

                                    // Check if all vK terminated.
                                    if (needUpdateVK())
                                    {
                                        break;
                                    }
                                }
                            } // if (!terminated.at(j))
                        } // for (int j=0;j<vK.count(); ++j)
                    } // if (parents.at(i) < 0)

                    // Update inputCount for consequent none v* points.
                    if(parents.at(i)>=0 && inputCount == i)
                    {
                        // Update start point of DAG search.
                        ++inputCount;
                    }

                    // Terminate loop early if we need to update vK.
                    if (needUpdateVK())
                        break;
                } // for (int i=inputCount; i<numPoints; ++i)

                // Minimize ISSED.
                minimizeISSED();

                // Force swap vK/vL no matter we need update vK.
                updateVK();

                // Omit viterbi decoding!!
            }

            // Decode the simplified data from back to front.
            QVector<int> temp;
            int idx = numPoints-1;
            int currentIndex = -1;
            if (!simplifiedIndex.empty())
                currentIndex = simplifiedIndex[simplifiedIndex.count()-1];
            while(idx>currentIndex)
            {
                temp.append(idx);
                idx = parents.at(idx);
            }
            for (int k=temp.count()-1; k>=0; --k)
                simplifiedIndex.append(temp.at(k));
        } // if (finished)
        else
        {
            while (true)
            {
                bool vKUpdated = false;
                if (inputCount>=numPoints)
                    break;

                // Move v* points to vL.
                for (int i=inputCount; i<numPoints; ++i)
                {
                    // Update vK if parent was not assigned yet.
                    if (parents.at(i) < 0)
                    {
                        for (int j=0;j<vK.count(); ++j)
                        {
                            int jIndex = vK.at(j);
                            if (!terminated.at(j))
                            {
                                double distance = getLSSD(jIndex,i);
                                if (distance < lssdTh)
                                {
                                    vL.append(i);
                                    issed[i] = issed[jIndex]+distance;
                                    parents[i] = jIndex;

                                    // Check if vL exceeds max size.
                                    if (needUpdateVK())
                                    {
                                        // Minimize ISSED.
                                        minimizeISSED();

                                        // Swap vK/vL
                                        updateVK();

                                        // Viterbi decoding.
                                        viterbiDecode();

                                        // Break search loop.
                                        vKUpdated = true;
                                    }
                                    break;
                                }
                                else if (distance > lssdUpperBound)
                                {
                                    terminated[j] = true;
                                    ++numTerminated;

                                    // Check if all vK terminated.
                                    if (needUpdateVK())
                                    {
                                        // Minimize ISSED.
                                        minimizeISSED();

                                        // Swap vK/vL
                                        updateVK();

                                        // Viterbi decoding.
                                        viterbiDecode();

                                        // Break search loop.
                                        vKUpdated = true;
                                        break;
                                    }
                                }
                            } // if (!terminated.at(j))
                        } // for (int j=0;j<vK.count(); ++j)
                    } // if (parents.at(i) < 0)

                    // Update inputCount for consequent none v* points.
                    if(parents.at(i)>=0 && inputCount == i)
                    {
                        // Update start point of DAG search.
                        ++inputCount;
                    }

                    if (vKUpdated)
                        break;
                } // for (int i=inputCount; i<numPoints; ++i)

                // Stop DAG search if no new vK generated.
                if (!vKUpdated)
                    break;
            }
        }
    }

    /**
     * @brief getLSSD calculates the LSSD (logogram of Local integral Square Synchronous Euclidean Distance ) between
     * two points indexed by fst and lst.
     * @param fst is index of the first point.
     * @param lst is index of the second point.
     * @return the LSSD distance.
     */
    inline double getLSSD(int fst, int lst)
    {
        fst = ptIndex.at(fst);
        lst = ptIndex.at(lst);
        if (fst+1>=lst)
            return 0;
        if (fst<0 || lst>=pT->count())
            DotsException(QString("Index out of bound error.")).raise();

        int plst = lst-1;
        double distance = 0.0;
        for (int a = 0; a < dim_; ++a)
        {
            const QVector<double>& xS  = pAXSum->at(a);
            const QVector<double>& x2  = pAX2Sum->at(a);
            const QVector<double>& xt  = pAXTSum->at(a);
            const QVector<double>& tR  = *pT;
            const QVector<double>& tS  = *pTSum;
            const QVector<double>& t2  = *pT2Sum;

            const double c1 = pC->at(fst).at(a)*tR.at(lst) - pC->at(lst).at(a)*tR.at(fst);
            const double c2 = c1*c1;
            const double c3 = tR.at(lst) - tR.at(fst);
            const double c4 = c3*c3;
            const double c5 = pC->at(lst).at(a) - pC->at(fst).at(a);
            const double c6 = c5*c5;

            distance += (plst-fst)*c2/c4
                    + c6/c4*(t2.at(plst)-t2.at(fst))
                    + (x2.at(plst)-x2.at(fst))
                    + 2*c1*c5/c4*(tS.at(plst)-tS.at(fst))
                    - 2*c1/c3*(xS.at(plst)-xS.at(fst))
                    - 2*c5/c3*(xt.at(plst)-xt.at(fst));
        }
        (void)0;
        return distance;
    }

    /**
     * @brief needUpdateVK Checks if we need to swap Vl and Vk sets.
     * @return true if a swap operation is necessary, false otherwise.
     */
    inline bool needUpdateVK()
    {
        return (vK.count() == numTerminated || vL.count()>= maxVkSize);
    }

    /**
     * @brief updateVK swaps Vl and Vk sets and updates the DAG paths from root node to each elements of current Vk set.
     */
    inline void updateVK()
    {
        // Update following-path.
        QVector< QVector<int> > newPath;
        for (int k=0; k<vL.count(); ++k)
        {
            int indexK = -1;
            int parentPos = parents.at(vL.at(k));
            for (int m=0; m<vK.count(); ++m)
            {
                if (vK.at(m) == parentPos)
                {
                    indexK = m;
                    break;
                }
            }

            QVector<int> p = pathK.at(indexK);
            p.append(vL.at(k));
            newPath.append(p);
        }
        pathK = newPath;

        // Update vK set.
        vK = vL;
        terminated.resize(vK.count());
        for (int k=0; k<terminated.count(); ++k)
            terminated[k] = false;
        numTerminated = 0;
        vL.clear();
    }

    /**
     * @brief minimizeISSED minimizes the total error from root node to each element of Vl set. The minimization is
     * done by choosing the best parents of Vl elements among Vk elements.
     */
    inline void minimizeISSED()
    {
        foreach (int i, vL) {
            double minDistance = issed.at(i);
            double minParent = parents.at(i);
            foreach (int j, vK) {
                double localDistance = getLSSD(j,i);
                double distance = issed.at(j) + localDistance;
                if (localDistance<lssdTh && distance<minDistance)
                {
                    minDistance = distance;
                    minParent = j;
                }
            }
            issed[i] = minDistance;
            parents[i] = minParent;
        }
    }

    /**
     * @brief viterbiDecode decodes output indices of simplified points in a viterbi-like manner.
     */
    inline void viterbiDecode()
    {
        int inputLayer = pathK.at(0).count();
        int startLayer = simplifiedIndex.count();
        while (startLayer < inputLayer) {
            bool equal = true;
            int reference = pathK.at(0).at(startLayer);
            foreach (QVector<int> path, pathK) {
                if (path.at(startLayer) != reference)
                {
                    equal = false;
                    break;
                }
            }
            if (!equal)
                break;

            simplifiedIndex.append(reference);
            ++startLayer;
        }
    }

protected:
    // DOTS settings.
    /**
     * @brief lssdTh
     */
    double lssdTh;
    double lssdUpperBound;
    int maxVkSize;

    // Input sequence.
    QVector< QVector<double> > ptc;   // one d-vector of coordinates per point
    QVector<double> ptt;
    int dim_;
    QVector<int> ptIndex;
    QVector< QVector<double> > *pC;
    QVector<double> *pT;

    // DOTS algorithm internal data.
    QVector< QVector<double> > axSum, ax2Sum, axTSum;   // prefix sums, one array per axis
    QVector<double> tSum, t2Sum;
    QVector< QVector<double> > *pAXSum, *pAX2Sum, *pAXTSum;
    QVector<double> *pTSum, *pT2Sum;
    QVector<double> vK,vL;
    QVector<bool> terminated;
    int numTerminated;
    QVector< QVector<int> > pathK;
    QVector<double> issed;
    QVector<int> parents;

    // Output sequence.
    QVector<int> simplifiedIndex;
    int inputCount;
    int outputCount;

    // Indicates if the input got EOF. No more data could be input after this flag was set.
    bool finished;

    bool isCascadeRoot;

signals:

public slots:
};

#endif // DOTSSIMPLIFIERD_H
