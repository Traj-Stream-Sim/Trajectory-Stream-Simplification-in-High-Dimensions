/* Definitions for the d-dimensional port of the author's DOTS.

   published/DotsSimplifier.cpp cannot be reused directly: it declares methods on
   DotsSimplifier and manipulates the 2D arrays ptx/pty/xSum/ySum/... Each function below
   is the published body with the coordinate handling replaced by the per-axis arrays.
   The logic is unchanged.
*/
#pragma once
#include "dots_published_d.h"

inline DotsSimplifierD::DotsSimplifierD(QObject *parent, DotsSimplifierD *cascadeRoot)
    : QObject(parent)
{
    lssdTh = 10000.0;
    lssdUpperBound = lssdTh * 2.0;
    maxVkSize = 1e6;
    dim_ = 2;
    resetInternalData();

    isCascadeRoot = (cascadeRoot == NULL);
    if (isCascadeRoot) cascadeRoot = this;
    pC     = &(cascadeRoot->ptc);
    pT     = &(cascadeRoot->ptt);
    pTSum  = &(cascadeRoot->tSum);
    pT2Sum = &(cascadeRoot->t2Sum);
    pAXSum  = &(cascadeRoot->axSum);
    pAX2Sum = &(cascadeRoot->ax2Sum);
    pAXTSum = &(cascadeRoot->axTSum);
}

inline void DotsSimplifierD::setParameters(double lssdTh_, double k, int maxVkSize_)
{
    this->lssdTh = lssdTh_;
    lssdUpperBound = lssdTh_ * k;
    this->maxVkSize = maxVkSize_;
}

inline void DotsSimplifierD::resetInternalData()
{
    ptc.clear();
    ptt.clear();
    ptIndex.clear();
    axSum.clear();
    ax2Sum.clear();
    axTSum.clear();
    tSum.clear();
    t2Sum.clear();

    vK.clear();
    vL.clear();
    terminated.clear();
    numTerminated = 0;
    pathK.clear();
    issed.clear();
    parents.clear();

    simplifiedIndex.clear();
    inputCount = 0;
    outputCount = 0;

    finished = false;
}

inline void DotsSimplifierD::finish()
{
    if (!finished)
    {
        finished = true;
        directedAcyclicGraphSearch();
    }
}

inline double DotsSimplifierD::getLssdThreshold() { return lssdTh; }

/* Convenience driver: feed all points, finish, return the kept indices. */
inline void DotsSimplifierD::batchDotsByIndexD(
        const QVector< QVector<double> > &coord, const QVector<double> &t,
        QVector<int> &simplifiedIndexOut, double lssdThreshold)
{
    DotsSimplifierD simplifier;
    simplifier.setParameters(lssdThreshold);
    simplifiedIndexOut.clear();
    for (int i = 0; i < coord.count(); ++i)
    {
        simplifier.feedData(coord[i], t[i]);
        int idx;
        if (simplifier.readOutputIndex(idx))
            simplifiedIndexOut.append(idx);
    }
    simplifier.finish();
    int idx;
    while (simplifier.readOutputIndex(idx))
        simplifiedIndexOut.append(idx);
}
