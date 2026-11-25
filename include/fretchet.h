#pragma once
#include "common.h"
#include <vector>
#include <algorithm>
#include <cmath>

// Continuous Frechet distance between two polylines.
//
// Decision problem solved with the queue-based free-space sweep of
// Buchin et al. (2017), "Efficient trajectory queries under the Frechet
// distance (GIS Cup)". Instead of visiting all m*n cells of the Alt-Godau
// free-space diagram, it maintains the reachable portions of each column and
// jumps over rows that cannot contribute. The distance itself is obtained by
// binary search over this decision oracle, with an equal-time upper bound
// used to shrink the initial interval.
//
// Decision: worst case O(m*n), but output-sensitive (skips unreachable rows);
//           O(m+n) transient memory.
// Distance: O(log(UB-LB)/tol) decision calls.

struct Interval {
    double lo, hi;
    bool empty() const { return lo > hi; }
};

static inline Interval emptyInterval() { return Interval{1.0, 0.0}; }

static inline double squaredDist(const Vertex& a, const Vertex& b) {
    double dx = a.x - b.x, dy = a.y - b.y;
    return dx * dx + dy * dy;
}

// Parameter interval u in [0,1] along segment a->b where dist(a + u*(b-a), c)
// is at most delta. Returns [1,0] (empty) if there is no such interval.
static inline Interval freeInterval(const Vertex& a, const Vertex& b,
                                    const Vertex& c, double delta) {
    double vx = b.x - a.x, vy = b.y - a.y;
    double wx = a.x - c.x, wy = a.y - c.y;
    double A = vx * vx + vy * vy;
    double B = 2.0 * (wx * vx + wy * vy);
    double C = wx * wx + wy * wy - delta * delta;

    if (A == 0.0) {            // degenerate segment (a == b)
        return (C <= 0.0) ? Interval{0.0, 1.0} : emptyInterval();
    }
    double D = B * B - 4.0 * A * C;
    if (D < 0.0) {
        double scale = std::max(1.0, std::max(std::abs(B * B), std::abs(4.0 * A * C)));
        if (D > -1e-12 * scale) D = 0.0;   // tangent case rounded below zero
        else return emptyInterval();
    }
    double s = std::sqrt(D);
    double u1 = (-B - s) / (2.0 * A);
    double u2 = (-B + s) / (2.0 * A);
    if (u1 > u2) std::swap(u1, u2);
    double lo = std::max(0.0, u1);
    double hi = std::min(1.0, u2);
    if (lo > hi) return emptyInterval();
    return Interval{lo, hi};
}

// Reachability recurrence for one free-space cell (Buchin et al. Sec. 6.1).
// L/R are P-parameter intervals, B/T are Q-parameter intervals.
static inline Interval propagateRight(const Interval& L, const Interval& B,
                                      const Interval& freeR) {
    if (freeR.empty()) return emptyInterval();
    if (!B.empty()) return freeR;
    if (!L.empty()) {
        double lo = std::max(freeR.lo, L.lo);
        return (lo <= freeR.hi) ? Interval{lo, freeR.hi} : emptyInterval();
    }
    return emptyInterval();
}

static inline Interval propagateTop(const Interval& L, const Interval& B,
                                    const Interval& freeT) {
    if (freeT.empty()) return emptyInterval();
    if (!L.empty()) return freeT;
    if (!B.empty()) {
        double lo = std::max(freeT.lo, B.lo);
        return (lo <= freeT.hi) ? Interval{lo, freeT.hi} : emptyInterval();
    }
    return emptyInterval();
}

// Index-parameter evaluation: C[floor(t)] + frac * (C[floor(t)+1] - C[floor(t)]).
static inline Vertex evalCurve(const std::vector<Vertex>& C, double t) {
    int m = static_cast<int>(C.size()) - 1;
    if (t <= 0.0) return C[0];
    if (t >= static_cast<double>(m)) return C[m];
    int i = static_cast<int>(std::floor(t));
    double f = t - static_cast<double>(i);
    return Vertex{C[i].x + f * (C[i + 1].x - C[i].x),
                  C[i].y + f * (C[i + 1].y - C[i].y)};
}

// Equal-time distance d_ET(P,Q) = max_t ||P(t) - Q(t*n/m)||, an admissible
// matching, hence an upper bound on the Frechet distance. O(m+n) time.
static inline double equalTimeDistance(const std::vector<Vertex>& P,
                                       const std::vector<Vertex>& Q) {
    int m = static_cast<int>(P.size()) - 1;
    int n = static_cast<int>(Q.size()) - 1;
    double best2 = squaredDist(P[0], Q[0]);

    int i = 0, j = 0;
    while (i < m || j < n) {
        double tP, tQ;
        if (i < m && j < n) {
            long long lhs = static_cast<long long>(i + 1) * n;
            long long rhs = static_cast<long long>(j + 1) * m;
            if (lhs < rhs)      { ++i; tP = static_cast<double>(i); tQ = tP * n / m; }
            else if (lhs > rhs) { ++j; tQ = static_cast<double>(j); tP = tQ * m / n; }
            else                { ++i; ++j; tP = static_cast<double>(i); tQ = static_cast<double>(j); }
        } else if (i < m) {
            ++i; tP = static_cast<double>(i); tQ = tP * n / m;
        } else {
            ++j; tQ = static_cast<double>(j); tP = tQ * m / n;
        }
        Vertex p = evalCurve(P, tP);
        Vertex q = evalCurve(Q, tQ);
        double d2 = squaredDist(p, q);
        if (d2 > best2) best2 = d2;
    }
    return std::sqrt(best2);
}

// Decision oracle: is d_F(P,Q) <= delta?
static inline bool frechetDecision(const std::vector<Vertex>& P,
                                   const std::vector<Vertex>& Q, double delta) {
    int m = static_cast<int>(P.size()) - 1;
    int n = static_cast<int>(Q.size()) - 1;
    double d2 = delta * delta;

    // Degenerate curves.
    if (m == 0 && n == 0) return squaredDist(P[0], Q[0]) <= d2;
    if (m == 0) { for (const auto& q : Q) if (squaredDist(P[0], q) > d2) return false; return true; }
    if (n == 0) { for (const auto& p : P) if (squaredDist(Q[0], p) > d2) return false; return true; }

    // Endpoint feasibility filter.
    if (squaredDist(P[0], Q[0]) > d2) return false;
    if (squaredDist(P[m], Q[n]) > d2) return false;

    // Symmetry: keep the curve with fewer edges as the rows.
    const std::vector<Vertex>* A = &P;
    const std::vector<Vertex>* B = &Q;
    int ma = m, nb = n;
    if (ma > nb) { std::swap(A, B); std::swap(ma, nb); }

    struct QueueItem { int row; Interval reach; };
    std::vector<QueueItem> queues[2];
    queues[0].reserve(ma + 1);
    queues[1].reserve(ma + 1);

    int prev = 0, next = 1;
    // Lambda_0 = {(row 1, reach [0,0])}: the start corner.
    queues[prev].push_back({1, Interval{0.0, 0.0}});

    const double tol = 1e-9;

    for (int j = 1; j <= nb; ++j) {
        queues[next].clear();
        const Vertex& qj = (*B)[j];
        const Vertex& qjm1 = (*B)[j - 1];

        size_t head = 0;
        while (head < queues[prev].size()) {
            int i = queues[prev][head].row;
            Interval Bc = emptyInterval();   // bottom interval, carried across rows
            while (true) {
                Interval L = emptyInterval();
                if (head < queues[prev].size() && queues[prev][head].row == i) {
                    L = queues[prev][head].reach;
                    ++head;
                }
                Interval freeR = freeInterval((*A)[i - 1], (*A)[i], qj, delta);
                Interval freeT = freeInterval(qjm1, qj, (*A)[i], delta);

                Interval R = propagateRight(L, Bc, freeR);
                Interval T = propagateTop(L, Bc, freeT);

                if (!R.empty()) queues[next].push_back({i, R});
                Bc = T;
                ++i;
                if (Bc.empty() || i > ma) break;
            }
        }

        if (queues[next].empty()) return false;  // no reachability can be restored
        std::swap(prev, next);
    }

    // Corner (P=m, Q=n) is reachable iff row ma reaches local coordinate 1.
    for (const auto& item : queues[prev]) {
        if (item.row == ma && item.reach.hi >= 1.0 - tol) return true;
    }
    return false;
}

// Continuous Frechet distance, obtained by binary search over the decision
// oracle. Returns a feasible (upper) approximation within tolerance.
inline double continuousFrechet(const std::vector<Vertex>& P,
                                const std::vector<Vertex>& Q) {
    int m = static_cast<int>(P.size()) - 1;
    int n = static_cast<int>(Q.size()) - 1;

    if (m == 0 && n == 0) return dist(P[0], Q[0]);
    if (m == 0) { double d = 0.0; for (const auto& q : Q) d = std::max(d, dist(P[0], q)); return d; }
    if (n == 0) { double d = 0.0; for (const auto& p : P) d = std::max(d, dist(Q[0], p)); return d; }

    double lower = std::max(dist(P[0], Q[0]), dist(P[m], Q[n]));
    double upper = equalTimeDistance(P, Q);
    if (upper <= lower + 1e-12) return upper;

    const double absTol = 1e-12;   // was 1e-6: too coarse to resolve near-zero errors
    const double relTol = 1e-9;
    const int maxIter = 100;

    for (int it = 0; it < maxIter; ++it) {
        double width = upper - lower;
        double stopTol = absTol + relTol * std::max(1.0, std::abs(upper));
        if (width <= stopTol) break;

        double mid = 0.5 * (lower + upper);
        if (mid <= lower || mid >= upper) break;

        if (frechetDecision(P, Q, mid)) upper = mid;
        else                            lower = mid;
    }
    return upper;
}
