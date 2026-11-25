#pragma once
/* Continuous Frechet distance for d-dimensional curves.

   A d-dimensional counterpart of include/fretchet.h, which operates on the 2D `Vertex`.
   The 2D file is not touched, so the 2D results stay reproducible from the existing
   binaries and this header can be adopted only where d > 2 needs it.

   ONLY THREE FUNCTIONS IN THE 2D IMPLEMENTATION ARE DIMENSION-SPECIFIC:

       squaredDist(a, b)              dx, dy
       freeInterval(a, b, c, delta)   vx/vy, wx/wy  ->  A, B, C
       evalCurve(C, t)                interpolates x and y

   Everything else -- the free-space propagation, the reachability queues, the binary
   search -- manipulates scalars and intervals only, so it is copied here unchanged.

   The free-interval derivation is the same in any dimension. With v = b - a and
   w = a - c, solving |a + u v - c|^2 = delta^2 for u gives

       A u^2 + B u + C = 0,   A = |v|^2,   B = 2 (w . v),   C = |w|^2 - delta^2

   so the generalisation is a dot product and a squared norm, both of which are sums of
   d component terms. At d = 2 the arithmetic is identical operation for operation.
*/
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include "metrics_d.h"

namespace dim {

struct Interval {
    double lo, hi;
    bool empty() const { return lo > hi; }
};
inline Interval emptyInterval() { return Interval{1.0, 0.0}; }

/* squared distance, summed over d components */
inline double sqDist(const VertexD& a, const VertexD& b) {
    double s = 0.0;
    for (int k = 0; k < a.dim(); ++k) {
        const double d = a.c[k] - b.c[k];
        s += d * d;
    }
    return s;
}
inline double distV(const VertexD& a, const VertexD& b) { return std::sqrt(sqDist(a, b)); }

/* Parameter interval u in [0,1] along segment a->b such that dist(a + u*(b-a), c) <= delta. */
inline Interval freeInterval(const VertexD& a, const VertexD& b,
                             const VertexD& c, double delta) {
    const int d = a.dim();
    double A = 0.0, B = 0.0, C = 0.0;
    for (int k = 0; k < d; ++k) {
        const double v = b.c[k] - a.c[k];      /* v = b - a */
        const double w = a.c[k] - c.c[k];      /* w = a - c */
        A += v * v;
        B += 2.0 * w * v;
        C += w * w;
    }
    C -= delta * delta;

    if (A == 0.0) {                            /* degenerate segment (a == b) */
        return (C <= 0.0) ? Interval{0.0, 1.0} : emptyInterval();
    }
    double D = B * B - 4.0 * A * C;
    if (D < 0.0) {
        const double scale = std::max(1.0, std::max(std::abs(B * B), std::abs(4.0 * A * C)));
        if (D > -1e-12 * scale) D = 0.0;       /* tangent rounded below zero */
        else return emptyInterval();
    }
    const double s = std::sqrt(D);
    double u1 = (-B - s) / (2.0 * A);
    double u2 = (-B + s) / (2.0 * A);
    if (u1 > u2) std::swap(u1, u2);
    const double lo = std::max(0.0, u1);
    const double hi = std::min(1.0, u2);
    if (lo > hi) return emptyInterval();
    return Interval{lo, hi};
}

/* Reachability recurrence for one free-space cell (Buchin et al. Sec. 6.1). */
inline Interval propagateRight(const Interval& L, const Interval& Bt, const Interval& freeR) {
    if (freeR.empty()) return emptyInterval();
    if (!Bt.empty()) return freeR;
    if (!L.empty()) {
        const double lo = std::max(freeR.lo, L.lo);
        return (lo <= freeR.hi) ? Interval{lo, freeR.hi} : emptyInterval();
    }
    return emptyInterval();
}
inline Interval propagateTop(const Interval& L, const Interval& Bt, const Interval& freeT) {
    if (freeT.empty()) return emptyInterval();
    if (!L.empty()) return freeT;
    if (!Bt.empty()) {
        const double lo = std::max(freeT.lo, Bt.lo);
        return (lo <= freeT.hi) ? Interval{lo, freeT.hi} : emptyInterval();
    }
    return emptyInterval();
}

/* Index-parameter evaluation, interpolating every coordinate. */
inline VertexD evalCurve(const std::vector<VertexD>& C, double t) {
    const int m = static_cast<int>(C.size()) - 1;
    if (t <= 0.0) return C[0];
    if (t >= static_cast<double>(m)) return C[m];
    const int i = static_cast<int>(std::floor(t));
    const double f = t - static_cast<double>(i);
    VertexD r(C[i].dim());
    for (int k = 0; k < C[i].dim(); ++k)
        r.c[k] = C[i].c[k] + f * (C[i + 1].c[k] - C[i].c[k]);
    return r;
}

/* Equal-time distance, an admissible matching and therefore an upper bound. O(m+n). */
inline double equalTimeDistance(const std::vector<VertexD>& P,
                                const std::vector<VertexD>& Q) {
    const int m = static_cast<int>(P.size()) - 1;
    const int n = static_cast<int>(Q.size()) - 1;
    double best2 = sqDist(P[0], Q[0]);

    int i = 0, j = 0;
    while (i < m || j < n) {
        double tP, tQ;
        if (i < m && j < n) {
            const long long lhs = static_cast<long long>(i + 1) * n;
            const long long rhs = static_cast<long long>(j + 1) * m;
            if (lhs < rhs)      { ++i; tP = static_cast<double>(i); tQ = tP * n / m; }
            else if (lhs > rhs) { ++j; tQ = static_cast<double>(j); tP = tQ * m / n; }
            else                { ++i; ++j; tP = static_cast<double>(i); tQ = static_cast<double>(j); }
        } else if (i < m) {
            ++i; tP = static_cast<double>(i); tQ = tP * n / m;
        } else {
            ++j; tQ = static_cast<double>(j); tP = tQ * m / n;
        }
        const double d2 = sqDist(evalCurve(P, tP), evalCurve(Q, tQ));
        if (d2 > best2) best2 = d2;
    }
    return std::sqrt(best2);
}

/* Decision oracle: is d_F(P, Q) <= delta? */
inline bool frechetDecision(const std::vector<VertexD>& P,
                            const std::vector<VertexD>& Q, double delta) {
    const int m = static_cast<int>(P.size()) - 1;
    const int n = static_cast<int>(Q.size()) - 1;
    const double d2 = delta * delta;

    if (m == 0 && n == 0) return sqDist(P[0], Q[0]) <= d2;
    if (m == 0) { for (const auto& q : Q) if (sqDist(P[0], q) > d2) return false; return true; }
    if (n == 0) { for (const auto& p : P) if (sqDist(Q[0], p) > d2) return false; return true; }

    if (sqDist(P[0], Q[0]) > d2) return false;
    if (sqDist(P[m], Q[n]) > d2) return false;

    const std::vector<VertexD>* A = &P;
    const std::vector<VertexD>* B = &Q;
    int ma = m, nb = n;
    if (ma > nb) { std::swap(A, B); std::swap(ma, nb); }

    struct QueueItem { int row; Interval reach; };
    std::vector<QueueItem> queues[2];
    queues[0].reserve(ma + 1);
    queues[1].reserve(ma + 1);

    int prev = 0, next = 1;
    queues[prev].push_back({1, Interval{0.0, 0.0}});

    const double tol = 1e-9;

    for (int j = 1; j <= nb; ++j) {
        queues[next].clear();
        const VertexD& qj = (*B)[j];
        const VertexD& qjm1 = (*B)[j - 1];

        std::size_t head = 0;
        while (head < queues[prev].size()) {
            int i = queues[prev][head].row;
            Interval Bc = emptyInterval();
            while (true) {
                Interval L = emptyInterval();
                if (head < queues[prev].size() && queues[prev][head].row == i) {
                    L = queues[prev][head].reach;
                    ++head;
                }
                const Interval freeR = freeInterval((*A)[i - 1], (*A)[i], qj, delta);
                const Interval freeT = freeInterval(qjm1, qj, (*A)[i], delta);

                const Interval R = propagateRight(L, Bc, freeR);
                const Interval T = propagateTop(L, Bc, freeT);

                if (!R.empty()) queues[next].push_back({i, R});
                Bc = T;
                ++i;
                if (Bc.empty() || i > ma) break;
            }
        }
        if (queues[next].empty()) return false;
        std::swap(prev, next);
    }

    for (const auto& item : queues[prev])
        if (item.row == ma && item.reach.hi >= 1.0 - tol) return true;
    return false;
}

/* Continuous Frechet distance by binary search over the decision oracle. */
inline double continuousFrechet(const std::vector<VertexD>& P,
                                const std::vector<VertexD>& Q) {
    const int m = static_cast<int>(P.size()) - 1;
    const int n = static_cast<int>(Q.size()) - 1;

    if (m == 0 && n == 0) return distV(P[0], Q[0]);
    if (m == 0) { double d = 0.0; for (const auto& q : Q) d = std::max(d, distV(P[0], q)); return d; }
    if (n == 0) { double d = 0.0; for (const auto& p : P) d = std::max(d, distV(Q[0], p)); return d; }

    double lower = std::max(distV(P[0], Q[0]), distV(P[m], Q[n]));
    double upper = equalTimeDistance(P, Q);
    if (upper <= lower + 1e-12) return upper;

    const double absTol = 1e-12;   // was 1e-6: too coarse to resolve near-zero errors
    const double relTol = 1e-9;
    const int maxIter = 100;

    for (int it = 0; it < maxIter; ++it) {
        const double width = upper - lower;
        const double stopTol = absTol + relTol * std::max(1.0, std::abs(upper));
        if (width <= stopTol) break;
        const double mid = 0.5 * (lower + upper);
        if (mid <= lower || mid >= upper) break;
        if (frechetDecision(P, Q, mid)) upper = mid;
        else                            lower = mid;
    }
    return upper;
}

} /* namespace dim */
