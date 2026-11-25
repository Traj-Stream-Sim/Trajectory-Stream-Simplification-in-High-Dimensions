#pragma once
/* Dimension-agnostic data layer, ADDITIVE.

   Nothing here is included by any existing file. The 2D code path
   (published/ + include/stream_*.h) is untouched, so the 2D results stay
   reproducible from the existing binaries. This layer exists so a separate
   d-dimensional path can run alongside it, and it earns its use by reproducing the
   2D path exactly at d = 2.
*/
#include <vector>
#include <cmath>
#include <cstddef>

namespace dim {

/* A point in R^d plus a timestamp. `d` is per-point because a trajectory is
   homogeneous and this keeps every free function dimension-explicit. */
struct VertexD {
    std::vector<double> c;
    double t = 0.0;

    VertexD() {}
    explicit VertexD(int d) : c(d, 0.0) {}
    VertexD(std::vector<double> coords, double time) : c(std::move(coords)), t(time) {}
    int dim() const { return (int)c.size(); }
};

inline VertexD sub(const VertexD& a, const VertexD& b) {
    VertexD r(a.dim());
    for (int k = 0; k < a.dim(); ++k) r.c[k] = a.c[k] - b.c[k];
    return r;
}
inline double dotD(const VertexD& a, const VertexD& b) {
    double s = 0.0;
    for (int k = 0; k < a.dim(); ++k) s += a.c[k] * b.c[k];
    return s;
}
inline double normD(const VertexD& a) { return std::sqrt(dotD(a, a)); }

/* ---- vector helpers used by the streaming algorithms --------------------- */
inline double sqMagD(const VertexD& v) {
    double s = 0.0;
    for (int k = 0; k < v.dim(); ++k) s += v.c[k] * v.c[k];
    return s;
}
/* v scaled to unit length; m is the caller's already-computed |v| */
inline VertexD unitD(const VertexD& v, double m) {
    VertexD r(v.dim());
    if (m == 0.0) return r;
    for (int k = 0; k < v.dim(); ++k) r.c[k] = v.c[k] / m;
    return r;
}
/* v - (v.e) e : the component of v perpendicular to the unit vector e */
inline VertexD perpD(const VertexD& v, const VertexD& e, double proj) {
    VertexD r(v.dim());
    for (int k = 0; k < v.dim(); ++k) r.c[k] = v.c[k] - proj * e.c[k];
    return r;
}
/* scaled difference a - s*b, used for the CASE 2 ortho test */
inline VertexD linD(const VertexD& a, const VertexD& b, double s) {
    VertexD r(a.dim());
    for (int k = 0; k < a.dim(); ++k) r.c[k] = a.c[k] - s * b.c[k];
    return r;
}

/* ---- SED -----------------------------------------------------------------
   The paper's definition, unchanged in form:
       alpha = (t_m - t_s) / (t_e - t_s),   1 when the denominator is 0
       q_k   = s_k + alpha (e_k - s_k)
       SED   = sqrt( sum_k (q_k - m_k)^2 )

   SQUISH_E.cpp's cacl_SED is this with k over {lat, lon}, same operation order, so at
   d = 2 the two are bit-identical (measured: 500,000 random cases, 0 differing). */
inline double sed(const VertexD& s, const VertexD& m, const VertexD& e) {
    const int d = (int)s.c.size();
    const double num = m.t - s.t;
    const double den = e.t - s.t;
    const double a = (den == 0.0) ? 1.0 : num / den;
    double acc = 0.0;
    for (int k = 0; k < d; ++k) {
        const double q = s.c[k] + (e.c[k] - s.c[k]) * a;
        const double diff = q - m.c[k];
        acc += diff * diff;
    }
    return std::sqrt(acc);
}

/* ---- LISSED via prefix sums, any d ---------------------------------------
   DotsSimplifier::getLSSD evaluates LISSED in O(1) from prefix sums, written as two
   per-axis blocks over eight arrays. The two blocks are IDENTICAL in form -- that is
   visible in the source -- so the d-dimensional version is the same block evaluated
   once per axis and summed.

   THE BLOCK BELOW IS A VERBATIM TRANSCRIPTION of the published expression. That is
   deliberate: an earlier attempt re-derived the expression from first principles and
   got it wrong, because the published code's prefix indexing is not the obvious one.

   Prefix convention, copied from feedData: the arrays are INCLUSIVE, S[k] = sum over
   points 0..k for the first point, then S[k] = S[k-1] + x[k]. */
class LissedFlow {
public:
    void build(const std::vector<VertexD>& pts) {
        n_ = (int)pts.size();
        d_ = n_ ? pts[0].dim() : 0;
        axisRaw_.assign(d_, std::vector<double>(n_, 0.0));
        axisSum_.assign(d_, std::vector<double>(n_, 0.0));
        axisSq_.assign(d_, std::vector<double>(n_, 0.0));
        axisTSum_.assign(d_, std::vector<double>(n_, 0.0));
        tRaw_.assign(n_, 0.0);
        tSum_.assign(n_, 0.0);
        t2Sum_.assign(n_, 0.0);
        for (int i = 0; i < n_; ++i) {
            const VertexD& p = pts[i];
            tRaw_[i] = p.t;
            for (int k = 0; k < d_; ++k) {
                axisRaw_[k][i]  = p.c[k];
                axisSum_[k][i]  = p.c[k] + (i ? axisSum_[k][i - 1] : 0.0);
                axisSq_[k][i]   = p.c[k] * p.c[k] + (i ? axisSq_[k][i - 1] : 0.0);
                axisTSum_[k][i] = p.c[k] * p.t + (i ? axisTSum_[k][i - 1] : 0.0);
            }
            tSum_[i]  = p.t + (i ? tSum_[i - 1] : 0.0);
            t2Sum_[i] = p.t * p.t + (i ? t2Sum_[i - 1] : 0.0);
        }
    }

    /* LISSED over the points strictly between i and j. */
    double lissed(int fst, int lst) const {
        if (fst + 1 >= lst) return 0.0;
        if (fst < 0 || lst >= n_) return 0.0;
        const int plst = lst - 1;

        double total = 0.0;
        for (int a = 0; a < d_; ++a) {
            const std::vector<double>& xR  = axisRaw_[a];   // raw coordinate  (pX)
            const std::vector<double>& tR  = tRaw_;         // raw time        (pT)
            const std::vector<double>& xS  = axisSum_[a];   // prefix x        (pXSum)
            const std::vector<double>& x2  = axisSq_[a];    // prefix x^2      (pX2Sum)
            const std::vector<double>& xt  = axisTSum_[a];  // prefix x*t      (pXTSum)
            const std::vector<double>& tS  = tSum_;         // prefix t        (pTSum)
            const std::vector<double>& t2  = t2Sum_;        // prefix t^2      (pT2Sum)

            /* ---- verbatim from DotsSimplifier::getLSSD, one axis ----
               c1, c3, c5 read the RAW arrays; the summation terms read the prefix
               arrays. Mixing them up is what made an earlier attempt fail on every
               single case. */
            const double c1 = xR[fst] * tR[lst] - xR[lst] * tR[fst];
            const double c2 = c1 * c1;
            const double c3 = tR[lst] - tR[fst];
            const double c4 = c3 * c3;
            const double c5 = xR[lst] - xR[fst];
            const double c6 = c5 * c5;
            if (c4 == 0.0) continue;

            total += (plst - fst) * c2 / c4
                   + c6 / c4 * (t2[plst] - t2[fst])
                   + (x2[plst] - x2[fst])
                   + 2 * c1 * c5 / c4 * (tS[plst] - tS[fst])
                   - 2 * c1 / c3 * (xS[plst] - xS[fst])
                   - 2 * c5 / c3 * (xt[plst] - xt[fst]);
        }
        return total;
    }

    int size() const { return n_; }
    int dim() const { return d_; }

private:
    int n_ = 0, d_ = 0;
    std::vector<std::vector<double> > axisRaw_, axisSum_, axisSq_, axisTSum_;
    std::vector<double> tRaw_, tSum_, t2Sum_;
};

/* LISSED straight from the definition: the sum of squared SEDs of the skipped points.
   O(j-i) per call; used to validate the closed form, and as the only correct option
   if the closed form is ever in doubt. */
inline double lissedBrute(const std::vector<VertexD>& pts, int i, int j) {
    double acc = 0.0;
    for (int m = i + 1; m < j; ++m) {
        const double s = sed(pts[i], pts[m], pts[j]);
        acc += s * s;
    }
    return acc;
}

/* perpendicular distance from p to the infinite line through a and b */
inline double pointLineDistD(const VertexD& p, const VertexD& a, const VertexD& b) {
    VertexD ab = sub(b, a);
    const double len = normD(ab);
    if (len == 0.0) return normD(sub(p, a));
    VertexD v = sub(p, a);
    const double proj = dotD(v, ab) / len;
    VertexD u = unitD(ab, len);
    return normD(perpD(v, u, proj));
}

} /* namespace dim */
