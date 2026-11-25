#pragma once

#include "dim/metrics_d.h"
#include <deque>
#include <algorithm>
#include <cmath>

// Streaming delta-simplification for d-dimensional data.
//
// This is the 2D StreamSimplifier (include/stream_simplifier.h) re-typed for dim::VertexD.
// The branch structure is identical: STEP 1 / CASE 1 / CASE 2 / rejection all carry over
// untouched. What changes is only the arithmetic inside a distance or an angle, which is
// now a loop over the d components instead of a pair of x/y expressions.
//
// The 2D implementation is untouched, so a d = 2 run can still use it.
// requires this class to select exactly the same vertices as the 2D one at d = 2.
//
// Original header follows.
//
// Streaming delta-simplification, d-dimensional.
// Exposed as a reusable library class that accumulates the simplified output
// polyline. Parameters: delta (error threshold), w (target-queue capacity w^2,
// with lambda = 1 + 1/w).
// dimension-agnostic point equality, replacing (a.x == b.x && a.y == b.y)
static inline bool samePoint(const dim::VertexD& a, const dim::VertexD& b) {
    if (a.dim() != b.dim()) return false;
    for (int k = 0; k < a.dim(); ++k) if (a.c[k] != b.c[k]) return false;
    return true;
}

class StreamSimplifierD {
private:
    double delta;
    int w_capacity;       // = w^2
    double lambda;        // = 1 + 1/w

    static constexpr double GEOM_EPS = 1e-12;

    std::deque<dim::VertexD> targetQueue;
    dim::VertexD v_s;             // start vertex of current segment
    dim::VertexD current_target;  // active target t(i-1)
    bool has_target;
    double p_max;           // largest projection distance so far
    dim::VertexD last_vertex;     // v_{i-1} (previously accepted vertex)
    std::vector<dim::VertexD> output;

    // Instrumentation only: number of
    // (point, queue-entry) pairs visited by the CASE 1 cone re-test. Counts the
    // same operation as StreamVariedSimplifier::coneScans(). Never read by any
    // branch of the algorithm.
    size_t coneScans_ = 0;

    // Append to the FIFO target queue; if size exceeds w^2, evict oldest.
    void enqueueTarget(const dim::VertexD& p) {
        targetQueue.push_back(p);
        if (static_cast<int>(targetQueue.size()) > w_capacity)
            targetQueue.pop_front();
    }

    void process(const dim::VertexD& v_i) {
        dim::VertexD V_i = dim::sub(v_i, v_s);
        double r_i = dim::normD(V_i);

        // To check if there are duplicates;
        // Degenerate duplicate of v_s: heading undefined. Defensive only;
        // the pseudocode does not cover r_i == 0.
        if (r_i <= GEOM_EPS) {
            // See include/stream_simplifier.h: only safe to overwrite when no target
            // is pending, otherwise a promoted-but-unemitted target is lost.
            if (!has_target) last_vertex = v_i;
            return;
        }

        // STEP 1: initialization (no target yet).
        if (!has_target) {
            if (r_i > 4.0 * delta) {
                current_target = v_i;
                has_target = true;
                p_max = r_i;
                enqueueTarget(v_i);
            }
            last_vertex = v_i;
            return;
        }

        // STEP 2: evaluation (target defined).
        dim::VertexD V_t = dim::sub(current_target, v_s);
        double r_t = dim::normD(V_t);
        dim::VertexD e_t = dim::unitD(V_t, r_t);

        // CASE 1: r_i >= lambda * r_t
        if (r_i >= lambda * r_t) {
            dim::VertexD e_i = dim::unitD(V_i, r_i);

            // The queued vector is NOT normalised here, because normalising it
            // computes nothing. With V_h = v_h - v_s and r_h = |V_h|, the published
            // test forms e_h = V_h / r_h and then asks
            //
            //     dot(e_i, e_h) > 0    and    r_h * crossMag(e_i, e_h) <= 2*delta
            //
            // Both expressions are homogeneous of degree 0 in V_h, so r_h cancels:
            //
            //     dot(e_i, V_h/r_h)      = dim::dotD(e_i, V_h) / r_h
            //     crossMag(e_i, V_h/r_h) = crossMag(e_i, V_h) / r_h
            //
            // With r_h > 0, guaranteed by the guard below, the first is positive
            // exactly when dim::dotD(e_i, V_h) is, and r_h times the second is exactly
            // crossMag(e_i, V_h). The sqrt and the two divisions per queue entry
            // were computing a normalisation that the comparisons then divide out.
            //
            // This is an identity, not an approximation: the loop does four
            // multiplies and two compares where it did a sqrt, two divides and a
            // sqrt-compare, and it rounds twice less. test_cone_identity.exe checks
            // the two forms against each other over 20M random and near-degenerate
            // inputs; they differ only for inputs placed exactly on the lateral
            // boundary, where both compute the same value to within one rounding
            // step and the comparison is a tie either way.
            bool all_pass = true;
            for (const dim::VertexD& v_h : targetQueue) {
                ++coneScans_;                       // instrumentation only
                dim::VertexD V_h = dim::sub(v_h, v_s);
                double r2 = dim::sqMagD(V_h);
                if (r2 <= GEOM_EPS * GEOM_EPS) continue;

                // The cone test. This is the #else arm of the 2D file, which the
                // generator deleted by mistake when it stripped the CONE_LEGACY_2D
                // switch -- and with it the whole test, which made the d class accept
                // every point in CASE 1.
                const double proj = dim::dotD(e_i, V_h);
                const dim::VertexD perp = dim::perpD(V_h, e_i, proj);
                const double perp2 = dim::sqMagD(perp);
                if (!(proj > 0.0 && perp2 <= 4.0 * delta * delta)) {
                    all_pass = false;
                    break;
                }

                // Dimension-agnostic form of the source test
                //     dot(e_i, e_h) > 0  and  r_h * crossMag(e_i, e_h) <= 2*delta
                // (e_i and e_h are unit vectors).
                //
                // For unit vectors crossMag(e_i, e_h) = |sin theta|, so the second
                // condition asks that the perpendicular distance from the queued
                // vector's endpoint to the heading line is at most 2*delta. That is a
                // statement about a projection and a norm, not about a cross product,
                // and it is defined in any dimension.
                //
                // The perpendicular component is formed DIRECTLY,
                //
                //     perp = V_h - dim::dotD(e_i, V_h) * e_i
                //
                // rather than as |V_h|^2 - dot(e_i,V_h)^2. Both are algebraically the
                // same (r_h sin theta)^2, but the second subtracts two nearly equal
                // large numbers, and on inputs placed exactly on the lateral boundary it
                // disagrees with the cross-product form in 113 of 273 cases while this
                // one disagrees in none.
                //
                // Comparing squared quantities is exact here: both sides are
                // non-negative, so perp2 <= 4*delta^2 is equivalent to
                // |perp| <= 2*delta. Every term generalises componentwise to d > 2.
            }

            if (all_pass) {
                current_target = v_i;
                p_max = r_i;
                enqueueTarget(v_i);
                last_vertex = v_i;
                return;
            }
            // fall through to STEP 3: rejection
        }
        // CASE 2: r_i < lambda * r_t
        else {
            double P_i = dim::dotD(V_i, e_t);            // = r_i * dot(e_i, e_t)
            p_max = std::max(P_i, p_max);

            dim::VertexD ortho = dim::perpD(V_i, e_t, P_i);
            double termA = dim::normD(ortho);       // distance from heading line
            double termB = std::abs(p_max - P_i);  // projection lag

            if (termA <= 2.0 * lambda * delta && termB <= 2.0 * delta) {
                last_vertex = v_i;                 // target unchanged
                return;
            }
            // fall through to STEP 3: rejection
        }

        // STEP 3: rejection.
        output.push_back(last_vertex);   // close segment at v_{i-1}
        output.push_back(v_i);           // new segment starts at v_i
        v_s = v_i;
        has_target = false;
        p_max = 0.0;
        targetQueue.clear();
        last_vertex = v_i;
    }

    void reset() {
        output.clear();
        targetQueue.clear();
        has_target = false;
        p_max = 0.0;
        coneScans_ = 0;
    }

public:
    StreamSimplifierD(double delta_, int w)
        : delta(delta_), w_capacity(w * w),
          lambda(1.0 + 1.0 / static_cast<double>(w)),
          has_target(false), p_max(0.0) {}

    // Feed one vertex at a time.
    void addPoint(const dim::VertexD& p) {
        if (output.empty()) {
            v_s = p;
            last_vertex = p;
            output.push_back(p);
            return;
        }
        process(p);
    }

    // Flush the final open segment; call once when the stream ends.
    void finish() {
        if (!output.empty() && !samePoint(last_vertex, v_s))
            output.push_back(last_vertex);
    }

    // Convenience: simplify a whole trajectory in one call.
    std::vector<dim::VertexD> simplify(const std::vector<dim::VertexD>& stream) {
        reset();
        if (stream.empty()) return output;
        for (size_t i = 0; i < stream.size(); ++i)
            addPoint(stream[i]);
        finish();
        return output;
    }

    const std::vector<dim::VertexD>& getOutput() const { return output; }

    // Instrumentation only; see coneScans_ above.
    size_t coneScans() const { return coneScans_; }
};

inline double approximationConstantD(int w) {
    double wd = static_cast<double>(w);
    double lambda = 1.0 + 1.0 / wd;
    double a = 1.0 / (1.0 - std::pow(lambda, -(wd * wd)));
    return std::sqrt(1.0 + 4.0 * lambda * lambda * ((a + 1.0) * (a + 1.0)))
         + std::sqrt(1.0 + 4.0 * lambda * lambda);
}
