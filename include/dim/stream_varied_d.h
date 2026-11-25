#pragma once

// Streaming delta-simplification with an explicit lambda, decoupled from the queue size w,
// and a commit rule that selects the segment endpoint by measured error.
//
// lambda defaults to 1 + 1/w^2, the smallest value in that family that keeps the growth gate
// strictly above 1; the original 1 + 1/w remains selectable through LambdaMode.
//
// The bound Frechet <= c(w)*delta is proved for lambda = 1 + 1/w. Choosing lambda
// independently invalidates that proof, so errors from this variant are measured, not
// guaranteed.
//
// The target queue is O(w^2). Scoring a candidate requires the points it covers, so the
// current segment's input points are buffered; segment length is the only growth term.

#include "dim/metrics_d.h"
#include "dim/frechet_d.h"
#include "fretchet.h"

#include <deque>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <algorithm>

// The policy enums and the settings bundle live at file scope: a nested class
// cannot be used as a default argument inside its own enclosing class.
namespace svd {

// How the segment endpoint is chosen from the candidate set at commit time.
//
// MEASURED (9 tracks, tuned to equal retention): the useful
// way to spend the error is as a CONSTRAINT, not as an objective. Minimising the
// error outright is degenerate -- the smallest error always belongs to the
// shortest segment -- and the error-per-point-saved form mis-spends the budget on
// tracks with a spike. BudgetLatest beat the current reference on BOTH metrics at
// 20/10/5% retention, so it is the default.
enum class Commit {
    Latest,            // baseline: the most recent accepted vertex (stream_simplifier behaviour)
    MinError,          // argmin error( v_s -> candidate )          [the literal "minimum error" rule]
    MinErrorPerPoint,  // argmin error / points-saved                [non-degenerate, but mis-spends]
    BudgetLatest       // longest candidate whose error <= budget    [DEFAULT: max compression inside an error budget]
};

// What "error" means for one candidate chord.
enum class Metric {
    ChordDev,   // max perpendicular deviation of the covered points from the chord
                // (v_s -> candidate). This is Douglas-Peucker's criterion and is
                // O(segment). Units: same as delta when the track is normalised.
    Frechet     // exact continuous Frechet distance between the covered
                // sub-trajectory and the chord. Faithful to the reported metric,
                // O(segment^2 log) per candidate -- study use, not a hot path.
};

// How lambda is derived. lambda is BOTH the target-growth gate (a new target must
// be at least lambda*r_t away) and the lateral allowance (2*lambda*delta).
//
//   OnePlusInvW   1 + 1/w    -- what stream_simplifier.h hardwires. With w = 8 the
//                              lateral allowance is 2.25*delta while the strip and
//                              projection-lag tests are both 2.0*delta, so the
//                              envelope is not uniform.
//   OnePlusInvW2  1 + 1/w^2  -- DEFAULT. One over the QUEUE CAPACITY (the queue is
//                              w^2 targets), and the smallest value in this family
//                              that stays strictly above 1, so the target can never
//                              regress. At w = 8 this is 1.015625, which makes the
//                              lateral allowance 2.03*delta: strip, lateral and lag
//                              then share one ~2*delta envelope.
//   Explicit                  use Settings::lambda verbatim (clamped to >= 1).
enum class LambdaMode { OnePlusInvW, OnePlusInvW2, Explicit };

struct Settings {
    int        w          = 8;                        // queue capacity = w*w
    LambdaMode lambdaMode = LambdaMode::OnePlusInvW2; // see the table above
    double     lambda     = 0.0;                      // > 0 forces Explicit
    Commit     commit     = Commit::BudgetLatest;     // measured best; see the note on Commit
    Metric     metric     = Metric::ChordDev;
    double     budgetK    = 2.0;                      // BudgetLatest: budget = budgetK * delta
    // Evaluate the BudgetLatest candidates longest-first and stop at the first one
    // inside the budget. That is the SAME candidate the ascending scan returns --
    // BudgetLatest asks for the LONGEST candidate within budget, so the first hit
    // walking down is the answer -- but it does not have to score the rest of the
    // candidate set, which at w = 8 is up to 65 chords of up to L points each.
    //
    // The candidate set is already in ascending index order (the queue is a FIFO of
    // increasing indices and lastAcc_ is the newest of them), so no sort and no
    // auxiliary vector are needed either.
    //
    // Output is IDENTICAL over 259 tracks x
    // 10 deltas x every commit rule, and compares the two paths point by point.
    bool       fastCommit = false;                    // false = the original evaluation order
};

} // namespace sv

// The configuration this project measures and reports.
//
// It is the header's design with the two changes from section 10.5 of
// Commit policy:
//
//   commit     = BudgetLatest      the header's measured-best rule
//   lambdaMode = OnePlusInvW       NOT the header's default. Dropping the
//                                  decoupling returns lambda to the original
//                                  algorithm's 1 + 1/w, which is 1.8-2.9x faster
//                                  across the useful band and equal or better on
//                                  both accuracy metrics at 20% retention.
//   fastCommit = true              early-exit candidate scan. Output is
//                                  provably identical; 1.3-1.9x faster.
//
// lambda = 1 + 1/w is also the value the published Frechet <= c(w)*delta proof
// assumes, so this configuration is the one that comes closest to having a
// guarantee as well as the best measured cost.
//
// Settings' own defaults remain exactly what the header specifies
// (BudgetLatest / OnePlusInvW2 / fastCommit = false) so that the file under test
// is unchanged; this is the configuration the harnesses select.
//
// THE SELECTED CONFIGURATION IS K3 ONLY.
//
// `Commit::Latest` returns lastAcc_, and on a rejection lastAcc_ is always idx-1 --
// every path through feed() that does not reject sets lastAcc_ = idx. So Latest
// selects exactly the vertex the ORIGINAL algorithm would have closed at, and
// `OnePlusInvW` is the original's own lambda. This configuration therefore differs
// from the original in exactly one respect: on rejection it re-anchors at the commit
// vertex and replays the tail instead of emitting the rejected point.
//
// The ablation in the report measures why: K3 alone captures essentially all of the
// accuracy gain over the original (98.7% of the best Fréchet reduction available to
// any configuration), while K1 and K2 together cost more runtime than they return on
// the primary metric. At equal output size K3-only is 14-52% cheaper than the
// previous BudgetLatest configuration and matches or beats it on Fréchet.
//
// The one thing given up is the max-SED gain from K2, which is real but modest:
// K3-only retains about 79% of the temporal improvement that the full configuration
// reached. If max SED is the objective rather than Fréchet, set
// `s.commit = svd::Commit::BudgetLatest`; nothing else needs to change.
inline svd::Settings measuredSettingsD(int w = 8) {
    svd::Settings s;
    s.w          = w;
    s.commit     = svd::Commit::Latest;
    s.lambdaMode = svd::LambdaMode::OnePlusInvW;
    s.fastCommit = true;      // inert for Latest; the fast path is BudgetLatest-only
    return s;
}

// dimension-agnostic point equality
static inline bool samePointD(const dim::VertexD& a, const dim::VertexD& b) {
    if (a.dim() != b.dim()) return false;
    for (int k = 0; k < a.dim(); ++k) if (a.c[k] != b.c[k]) return false;
    return true;
}

class StreamVariedSimplifierD {
public:
    using Commit     = svd::Commit;
    using Metric     = svd::Metric;
    using Settings   = svd::Settings;
    using LambdaMode = svd::LambdaMode;

private:
    static constexpr double GEOM_EPS = 1e-12;

    double delta_;
    int    w_;
    double lambda_;
    Commit commit_;
    Metric metric_;
    double budget_;
    bool   fastCommit_;

    // Current segment state. buf_[0] is always v_s, so buf_ indices ARE the
    // original indices of this segment's points.
    std::vector<dim::VertexD> buf_;
    std::deque<size_t>  queue_;        // candidate indices into buf_ (FIFO, capacity w^2)

    dim::VertexD v_s_;
    dim::VertexD target_;
    bool   hasTarget_ = false;
    double pMax_      = 0.0;
    size_t lastAcc_   = 0;             // buf_ index of the most recent accepted vertex

    std::vector<dim::VertexD> out_;
    size_t commits_   = 0;

    // --- instrumentation ------------------------------------------------------
    size_t maxBuf_      = 0;   // peak points buffered for the current segment
    size_t maxQueue_    = 0;   // peak candidate-queue length (<= w^2 by construction)
    size_t replaySteps_ = 0;   // points re-fed after a commit (the recomputation cost)
    mutable size_t errScans_ = 0;  // points visited by candidateError() in total
    // Scratch space for the fast commit path, reused across calls so that the path
    // allocates nothing. Holds (candidate index, its chord error) in descending index
    // order. See chooseCommit().
    mutable std::vector<size_t> candJ_;
    mutable std::vector<double> candE_;
    // Added for this project's complexity harness; read-only, not used by any
    // branch of the algorithm.
    size_t coneScans_   = 0;   // (point, queue-entry) pairs visited by the CASE 1 re-test
    size_t feedCalls_   = 0;   // feed() invocations = n_in + replaySteps_

    int capacity() const { return w_ * w_; }

    void enqueueTarget(size_t idx) {
        queue_.push_back(idx);
        if (queue_.size() > maxQueue_) maxQueue_ = queue_.size();
        while (static_cast<int>(queue_.size()) > capacity()) queue_.pop_front();
    }

    // --- the candidate error -------------------------------------------------
    // Error of closing the current segment at buf_[j]: how badly the original
    // points strictly between v_s and buf_[j] deviate from the straight chord
    // v_s -> buf_[j]. j <= 1 means "nothing covered", error 0.
    double candidateError(size_t j) const {
        if (j <= 1) return 0.0;
        const dim::VertexD& a = buf_[0];
        const dim::VertexD& b = buf_[j];

        if (metric_ == Metric::ChordDev) {
            double worst = 0.0;
            errScans_ += j - 1;                 // instrumentation: O(j) work
            for (size_t k = 1; k < j; ++k)
                worst = std::max(worst, dim::pointLineDistD(buf_[k], a, b));
            return worst;
        }
        // Exact continuous Frechet of the covered sub-trajectory vs the chord.
        errScans_ += j - 1;
        std::vector<dim::VertexD> sub(buf_.begin(), buf_.begin() + static_cast<long>(j) + 1);
        std::vector<dim::VertexD> chord{a, b};
        return dim::continuousFrechet(sub, chord);
    }

    // --- the commit rule -----------------------------------------------------
    size_t chooseCommit(size_t upto) const {
        // FAST PATH for the default rule. One descending pass over the candidate set,
        // returning the first candidate inside the budget -- which, because the set is
        // in descending index order, is the longest one inside it, exactly what
        // BudgetLatest returns.
        //
        // The errors computed on the way are cached, so that when NOTHING is inside the
        // budget the fallback can be evaluated from the cache instead of re-scanning.
        // That matters: the fallback is a full pass, and at lambda = 1 + 1/w it fires
        // often enough that re-scanning makes this path slower than not having it.
        //
        // The fallback's tie-break is order-dependent in the original code -- it starts
        // from the largest candidate and takes STRICT improvements while scanning
        // ascending, so an all-equal set returns the largest, and otherwise the first
        // index attaining the minimum. The cache is in descending order, so replaying
        // it backwards reproduces that exactly; a plain "minimum, ties to the longest"
        // rule does not, and produced 3 differing outputs out of 20,720 checks before
        // this was written this way.
        if (fastCommit_ && commit_ == Commit::BudgetLatest) {
            candJ_.clear();
            candE_.clear();
            if (lastAcc_ >= 1 && lastAcc_ <= upto) {
                candJ_.push_back(lastAcc_);
                candE_.push_back(candidateError(lastAcc_));
                if (candE_.back() <= budget_) return lastAcc_;
            }
            for (auto it = queue_.rbegin(); it != queue_.rend(); ++it) {
                size_t j = *it;
                if (j < 1 || j > upto) continue;
                if (!candJ_.empty() && candJ_.back() == j) continue;   // adjacent duplicate
                candJ_.push_back(j);
                candE_.push_back(candidateError(j));
                if (candE_.back() <= budget_) return j;
            }
            const std::size_t n = candJ_.size();
            if (n == 0) return upto;
            // Nothing fitted. Reproduce the fallback on the cached ascending order.
            size_t fallback = candJ_[0];          // the largest candidate
            double bestErr  = candE_[0];
            for (std::size_t k = n; k-- > 0; ) {  // ascending
                if (candE_[k] < bestErr) { bestErr = candE_[k]; fallback = candJ_[k]; }
            }
            return fallback;
        }

        // Latest reads none of the candidate set, so build it only for the rules that
        // do. The construction below is a heap allocation, a sort of up to w^2
        // elements and a unique, on every commit; for Latest that is pure waste. The
        // early return is output-identical: if the set would have been empty then
        // `lastAcc_` is outside [1, upto] and Latest returns `upto` either way.
        if (commit_ == Commit::Latest)
            return (lastAcc_ >= 1 && lastAcc_ <= upto) ? lastAcc_ : upto;

        // Candidate set: every live queued target plus the latest accepted vertex.
        std::vector<size_t> cands;
        cands.reserve(queue_.size() + 1);
        for (size_t h : queue_) if (h >= 1 && h <= upto) cands.push_back(h);
        if (lastAcc_ >= 1 && lastAcc_ <= upto) cands.push_back(lastAcc_);
        std::sort(cands.begin(), cands.end());
        cands.erase(std::unique(cands.begin(), cands.end()), cands.end());
        if (cands.empty()) return upto;

        switch (commit_) {
        case Commit::Latest: {
            return (lastAcc_ >= 1 && lastAcc_ <= upto) ? lastAcc_ : upto;
        }

        case Commit::MinError: {
            // Walk from the longest candidate down; only a STRICTLY smaller error
            // replaces it, so ties resolve toward more compression.
            //
            // The tolerance is RELATIVE to the incumbent. It used to read
            // `1e-12 * std::max(1.0, std::abs(bestErr))`, and that max() made it an
            // absolute floor of 1e-12 whenever the error was below 1 -- which is every
            // track at these scales, since coordinates are normalised to unit span.
            // An absolute tie-break does not scale with the geometry, so the same
            // track in different units could resolve a tie differently and, through
            // the commit, produce a different segmentation downstream.
            // test_scale_invariance.exe caught it: 156 of 240 comparisons on GeoLife
            // and 224 of 240 on taxi changed under rescaling, where every other
            // algorithm held. Purely relative, the comparison is scale-free.
            size_t best = cands.back();
            double bestErr = candidateError(best);
            for (auto it = cands.rbegin(); it != cands.rend(); ++it) {
                double e = candidateError(*it);
                double tol = 1e-12 * std::abs(bestErr);
                if (e < bestErr - tol) { bestErr = e; best = *it; }
            }
            return best;
        }

        case Commit::MinErrorPerPoint: {
            // score = error per point saved. j == 1 saves nothing -> +inf.
            // Ties (including an all-zero-error straight segment) resolve to the
            // LONGEST candidate, which is what keeps this rule from collapsing
            // into "emit every point".
            double bestScore = 0.0; size_t best = 0; bool have = false;
            for (size_t j : cands) {
                if (j <= 1) continue;
                double score = candidateError(j) / static_cast<double>(j - 1);
                if (!have || score <= bestScore) { bestScore = score; best = j; have = true; }
            }
            if (!have) return (lastAcc_ >= 2) ? lastAcc_ : upto;
            return best;
        }

        case Commit::BudgetLatest: {
            // Longest candidate that respects the error budget; if none does,
            // fall back to the least-error candidate.
            size_t best = 0;
            for (size_t j : cands)
                if (candidateError(j) <= budget_) best = j;      // cands ascending -> last wins
            if (best >= 1) return best;
            size_t fallback = cands.back();
            double bestErr = candidateError(fallback);
            for (size_t j : cands) {
                double e = candidateError(j);
                if (e < bestErr) { bestErr = e; fallback = j; }
            }
            return fallback;
        }
        }
        return upto;
    }

    // Close the current segment at buf_[j] and restart the segment there.
    void commitAt(size_t j) {
        if (j < 1) j = 1;                       // always make progress
        if (j >= buf_.size()) j = buf_.size() - 1;
        out_.push_back(buf_[j]);
        ++commits_;

        buf_.erase(buf_.begin(), buf_.begin() + static_cast<long>(j));
        v_s_       = buf_[0];
        lastAcc_   = 0;
        hasTarget_ = false;
        pMax_      = 0.0;
        queue_.clear();
    }

    // Feed buf_[idx] into the state machine. Returns true if a commit happened,
    // in which case the caller must replay the (shrunken) buffer from index 1.
    bool feed(size_t idx) {
        ++feedCalls_;                       // instrumentation only
        const dim::VertexD& v_i = buf_[idx];
        dim::VertexD V_i = dim::sub(v_i, v_s_);
        double r_i = dim::normD(V_i);

        // See include/stream_varied.h.
        if (r_i <= GEOM_EPS) { if (!hasTarget_) lastAcc_ = idx; return false; }

        // STEP 1: initialization (no target yet).
        if (!hasTarget_) {
            if (r_i > 4.0 * delta_) {
                target_ = v_i; hasTarget_ = true;
                pMax_ = r_i;
                enqueueTarget(idx);
            }
            lastAcc_ = idx;
            return false;
        }

        // STEP 2: evaluation (target defined).
        dim::VertexD V_t = dim::sub(target_, v_s_);
        double r_t = dim::normD(V_t);
        dim::VertexD e_t = dim::unitD(V_t, r_t);

        if (r_i >= lambda_ * r_t) {
            // CASE 1: the candidate is at least lambda*r_t away -- re-test the
            // whole cone against the NEW heading.
            dim::VertexD e_i = dim::unitD(V_i, r_i);
            // Same identity as in stream_simplifier.h: the cone test is homogeneous
            // of degree 0 in V_h, so normalising the queued vector and then
            // multiplying the lateral term by its length again cancels. Written
            // un-normalised the body is four multiplies and two compares instead of
            // a sqrt, two divides and a sqrt-compare. See the derivation there.
            bool allPass = true;
            for (size_t h : queue_) {
                ++coneScans_;               // instrumentation only
                dim::VertexD V_h = dim::sub(buf_[h], v_s_);
                double r2 = dim::sqMagD(V_h);
                if (r2 <= GEOM_EPS * GEOM_EPS) continue;
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
                //     perp = V_h - dot(e_i, V_h) * e_i
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

                const double proj = dim::dotD(e_i, V_h);
                const dim::VertexD perp = dim::perpD(V_h, e_i, proj);
                
                const double perp2 = dim::sqMagD(perp);
                if (!(proj > 0.0 && perp2 <= 4.0 * delta_ * delta_)) {
                    allPass = false;
                    break;
                }

            }
            if (allPass) {
                target_ = v_i; pMax_ = r_i;
                enqueueTarget(idx);
                lastAcc_ = idx;
                return false;
            }
            // fall through: rejection
        } else {
            // CASE 2: absorb, if the point stays inside the heading envelope.
            double P_i = dim::dotD(V_i, e_t);
            pMax_ = std::max(P_i, pMax_);
            dim::VertexD ortho = dim::perpD(V_i, e_t, P_i);
            double termA = dim::normD(ortho);
            double termB = std::abs(pMax_ - P_i);
            if (termA <= 2.0 * lambda_ * delta_ && termB <= 2.0 * delta_) {
                lastAcc_ = idx;
                return false;
            }
            // fall through: rejection
        }

        // STEP 3: rejection -> commit at the candidate the RULE picks.
        commitAt(chooseCommit(idx));
        return true;
    }

    // Close every remaining open segment at end of stream.
    void closeOut() {
        while (buf_.size() >= 2)
            commitAt(chooseCommit(buf_.size() - 1));
        if (!buf_.empty()) {
            const dim::VertexD& last = buf_[0];
            if (out_.empty() || !samePointD(out_.back(), last))
                out_.push_back(last);
        }
    }

public:
    StreamVariedSimplifierD(double delta, const Settings& s = Settings())
        : delta_(delta), w_(s.w),
          lambda_(resolveLambda(s)),
          commit_(s.commit), metric_(s.metric),
          budget_(s.budgetK * delta), fastCommit_(s.fastCommit) {}

    // lambda = 1 + 1/w  (original coupling) | 1 + 1/w^2 (default) | explicit value.
    static double resolveLambda(const Settings& s) {
        double lam;
        if (s.lambda > 0.0)                       lam = s.lambda;
        else if (s.lambdaMode == LambdaMode::OnePlusInvW)
                                                  lam = 1.0 + 1.0 / (double)s.w;
        else                                      lam = 1.0 + 1.0 / ((double)s.w * (double)s.w);
        return lam < 1.0 ? 1.0 : lam;             // < 1 would let the target regress
    }

    const char* lambdaModeName() const {
        if (lambda_ == 1.0 + 1.0 / ((double)w_ * (double)w_)) return "1+1/w^2";
        if (lambda_ == 1.0 + 1.0 / (double)w_)                return "1+1/w";
        return "explicit";
    }

    std::vector<dim::VertexD> simplify(const std::vector<dim::VertexD>& in) {
        buf_.clear(); queue_.clear(); out_.clear();
        hasTarget_ = false; pMax_ = 0.0; lastAcc_ = 0; commits_ = 0;
        maxBuf_ = maxQueue_ = replaySteps_ = errScans_ = 0;
        coneScans_ = feedCalls_ = 0;
        if (in.empty()) return out_;

        buf_.push_back(in[0]);
        v_s_ = in[0];
        out_.push_back(in[0]);

        size_t nextIn = 1;
        size_t cursor = 1;
        for (;;) {
            if (cursor < buf_.size()) {
                if (buf_.size() > maxBuf_) maxBuf_ = buf_.size();
                bool committed = feed(cursor);
                if (committed) {
                    replaySteps_ += buf_.size() - 1;   // recomputation cost
                    cursor = 1;                        // replay after a commit
                } else {
                    cursor = cursor + 1;
                }
                continue;
            }
            if (nextIn >= in.size()) break;
            buf_.push_back(in[nextIn++]);
            cursor = buf_.size() - 1;
        }
        closeOut();
        return out_;
    }

    size_t commits() const { return commits_; }
    double lambda()  const { return lambda_; }

    // --- instrumentation accessors -------------------------------------------
    size_t maxBuffered()   const { return maxBuf_; }       // peak O(segment length)
    size_t maxQueueLen()   const { return maxQueue_; }     // peak <= w^2
    size_t replaySteps()   const { return replaySteps_; }  // re-fed points after commits
    size_t errorScans()    const { return errScans_; }     // total candidateError() work
    size_t queueCapacity() const { return (size_t)capacity(); }
    size_t coneScans()     const { return coneScans_; }    // CASE 1 queue-entry visits
    size_t feedCalls()     const { return feedCalls_; }    // n_in + replaySteps_

    const char* commitName() const {
        switch (commit_) {
        case Commit::Latest:           return "Latest";
        case Commit::MinError:         return "MinError";
        case Commit::MinErrorPerPoint: return "MinErrorPerPoint";
        case Commit::BudgetLatest:     return "BudgetLatest";
        }
        return "?";
    }
    const char* metricName() const {
        return metric_ == Metric::ChordDev ? "ChordDev" : "Frechet";
    }
    std::string label() const {
        char b[192];
        std::snprintf(b, sizeof(b), "%s/%s/l=%.6f(%s)/w=%d",
                      commitName(), metricName(), lambda_, lambdaModeName(), w_);
        return std::string(b);
    }
};

// Bounded-memory alternative (not implemented here).
//
// candidateError() is the only reason the current segment's points are buffered.
// If O(w^2) memory must be a hard guarantee, score each queued candidate
// INCREMENTALLY instead: a candidate's chord (v_s -> v_h) is fixed the moment it
// is enqueued, so its error can only grow, and every new point can update all
// queued candidates in O(w^2) -- the same cost as the cone test that already runs
// per point. The catch is that points which arrived BEFORE the candidate was
// enqueued are then missed, so the score becomes "how well does this candidate
// fit the data that followed it" rather than the error of the whole segment.
