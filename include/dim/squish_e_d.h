#pragma once
/* SQUISH-E for d-dimensional points.

   A line-for-line counterpart of SQUISH_E.cpp's SQUISH_E / adjust_priority / reduce /
   find_mini_priority. The ONE dimension-dependent call in the whole algorithm is
   cacl_SED, which here becomes dim::sed (already verified bit-identical to cacl_SED at
   d = 2 over 500,000 random cases).

   Everything else -- the growing capacity, the neighbour priority recomputation, the
   accumulation of `pi`, the removal of the minimum -- is scalar bookkeeping and is copied
   unchanged. src/check_d_published_equiv.cpp requires this to select the same points as
   the published routine at d = 2.
*/
#include <vector>
#include <limits>
#include <algorithm>
#include "metrics_d.h"

namespace dim {

struct SquishEItem {
    double priority;
    double pi;
    VertexD point;
    int    index;
};

class SquishE {
public:
    /* Returns the indices of the points kept, in order. */
    std::vector<int> run(const std::vector<VertexD>& pts, int cmp_ratio, double sed_error) {
        Q_.clear();
        capacity_ = 4;
        const int n = (int)pts.size();

        for (int i = 0; i < n; ++i) {
            if (i / cmp_ratio >= capacity_) capacity_ += 1;

            SquishEItem it;
            it.priority = std::numeric_limits<double>::max();
            it.pi       = 0.0;
            it.point    = pts[i];
            it.index    = i;
            Q_.push_back(it);

            if (i > 0) adjustPriority((int)Q_.size() - 3, (int)Q_.size() - 2, (int)Q_.size() - 1);

            if ((int)Q_.size() == capacity_) {
                const int mi = findMinPriority();
                const double mp = Q_[mi].priority;
                reduce(mi, mp);
            }
        }

        int mi = findMinPriority();
        double mp = Q_[mi].priority;
        while (mp <= sed_error) {
            reduce(mi, mp);
            mi = findMinPriority();
            mp = Q_[mi].priority;
        }

        std::vector<int> kept;
        kept.reserve(Q_.size());
        for (const SquishEItem& q : Q_) kept.push_back(q.index);
        return kept;
    }

private:
    std::vector<SquishEItem> Q_;
    int capacity_ = 4;

    void adjustPriority(int pre, int cur, int succ) {
        if (cur == (int)Q_.size() - 1 || cur == 0) return;
        const double p = Q_[cur].pi + sed(Q_[pre].point, Q_[cur].point, Q_[succ].point);
        Q_[cur].priority = p;
    }

    void reduce(int min_index, double min_p) {
        Q_[min_index - 1].pi = std::max(min_p, Q_[min_index - 1].pi);
        Q_[min_index + 1].pi = std::max(min_p, Q_[min_index + 1].pi);
        adjustPriority(min_index - 2, min_index - 1, min_index + 1);
        adjustPriority(min_index - 1, min_index + 1, min_index + 2);
        Q_.erase(Q_.begin() + min_index);
    }

    int findMinPriority() {
        int min_index = 1;
        int k = 1;
        const int sz = (int)Q_.size();
        for (int i = 1; i < sz - 1; ++i) {
            if (Q_[i].priority < Q_[min_index].priority) min_index = i;
            ++k;
        }
        return min_index;
    }
};

} /* namespace dim */
