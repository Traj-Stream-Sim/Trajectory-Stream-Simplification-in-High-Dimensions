//  run_bisect.cpp -- exact table cells by bisection instead of interpolation.
//
//  For each (dataset, algorithm, target) the parameter is bisected in log space until the
//  measured value matches the target within 0.5%. Every reported cell is therefore a real
//  run of the algorithm, not a reconstruction between two runs.
//
//  Two modes:
//     --mode=error   target is a Frechet error; solve for the parameter that achieves it,
//                    and report the output size needed.
//     --mode=size    target is a retention fraction; solve for the parameter that achieves
//                    it, and report the Frechet error obtained.
//
//  Preparation, algorithms and measurement are identical to run_unified.cpp. The corpus is
//  loaded and standardised ONCE into memory, which is the main reason this is affordable:
//  run_unified re-reads every file for every parameter setting.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include "qtcompat.h"
#include "dim/metrics_d.h"
#include "dim/frechet_d.h"
#include "dim/stream_d.h"
#include "dim/stream_varied_d.h"
#include "dim/squish_e_d.h"
#include "dim/dots_published_d_impl.h"

namespace fs = std::filesystem;

enum class Layout { TimeFirst, TimeLast };
struct Track { std::vector<dim::VertexD> pts; int d = 0; };

static bool load(const std::string& path, Layout lay, Track& tr) {
    tr.pts.clear();
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<double> v;
        const char* p = line.c_str();
        while (*p) {
            char* e = nullptr;
            const double x = std::strtod(p, &e);
            if (e == p) break;
            v.push_back(x);
            p = e;
            while (*p == ' ' || *p == ',' || *p == '\t') ++p;
        }
        if (v.size() < 2) continue;
        double t;
        if (lay == Layout::TimeFirst) {
            if (v.size() < 3) continue;
            t = v[0]; v.erase(v.begin());
        } else {
            t = v.back(); v.pop_back();
        }
        if (v.empty()) continue;
        tr.pts.push_back(dim::VertexD(v, t));
    }
    tr.d = tr.pts.empty() ? 0 : tr.pts[0].dim();
    return tr.pts.size() > 8;
}

static std::vector<dim::VertexD> standardise(const Track& tr) {
    const int d = tr.d;
    std::vector<double> lo(d), hi(d);
    for (int k = 0; k < d; ++k) lo[k] = hi[k] = tr.pts[0].c[k];
    for (const auto& p : tr.pts)
        for (int k = 0; k < d; ++k) {
            if (p.c[k] < lo[k]) lo[k] = p.c[k];
            if (p.c[k] > hi[k]) hi[k] = p.c[k];
        }
    double span = 0.0;
    for (int k = 0; k < d; ++k) span = std::max(span, hi[k] - lo[k]);
    if (span <= 0.0) span = 1.0;
    std::vector<dim::VertexD> out;
    out.reserve(tr.pts.size());
    for (const auto& p : tr.pts) {
        std::vector<double> c(d);
        for (int k = 0; k < d; ++k) c[k] = (p.c[k] - lo[k]) / span;
        out.push_back(dim::VertexD(c, p.t));
    }
    return out;
}

struct Corpus {
    std::vector<std::vector<dim::VertexD>> tracks;
    double meanIn = 0.0;
    int d = 0;
};

struct Meas { double out = 0, err = 0, us = 0; bool ok = false; };

static Meas evaluate(const Corpus& C, int alg, double par, int w) {
    Meas m;
    long long ntr = 0;
    double sout = 0, serr = 0, cpu = 0;
    double cores = 0;
    for (const auto& in : C.tracks) {
        std::vector<dim::VertexD> out;
        const auto t0 = std::chrono::steady_clock::now();
        if (alg == 0) {
            StreamSimplifierD s(par, w);
            out = s.simplify(in);
        } else if (alg == 1) {
            svd::Settings st = measuredSettingsD(w);
            StreamVariedSimplifierD s(par, st);
            out = s.simplify(in);
        } else if (alg == 2) {
            dim::SquishE se;
            std::vector<int> keep = se.run(in, 1, par);
            out.reserve(keep.size());
            for (int i : keep) out.push_back(in[i]);
        } else {
            QVector< QVector<double> > coord;
            QVector<double> tt;
            for (const auto& v : in) {
                QVector<double> c;
                for (int k = 0; k < v.dim(); ++k) c.append(v.c[k]);
                coord.append(c); tt.append(v.t);
            }
            QVector<int> idx;
            DotsSimplifierD::batchDotsByIndexD(coord, tt, idx, par);
            out.reserve(idx.count());
            for (int i = 0; i < idx.count(); ++i) out.push_back(in[idx[i]]);
        }
        const auto t1 = std::chrono::steady_clock::now();
        cpu += std::chrono::duration<double>(t1 - t0).count();
        sout += (double)out.size();
        if (out.size() >= 2) serr += dim::continuousFrechet(in, out);
        cores += (double)in.size();
        ++ntr;
    }
    if (!ntr) return m;
    m.ok = true;
    m.out = sout / (double)ntr;
    m.err = serr / (double)ntr;
    m.us = cpu * 1e6 / cores;
    return m;
}

/* Bisect log(par) for the parameter at which `quant` equals `target`.
   quant is monotonically decreasing in par for every algorithm here. */
template <typename Q>
static bool bisect(const Corpus& C, int alg, int w, Q quant, double target,
                   double pLo, double pHi, Meas& result, int maxIter) {
    Meas a = evaluate(C, alg, pLo, w);
    Meas b = evaluate(C, alg, pHi, w);
    if (!a.ok || !b.ok) return false;
    const double va = quant(a), vb = quant(b);
    const bool increasing = (vb > va);                 // error rises with par, size falls
    if (increasing) { if (!(va <= target && target <= vb)) return false; }
    else            { if (!(va >= target && target >= vb)) return false; }

    double lo = pLo, hi = pHi;
    Meas best = b; double bestGap = 1e30;
    for (int it = 0; it < maxIter && hi / lo > 1.0005; ++it) {
        const double mid = std::sqrt(lo * hi);
        Meas m = evaluate(C, alg, mid, w);
        if (!m.ok) break;
        const double v = quant(m);
        const double gap = std::fabs(v - target) / std::fabs(target);
        if (gap < bestGap) { bestGap = gap; best = m; }
        if (gap <= 0.005) { result = m; return true; }
        if (increasing) { if (v < target) lo = mid; else hi = mid; }
        else            { if (v < target) hi = mid; else lo = mid; }
    }
    if (bestGap <= 0.02) { result = best; return true; }
    return false;
}

int main(int argc, char** argv) {
    std::string dir;
    int maxFiles = 60, w = 8;
    Layout lay = Layout::TimeLast;
    bool errorMode = true;
    for (int i = 1; i < argc; ++i) {
        if (!std::strncmp(argv[i], "--max-files=", 12)) maxFiles = std::atoi(argv[i] + 12);
        else if (!std::strncmp(argv[i], "--w=", 4)) w = std::atoi(argv[i] + 4);
        else if (!std::strncmp(argv[i], "--layout=tlast", 14)) lay = Layout::TimeLast;
        else if (!std::strncmp(argv[i], "--layout=tfirst", 15)) lay = Layout::TimeFirst;
        else if (!std::strncmp(argv[i], "--mode=", 7)) errorMode = std::strncmp(argv[i]+7, "size", 4) != 0;
        else if (std::strncmp(argv[i], "--", 2)) dir = argv[i];
    }
    if (dir.empty()) { std::printf("usage: run_bisect <dir> [--layout=] [--max-files=] [--w=] [--mode=error|size]\n"); return 2; }

    std::vector<std::string> files;
    std::error_code ec;
    for (auto& e : fs::recursive_directory_iterator(dir, ec)) {
        if (!e.is_regular_file()) continue;
        const std::string ext = e.path().extension().string();
        if (ext == ".txt" || ext == ".csv") files.push_back(e.path().string());
    }
    std::sort(files.begin(), files.end());

    Corpus C;
    for (const std::string& f : files) {
        if ((int)C.tracks.size() >= maxFiles) break;
        Track tr;
        if (!load(f, lay, tr)) continue;
        C.tracks.push_back(standardise(tr));
        C.d = tr.d;
    }
    if (C.tracks.empty()) { std::printf("no tracks\n"); return 1; }
    double tot = 0;
    for (const auto& t : C.tracks) tot += (double)t.size();
    C.meanIn = tot / (double)C.tracks.size();

    const double ERR_T[] = {1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1};
    const double RET_T[] = {0.60, 0.40, 0.30, 0.20, 0.10, 0.05};
    const char* names[4] = {"stream", "svaried", "squishe", "dots"};
    const double P_LO = 1e-12, P_HI = 1e4;

    std::printf("dataset: %s\n", dir.c_str());
    std::printf("tracks: %zu  dimension: %d  mean_in: %.1f  mode: %s\n\n",
                C.tracks.size(), C.d, C.meanIn, errorMode ? "error" : "size");

    if (errorMode) {
        std::printf(" %-9s %-10s %-12s %-10s %-13s %-9s %s\n",
                    "algorithm", "target err", "achieved err", "out size", "retained", "us/pt", "status");
        std::printf(" %s\n", std::string(84, '-').c_str());
    } else {
        std::printf(" %-9s %-10s %-12s %-10s %-13s %-9s %s\n",
                    "algorithm", "target ret", "achieved ret", "out size", "frechet", "us/pt", "status");
        std::printf(" %s\n", std::string(84, '-').c_str());
    }

    for (int a = 0; a < 4; ++a) {
        const int nT = errorMode ? 7 : 6;
        for (int ti = 0; ti < nT; ++ti) {
            const double target = errorMode ? ERR_T[ti] : RET_T[ti] * C.meanIn;
            Meas r;
            bool ok;
            if (errorMode) {
                ok = bisect(C, a, w, [](const Meas& m) { return m.err; }, target, P_LO, P_HI, r, 40);
                std::printf(" %-9s %-10.3g %-12.3e %-10.1f %-13.1f%% %-9.3f %s\n",
                            names[a], target, r.err, r.out, 100.0 * r.out / C.meanIn, r.us,
                            ok ? "ok" : "UNREACHABLE");
            } else {
                ok = bisect(C, a, w, [](const Meas& m) { return m.out; }, target, P_LO, P_HI, r, 40);
                std::printf(" %-9s %-9.0f%% %-11.1f%% %-10.1f %-13.6e %-9.3f %s\n",
                            names[a], 100.0 * target / C.meanIn, 100.0 * r.out / C.meanIn, r.out,
                            r.err, r.us, ok ? "ok" : "UNREACHABLE");
            }
        }
        std::printf("\n");
    }
    return 0;
}
