//  run_unified.cpp -- ONE benchmark driver for all datasets and all algorithms.
//
//  AUDIT FINDING. The previous GeoLife/Rome pipeline prepared data differently per
//  algorithm: `stream` and `svaried` received normalizeTrajectory(), while `squishe` and
//  `dots` received projectToMetres(). That is algorithm-specific preprocessing, and
//  projectToMetres is not even isotropic -- it scales longitude by 111320*cos(lat) against
//  latitude's 110540, which compresses one axis relative to the other by up to 30% at high
//  latitudes. The published algorithms were therefore solving a geometrically different
//  problem from ours.
//
//  STANDARDISED METHODOLOGY, applied identically to every dataset and every algorithm:
//
//       span = max over axes of (max_c - min_c)          one common factor
//       c'[k] = (c[k] - min[k]) / span                    all axes divided by the same number
//       t'   = t                                          untouched
//
//  Time is left alone deliberately: SED's interpolation parameter is
//  (t_m - t_s)/(t_e - t_s), a ratio, so any uniform scaling of the clock cancels exactly and
//  cannot affect any result. Introducing a time scale would be a no-op dressed as care.
//
//  Consequences:
//    * no per-axis rescaling of any kind, on any dataset;
//    * no latitude-dependent projection;
//    * every parameter is relative to the trajectory's own extent, so the same grid means
//      the same thing for every algorithm and every dataset;
//    * the Frechet error is in the same normalised units everywhere.

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

struct Track {
    std::vector<dim::VertexD> pts;
    int d = 0;
};

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
            t = v[0];
            v.erase(v.begin());
        } else {
            t = v.back();
            v.pop_back();
        }
        if (v.empty()) continue;
        tr.pts.push_back(dim::VertexD(v, t));
    }
    tr.d = tr.pts.empty() ? 0 : tr.pts[0].dim();
    return tr.pts.size() > 8;
}

/* The one and only preparation. One common factor across every axis, on every dataset,
   for every algorithm. */
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
        for (int k = 0; k < d; ++k) c[k] = (p.c[k] - lo[k]) / span;   // SAME divisor, all axes
        out.push_back(dim::VertexD(c, p.t));
    }
    return out;
}

int main(int argc, char** argv) {
    std::string dir;
    int maxFiles = 60, w = 8;
    Layout lay = Layout::TimeLast;
    std::vector<double> grid;
    for (int i = 1; i < argc; ++i) {
        if (!std::strncmp(argv[i], "--max-files=", 12)) maxFiles = std::atoi(argv[i] + 12);
        else if (!std::strncmp(argv[i], "--w=", 4))    w = std::atoi(argv[i] + 4);
        else if (!std::strncmp(argv[i], "--layout=tlast", 14))  lay = Layout::TimeLast;
        else if (!std::strncmp(argv[i], "--layout=tfirst", 15)) lay = Layout::TimeFirst;
        else if (!std::strncmp(argv[i], "--grid=", 7)) {
            std::string s(argv[i] + 7); std::size_t p = 0;
            while (p < s.size()) {
                std::size_t q = s.find(',', p);
                if (q == std::string::npos) q = s.size();
                grid.push_back(std::atof(s.substr(p, q - p).c_str()));
                p = q + 1;
            }
        } else dir = argv[i];
    }
    if (grid.empty()) grid = {1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1};
    if (dir.empty()) { std::printf("usage: run_unified <dir> [--layout=] [--max-files=] [--grid=]\n"); return 2; }

    std::vector<std::string> files;
    std::error_code ec;
    for (auto& e : fs::recursive_directory_iterator(dir, ec)) {
        if (!e.is_regular_file()) continue;
        const std::string ext = e.path().extension().string();
        if (ext == ".txt" || ext == ".csv") files.push_back(e.path().string());
    }
    std::sort(files.begin(), files.end());
    if (files.empty()) { std::printf("no .txt files under %s\n", dir.c_str()); return 1; }

    Track probe;
    bool probed = load(files[0], lay, probe);
    std::printf("dataset: %s\n", dir.c_str());
    std::printf("files: %d  dimension: %d  layout: %s\n\n",
                (int)files.size(), probed ? probe.d : 0,
                lay == Layout::TimeFirst ? "time first" : "time last");

    const char* names[4] = {"stream", "svaried", "squishe", "dots"};
    for (int a = 0; a < 4; ++a) {
        std::printf("===========================================================================\n");
        std::printf(" %s\n", names[a]);
        std::printf("===========================================================================\n");
        std::printf(" %-10s %-7s %-9s %-10s %-13s %s\n",
                    "param", "tracks", "mean_in", "mean_out", "frechet", "us/point");
        std::printf(" %s\n", std::string(74, '-').c_str());

        for (double par : grid) {
            long long ntr = 0; double sin_ = 0, sout = 0, serr = 0, core = 0, cpu = 0;
            for (const std::string& f : files) {
                if ((int)ntr >= maxFiles) break;
                Track tr;
                if (!load(f, lay, tr)) continue;
                ++ntr;
                std::vector<dim::VertexD> in = standardise(tr);

                std::vector<dim::VertexD> out;
                const auto t0 = std::chrono::steady_clock::now();
                if (a == 0) {
                    StreamSimplifierD s(par, w);
                    out = s.simplify(in);
                } else if (a == 1) {
                    svd::Settings st = measuredSettingsD(w);
                    StreamVariedSimplifierD s(par, st);
                    out = s.simplify(in);
                } else if (a == 2) {
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
                        coord.append(c);
                        tt.append(v.t);
                    }
                    QVector<int> idx;
                    DotsSimplifierD::batchDotsByIndexD(coord, tt, idx, par);
                    out.reserve(idx.count());
                    for (int i = 0; i < idx.count(); ++i) out.push_back(in[idx[i]]);
                }
                const auto t1 = std::chrono::steady_clock::now();
                cpu += std::chrono::duration<double>(t1 - t0).count();

                sin_ += (double)in.size();
                sout += (double)out.size();
                core += (double)in.size();
                if (out.size() >= 2) serr += dim::continuousFrechet(in, out);
            }
            if (ntr == 0) continue;
            std::printf(" %-10g %-7lld %-9.0f %-9.1f %-13.6g %.3f\n",
                        par, ntr, sin_ / ntr, sout / ntr, serr / ntr,
                        core > 0 ? (cpu * 1e6 / core) : 0.0);
        }
        std::printf("\n");
    }
    return 0;
}
