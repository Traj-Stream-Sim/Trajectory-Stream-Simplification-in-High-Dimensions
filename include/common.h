#pragma once
#include <cmath>
#include <vector>

// A trajectory point.
//
// `t` is an optional timestamp (seconds, relative to the trajectory start).
// It is only required by the time-aware algorithms (SQUISH, DOT), which score
// candidates with the Synchronous Euclidean Distance; the geometric algorithms
// (streaming delta-simplification, Douglas-Peucker) and the continuous Frechet
// metric ignore it entirely. Loaders that read x,y-only files leave it at 0.
struct Vertex {
    double x, y;
    double t;

    Vertex() : x(0.0), y(0.0), t(0.0) {}
    Vertex(double x_, double y_) : x(x_), y(y_), t(0.0) {}
    Vertex(double x_, double y_, double t_) : x(x_), y(y_), t(t_) {}
};

struct Vec2D {
    double x, y;
};

inline Vec2D operator-(const Vertex& a, const Vertex& b) {
    return {a.x - b.x, a.y - b.y};
}

inline double dot(const Vec2D& a, const Vec2D& b) {
    return a.x * b.x + a.y * b.y;
}

inline double crossMag(const Vec2D& a, const Vec2D& b) {
    // |a x b| == |a| * |b| * sin(angle)
    return std::abs(a.x * b.y - a.y * b.x);
}

inline double magnitude(const Vec2D& v) {
    return std::sqrt(v.x * v.x + v.y * v.y);
}

inline double dist(const Vertex& a, const Vertex& b) {
    double dx = a.x - b.x, dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Perpendicular distance from point p to the infinite line through a-b.
inline double pointLineDist(const Vertex& p, const Vertex& a, const Vertex& b) {
    Vec2D ab = b - a;
    double len = magnitude(ab);
    if (len == 0.0) return dist(p, a);
    return crossMag(p - a, ab) / len;
}

// Characteristic scale of a trajectory: max(x-span, y-span) of its bounding box.
inline double trajectoryScale(const std::vector<Vertex>& pts) {
    if (pts.empty()) return 1.0;
    double minx = pts[0].x, maxx = pts[0].x;
    double miny = pts[0].y, maxy = pts[0].y;
    for (const Vertex& p : pts) {
        if (p.x < minx) minx = p.x;
        if (p.x > maxx) maxx = p.x;
        if (p.y < miny) miny = p.y;
        if (p.y > maxy) maxy = p.y;
    }
    double sx = maxx - minx, sy = maxy - miny;
    double s = (sx > sy) ? sx : sy;
    return (s > 0.0) ? s : 1.0;
}

// Isotropically rescale + translate a trajectory so max(x-span, y-span) == 1.
inline std::vector<Vertex> normalizeTrajectory(const std::vector<Vertex>& pts) {
    if (pts.empty()) return pts;
    double minx = pts[0].x, miny = pts[0].y;
    double scale = trajectoryScale(pts);
    std::vector<Vertex> out;
    out.reserve(pts.size());
    for (const Vertex& p : pts)
        out.push_back(Vertex{(p.x - minx) / scale, (p.y - miny) / scale, p.t});
    return out;
}

// --- metric projection for the time-aware algorithms -----------------------
//
// SQUISH and DOTS measure the Synchronous Euclidean Distance, which compares
// "where the point is" against "where a constant-velocity model says it should
// be at that timestamp". That comparison is only meaningful when space and time
// are in physically commensurate units.
inline double metresPerDegreeLon(double latDeg) {
    return 111320.0 * std::cos(latDeg * 3.14159265358979323846 / 180.0);
}
inline double metresPerDegreeLat() {
    return 110540.0;   // mean meridional degree
}

// Metres per world unit, for input that is already planar: 1 Source/CS2 Hammer
// unit is 1.905 cm. Only a uniform factor, and every relative result below is
// invariant to it, but it is the honest one for game data.
inline double metresPerUnit() {
    return 0.01905;
}

// Is this coordinate pair outside anything a lon/lat pair could be?
//
// Latitude cannot exceed +-90 and longitude +-180, so a point beyond those bounds
// is not a degree pair. The CS2 tracks in data/cs2 are Hammer world units (|x|,|y|
// in the thousands), and feeding those through cos(y) below is not merely useless
// but wrong: y is then tens of radians, so the x axis gets multiplied by an
// arbitrary factor in [-1,1]. Measured on the 630 CS2 tracks, the resulting x:y
// metric scale is 0.77x at the median, below 0.5x on 31% of tracks, and below
// 0.05x on 10 of them -- on those, SED would effectively ignore horizontal motion
// and SQUISH/DOTS would be scored in a space that is not the one the player moved in.
inline bool isPlanarInput(const std::vector<Vertex>& pts) {
    for (const Vertex& p : pts)
        if (std::abs(p.y) > 90.0 || std::abs(p.x) > 180.0) return true;
    return false;
}

// Project to metres relative to the first point. Timestamps are untouched.
inline std::vector<Vertex> projectToMetres(const std::vector<Vertex>& pts) {
    if (pts.empty()) return pts;
    const double x0 = pts[0].x, y0 = pts[0].y;
    std::vector<Vertex> out;
    out.reserve(pts.size());

    // Already a planar system (game world units): rescale uniformly, do NOT
    // interpret y as latitude. Uniform scaling leaves the SED ordering, the
    // relative parameter grids and every span-relative error identical.
    if (isPlanarInput(pts)) {
        const double m = metresPerUnit();
        for (const Vertex& p : pts)
            out.push_back(Vertex{(p.x - x0) * m, (p.y - y0) * m, p.t});
        return out;
    }

    const double mlon = metresPerDegreeLon(y0);
    const double mlat = metresPerDegreeLat();
    for (const Vertex& p : pts)
        out.push_back(Vertex{(p.x - x0) * mlon, (p.y - y0) * mlat, p.t});
    return out;
}

// True when the trajectory carries real (non-degenerate) timestamps.
inline bool hasTime(const std::vector<Vertex>& pts) {
    if (pts.size() < 2) return false;
    for (std::size_t i = 1; i < pts.size(); ++i)
        if (pts[i].t != pts[i - 1].t) return true;
    return false;
}
