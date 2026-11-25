/* Minimal Qt compatibility layer.

   The published baselines are Qt sources: QVector, QString, QObject, QtMath and the
   foreach/qSort macros. No Qt is installed on this toolchain, and none of the Qt
   machinery they use is actually needed -- QVector is a growable array, QString is
   only ever a filename, QObject is a NULL parent, and QtMath is a thin wrapper over
   <cmath>. This header maps each used name onto its standard-library equivalent so
   the published sources compile unmodified under g++ -O2 -std=c++17.

   Only the API the provided files actually touch is implemented; anything missing is
   a compile error rather than a silent behaviour change.
*/
#ifndef QTCOMPAT_H
#define QTCOMPAT_H

#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- Qt macros ---------------------------------------------------------- */
#define Q_OBJECT
/* The sources write `signals:` and `public slots:` -- Qt's macros expand to the access
   specifier, and the colon is in the source. */
#define signals public
#define slots
#define emit
#define Q_UNUSED(x) (void)(x)
/* Qt's foreach(variable, container). Qt copies the container; a range-for over a
   const reference does not, which is safe here because none of the loops mutate. */
#define foreach(variable, container) for (variable : container)
#define qSort(v) std::sort((v).begin(), (v).end())
#define qDebug(...) do { } while (0)

/* ---- QObject ------------------------------------------------------------ */
class QObject {
public:
    explicit QObject(QObject* = nullptr) {}
    virtual ~QObject() {}
};

/* ---- QException --------------------------------------------------------- */
class QException : public std::exception {
public:
    virtual void raise() const = 0;
    virtual QException* clone() const = 0;
    const char* what() const noexcept override { return "QException"; }
};

/* ---- QString ------------------------------------------------------------ */
/* Only used for filenames and exception text. QString's `arg` performs %1, %2 ...
   substitution, which the sources rely on for two messages; implemented so those
   messages read correctly rather than becoming placeholder soup. */
class QString : public std::string {
public:
    QString() {}
    QString(const char* s) : std::string(s ? s : "") {}
    QString(const std::string& s) : std::string(s) {}

    QString arg(const QString& a) const
    {
        /* `slots` is a macro in this header, so the table cannot be called that. */
        static const char* ph[9] = { "%1", "%2", "%3", "%4", "%5", "%6", "%7", "%8", "%9" };
        std::string out = *this;
        for (int k = 0; k < 9; ++k) {
            std::size_t p = out.find(ph[k]);
            if (p != std::string::npos) { out.replace(p, 2, a); break; }
        }
        return QString(out);
    }
    template <class T>
    QString arg(T v) const { return arg(QString(std::to_string(v))); }
};
typedef std::string QByteArray;

/* ---- QVector ------------------------------------------------------------ */
/* A std::vector with the Qt spelling of the methods the sources call. Inheriting
   rather than wrapping keeps pointer/iterator/reference semantics identical. */
template <class T>
class QVector : public std::vector<T> {
public:
    typedef std::vector<T> Base;
    QVector() {}
    QVector(int n) : Base((std::size_t)(n > 0 ? n : 0)) {}

    void append(const T& v) { this->push_back(v); }
    int  count() const { return (int)this->Base::size(); }
    bool isEmpty() const { return this->empty(); }
    const T& at(int i) const { return this->Base::operator[]((std::size_t)i); }
    T&       at(int i)       { return this->Base::operator[]((std::size_t)i); }
    const T& first() const { return this->front(); }
    T&       first()       { return this->front(); }
    const T& last() const  { return this->back(); }
    T&       last()        { return this->back(); }
    void removeLast() { this->pop_back(); }
    void clear() { this->Base::clear(); }
};

/* std::vector<bool> is the bitset specialization, so its operator[] returns a proxy
   rather than a bool&. The sources use `terminated[j] = true` (works through the
   proxy) and `terminated.at(j)` in a condition (needs a value, not a reference), so
   only at() has to be overridden. */
template <>
class QVector<bool> : public std::vector<bool> {
public:
    QVector() {}
    QVector(int n) : std::vector<bool>((std::size_t)(n > 0 ? n : 0)) {}
    void append(bool v) { this->push_back(v); }
    int  count() const { return (int)this->size(); }
    bool isEmpty() const { return this->empty(); }
    bool at(int i) const { return this->operator[]((std::size_t)i); }
    void clear() { this->std::vector<bool>::clear(); }
};

/* ---- QtMath ------------------------------------------------------------- */
inline double qSqrt(double x)  { return std::sqrt(x); }
inline double qAbs(double x)   { return std::fabs(x); }
inline int    qAbs(int x)      { return x < 0 ? -x : x; }
inline double qFabs(double x)  { return std::fabs(x); }
inline double qFloor(double x) { return std::floor(x); }
inline double qCeil(double x)  { return std::ceil(x); }
inline double qLn(double x)    { return std::log(x); }
inline double qPow(double a, double b) { return std::pow(a, b); }
inline double qCos(double x)   { return std::cos(x); }
inline double qSin(double x)   { return std::sin(x); }
inline double qTan(double x)   { return std::tan(x); }
inline double qAtan2(double y, double x) { return std::atan2(y, x); }
inline double qDegreesToRadians(double d) { return d * M_PI / 180.0; }
inline double qRadiansToDegrees(double r) { return r * 180.0 / M_PI; }
inline int    qRound(double x) { return (int)std::floor(x + 0.5); }
template <class T> inline T qMax(T a, T b) { return a > b ? a : b; }
template <class T> inline T qMin(T a, T b) { return a < b ? a : b; }

#endif /* QTCOMPAT_H */
