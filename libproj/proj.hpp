#ifndef PROJ_H
#define PROJ_H

#include <cmath>
#include <string>

namespace Projection {

enum Projection {
    ProjLambertAzimuthalEqualArea,
    ProjStereographic,
    ProjBreusingHarmonicMean,
    ProjMixLambertStereographic,
    ProjEquidistantAzimuthal,
    ProjPerspective
};

enum Center {
    CenterNorthPole,
    CenterSouthPole
};

enum Layout {
    LayoutHemisphere,
    LayoutSphere,
    LayoutQuincuncial,
};

enum SquareMethod {
    SquareMethodNone,
    SquareMethodStretch,
    SquareMethodShirley,
    SquareMethodSquircle,
    SquareMethodElliptical,
    SquareMethodConformal
};

class Implementation;

class Context
{
private:
    bool _valid;
    Implementation* _imp;
    void (*_forward)(const Context* ctx, const Implementation* imp, double lat, double lon, double& x, double& y);
    void (*_inverse)(const Context* ctx, const Implementation* imp, double x, double y, double& lat, double& lon);
    void (*_disk_to_square)(double u, double v, double& x, double& y);
    void (*_square_to_disk)(double x, double y, double& u, double& v);

public:

    /* Information about the Ellipsoid */
    const double semi_major_axis;
    const double semi_minor_axis;
    const double eccentricity_squared;
    const bool is_sphere;

    /* Generic parameters of the projection */
    const enum Projection projection;
    const enum Center center;
    const enum Layout layout;
    const enum SquareMethod square_method;

    /* User interface */
    Context(enum Projection P,
            const std::string& params,
            enum Center center,
            enum Layout layout,
            enum SquareMethod square_method,
            double semi_major_axis = 1, double semi_minor_axis = 1);

    bool is_valid() const { return _valid; }

    void forward(double lat, double lon, double& x, double& y) const;
    void inverse(double x, double y, double& lat, double& lon) const;

    void derivatives(double lat, double lon,
            double& dx_dlat, double& dx_dlon, double& dy_dlat, double& dy_dlon, double delta = 1e-5) const;

    void analysis(double lat, double lon,
            double& h, double& k, double& theta_prime,
            double& a, double& b,
            double& omega, double& s,
            double delta = 1e-5) const;

    /* Implementation interface */
    double geodetic_to_geocentric(double lat) const;
    double geocentric_to_geodetic(double lat) const;
};

using std::sin;
using std::asin;
using std::cos;
using std::acos;
using std::tan;
using std::atan;
using std::atan2;
using std::fabs;
using std::hypot;
using std::sqrt;
using std::isfinite;
using std::min;
using std::max;
using std::floor;
using std::ceil;
using std::hypot;
using std::log;
using std::exp;

constexpr double sign(double x) { return (x < 0 ? -1 : x > 0 ? +1 : 0); }
constexpr double sign_not_zero(double x) { return (x < 0 ? -1 : +1); }

constexpr double degrees(double x) { return 180 / M_PI * x; }
constexpr double radians(double x) { return M_PI / 180 * x; }

constexpr bool in_hemisphere(Center c, double lat, double /* lon */)
{
    return (c == CenterNorthPole ? (lat >= 0) : (lat <= 0));
}

class Implementation
{
public:
    static void forward(const Context* /* ctx */, const Implementation* /* imp */, double /* lat */, double /* lon */, double& /* x */, double& /* y */) {}
    static void inverse(const Context* /* ctx */, const Implementation* /* imp */, double /* x */, double /* y */, double& /* lat */, double& /* lon */) {}
};

}

#endif
