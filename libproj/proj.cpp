#include "proj.hpp"

#include "sphere-lambert-azimuthal-equal-area.hpp"
#include "sphere-breusing-harmonic.hpp"
#include "sphere-equidistant-azimuthal.hpp"
#include "hemisphere-stereographic.hpp"
#include "hemisphere-mix-lambert-stereographic.hpp"
#include "hemisphere-perspective.hpp"

#include "disk-square.hpp"

namespace Projection {

Context::Context(enum Projection P,
        const std::string& params,
        enum Center center,
        enum Layout layout,
        enum SquareMethod square_method,
        double semi_major_axis, double semi_minor_axis) :
    _valid(false),
    _imp(nullptr), _forward(nullptr), _inverse(nullptr),
    _disk_to_square(nullptr), _square_to_disk(nullptr),
    semi_major_axis(semi_major_axis), semi_minor_axis(semi_minor_axis),
    eccentricity_squared(1 - (semi_minor_axis * semi_minor_axis) / (semi_major_axis * semi_major_axis)),
    is_sphere(semi_major_axis >= semi_minor_axis && semi_major_axis <= semi_minor_axis),
    projection(P), center(center), layout(layout), square_method(square_method)
{
    switch (P) {
    case ProjLambertAzimuthalEqualArea:
        if (center == CenterNorthPole || center == CenterSouthPole) {
            _imp = new LambertAzimuthalEqualArea();
            _forward = LambertAzimuthalEqualArea::forward;
            _inverse = LambertAzimuthalEqualArea::inverse;
            _valid = true;
        }
        break;
    case ProjStereographic:
        if (layout != LayoutSphere && (center == CenterNorthPole || center == CenterSouthPole)) {
            _imp = new Stereographic();
            _forward = Stereographic::forward;
            _inverse = Stereographic::inverse;
            _valid = true;
        }
        break;
    case ProjBreusingHarmonicMean:
        if (center == CenterNorthPole || center == CenterSouthPole) {
            _imp = new BreusingHarmonicMean();
            _forward = BreusingHarmonicMean::forward;
            _inverse = BreusingHarmonicMean::inverse;
            _valid = true;
        }
        break;
    case ProjMixLambertStereographic:
        if (layout != LayoutSphere && (center == CenterNorthPole || center == CenterSouthPole)) {
            double beta;
            if (sscanf(params.c_str(), "beta=%lf", &beta) == 1) {
                _imp = new MixLambertStereographic(beta);
                _forward = MixLambertStereographic::forward;
                _inverse = MixLambertStereographic::inverse;
                _valid = true;
            }
        }
        break;
    case ProjEquidistantAzimuthal:
        if (center == CenterNorthPole || center == CenterSouthPole) {
            _imp = new EquidistantAzimuthal();
            _forward = EquidistantAzimuthal::forward;
            _inverse = EquidistantAzimuthal::inverse;
            _valid = true;
        }
        break;
    case ProjPerspective:
        if (layout != LayoutSphere && (center == CenterNorthPole || center == CenterSouthPole)) {
            _imp = new Perspective(-1e10,
                    50.90592 / 180 * M_PI, 8.0285  / 180 * M_PI); // my office ;)
            _forward = Perspective::forward;
            _inverse = Perspective::inverse;
            _valid = true;
        }
        break;
    }
    if (layout == LayoutQuincuncial && square_method == SquareMethodNone) {
        _valid = false;
    } else {
        switch (square_method) {
        case SquareMethodNone:
            break;
        case SquareMethodStretch:
            _disk_to_square = disk_to_square_stretch;
            _square_to_disk = square_to_disk_stretch;
            break;
        case SquareMethodShirley:
            _disk_to_square = disk_to_square_shirley;
            _square_to_disk = square_to_disk_shirley;
            break;
        case SquareMethodSquircle:
            _disk_to_square = disk_to_square_squircle;
            _square_to_disk = square_to_disk_squircle;
            break;
        case SquareMethodElliptical:
            _disk_to_square = disk_to_square_elliptical;
            _square_to_disk = square_to_disk_elliptical;
            break;
        case SquareMethodConformal:
            _disk_to_square = disk_to_square_conformal;
            _square_to_disk = square_to_disk_conformal;
            break;
        }
    }
}

void Context::forward(double lat, double lon, double& x, double& y) const
{
    if (projection != ProjPerspective && layout == LayoutHemisphere && !in_hemisphere(center, lat, lon)) {
        x = y = NAN;
    } else if (layout == LayoutQuincuncial) {
        bool center_hemisphere =
            (center == CenterNorthPole && lat >= 0)
            || (center == CenterSouthPole && lat <= 0);
        _forward(this, _imp, center_hemisphere ? lat : -lat, lon + M_PI_4, x, y);
        double sx, sy, t;
        // Transform unit disk to [-1,+1]^2 square
        _disk_to_square(x, y, sx, sy);
        // rotate by -45 degrees while keeping square in [-1,+1]^2
        x =  sx / 2 + sy / 2;
        y = -sx / 2 + sy / 2;
        // place the southern hemisphere into the four free corners
        if (!center_hemisphere) {
            if (x * y >= 0) {
                // I or III
                if (x >= 0) {
                    t = -y + 1;
                    y = -x + 1;
                    x = t;
                } else {
                    t = -y - 1;
                    y = -x - 1;
                    x = t;
                }
            } else {
                // II or IV
                if (x < 0) {
                    t = y - 1;
                    y = x + 1;
                    x = t;
                } else {
                    t = y + 1;
                    y = x - 1;
                    x = t;
                }
            }
        }
    } else {
        _forward(this, _imp, lat, lon, x, y);
        if (square_method != SquareMethodNone)
            _disk_to_square(x, y, x, y);
    }
}

void Context::inverse(double x, double y, double& lat, double& lon) const
{
    if (layout == LayoutQuincuncial) {
        bool center_hemisphere = (fabs(x) + fabs(y) <= 1);
        // remap the southern hemisphere from the corners to a center square
        if (!center_hemisphere) {
            double t;
            if (x * y >= 0) {
                // I or III
                if (x >= 0) {
                    t = -y + 1;
                    y = -x + 1;
                    x = t;
                } else {
                    t = -y - 1;
                    y = -x - 1;
                    x = t;
                }
            } else {
                // II or IV
                if (x < 0) {
                    t = y - 1;
                    y = x + 1;
                    x = t;
                } else {
                    t = y + 1;
                    y = x - 1;
                    x = t;
                }
            }
        }
        // rotate by 45 degrees while keeping square in [-1,+1]^2
        double sx = x - y;
        double sy = x + y;
        // transform square to disk
        _square_to_disk(sx, sy, x, y);
        // inverse projection
        _inverse(this, _imp, x, y, lat, lon);
        if (!center_hemisphere)
            lat = -lat;
        lon -= M_PI_4;
        if (lon < -M_PI)
            lon += 2 * M_PI;
    } else {
        if (square_method != SquareMethodNone)
            _square_to_disk(x, y, x, y);
        _inverse(this, _imp, x, y, lat, lon);
    }
}

void Context::derivatives(double lat, double lon,
        double& dx_dlat, double& dx_dlon, double& dy_dlat, double& dy_dlon, double delta) const
{
    double mlat = lat - delta;
    double plat = lat + delta;
    double mlon = lon - delta;
    double plon = lon + delta;

    if (mlat < -M_PI_2 || plat > M_PI_2)
        goto fail;
    if (mlon < -M_PI)
        mlon += 2 * M_PI;
    if (plon > M_PI)
        plon -= 2 * M_PI;

    double x_mlat_mlon, x_mlat_plon, x_plat_mlon, x_plat_plon;
    double y_mlat_mlon, y_mlat_plon, y_plat_mlon, y_plat_plon;

    forward(mlat, mlon, x_mlat_mlon, y_mlat_mlon);
    forward(mlat, plon, x_mlat_plon, y_mlat_plon);
    forward(plat, mlon, x_plat_mlon, y_plat_mlon);
    forward(plat, plon, x_plat_plon, y_plat_plon);

    if (!isfinite(x_mlat_mlon) || !isfinite(x_mlat_plon) || !isfinite(x_plat_mlon) || !isfinite(x_plat_plon)
            || !isfinite(y_mlat_mlon) || !isfinite(y_mlat_plon) || !isfinite(y_plat_mlon) || !isfinite(y_plat_plon))
        goto fail;

    dx_dlat = (- x_plat_plon + x_mlat_plon + x_mlat_mlon - x_plat_mlon) / (4 * delta);
    dx_dlon = (+ x_plat_plon + x_mlat_plon - x_mlat_mlon - x_plat_mlon) / (4 * delta);
    dy_dlat = (+ y_plat_plon - y_mlat_plon - y_mlat_mlon + y_plat_mlon) / (4 * delta);
    dy_dlon = (- y_plat_plon - y_mlat_plon + y_mlat_mlon + y_plat_mlon) / (4 * delta);
    return;

fail:
    dx_dlat = dx_dlon = dy_dlat = dy_dlon = NAN;
}

void Context::analysis(double lat, double lon,
        double& h, double& k,
        double& theta_prime,
        double& a, double& b,
        double& omega, double& s,
        double delta) const
{
    double dx_dlat, dx_dlon, dy_dlat, dy_dlon;

    derivatives(lat, lon, dx_dlat, dx_dlon, dy_dlat, dy_dlon, delta);
    if (!isfinite(dx_dlat) || !isfinite(dx_dlon) || !isfinite(dy_dlat) || !isfinite(dy_dlon)) {
        h = k = theta_prime = a = b = omega = s = NAN;
        return;
    }

    double coslat = cos(lat);
    h = hypot(dx_dlat, dy_dlat);
    k = hypot(dx_dlon, dy_dlon) / coslat;
    double sin_theta_prime = (dy_dlat * dx_dlon - dx_dlat * dy_dlon) / (h * k * coslat);
    s = h * k * sin_theta_prime; // equivalent to s = a * b
    theta_prime = asin(min(max(sin_theta_prime, -1.0), +1.0));
    double a_prime = sqrt(max(h * h + k * k + 2 * s, 0.0));
    double b_prime = sqrt(max(h * h + k * k - 2 * s, 0.0));
    a = (a_prime + b_prime) / 2;
    b = (a_prime - b_prime) / 2;
    omega = 2 * asin(min(max(b_prime / a_prime, -1.0), +1.0));
}

double Context::geodetic_to_geocentric(double lat) const
{
    return is_sphere ? lat : atan(tan(lat) * (1 - eccentricity_squared));
}

double Context::geocentric_to_geodetic(double lat) const
{
    return is_sphere ? lat : atan(tan(lat) / (1 - eccentricity_squared));
}

}
