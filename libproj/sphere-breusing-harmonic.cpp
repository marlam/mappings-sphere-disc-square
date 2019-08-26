#include "sphere-breusing-harmonic.hpp"

namespace Projection {

BreusingHarmonicMean::BreusingHarmonicMean()
{
}

void BreusingHarmonicMean::forward(const Context* ctx, const Implementation* /* imp */, double lat, double lon, double& x, double& y)
{
    if (!ctx->is_sphere)
        lat = ctx->geodetic_to_geocentric(lat);

    if (ctx->center == CenterSouthPole)
        lat = -lat;

    // Lambert az. eq-ar: r=sin((pi/2-lat)/2)
    // Stereographic: r=tan((pi/2-lat)/2)
    // Harmonic mean: r = 2*tan((pi/2-lat)/4)
    // Geometric mean: r = 2*sqrt(sin((pi/2-lat)/2)*tan((pi/2-lat)/2))
    // Arithm mean: r = (sin((pi/2-lat)/2)+tan((pi/2-lat)/2))/2
    double r = tan((M_PI_2 - lat) / 4); // harmonic mean
    double phi = lon - M_PI_2;

    if (ctx->layout != LayoutSphere)
        r /= (M_SQRT2 - 1);

    x = r * cos(phi);
    y = r * sin(phi);
}

void BreusingHarmonicMean::inverse(const Context* ctx, const Implementation* /* imp */, double x, double y, double& lat, double& lon)
{
    double r = hypot(x, y);
    double phi = atan2(y, x);

    if (ctx->layout != LayoutSphere)
        r *= (M_SQRT2 - 1);

    lat = M_PI_2 - 4 * atan(r);
    lon = phi + M_PI_2;

    if (ctx->center == CenterSouthPole)
        lat = -lat;

    if (!ctx->is_sphere)
        lat = ctx->geocentric_to_geodetic(lat);
}

}
