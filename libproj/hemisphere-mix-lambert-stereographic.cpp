#include "hemisphere-mix-lambert-stereographic.hpp"

namespace Projection {

MixLambertStereographic::MixLambertStereographic(double beta) : beta(beta)
{
}

void MixLambertStereographic::forward(const Context* ctx, const Implementation* imp, double lat, double lon, double& x, double& y)
{
    const MixLambertStereographic* me = reinterpret_cast<const MixLambertStereographic*>(imp);

    if (!ctx->is_sphere)
        lat = ctx->geodetic_to_geocentric(lat);

    if (ctx->center == CenterSouthPole)
        lat = -lat;

    double t = tan((M_PI_2 - lat) / 2);
    double r = sqrt(1 + me->beta) * t / sqrt(1 + me->beta * t * t);
    double phi = lon - M_PI_2;

    x = r * cos(phi);
    y = r * sin(phi);
}

void MixLambertStereographic::inverse(const Context* ctx, const Implementation* imp, double x, double y, double& lat, double& lon)
{
    const MixLambertStereographic* me = reinterpret_cast<const MixLambertStereographic*>(imp);

    double r = hypot(x, y);
    double phi = atan2(y, x);

    lat = M_PI_2 - 2 * atan(r / sqrt(1 + me->beta * (1 - r * r)));
    lon = phi + M_PI_2;

    if (ctx->center == CenterSouthPole)
        lat = -lat;

    if (!ctx->is_sphere)
        lat = ctx->geocentric_to_geodetic(lat);
}

}
