#include <cmath>
#include <cstdio>

#include "proj.hpp"

#include "sphere-map-common.hpp"


void latlon_to_sphere(double lat, double lon, double& u, double& v, double& w)
{
    u = std::cos(lat) * std::cos(lon);
    v = std::cos(lat) * std::sin(lon);
    w = std::sin(lat);
}

double dist_rad_sphere(double lat0, double lon0, double lat1, double lon1)
{
    /* angle between vectors: not precise enough!
    double u0, v0, w0, u1, v1, w1;
    latlon_to_sphere(lat0, lon0, u0, v0, w0);
    latlon_to_sphere(lat1, lon1, u1, v1, w1);
    double angle = std::acos(u0 * u1 + v0 * v1 + w0 * w1);
    return angle;
    */
    double slat2 = std::sin(std::fabs(lat0 - lat1) / 2.0);
    slat2 *= slat2;
    double slon2 = std::sin(std::fabs(lon0 - lon1) / 2.0);
    slon2 *= slon2;
    double t = slat2 + std::cos(lat0) * std::cos(lat1) * slon2;
    return 2.0 * std::asin(std::sqrt(t));
}

double check_precision(Context& ctx)
{
    const int xsteps = 2000;
    const int ysteps = 2000;

    double max_error = -1.0;
    int max_considered_points = 0;
    int considered_points = 0;
    for (int iy = 0; iy < ysteps; iy++) {
        double y = (iy / (ysteps - 1.0)) * 2.0 - 1.0;
        for (int ix = 0; ix < xsteps; ix++) {
            double x = (ix / (xsteps - 1.0)) * 2.0 - 1.0;
            double lat, lon;
            if ((ctx.layout == LayoutHemisphere || ctx.layout == LayoutSphere)
                    && x * x + y * y > 1) {
                continue;
            }
            max_considered_points++;
            ctx.inverse(x, y, lat, lon);
            if (!std::isfinite(lat) || !std::isfinite(lon)) {
                //fprintf(stderr, "!cannot inverse map %g,%g! ", x, y);
                return NAN;
            }
            double fx, fy, flat, flon;
            ctx.forward(lat, lon, fx, fy);
            if (!std::isfinite(fx) || !std::isfinite(fy)) {
                // we cannot compute an error for this point
                //fprintf(stderr, "!cannot map %.19g,%.19g (comes from %g,%g)! ", degrees(lat), degrees(lon), x, y);
                continue;
            }
            ctx.inverse(fx, fy, flat, flon);
            if (!std::isfinite(flat) || !std::isfinite(flon)) {
                //fprintf(stderr, "!cannot inverse map %g,%g! ", fx, fy);
                return NAN;
            }
            considered_points++;
            double error = dist_rad_sphere(lat, lon, flat, flon);
            /*
            if (error > 1e-12) {
                fprintf(stderr, "suspiciously large error at "
                        "x=%.17g y=%.17g lat=%.17g lon=%.17g "
                        "fx=%.17g fy=%.17g flat=%.17g flon=%.17g\n",
                        x, y, degrees(lat), degrees(lon),
                        fx, fy, degrees(flat), degrees(flon));
            }
            */
            if (error > max_error) {
                max_error = error;
            }
        }
    }
    if (considered_points < max_considered_points - max_considered_points / 100) {
        fprintf(stderr, "!only %d out of %d points could be considered! ", considered_points, max_considered_points);
        return NAN;
    }
    return max_error;
}

double in_earth_millimeters(double rad_error)
{
    return rad_error * 40075.017 * 1e3 * 1e3;
}

/*
 * Main function.
 */

int main(void)
{
    // Loop over the projection methods
    for (int method = 0; method < n_projections; method++) {
        Context ctx(projections[method].proj,
                projections[method].params,
                projections[method].center,
                projections[method].layout,
                projections[method].square_method);
        double max_error = check_precision(ctx);
        fprintf(stderr, "%s: %g (ca. %g mm on Earth)\n", projections[method].name, max_error, in_earth_millimeters(max_error));
    }
    return 0;
}
