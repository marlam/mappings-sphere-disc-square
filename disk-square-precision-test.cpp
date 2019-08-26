#include <cmath>
#include <cstdio>

#include "disk-square-common.hpp"


double check_precision(int method)
{
    const int xsteps = 2000;
    const int ysteps = 2000;

    double max_error = -1.0;
    for (int iy = 0; iy < ysteps; iy++) {
        double y = (iy / (ysteps - 1.0)) * 2.0 - 1.0;
        for (int ix = 0; ix < xsteps; ix++) {
            double x = (ix / (xsteps - 1.0)) * 2.0 - 1.0;
            double u, v, fx, fy, fu, fv;
            projections[method].square_to_disk(x, y, u, v);
            projections[method].disk_to_square(u, v, fx, fy);
            projections[method].square_to_disk(fx, fy, fu, fv);
            double error = hypot(u - fu, v - fv);
            if (!std::isfinite(u) || !std::isfinite(v)
                    || !std::isfinite(fx) || !std::isfinite(fy)
                    || !std::isfinite(fu) || !std::isfinite(fv)) {
                fprintf(stderr, "projection failure "
                        "ix=%d iy=%d x=%g y=%g u=%g v=%g fx=%g fy=%g fu=%g fv=%g\n",
                        ix, iy, x, y, u, v, fx, fy, fu, fv);
            }
            if (error > max_error) {
                max_error = error;
            }
        }
    }
    return max_error;
}

double in_earth_millimeters(double error)
{
    return error * 40075.017 * 1e3 * 1e3;
}

/*
 * Main function.
 */

int main(void)
{
    // Loop over the projection methods
    for (int method = 0; method < n_projections; method++) {
        double max_error = check_precision(method);
        fprintf(stderr, "%s: %g (ca. %g mm on Earth equatorial disk)\n", projections[method].name, max_error, in_earth_millimeters(max_error));
    }
    return 0;
}
