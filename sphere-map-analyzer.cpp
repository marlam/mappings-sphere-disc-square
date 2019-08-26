#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

#include "proj.hpp"

#include "sphere-map-common.hpp"
#include "pngvis.hpp"


void analyze(const Context& ctx, int w, int h, double* dist_area, double* dist_isot)
{
    const double area_unit_sphere = 4 * M_PI;
    const double area_unit_disk = M_PI;
    const double area_unit_square = 4;
    double R;
    switch (ctx.layout) {
    case LayoutHemisphere:
        R = (ctx.square_method == SquareMethodNone ? area_unit_disk : area_unit_square) / (0.5 * area_unit_sphere);
        break;
    case LayoutSphere:
        R = (ctx.square_method == SquareMethodNone ? area_unit_disk : area_unit_square) / area_unit_sphere;
        break;
    case LayoutQuincuncial:
        R = area_unit_square / area_unit_sphere;
        break;
    }
    for (int iy = 0; iy < h; iy++) {
        double y = ((iy + 0.5) / h) * 2 - 1;
        for (int ix = 0; ix < w; ix++) {
            double x = ((ix + 0.5) / h) * 2 - 1;
            double lat, lon;
            ctx.inverse(x, y, lat, lon);
            if (std::isfinite(lat) && std::isfinite(lon)) {
                double h, k, theta_prime, a, b, omega, s;
                ctx.analysis(lat, lon, h, k, theta_prime, a, b, omega, s);
                dist_area[iy * w + ix] = s / R;
                dist_isot [iy * w + ix] = a / b;
            } else {
                dist_area[iy * w + ix] = NAN;
                dist_isot [iy * w + ix] = NAN;
            }
        }
    }
}

void scan_array(const double* data, int size,
        double& min, double& min_1percent,
        double& max, double& max_1percent,
        double& mean, double& median)
{
    std::vector<double> valid_data;
    valid_data.reserve(size);
    min = max = NAN;
    mean = 0;
    for (int i = 0; i < size; i++) {
        if (std::isfinite(data[i])) {
            if (!std::isfinite(min) || data[i] < min)
                min = data[i];
            if (!std::isfinite(max) || data[i] > max)
                max = data[i];
            mean += data[i];
            valid_data.push_back(data[i]);
        }
    }
    if (valid_data.size() == 0) {
        mean = NAN;
        min_1percent = NAN;
        max_1percent = NAN;
        median = NAN;
    } else {
        mean /= valid_data.size();
        std::sort(valid_data.begin(), valid_data.end());
        if (valid_data.size() % 2 == 0) {
            median = (valid_data[valid_data.size() / 2 - 1]
                    + valid_data[valid_data.size() / 2]) / 2;
        } else {
            median = valid_data[valid_data.size() / 2];
        }
        min_1percent = valid_data[valid_data.size() / 100];
        max_1percent = valid_data[valid_data.size() - 1 - valid_data.size() / 100];
    }
}

/*
 * Main function.
 */

int main(void)
{
    const int w = 512;
    const int h = 512;

    double *dist_area = new double[w * h];
    double *dist_isot = new double[w * h];
    double min, max, mean;
    double min_1percent, max_1percent, median;

    // Loop over the projection methods
    for (int method = 0; method < n_projections; method++) {
        Context ctx(projections[method].proj,
                projections[method].params,
                projections[method].center,
                projections[method].layout,
                projections[method].square_method);
        fprintf(stderr, "%s\n", projections[method].name);
        analyze(ctx, w, h, dist_area, dist_isot);

        scan_array(dist_area, w * h, min, min_1percent, max, max_1percent, mean, median);
        fprintf(stderr, "  area dist min=%g min_1percent=%g max=%g max_1percent=%g mean=%g median=%g\n", min, min_1percent, max, max_1percent, mean, median);
        png_vis(std::string(projections[method].name) + "-dist-area", w, h, dist_area, true, 2.5);

        scan_array(dist_isot , w * h, min, min_1percent, max, max_1percent, mean, median);
        fprintf(stderr, "  isot dist min=%g min_1percent=%g max=%g max_1percent=%g mean=%g median=%g\n", min, min_1percent, max, max_1percent, mean, median);
        png_vis(std::string(projections[method].name) + "-dist-isot", w, h, dist_isot, true, 2.5);
    }
    return 0;
}
