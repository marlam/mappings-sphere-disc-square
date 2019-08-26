#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

#include "proj.hpp"

#include "disk-square-common.hpp"
#include "pngvis.hpp"


using std::isfinite;
using std::min;
using std::max;
using Projection::degrees;
using Projection::radians;

void derivatives(int method, double u, double v,
        double& dx_dv, double& dx_du, double& dy_dv, double& dy_du, double delta = 1e-5)
{
    double mv = v - delta;
    double pv = v + delta;
    double mu = u - delta;
    double pu = u + delta;

    if (mv * mv + mu * mu > 1 || mv * mv + pu * pu > 1
            || pv * pv + mu * mu > 1 || pv * pv + pu * pu > 1)
        goto fail;

    double x_mv_mu, x_mv_pu, x_pv_mu, x_pv_pu;
    double y_mv_mu, y_mv_pu, y_pv_mu, y_pv_pu;

    projections[method].disk_to_square(mu, mv, x_mv_mu, y_mv_mu);
    projections[method].disk_to_square(pu, mv, x_mv_pu, y_mv_pu);
    projections[method].disk_to_square(mu, pv, x_pv_mu, y_pv_mu);
    projections[method].disk_to_square(pu, pv, x_pv_pu, y_pv_pu);

    if (!isfinite(x_mv_mu) || !isfinite(x_mv_pu) || !isfinite(x_pv_mu) || !isfinite(x_pv_pu)
            || !isfinite(y_mv_mu) || !isfinite(y_mv_pu) || !isfinite(y_pv_mu) || !isfinite(y_pv_pu))
        goto fail;

    dx_dv = (- x_pv_pu + x_mv_pu + x_mv_mu - x_pv_mu) / (4 * delta);
    dx_du = (+ x_pv_pu + x_mv_pu - x_mv_mu - x_pv_mu) / (4 * delta);
    dy_dv = (+ y_pv_pu - y_mv_pu - y_mv_mu + y_pv_mu) / (4 * delta);
    dy_du = (- y_pv_pu - y_mv_pu + y_mv_mu + y_pv_mu) / (4 * delta);
    return;

fail:
    dx_dv = dx_du = dy_dv = dy_du = NAN;
}

void analysis(int method, double u, double v,
        double& h, double& k,
        double& theta_prime,
        double& a, double& b,
        double& omega, double& s,
        double delta = 1e-5)
{
    double dx_dv, dx_du, dy_dv, dy_du;

    derivatives(method, u, v, dx_dv, dx_du, dy_dv, dy_du, delta);
    if (!isfinite(dx_dv) || !isfinite(dx_du) || !isfinite(dy_dv) || !isfinite(dy_du)) {
        h = k = theta_prime = a = b = omega = s = NAN;
        return;
    }

    h = hypot(dx_dv, dy_dv);
    k = hypot(dx_du, dy_du);
    double sin_theta_prime = (dy_dv * dx_du - dx_dv * dy_du) / (h * k);
    s = h * k * sin_theta_prime; // equivalent to s = a * b
    theta_prime = asin(min(max(sin_theta_prime, -1.0), +1.0));
    double a_prime = sqrt(max(h * h + k * k + 2 * s, 0.0));
    double b_prime = sqrt(max(h * h + k * k - 2 * s, 0.0));
    a = (a_prime + b_prime) / 2;
    b = (a_prime - b_prime) / 2;
    omega = 2 * asin(min(max(b_prime / a_prime, -1.0), +1.0));
}

void analyze(int method, int w, int h, double* dist_area, double* dist_isot)
{
    for (int iy = 0; iy < h; iy++) {
        double y = ((iy + 0.5) / h) * 2 - 1;
        for (int ix = 0; ix < w; ix++) {
            double x = ((ix + 0.5) / h) * 2 - 1;
            double u, v;
            projections[method].square_to_disk(x, y, u, v);
            if (std::isfinite(u) && std::isfinite(v)) {
                double h, k, theta_prime, a, b, omega, s;
                analysis(method, u, v, h, k, theta_prime, a, b, omega, s);
                dist_area[iy * w + ix] = s / (4 / M_PI);
                dist_isot[iy * w + ix] = a / b;
            } else {
                dist_area[iy * w + ix] = NAN;
                dist_isot[iy * w + ix] = NAN;
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
        fprintf(stderr, "%s\n", projections[method].name);
        analyze(method, w, h, dist_area, dist_isot);

        scan_array(dist_area, w * h, min, min_1percent, max, max_1percent, mean, median);
        fprintf(stderr, "  area dist min=%g min_1percent=%g max=%g max_1percent=%g mean=%g median=%g\n", min, min_1percent, max, max_1percent, mean, median);
        png_vis(std::string(projections[method].name) + "-dist-area", w, h, dist_area, true, 2.5);

        scan_array(dist_isot , w * h, min, min_1percent, max, max_1percent, mean, median);
        fprintf(stderr, "  isot dist min=%g min_1percent=%g max=%g max_1percent=%g mean=%g median=%g\n", min, min_1percent, max, max_1percent, mean, median);
        png_vis(std::string(projections[method].name) + "-dist-isot", w, h, dist_isot, true, 2.5);
    }
    return 0;
}
