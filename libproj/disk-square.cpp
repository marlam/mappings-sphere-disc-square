#include <algorithm>
#include <cmath>

#include "disk-square.hpp"

// Define exactly one of the following to 1 to choose which numerical methods to use
// for the conformal mapping.
// The implementation of Stark's method still has bugs, see comments there
#define WITH_MYOWN 1
#define WITH_STARK 0

#if WITH_MYOWN
# include <cfloat>
#elif WITH_STARK
# include <complex>
  using std::complex;
  using std::polar;
#endif

namespace Projection {

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


/* Radius Stretching */

void disk_to_square_stretch(double u, double v, double& x, double& y)
{
    double r = hypot(u, v);
    if (r <= 0.0) {
        x = 0.0;
        y = 0.0;
    } else {
        x = sign(u) * r;
        y = sign(v) * r;
        if (u * u >= v * v)
            y *= v / u;
        else
            x *= u / v;
    }
}

void square_to_disk_stretch(double x, double y, double& u, double& v)
{
    double t = hypot(x, y);
    if (t <= 0.0) {
        u = 0.0;
        v = 0.0;
    } else {
        double c = (x * x >= y * y ? x : y);
        u = sign(x) * c * x / t;
        v = sign(y) * c * y / t;
    }
}

/* Shirley's equal-area method, described in Shirley, P., & Chiu, K. (1997).
 * A low distortion map between disk and square. Journal of graphics tools,
 * 2(3), 45-52.
 *
 * Note that this is equivalent to the independently derived mapping described
 * in Roşca, D. "New uniform grids on the sphere." Astronomy & Astrophysics
 * 520 (2010): A63.
 *
 * For the square-to-disk functino, we use a trick to simplify the original
 * equations, first introduced by Dave Cline here:
 * http://psgraphics.blogspot.de/2011/01/improved-code-for-concentric-map.html
 */

void disk_to_square_shirley(double u, double v, double& x, double& y)
{
    double r = hypot(u, v);
    double phi = atan2(v, u);
    if (phi < -M_PI_4)
        phi += 2 * M_PI;
    if (phi < M_PI_4) {
        x = r;
        y = phi * x / M_PI_4;
    } else if (phi < M_PI_2 + M_PI_4) {
        y = r;
        x = -(phi - M_PI_2) * y / M_PI_4;
    } else if (phi < M_PI + M_PI_4) {
        x = -r;
        y = (phi - M_PI) * x / M_PI_4;
    } else {
        y = -r;
        x = -(phi - (M_PI + M_PI_2)) * y / M_PI_4;
    }
}

void square_to_disk_shirley(double x, double y, double& u, double& v)
{
    double r, phi;
    if (x * x > y * y) {
        r = x;
        phi = M_PI_4 * y / x;
    } else {
        r = y;
        if (fabs(y) > 0.0)
            phi = M_PI_2 - M_PI_4 * x / y;
        else
            phi = 0.0;
    }
    u = r * cos(phi);
    v = r * sin(phi);
}

/* Method using Guasti's Squircle: map concentric circles to concentric
 * squircles that become more square-like with growing size.
 * Details of the mapping are described in Fong, C. (2014). An Indoor
 * Alternative to Stereographic Spherical Panoramas. In Proceedings of
 * Bridges 2014: Mathematics, Music, Art, Architecture, Culture
 * (pp. 103-110).
 *
 * This mapping is neither conformal nor equal-area. */

void disk_to_square_squircle(double u, double v, double& x, double& y)
{
    double u2 = u * u;
    double v2 = v * v;
    double w = sign(u * v) / M_SQRT2 * sqrt(max(u2 + v2 - sqrt((u2 + v2) * (u2 + v2 - 4 * u2 * v2)), 0.0));
    if (fabs(w) <= 0) {
        x = u;
        y = v;
    } else {
        x = w / v;
        y = w / u;
    }
}

void square_to_disk_squircle(double x, double y, double& u, double& v)
{
    double x2 = x * x;
    double y2 = y * y;
    double t = x2 + y2;
    if (t <= 0) {
        u = v = 0;
    } else {
        double r = sqrt(max(t - x2 * y2, 0.0));
        double q = r / sqrt(t);
        u = x * q;
        v = y * q;
    }
}

/* Method using elliptical arc mapping.
 * This mapping is neither conformal nor equal-area. */

void disk_to_square_elliptical(double u, double v, double& x, double& y)
{
    x =   0.5 * sqrt(max(2 + u * u - v * v + u * 2 * M_SQRT2, 0.0))
        - 0.5 * sqrt(max(2 + u * u - v * v - u * 2 * M_SQRT2, 0.0));
    y =   0.5 * sqrt(max(2 - u * u + v * v + v * 2 * M_SQRT2, 0.0))
        - 0.5 * sqrt(max(2 - u * u + v * v - v * 2 * M_SQRT2, 0.0));
}

void square_to_disk_elliptical(double x, double y, double& u, double& v)
{
    u = x * sqrt(1 - y * y / 2);
    v = y * sqrt(1 - x * x / 2);
}

/* Conformal mapping based on Schwarz-Christoffel transformation. */

static const double K = 1.854074677301371918433850347195260046217598823521766905586;

#if WITH_STARK

/* This optimized implementation is based on Stark, M. M. (2009). Fast and
 * Stable Conformal Mapping Between a Disc and a Square. Journal of Graphics,
 * GPU, and Game Tools, 14(2), 1-23.
 *
 * XXX: this has errors in the disk-to-square function that I was unable to find!
 * Regions 2 and 3 are affected.
 */

static complex<double> stark_sum_series(complex<double> z, int n)
{
    complex<double> s = z;
    complex<double> p = z;
    for (int k = 1; k <= n; k++) {
        p = p * z * z * z * z * (2 * k - 1.0) / (2.0 * k);
        s = s + p / (4 * k + 1.0);
    }
    return s;
}

static double stark_F(double t)
{
    const int n = 4; // use n=3 for single precision
    double tt = t * t;
    double ttt = tt * t;
    double F0 = t - 0.193 * ttt + 0.026 * ttt * tt;
    double a = 1;
    double b = M_SQRT1_2;
    for (int k = 1; k <= n; k++) {
        t = (a + b) * t / (a - b * tt);
        a = (a + b) / 2;
        b = sqrt(a * b);
    }
    double phi_n = atan(t);
    return phi_n + M_PI * floor(0.5 + ((1 << n) * a * F0 - phi_n) / M_PI) / ((1 << n) * a);
}

void disk_to_square_conformal(double u, double v, double& x, double& y)
{
    double r = hypot(u, v);
    double phi = atan2(v, u);
    if (phi <= -M_PI)
        phi = M_PI;

    const double r0 = 0.6436; // sqrt(sqrt(2.0) - 1.0) + epsilon
    complex<double> z, w;
    if (r < r0) {
        z = polar(r, phi + M_PI_4);
        w = 2.0 * sqrt(complex<double>(0, -1)) / K * stark_sum_series(z, floor(6.5 / (1 - r) - 2));
    } else {
        double Q = floor(phi / M_PI_2) + 0.5;
        z = polar(r, phi - Q * M_PI_2);
        complex<double> xi = sqrt((1.0 - z * z) / (1.0 + z * z));
        if (abs(xi) < r0) {
            w = K / M_SQRT2 - stark_sum_series(xi, floor(6.5 / (1 - abs(xi)) - 2));
        } else {
//fprintf(stderr, "X4: Q=%g r=%g phi=%g |z|=%g re(z)=%g im(z)=%g\n", Q, r, phi * RAD_TO_DEG, cabs(z), creal(z), creal(z));
            double A = hypot(z.real() + 1, z.imag());
            double B = hypot(z.real() - 1, z.imag());
//fprintf(stderr, "    A=%g B=%g\n", A, B);
            double alpha = (A + B) / 2;
            double beta = (A - B) / 2;
//fprintf(stderr, "    alpha=%g beta=%g\n", alpha, beta);
            double gamma = alpha + sqrt(alpha * alpha - 1);
            double Tphi = (1 - beta * beta) / (beta * beta);
            double b = 2 / Tphi + (gamma - 1 / gamma) * (gamma - 1 / gamma) / (4 * (1 - beta * beta)) - 1;
//fprintf(stderr, "    gamma=%g Tphi=%g b=%g\n", gamma, Tphi, b);
            double Clambda = (b + sqrt(b * b + 8 / Tphi)) / 4;
            double Tmu = 2 * (Tphi * Clambda - 1);
//fprintf(stderr, "    Clambda=%g Tmu=%g\n", Clambda, Tmu);
            w = complex<double>(K - stark_F(sqrt(1 / Clambda)), -stark_F(sqrt(Tmu)));
            if (z.imag() > 0)
                w = conj(w);
//fprintf(stderr, "    |w|=%g re(w)=%g im(w)=%g\n", cabs(w), creal(w), cimag(w));
        }
        //w = polar(abs(w), arg(w) + Q * M_PI_2) * (2 / K);
        w *= polar(1.0, Q * M_PI_2) * (2 / K);
//fprintf(stderr, "    |w|=%g re(w)=%g im(w)=%g\n", cabs(w), creal(w), cimag(w));
    }
    x = w.real();
    y = w.imag();
}

static complex<double> stark_inv_int(complex<double> w)
{
    const double C0 = 13.750371636040937;
    const double c1 =  -1.6360491363469976e-1;
    const double c5 =  -1.5316508620083077e-3;
    const double c9 =  +5.9455890307966153e-7;
    const double c13 = +1.7520282395125552e-8;
    const double c17 = +2.8997255626121623e-11;
    const double c21 = -1.4786423004927015e-13;
    const double c25 = -6.1034023099548599e-16;
    const double c29 = +5.3527850055041484e-19;
    const double c33 = +7.9773298274614004e-21;
    const double c37 = +1.1683926152311516e-23;
    const complex<double> w2 = w * w;
    const complex<double> w4 = w2 * w2;
    const complex<double> w8 = w4 * w4;
    const complex<double> w12 = w8 * w4;
    const complex<double> w16 = w8 * w8;
    const complex<double> w20 = w16 * w4;
    const complex<double> w24 = w16 * w8;
    const complex<double> w28 = w16 * w12;
    const complex<double> w32 = w16 * w16;
    const complex<double> w36 = w32 * w4;
    return C0 * w / (K * K * K * K + w4)
        + w * (c1 + c5 * w4 + c9 * w8 + c13 * w12 + c17 * w16
                + c21 * w20 + c25 * w24 + c29 * w28 + c33 * w32 + c37 * w36);
}

void square_to_disk_conformal(double x, double y, double& u, double& v)
{
    complex<double> w = complex<double>(x, y);
    complex<double> z = 1.0 / sqrt(complex<double>(0, 1))
        * stark_inv_int((K * w) / (2.0 * sqrt(complex<double>(0, -1))));
    u = z.real();
    v = z.imag();
}

#else

static const long long myown_max_iter = 63;
// compute jacobi elliptic functions sn, cn, dn for fixed modulus k=1/sqrt(2)
// via Arithmetic-Geometric Mean (AGM), Abramowitz & Stegun 16.4
static void myown_jacobi_sn_cn_dn(double u, double* sn, double* cn, double* dn)
{
    static const long long myown_max_iter = 64;
    double a[myown_max_iter+1], g[myown_max_iter+1], c[myown_max_iter+1];
    a[0] = 1;
    g[0] = M_SQRT1_2;
    c[0] = M_SQRT1_2;
    long long i = 0;
    do {
        a[i+1] = 0.5 * (a[i] + g[i]);
        g[i+1] = sqrt(a[i] * g[i]);
        c[i+1] = 0.5 * (a[i] - g[i]);
        i++;
    } while (i < myown_max_iter && fabs(a[i]-g[i]) > DBL_EPSILON);
    if (i == myown_max_iter)
        fprintf(stderr, "WARNING: max iterations in myown_jacobi_sn_cn_dn()\n");

    double phi = (1LL << i) * a[i] * u;
    for (; i > 0; i--)
        phi = 0.5 * (phi + asin(c[i] * sin(phi) / a[i]));

    *sn = sin(phi);
    *cn = cos(phi);
    *dn = sqrt(1 - 0.5 * (*sn * *sn));
}
// compute incomplete elliptic integral F for fixed modulus k=1/sqrt(2)
// via Landen transformation, Abramowitz & Stegun 17.5
static double myown_ellint_F(double phi)
{
    double a = 1;
    double g = M_SQRT1_2;
    double last_a, last_g;
    double tan_2n_phi;

    long long i = 0;
    do {
        tan_2n_phi = tan(phi);
        i++;
        phi += phi - atan((a - g) * tan_2n_phi / (a + g * tan_2n_phi * tan_2n_phi));
        last_a = a;
        last_g = g;
        a = 0.5 * (last_a + last_g);
        g = sqrt(last_a * last_g);
    }
    while (i < myown_max_iter && fabs(a - g) > DBL_EPSILON);

    if (i == myown_max_iter)
        fprintf(stderr, "WARNING: max iterations in myown_ellint_F()\n");

    phi /= (1LL << i);
    return phi / g;
}

static double F(double phi)
{
    return myown_ellint_F(phi);
}

static void ccn(double re, double im, double* ret_re, double* ret_im)
{
    double sn_re, cn_re, dn_re;
    double sn_im, cn_im, dn_im;
    myown_jacobi_sn_cn_dn(re, &sn_re, &cn_re, &dn_re);
    myown_jacobi_sn_cn_dn(im, &sn_im, &cn_im, &dn_im);
    double t = 1 - dn_re * dn_re * sn_im * sn_im;
    *ret_re = cn_re * cn_im / t;
    *ret_im = - sn_re * dn_re * sn_im * dn_im / t;
}

void disk_to_square_conformal(double u, double v, double& x, double& y)
{
    // Rotate (u,v) by pi/4;
    double ru = (u - v) * M_SQRT1_2;
    double rv = (u + v) * M_SQRT1_2;
    // Map (ru,rv) to (rx,ry) on a rotated square.
    double A = ru * ru + rv * rv;
    double B = ru * ru - rv * rv;
    double U = 1 + 2 * B - A * A;
    double T = sqrt((1 + A * A) * (1 + A * A) - 4 * B  * B);
    double cos_a = (2 * A - T) / U;
    double cos_b = U / (2 * A + T);
    double a = acos(min(max(cos_a, -1.0), 1.0));
    double b = acos(min(max(cos_b, -1.0), 1.0));
    double rx = sign(ru) * (1 - F(a) / (2 * K));
    double ry = sign(rv) * (    F(b) / (2 * K));
    // Rotate square by -45 deg and make it fill [-1,+1]^2
    x = rx + ry;
    y = ry - rx;
}

void square_to_disk_conformal(double x, double y, double& u, double& v)
{
    // Rotate by 45 deg and keep the square inside [-1,+1]^2
    double z_re = x / 2 - y / 2;
    double z_im = x / 2 + y / 2;
    // Map z to unit disk
    double ru, rv;
    ccn(K * (1 - z_re), - K * z_im, &ru, &rv);
    // Rotate (u,v) by pi/4;
    u = (ru + rv) * M_SQRT1_2;
    v = (rv - ru) * M_SQRT1_2;
}

#endif
}
