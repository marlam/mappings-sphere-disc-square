#ifndef COMMON_2D_HPP
#define COMMON_2D_HPP

#include "disk-square.hpp"


static void disk_to_square_noop(double u, double v, double& x, double& y)
{
    x = u;
    y = v;
}

static void square_to_disk_noop(double x, double y, double& u, double& v)
{
    u = x;
    v = y;
}

typedef struct {
    const char* name;
    void (*disk_to_square)(double u, double v, double& x, double& y);
    void (*square_to_disk)(double x, double y, double& u, double& v);
} projection_t;

projection_t projections[] = {
    { "original",  disk_to_square_noop,                   square_to_disk_noop                   },
    { "stretch",   Projection::disk_to_square_stretch,    Projection::square_to_disk_stretch    },
    { "shirley",   Projection::disk_to_square_shirley,    Projection::square_to_disk_shirley    },
    { "squircle",  Projection::disk_to_square_squircle,   Projection::square_to_disk_squircle   },
    { "elliptic",  Projection::disk_to_square_elliptical, Projection::square_to_disk_elliptical },
    { "conformal", Projection::disk_to_square_conformal,  Projection::square_to_disk_conformal  },
};

int n_projections = sizeof(projections) / sizeof(projections[0]);

#endif
