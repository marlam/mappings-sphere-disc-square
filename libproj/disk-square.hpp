#ifndef PROJ_DISC_SQUARE_H
#define PROJ_DISC_SQUARE_H

namespace Projection {

void disk_to_square_stretch(double u, double v, double& x, double& y);
void square_to_disk_stretch(double x, double y, double& u, double& v);

void disk_to_square_shirley(double u, double v, double& x, double& y);
void square_to_disk_shirley(double x, double y, double& u, double& v);

void disk_to_square_squircle(double u, double v, double& x, double& y);
void square_to_disk_squircle(double x, double y, double& u, double& v);

void disk_to_square_elliptical(double u, double v, double& x, double& y);
void square_to_disk_elliptical(double x, double y, double& u, double& v);

void disk_to_square_conformal(double u, double v, double& x, double& y);
void square_to_disk_conformal(double x, double y, double& u, double& v);

}

#endif
