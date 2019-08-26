#include <vector>
#include <cmath>

#include <QCoreApplication>
#include <QPainter>

#include "disk-square-common.hpp"
#include "pdfplot.hpp"

/*
 * Get polygon data
 */

// The teapot outlines as polygons were obtained by rendering a teapot into a
// black/white bitmap and then running
// $ potrace -b geojson -a 0 --flat teapot.pbm -o teapot.geojson
const std::vector<std::vector<std::pair<double, double>>>& get_teapot_polygons()
{
    static std::vector<std::vector<std::pair<double, double>>> polygons;
    if (polygons.size() == 0) {
        static const float p0[] = {
            413.5, 478.4, 409.5, 478.0, 400.0, 476.9, 390.5, 475.9, 386.5,
            474.9, 382.5, 473.8, 377.8, 470.7, 373.1, 467.5, 373.0, 463.1,
            373.0, 458.7, 375.2, 455.1, 377.4, 451.5, 385.0, 443.0, 392.5,
            434.5, 394.9, 431.2, 397.2, 427.8, 398.1, 423.9, 399.0, 419.9,
            397.3, 416.8, 395.5, 413.8, 389.9, 411.0, 384.4, 408.2, 378.4,
            406.5, 372.5, 404.9, 348.5, 404.0, 324.5, 403.1, 314.0, 402.5,
            303.5, 401.9, 286.0, 400.4, 268.5, 399.0, 257.0, 397.4, 245.5,
            395.8, 237.5, 394.0, 229.5, 392.1, 226.0, 390.6, 222.5, 389.2,
            220.0, 386.8, 217.4, 384.5, 214.1, 379.6, 210.8, 374.8, 204.9,
            362.5, 199.1, 350.3, 155.3, 349.7, 111.5, 349.1, 95.0, 347.6,
            78.6, 346.1, 68.0, 344.5, 57.5, 342.9, 47.0, 339.8, 36.5,
            336.8, 28.5, 332.4, 20.5, 328.1, 15.8, 323.3, 11.1, 318.5, 8.0,
            311.7, 5.0, 305.0, 3.9, 297.2, 2.8, 289.5, 3.4, 279.6, 4.1,
            269.7, 6.1, 261.1, 8.1, 252.5, 12.1, 242.7, 16.1, 232.9, 21.7,
            223.7, 27.2, 214.5, 34.1, 206.0, 40.9, 197.5, 49.7, 188.9,
            58.5, 180.3, 67.0, 173.4, 75.5, 166.4, 86.7, 158.8, 97.9,
            151.1, 111.7, 143.2, 125.5, 135.3, 133.7, 131.2, 141.8, 127.1,
            143.5, 121.9, 145.2, 116.8, 149.3, 109.6, 153.4, 102.5, 159.1,
            95.5, 164.8, 88.5, 172.6, 80.9, 180.5, 73.3, 195.3, 61.6,
            210.1, 49.8, 212.4, 44.6, 214.7, 39.5, 218.7, 36.4, 222.6,
            33.2, 230.0, 29.5, 237.3, 25.8, 247.6, 22.4, 257.9, 19.1,
            272.7, 16.0, 287.5, 12.9, 300.2, 10.9, 312.9, 9.0, 327.2, 7.5,
            341.5, 5.9, 365.0, 4.5, 388.5, 3.0, 421.5, 3.0, 454.5, 3.0,
            478.0, 4.5, 501.5, 6.0, 515.5, 7.4, 529.5, 8.9, 540.0, 10.5,
            550.5, 12.0, 566.0, 15.1, 581.5, 18.2, 592.0, 21.3, 602.5,
            24.3, 611.0, 28.2, 619.5, 32.1, 624.6, 36.1, 629.8, 40.1,
            631.8, 44.8, 633.9, 49.5, 637.2, 52.6, 640.5, 55.7, 651.2,
            63.8, 661.8, 72.0, 670.0, 79.7, 678.2, 87.4, 683.9, 94.4,
            689.7, 101.5, 694.0, 109.0, 698.4, 116.5, 701.1, 124.0, 703.7,
            131.5, 709.3, 133.1, 714.9, 134.8, 732.5, 144.6, 750.2, 154.4,
            762.2, 167.2, 774.1, 179.9, 781.9, 194.2, 789.7, 208.5, 795.3,
            225.5, 800.9, 242.5, 803.5, 253.0, 806.1, 263.5, 811.6, 282.3,
            817.2, 301.0, 824.0, 315.3, 830.9, 329.5, 841.3, 341.5, 851.7,
            353.5, 866.1, 362.1, 880.5, 370.8, 889.1, 374.8, 897.7, 378.9,
            898.9, 380.4, 900.2, 381.9, 897.5, 382.6, 894.9, 383.2, 865.0,
            382.1, 835.2, 381.0, 823.8, 379.9, 812.5, 378.9, 806.1, 375.9,
            799.7, 373.0, 791.2, 364.7, 782.7, 356.5, 777.5, 348.5, 772.3,
            340.4, 768.2, 331.0, 764.1, 321.5, 759.0, 307.5, 753.9, 293.5,
            748.2, 282.1, 742.5, 270.8, 736.0, 263.3, 729.5, 255.8, 721.1,
            250.3, 712.6, 244.8, 703.2, 241.9, 693.8, 239.0, 693.1, 239.0,
            692.4, 239.0, 690.0, 245.8, 687.7, 252.5, 682.0, 267.5, 676.2,
            282.5, 669.2, 298.5, 662.1, 314.5, 652.7, 334.5, 643.3, 354.5,
            637.0, 367.2, 630.6, 380.0, 625.9, 384.8, 621.1, 389.7, 614.3,
            391.9, 607.5, 394.0, 596.0, 395.9, 584.5, 397.9, 568.5, 399.4,
            552.5, 400.9, 536.0, 402.0, 519.5, 403.1, 495.5, 404.0, 471.5,
            404.8, 465.0, 406.8, 458.4, 408.7, 454.3, 410.8, 450.2, 413.0,
            448.0, 416.0, 445.8, 419.1, 446.4, 423.3, 447.1, 427.4, 450.9,
            432.5, 454.7, 437.5, 460.9, 444.3, 467.0, 451.0, 469.5, 455.0,
            472.0, 459.0, 472.0, 462.9, 472.0, 466.8, 468.5, 469.4, 465.0,
            472.1, 461.5, 473.6, 457.9, 475.1, 447.7, 476.5, 437.5, 477.9,
            427.5, 478.3, 417.5, 478.8, 413.5, 478.4 ,
        };
        static const float p1[] = {
            184.0, 318.7, 184.0, 318.4, 178.1, 305.5, 172.3, 292.5, 165.2,
            275.0, 158.2, 257.5, 153.2, 242.5, 148.2, 227.5, 145.0, 214.5,
            141.8, 201.5, 139.9, 190.0, 138.0, 178.5, 138.0, 173.2, 138.0,
            168.0, 137.5, 168.0, 136.9, 168.0, 132.7, 170.1, 128.5, 172.2,
            119.0, 178.8, 109.5, 185.3, 102.7, 191.3, 95.9, 197.3, 89.4,
            204.0, 82.9, 210.7, 77.1, 218.6, 71.2, 226.5, 66.1, 235.4,
            61.0, 244.3, 57.9, 251.8, 54.8, 259.4, 52.3, 268.9, 49.9,
            278.5, 50.2, 287.1, 50.6, 295.8, 53.3, 299.4, 56.1, 302.9,
            60.3, 305.8, 64.5, 308.6, 70.5, 310.8, 76.5, 312.9, 84.5,
            314.4, 92.5, 315.9, 107.0, 317.0, 121.5, 318.1, 130.0, 318.3,
            138.5, 318.6, 161.2, 318.8, 184.0, 319.0, 184.0, 318.7
        };
        static const float* p[] = { p0, p1 };
        static const size_t n[] = { sizeof(p0) / sizeof(float), sizeof(p1) / sizeof(float) };
        for (size_t i = 0; i < sizeof(p) / sizeof(float*); i++) {
            std::vector<std::pair<double, double>> polygon;
            for (size_t j = 0; j < n[i] / 2; j++) {
                float x = p[i][2 * j];
                float y = p[i][2 * j + 1];
                x = ((x - 450) / 500);
                y = ((y - 250) / 500);
                polygon.push_back(std::pair<double, double>(x, y));
            }
            polygons.push_back(polygon);
        }
    }
    return polygons;
}

/*
 * Main function.
 */

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    QPdfWriter* writer;
    QPainter* painter;
    QPen pen;
    QBrush brush;

    const std::vector<std::vector<std::pair<double, double>>>& polygons = get_teapot_polygons();

    const int r_lines = 12;
    const int phi_lines = 24;
    const int x_lines = 18;
    const int y_lines = 18;
    const int steps = 1000;
    const double map_sub_point_dist = 1e-3;
    std::vector<std::pair<float, float>> point_list;
    double r, phi, u, v, x, y;

    pen.setCapStyle(Qt::FlatCap);
    pen.setStyle(Qt::SolidLine);
    pen.setWidth(4);
    brush.setColor(QColor(224, 224, 224));
    brush.setStyle(Qt::NoBrush);

    for (int method = 0; method < n_projections; method++) {

        /* Disk -> Square */

        std::string disk_square_id = std::string("plot-disk-to-square-") + projections[method].name;
        fprintf(stderr, "%s\n", disk_square_id.c_str());
        make_pdf_painter(disk_square_id, false, &writer, &painter);

        // Polygons
        pen.setStyle(Qt::NoPen);
        painter->setPen(pen);
        brush.setStyle(Qt::SolidPattern);
        painter->setBrush(brush);
        for (size_t i = 0; i < polygons.size(); i++) {
            const std::vector<std::pair<double, double>>& polygon = polygons[i];
            point_list.clear();
            if (i % 2 == 0)
                brush.setColor(QColor(128, 128, 128));
            else
                brush.setColor(QColor(255, 255, 255));
            painter->setBrush(brush);
            // first valid point
            u = polygon[0].first;
            v = polygon[0].second;
            projections[method].disk_to_square(u, v, x, y);
            point_list.push_back(std::pair<float, float>(x, y));
            double last_u = u;
            double last_v = v;
            double last_x = x;
            double last_y = y;
            for (size_t j = 1; j <= polygon.size(); j++) {
                size_t jj = (j < polygon.size() ? j : 0);
                u = polygon[jj].first;
                v = polygon[jj].second;
                projections[method].disk_to_square(u, v, x, y);
                double dist = hypot(x - last_x, y - last_y);
                if (dist < map_sub_point_dist)
                    continue;
                int substeps = steps / 2 * dist;
                for (int k = 1; k < substeps; k++) {
                    double alpha = static_cast<double>(k) / substeps;
                    double subu = alpha * u + (1.0 - alpha) * last_u;
                    double subv = alpha * v + (1.0 - alpha) * last_v;
                    double subx, suby;
                    projections[method].disk_to_square(subu, subv, subx, suby);
                    point_list.push_back(std::pair<float, float>(subx, suby));
                }
                if (jj != 0)
                    point_list.push_back(std::pair<float, float>(x, y));
                last_u = u;
                last_v = v;
                last_x = x;
                last_y = y;
            }
            plot_points(painter, false, point_list);
        }

        pen.setStyle(Qt::SolidLine);
        painter->setPen(pen);
        brush.setStyle(Qt::NoBrush);
        painter->setBrush(brush);

        // Lines of constant radius
        for (int i = 0; i < r_lines; i++) {
            point_list.clear();
            r = (i + 1.0) / r_lines;
            for (int j = 0; j <= steps; j++) {
                phi = -M_PI + j * 2 * M_PI / steps;
                u = r * cos(phi);
                v = r * sin(phi);
                projections[method].disk_to_square(u, v, x, y);
                if (point_list.size() > 0) {
                    float dist = hypot(point_list.back().first - x, point_list.back().second - y);
                    if (dist < map_sub_point_dist)
                        continue;
                }
                point_list.push_back(std::pair<float, float>(x, y));
            }
            plot_points(painter, true, point_list);
        }

        // Lines of constant phi
        for (int i = 0; i < phi_lines; i++) {
            point_list.clear();
            phi = -M_PI + i * 2 * M_PI / phi_lines;
            for (int j = 0; j <= steps; j++) {
                r = static_cast<double>(j) / steps;
                u = r * cos(phi);
                v = r * sin(phi);
                projections[method].disk_to_square(u, v, x, y);
                if (point_list.size() > 0) {
                    float dist = hypot(point_list.back().first - x, point_list.back().second - y);
                    if (dist < map_sub_point_dist)
                        continue;
                }
                point_list.push_back(std::pair<float, float>(x, y));
            }
            plot_points(painter, true, point_list);
        }
        finish_pdf_painter(writer, painter);

        /* Square -> Disk */

        std::string square_disk_id = std::string("plot-square-to-disk-") + projections[method].name;
        fprintf(stderr, "%s\n", square_disk_id.c_str());
        make_pdf_painter(square_disk_id, false, &writer, &painter);

        // Polygons
        pen.setStyle(Qt::NoPen);
        painter->setPen(pen);
        brush.setStyle(Qt::SolidPattern);
        painter->setBrush(brush);
        for (size_t i = 0; i < polygons.size(); i++) {
            const std::vector<std::pair<double, double>>& polygon = polygons[i];
            point_list.clear();
            if (i % 2 == 0)
                brush.setColor(QColor(128, 128, 128));
            else
                brush.setColor(QColor(255, 255, 255));
            painter->setBrush(brush);
            // first valid point
            x = polygon[0].first;
            y = polygon[0].second;
            projections[method].square_to_disk(x, y, u, v);
            point_list.push_back(std::pair<float, float>(u, v));
            double last_x = x;
            double last_y = y;
            double last_u = u;
            double last_v = v;
            for (size_t j = 1; j <= polygon.size(); j++) {
                size_t jj = (j < polygon.size() ? j : 0);
                x = polygon[jj].first;
                y = polygon[jj].second;
                projections[method].square_to_disk(x, y, u, v);
                double dist = hypot(u - last_u, v - last_v);
                if (dist < map_sub_point_dist)
                    continue;
                int substeps = steps / 2 * dist;
                for (int k = 1; k < substeps; k++) {
                    double alpha = static_cast<double>(k) / substeps;
                    double subx = alpha * x + (1.0 - alpha) * last_x;
                    double suby = alpha * y + (1.0 - alpha) * last_y;
                    double subu, subv;
                    projections[method].square_to_disk(subx, suby, subu, subv);
                    point_list.push_back(std::pair<float, float>(subu, subv));
                }
                if (jj != 0)
                    point_list.push_back(std::pair<float, float>(u, v));
                last_x = x;
                last_y = y;
                last_u = u;
                last_v = v;
            }
            plot_points(painter, false, point_list);
        }

        pen.setStyle(Qt::SolidLine);
        painter->setPen(pen);
        brush.setStyle(Qt::NoBrush);
        painter->setBrush(brush);

        // Lines of constant x
        for (int i = 0; i <= x_lines; i++) {
            point_list.clear();
            x = -1 + 2.0 * i / x_lines;
            for (int j = 0; j <= steps; j++) {
                y = -1 + 2.0 * j / steps;
                projections[method].square_to_disk(x, y, u, v);
                if (point_list.size() > 0) {
                    float dist = hypot(point_list.back().first - u, point_list.back().second - v);
                    if (dist < map_sub_point_dist)
                        continue;
                }
                point_list.push_back(std::pair<float, float>(u, v));
            }
            plot_points(painter, true, point_list);
        }

        // Lines of constant y
        for (int i = 0; i <= y_lines; i++) {
            point_list.clear();
            y = -1 + 2.0 * i / y_lines;
            for (int j = 0; j <= steps; j++) {
                x = -1 + 2.0 * j / steps;
                projections[method].square_to_disk(x, y, u, v);
                if (point_list.size() > 0) {
                    float dist = hypot(point_list.back().first - u, point_list.back().second - v);
                    if (dist < map_sub_point_dist)
                        continue;
                }
                point_list.push_back(std::pair<float, float>(u, v));
            }
            plot_points(painter, true, point_list);
        }
        finish_pdf_painter(writer, painter);
    }

    return 0;
}
