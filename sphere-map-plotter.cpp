#include <vector>
#include <cmath>

#include <QCoreApplication>
#include <QPainter>
#include <QFile>
#include <QTextStream>

#include "proj.hpp"

#include "sphere-map-common.hpp"
#include "pdfplot.hpp"


/*
 * Get land mass polygon data from buildin CSV file.
 */

const std::vector<std::vector<std::pair<double, double>>>& get_landmass_polygons()
{
    static std::vector<std::vector<std::pair<double, double>>> polygons;
    if (polygons.size() == 0) {
        /* The land mass polygon data (actually country border polygons) was
         * downloaded from http://thematicmapping.org/downloads/world_borders.php,
         * unzipped in a directory named 'tm', and converted to CSV with
         * ogr2ogr -f CSV tm-world-borders-0.3.csv tm -lco GEOMETRY=AS_WKT
         */
        QFile f(":tm-world-borders-0.3.csv");
        if (!f.open(QIODevice::ReadOnly)) {
            fprintf(stderr, "cannot load land polygon data\n");
            std::exit(1);
        }
        QTextStream in(&f);
        in.readLine(); // ignore first line 'WKT,...'
        while (!in.atEnd()) {
            std::vector<QStringList> line_polygons_as_pair_lists;
            QString line = in.readLine();
            if (line.startsWith("\"POLYGON")) {
                line.remove(0, 11); // ignore '"POLYGON (('
                line.remove(line.indexOf(')'), (1 << 30)); // ignore ')",...' THIS ALSO THROWS AWAY INTERIOR POLYGONS!
                QStringList pair_list = line.split(',');
                line_polygons_as_pair_lists.push_back(pair_list);
            } else if (line.startsWith("\"MULTIPOLYGON")) {
                line.remove(0, 17); // ignore '"MULTIPOLYGON ((('
                line.remove(line.indexOf(")))"), 256); // ignore ')))",<gedoense>'
                QStringList polygon_list = line.split(")),((");
                for (int i = 0; i < polygon_list.size(); i++) {
                    polygon_list[i].remove(polygon_list[i].indexOf(')'), (1 << 30)); // THIS THROWS AWAY INTERIOR POLYGONS!
                    QStringList pair_list = polygon_list[i].split(',');
                    line_polygons_as_pair_lists.push_back(pair_list);
                }
            }
            for (size_t i = 0; i < line_polygons_as_pair_lists.size(); i++) {
                std::vector<std::pair<double, double>> polygon;
                polygon.reserve(line_polygons_as_pair_lists[i].size());
                for (int j = 0; j < line_polygons_as_pair_lists[i].size(); j++) {
                    QStringList pair = line_polygons_as_pair_lists[i][j].split(' ');
                    double lon = radians(pair[0].toDouble());
                    double lat = radians(pair[1].toDouble());
                    polygon.push_back(std::pair<double, double>(lat, lon));
                }
                polygons.push_back(polygon);
            }
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

    const std::vector<std::vector<std::pair<double, double>>>& landmass_polygons = get_landmass_polygons();

    const int lat_lines = 18;
    const int lon_lines = 36;
    const double map_sub_point_dist = 1e-3;
    const int line_steps = 5000;

    std::vector<std::pair<float, float>> point_list;
    double lat, lon;
    double x, y;

    pen.setCapStyle(Qt::FlatCap);

    // Loop over the map methods
    for (int method = 0; method < n_projections; method++) {
        // Set up the PDF painter and the projection context
        std::string id = std::string("plot-") + projections[method].name;
        fprintf(stderr, "%s\n", id.c_str());
        Context ctx(projections[method].proj,
                projections[method].params,
                projections[method].center,
                projections[method].layout,
                projections[method].square_method);
        make_pdf_painter(id,
                (ctx.layout == LayoutHemisphere || ctx.layout == LayoutSphere) && ctx.square_method == SquareMethodNone,
                &writer, &painter);
        // Plot the landmass map
        pen.setColor(QColor(128, 128, 128));
        pen.setStyle(Qt::SolidLine);
        pen.setWidth(1);
        painter->setPen(pen);
        brush.setColor(QColor(128, 128, 128));
        brush.setStyle(Qt::SolidPattern);
        painter->setBrush(brush);
        for (size_t i = 0; i < landmass_polygons.size(); i++) {
            const std::vector<std::pair<double, double>>& polygon = landmass_polygons[i];
            // Check if this needs special handling of the south pole in the
            // quincuncial layout.
            // This is the case if the polygon intersects more than one of the four
            // parts that represent the southern hemisphere.
            bool special_handling = false;
            if (ctx.layout == LayoutQuincuncial) {
                unsigned int corner_seen = 0;
                for (size_t j = 0; j < polygon.size(); j++) {
                    lat = polygon[j].first;
                    lon = polygon[j].second;
                    ctx.forward(lat, lon, x, y);
                    if (std::fabs(x) + fabs(y) <= 1.0)
                        break;
                    unsigned int corner;
                    if (x > 0)
                        corner = (y > 0 ? (1 << 0) : (1 << 3));
                    else
                        corner = (y > 0 ? (1 << 1) : (1 << 2));
                    if (corner_seen & ~corner) {
                        special_handling = true;
                        break;
                    }
                    corner_seen = corner;
                }
            }
            point_list.clear();
            // find first valid point
            size_t j;
            for (j = 0; j < polygon.size(); j++) {
                lat = polygon[j].first;
                lon = polygon[j].second;
                ctx.forward(lat, lon, x, y);
                if (point_is_valid(x, y)) {
                    point_list.push_back(std::pair<float, float>(x, y));
                    break;
                }
            }
            size_t first_valid = j;
            // find next valid point
            double last_lat = lat;
            double last_lon = lon;
            double last_x = x;
            double last_y = y;
            for (; j <= polygon.size(); j++) {
                size_t jj = (j < polygon.size() ? j : first_valid);
                lat = polygon[jj].first;
                lon = polygon[jj].second;
                ctx.forward(lat, lon, x, y);
                if (point_is_valid(x, y)) {
                    // plot the line from last point to current point
                    double dist = hypot(x - last_x, y - last_y);
                    if (dist < map_sub_point_dist)
                        continue;
                    int substeps = dist / map_sub_point_dist;
                    for (int k = 1; k < substeps; k++) {
                        double alpha = static_cast<double>(k) / substeps;
                        double sublat = alpha * lat + (1.0 - alpha) * last_lat;
                        double sublon = alpha * lon + (1.0 - alpha) * last_lon;
                        double subx, suby;
                        ctx.forward(sublat, sublon, subx, suby);
                        if (point_is_valid(subx, suby))
                            point_list.push_back(std::pair<float, float>(subx, suby));
                    }
                    if (jj != first_valid)
                        point_list.push_back(std::pair<float, float>(x, y));
                    last_lat = lat;
                    last_lon = lon;
                    last_x = x;
                    last_y = y;
                }
            }
            if (special_handling) {
                for (unsigned int corner = 0; corner < 4; corner++) {
                    std::vector<std::pair<float, float>> special_point_list;
                    for (size_t j = 0; j < point_list.size(); j++) {
                        x = point_list[j].first;
                        y = point_list[j].second;
                        double t;
                        // first, map corners back to center
                        // see inverse quincuncial mapping
                        if (x * y >= 0) {
                            if (x >= 0 && y >= 0) {
                                t = -y + 1;
                                y = -x + 1;
                                x = t;
                            } else {
                                t = -y - 1;
                                y = -x - 1;
                                x = t;
                            }
                        } else {
                            if (x < 0) {
                                t = y - 1;
                                y = x + 1;
                                x = t;
                            } else {
                                t = y + 1;
                                y = x - 1;
                                x = t;
                            }
                        }
                        // then map the whole polygon into the single current corner
                        // (see forward quincuncial mapping)
                        // let the plotting library do the appropriate clipping
                        if (corner == 0) {
                            t = -y + 1;
                            y = -x + 1;
                            x = t;
                        } else if (corner == 1) {
                            t = y - 1;
                            y = x + 1;
                            x = t;
                        } else if (corner == 2) {
                            t = -y - 1;
                            y = -x - 1;
                            x = t;
                        } else {
                            t = y + 1;
                            y = x - 1;
                            x = t;
                        }
                        special_point_list.push_back(std::pair<float, float>(x, y));
                    }
                    plot_points(painter, false, special_point_list);
                }
            } else {
                plot_points(painter, false, point_list);
            }
        }
        // Plot the lat/lon grid
        pen.setColor(QColor(0, 0, 0));
        pen.setStyle(Qt::SolidLine);
        painter->setPen(pen);
        brush.setStyle(Qt::NoBrush);
        painter->setBrush(brush);
        for (int ilat = 0; ilat <= lat_lines; ilat++) {
            lat = M_PI_2 - ilat * M_PI / lat_lines;
            point_list.clear();
            for (int i = 0; i <= line_steps; i++) {
                lon = i * 2.0 * M_PI / line_steps;
                lon -= M_PI;
                ctx.forward(lat, lon, x, y);
                if (point_is_valid(x, y)) {
                    if (point_list.size() > 0) {
                        float dist = hypot(point_list.back().first - x, point_list.back().second - y);
                        if (dist < map_sub_point_dist)
                            continue;
                    }
                    point_list.push_back(std::pair<float, float>(x, y));
                }
            }
            pen.setWidth((ilat % 3 == 0) ? 8 : 2);
            painter->setPen(pen);
            plot_points(painter, true, point_list);
        }
        for (int ilon = 0; ilon <= lon_lines; ilon++) {
            lon = -M_PI + ilon * 2.0 * M_PI / lon_lines;
            point_list.clear();
            for (int i = 0; i <= line_steps; i++) {
                if (i == 0 || i == line_steps)
                    continue;
#if 0
                if (ctx.layout == LayoutSphere) {
                    // latitude range starts 5 deg from center pole and ends at outer pole
                    if (ctx.center == CenterNorthPole)
                        lat = -M_PI_2 + i * (M_PI - M_PI / 36) / line_steps;
                    else
                        lat = -M_PI_2 + M_PI / 36 + i * (M_PI - M_PI / 36) / line_steps;
                } else {
                    // latitude range starts/stops 5 deg from each pole
                    lat = -M_PI_2 + M_PI / 36 + i * (M_PI - M_PI / 18) / line_steps;
                }
#endif
                lat = -M_PI_2 + i * M_PI / line_steps;
                ctx.forward(lat, lon, x, y);
                if (point_is_valid(x, y)) {
                    if (point_list.size() > 0) {
                        float dist = hypot(point_list.back().first - x, point_list.back().second - y);
                        if (dist < map_sub_point_dist)
                            continue;
                    }
                    point_list.push_back(std::pair<float, float>(x, y));
                }
            }
            pen.setWidth((ilon % 3 == 0) ? 8 : 2);
            painter->setPen(pen);
            plot_points(painter, true, point_list);
        }
        finish_pdf_painter(writer, painter);
    }

    return 0;
}
