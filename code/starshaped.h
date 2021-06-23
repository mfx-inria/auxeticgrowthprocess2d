#ifndef STARSHAPED_H
#define STARSHAPED_H

#include <vector>
#include <limits>
#include <algorithm>

enum InterpolationRadialSpanType{PolarPiecewise, PolarLinear, PolarCubic,  Polygonal};
enum SymmetryType{NoSymmetry, RotationalSymmetry, ReflectionalSymmetry};

struct starshaped_struct {
    std::vector<double> radial_spans;
    InterpolationRadialSpanType interpolation_type;
    unsigned int symmetry_degree;
    double min_radial_span, max_radial_span;
};

void line_through_two_points(const Vec2d& p, const Vec2d& q, double& a, double& b, double& c) {
    a = q.y - p.y;
    b = (p.x - q.x);
    c = -(a * p.x + b * p.y);
}

double get_interpolated_radial_span(Vec2d& direction, starshaped_struct& starshaped) {
    if (direction.x == 0.0 && direction.y == 0.0) { return 0.0;}

    double angle = std::atan2(direction.y, direction.x);
    if (angle < 0.0) angle += 2.0 * M_PI;  // atan2 returns from [-pi, pi], for convenience convert to [0, 2*pi] range
    double angle_pos = angle * static_cast<double>(starshaped.radial_spans.size()) / (2.0 * M_PI);
    int li = static_cast<int>(std::floor(angle_pos));
    int ui = (li + 1) % starshaped.radial_spans.size();
    if (starshaped.interpolation_type == Polygonal) {
        // polygonal interpolation
        double ali = static_cast<double>(li) * 2.0 * M_PI / static_cast<double>(starshaped.radial_spans.size());
        double aui = static_cast<double>(ui) * 2.0 * M_PI / static_cast<double>(starshaped.radial_spans.size());
        Vec2d p1(cos(ali) * starshaped.radial_spans[li], sin(ali) * starshaped.radial_spans[li]);
        Vec2d p2(cos(aui) * starshaped.radial_spans[ui], sin(aui) * starshaped.radial_spans[ui]);
        double la, lb, lc;
        line_through_two_points(p1, p2, la, lb, lc);  // line equation with form ax + by + c = 0
        double rt = (-lc / (Vec2d(la, lb).dotProduct(direction)));  // intersection point of ray and line
        return (direction * rt).length();
    } else {
        // polar interpolation
        double t = angle_pos - static_cast<double>(li);
        if (starshaped.interpolation_type == PolarPiecewise) {
            if (t < 0.5) return starshaped.radial_spans[li];
            return starshaped.radial_spans[ui];
        } else if (starshaped.interpolation_type == PolarLinear) {
            return starshaped.radial_spans[li] + t * (starshaped.radial_spans[ui] - starshaped.radial_spans[li]);
        } else if (starshaped.interpolation_type == PolarCubic) {
            // polar cubic hermite interpolation with tangent = 0  (no overshoot) https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
            return (2.0 * t * t * t - 3.0 * t * t + 1.0) * starshaped.radial_spans[li]  + (-2.0 * t * t * t + 3.0 * t * t) * starshaped.radial_spans[ui];
        }
    }
    return 0.0;
}

double get_distance(const Vec2d& coordinate, const Vec2d& site, starshaped_struct& starshaped) {
    Vec2d direction = coordinate - site;
    return direction.length() / get_interpolated_radial_span(direction, starshaped);
}

double get_distance_squared(const Vec2d& coordinate, const Vec2d& site, starshaped_struct& starshaped) {
    Vec2d direction = coordinate - site;
    double interpolated_radial_span = get_interpolated_radial_span(direction, starshaped);
    return direction.length_squared() / (interpolated_radial_span * interpolated_radial_span);
}

void update_min_max_radial_spans(starshaped_struct& starshaped) {
    starshaped.min_radial_span = std::numeric_limits<double>::max();
    starshaped.max_radial_span = -std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < starshaped.radial_spans.size(); i++) {
        double radial_span = starshaped.radial_spans[i];
        starshaped.min_radial_span = std::min(starshaped.min_radial_span, radial_span);
        starshaped.max_radial_span = std::max(starshaped.max_radial_span, radial_span);
    }
}

void symmetrize_radial_spans(starshaped_struct& starshaped, SymmetryType& symmetry_type, unsigned int symmetry_degree) {
    if (symmetry_type == NoSymmetry) return;
    int ori_radial_spans_size = starshaped.radial_spans.size();
    for (unsigned int d = 1; d < symmetry_degree; d++) {
        for (int j = 0; j < ori_radial_spans_size; j++) {
            if (symmetry_type == RotationalSymmetry || (symmetry_type == ReflectionalSymmetry && (j == 0 || (d % 2 == 0) ))) {
                starshaped.radial_spans.push_back(starshaped.radial_spans[j]);
            } else {
                starshaped.radial_spans.push_back(starshaped.radial_spans[ori_radial_spans_size - j]);
            }
        }
    }
}

#endif
