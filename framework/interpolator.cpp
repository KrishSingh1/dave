#include <functional>

#include "interpolator.h"
#include "gutil.h"
#include "shell.h"

// Build a interpolator out of the vectors 'x' and 'y'.
void Interpolator::init(const std::vector<double>& x_,
			            const std::vector<double>& y_) {
    if (y_.size() != x_.size()) {
        xtring buf;
        buf.printf("Cannot initialise interpolator: x length (%d) is not equal to y length (%d)", x.size(), y.size());
        fail(buf);
    }
    if (x_.size() < 3) {
        fail("insufficient number of points");
    }
    x = x_;
    y = y_;
    initialise();
    set_limits();
}

void Interpolator::set_limits() {
    _size = x.size();
    _min = _size > 0 ? x.front() :  1e36;
    _max = _size > 0 ? x.back()  : -1e36;
}

// Compute the interpolated function from the points contained in 'x' and 'y'.
void Interpolator::initialise() {
    // https://stackoverflow.com/questions/17769114/stdis-sorted-and-strictly-less-comparison
    if (!std::is_sorted(x.begin(), x.end(), std::less_equal<double>())) {
        fail("spline control points must be unique and in ascending order");
    }
    if (x.size() > 0) {
        tk_spline.set_points(x, y);
        active = true;
    }
}

// Support for adding points in turn (assumes monotonic increasing in
// 'x', unchecked).
void Interpolator::add_point(double xi, double yi) {
    x.push_back(xi);
    y.push_back(yi);

    set_limits();
}

// adds point in sorted position (slower than above)
void Interpolator::add_point_sorted(double xi, double yi) {
    auto x_upper = std::upper_bound(x.begin(), x.end(), xi); // find smallest number larger than xi
    x.insert(x_upper, xi); // add xi below that number
    auto y_upper = std::upper_bound(y.begin(), y.end(), yi);
    y.insert(y_upper, yi);

    set_limits();
}

// Remove all the contents, being ready to be refilled.
void Interpolator::clear() {
    x.clear();
    y.clear();
    active = false;

    set_limits();
}

bool Interpolator::try_eval(double u, double &y) const {
    check_active();
    if (!extrapolate && (u < _min || u > _max)) {
        return false;
    }
    y = tk_spline(u);
    return true;
}

// Compute the value of the interpolated function at point `x=u`
double Interpolator::eval(double u) const {
    check_active();
    if (!extrapolate && (u < _min || u > _max)) {
        fail("Extrapolation disabled and evaluation point outside of interpolated domain.");
    }
    return tk_spline(u);
}

// faster version of above
double Interpolator::operator()(double u) const {
    return tk_spline(u);
}

void Interpolator::set_extrapolate(bool e) {
    extrapolate = e;
}

std::vector<double> Interpolator::get_x() const {
    return x;
}

std::vector<double> Interpolator::get_y() const {
    return y;
}

// // Get the (x,y) pairs in the Interpolator as a two-column matrix
// SEXP Interpolator::r_get_xy() const {
//     std::vector< std::vector<double> > xy;
//     xy.push_back(x);
//     xy.push_back(y);
//     return Rcpp::wrap(util::to_rcpp_matrix(xy));
// }

// Compute the value of the interpolated function at a vector of
// points `x=u`, returning a vector of the same length.
std::vector<double> Interpolator::r_eval(std::vector<double> u) const {
    check_active();
    auto ret = std::vector<double>();
    ret.reserve(u.size()); // fast to do this once rather than multiple times with push_back
    for (auto const& x : u) {
        ret.push_back(eval(x));
    }
    return ret;
}

void Interpolator::check_active() const {
    if (!active) {
        fail("Interpolator not initialised -- cannot evaluate");
    }
}
