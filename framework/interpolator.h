// -*-c++-*-
#ifndef LPJ_GUESS_INTERPOLATOR_H_
#define LPJ_GUESS_INTERPOLATOR_H_

#include <vector>
#include "spline.h"

class Interpolator {
public:
    void init(const std::vector<double>& x_,
              const std::vector<double>& y_);
    void initialise();

    void add_point(double xi, double yi);
    void add_point_sorted(double xi, double yi);
    void clear();

    double eval(double u) const;
    double operator()(double u) const;

    void set_extrapolate(bool e);

    std::vector<double> get_x() const;
    std::vector<double> get_y() const;
    
    /*
    // * R interface
    SEXP r_get_xy() const;
    */
    // change to const& vec?
    std::vector<double> r_eval(std::vector<double> u) const;

    // Compute the value of the interpolated function at point `x=u` The result
    // will be stored in y if the point lies within the defined range, in which
    // case true will be returned. If the point u lies outside the defined
    // range, false will be returned, and the value of y is undefined.
    bool try_eval(double u, double &y) const;
private:
    void check_active() const;
    std::vector<double> x, y;
    spline tk_spline;
    bool active = false;
    bool extrapolate = true;

    void set_limits();

    // The number of (x,y) pairs contained in the Interpolator. This is stored
    // as a member variable as a micro-optimisation.
    size_t _size;

    /*
    Minimum and maximum x values.

    These are chosen so that if an Interpolator is empty, functions looking to
    see if they will fall outside of the covered range will always find they do.
    This is the same principle as R's range(numeric(0)) -> c(Inf, -Inf)
    */
    double _min;
    double _max;
};

#endif // LPJ_GUESS_INTERPOLATOR_H_
