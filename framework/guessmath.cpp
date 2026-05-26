#include "guessmath.h"
#include "guess.h"

bool solve_quadratic(double a, double b, double c, double &root, bool lower) {

	double b_sq = b * b;
	double four_ac = 4.0 * a * c;
	double discriminant = b_sq - four_ac;

	// Variable tolerance based on magnitude of inputs.
	double discriminant_tol =
		64.0 * DBL_EPSILON * std::max(fabs(b_sq), fabs(four_ac));

    // Roundoff error can lead to a very slightly negative discriminant. In this
	// case, we set the discriminant to zero.
	if (discriminant < 0.0 && discriminant >= -discriminant_tol) {
		discriminant = 0.0;
	}

	if (discriminant < 0) {
		return false;
	}

	if (is_zero(a)) {
		// Equation degenerates to linear.
		if (is_zero(b)) {
			// No solution - not a function.
			return false;
		}
		root = -c / b;
		return true;
	}

	int sign = lower ? -1 : 1;
	root = (-b + sign * sqrt(discriminant)) / (2 * a);
	return true;
}

double sum(const Vegetation &vegetation, double (*func)(const Individual &indiv)) {
	double result = 0.0;
	for (int i = 0; i < vegetation.nobj; i++) {
		result += func(vegetation[i]);
	}
	return result;
}

double max(const Vegetation &vegetation, double (*func)(const Individual &indiv)) {
	double result = 0;
	for (int i = 0; i < vegetation.nobj; i++) {
		double value = func(vegetation[i]);
		result = std::max(result, value);
	}
	return result;
}
