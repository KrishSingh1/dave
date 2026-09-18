#include "equationsolvers.h"
#include <valarray>
#include <cmath>


IEquationSolver::IEquationSolver(int max_steps, double eps):max_steps(max_steps), eps(eps) {
}

IEquationSolver::~IEquationSolver() {
}

SecantSolver::SecantSolver(int max_steps, double eps) : IEquationSolver(max_steps, eps) {
}

double SecantSolver::FindRoot(double a, double b) {

	c = 0;
	xn = a;
	xn1 = b;


	do
	{
		fn = f(xn);
		fn1 = f(xn1);

		xn2 = xn1 - (xn1 - xn) / (fn1 - fn)*fn1;
		xn = xn1;
		xn1 = xn2;

		c++;
	} while (fabs(fn1 - fn) > eps || c < max_steps);


	return xn2;
}

NewtonSolver::NewtonSolver(int max_steps, double eps) : IEquationSolver(max_steps, eps) {
}

double NewtonSolver::FindRoot(double x0) {
	c = 0;
	xn = x0; 

	// Initialize with a value to get at least one time into the loop
	xn1 = x0 + 2 * eps;

	while (fabs(xn1 - xn) > eps && c < max_steps)
	{
		xn = xn1;
		xn1 = xn - f(xn) / df(xn);
		c++;
	}

	return xn1;
}


BisectionSolver::BisectionSolver(int max_steps, double eps) : IEquationSolver(max_steps, eps) {
}

double BisectionSolver::FindRoot(double a, double b) {
	c = 0;
	double mid = (b - a) / 2.0;
	xn = a;
	xn1 = b;

	while (fabs(xn1 - xn) > eps && c < max_steps)
	{		
		// Find middle point
		mid = (xn + xn1) / 2;

		if (f(mid) == 0.0)
			xn1 = xn;

		else if (f(mid) * f(xn) < 0.0)
			xn1 = mid;
		else
			xn = mid;

		c++;
	}

	return mid;
}

int BisectionSolver::steps_used() {
	return c;
}

bool BisectionSolver::solved() {
	return c < max_steps;
}

AllocationSolver::AllocationSolver(int maxsteps, double eps) :  BisectionSolver(maxsteps, eps) {
}

AllocationSolver::~AllocationSolver() {
}

double AllocationSolver::f(const double x) {
	switch (solver_state) {	
	case 1:
		return f1(x) - h;
	case 2:
		return f2(x) - h;
	default:
		return f1(x) - f2(x);
	}
}

double AllocationSolver::f1(double cmass_leaf_inc) {
	return h + (( k1 - cmass_leaf_inc * (1.0 + 1.0 / ltor_g) )/k2);
}

double AllocationSolver::f2(double cmass_leaf_inc) {
	return k3 * ((b - cmass_leaf_inc * (1.0 + 1.0 / ltor_g)) / (cmass_leaf_g + cmass_leaf_inc));
}

void AllocationSolver::set_solver_state(int new_state) {
	if(new_state > 0 && new_state < 3) {
		solver_state = new_state;
	} else {
		solver_state = 0;
	}
}

double AllocationSolver::get_current_state(double cmass_leaf_inc) {
	return f(cmass_leaf_inc);
}
/*
NormalAllocationSolver::NormalAllocationSolver(int maxsteps, double eps) :  AllocationSolver(maxsteps, eps) {
}

NormalAllocationSolver::~NormalAllocationSolver() {
}

AbnormalAllocationSolver::AbnormalAllocationSolver(int maxsteps, double eps) :  AllocationSolver(maxsteps, eps) {
}

AbnormalAllocationSolver::~AbnormalAllocationSolver() {
}
*/