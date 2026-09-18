///////////////////////////////////////////////////////////////////////////////////////
/// \file EquationSolvers.h
/// \brief Auxiliary root finding equation solver class
///
/// \author Phillip Papstefanou
/// $Date: 2020-04-23 15:04:09 +0200 (Wed, 23 Apr 2020) $
///
///////////////////////////////////////////////////////////////////////////////////////
#pragma once
class IEquationSolver {
public:
	IEquationSolver(int max_steps, double eps);
	virtual ~IEquationSolver();

protected:
	virtual double f(const double x) = 0;

	/// Maximum number of iterations
	int max_steps;

	/// Desired precision of root
	double eps;

	/// Running counter
	int c;
};


class SecantSolver : public IEquationSolver {
public:
	SecantSolver(int max_steps, double eps);

	double FindRoot(double xa, double xb);

private:
	double xn;
	double xn1;
	double xn2;

	double fn;
	double fn1;
};


class NewtonSolver : public IEquationSolver {
public:
	NewtonSolver(int max_steps, double eps);

	double FindRoot(double x0);

protected:
	virtual double df(double x) = 0;

private:
	double xn;
	double xn1;

	double fn;
	double dfn;
};


class BisectionSolver : public IEquationSolver {
public:
	BisectionSolver(int max_steps, double eps);

	double FindRoot(double a, double b);
	int steps_used();
	bool solved();
private:
	double xn;
	double xn1;

	double fn;
	double dfn;
};


class AllocationSolver : public BisectionSolver {

public:
	AllocationSolver(int maxsteps, double eps);
	~AllocationSolver();

	double h;
	double k1;
	double ltor_g;
	double k2;
	double b;
	double k3;
	double cmass_leaf_g;
	int solver_state = 0;
	double current_state = -1.0;

	void set_solver_state(int new_state);

	double get_current_state(double x);
	
protected:
	double f(const double x);/*{
		switch (solver_state) {	
		case 1:
			return f1(x) - h;
		case 2:
			return f2(x) - h;
		default:
			return f1(x) - f2(x);
		}
	}*/
	double f1(double cmass_leaf_inc);
	double f2(double cmass_leaf_inc);
	
};

/*
class NormalAllocationSolver : public AllocationSolver{
public:

NormalAllocationSolver(int maxsteps, double eps);
~NormalAllocationSolver();

private:

double f(const double x) override
{
	return f1(x)-f2(x);
}


};


class AbnormalAllocationSolver : public AllocationSolver{
public:

AbnormalAllocationSolver(int maxsteps, double eps);
~AbnormalAllocationSolver();

private:

double f(const double x) override
{
	return f2(x)-h;
}


};
*/

