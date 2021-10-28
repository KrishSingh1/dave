#ifndef LPJ_GUESS_LOOKUP_H
#define LPJ_GUESS_LOOKUP_H

#include <math.h>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <iostream>

class Lookup {

private:
	
	double _xmin;
	double _xmax;
	double _y_xmax;
	double _precision;
	double _inverse_precision;
	double _inverse_dx_last;
	std::vector<double> _data;

public:

	Lookup() {

		_xmin = -1;
		_xmax = -1;
		_y_xmax = -1;
		_precision = 1.;
		_inverse_precision = 1.;
		_inverse_dx_last = 0.;
		_data.resize(0);
	}

	void initialize(double xmin, double xmax, double precision, std::vector<double>& params) {

		assert ((xmax-xmin)>0.);
		assert (precision > 0.);
		assert ((xmax - xmin) > precision);
		_xmin = xmin;
		_xmax = xmax;
		_precision = precision;
		_inverse_precision = 1./precision;
		// I do this to avoid roundoff problems
		int nelem = int((xmax-xmin)/precision + 2);
		_data.resize(nelem);

		_y_xmax = fn(xmax, params);

		for (int i=0; i<_data.size(); i++) {
			double x = xmin + i*precision;
			_data[i] = fn(x, params);
			if (fabs(x-xmax + precision) < precision/2.) {
				_data.resize(i+2);
				//_data[i+1] = fn(x+precision, params);
				_data[i+1] = fn(xmax, params);
				break;
			}
		}

		_inverse_dx_last = 1./(xmax - (_data.size()-2)*_precision);

		/*
		std::cout << "Lookup members: " << "\n" <<
			"_xmin = " << _xmin << "; " <<
			"_xmax = " << _xmax << "; " <<
			"_precision = " << _precision << "; " <<
			"_inverse_precision = " << _inverse_precision << "; " <<
			"_inverse_dx_last = " << _inverse_dx_last << "; " <<
			"_y_xmax = " << _y_xmax << "; " <<
			"_data.size() = " << _data.size() << "; " <<
			"\n";
		*/
	}

	virtual double fn(double x, const std::vector<double>& params) = 0;

	double operator() (double x) {

		if (x < _xmin)
			x = _xmin;
		
		if (x > _xmax)
			x = _xmax;
		
		assert (x >= _xmin && x <= _xmax);
		if (x == _xmin) {
			return _data[0];
		} else if (fabs(x - _xmax) < 0.001*_precision) {
			return _y_xmax;
		} else {
			double idx;
			if (x > (_data.size() - 2)*_precision) {
				idx = _inverse_dx_last;
			} else {
				idx = _inverse_precision;
			}
			int i = int((x - _xmin)*_inverse_precision);
			double y = _data[i] + (_data[i+1] - _data[i])*idx*(x - (_xmin + i*_precision) );
			return y;
		}
	}

	~Lookup() {}

};

class NatExp_LUT : public Lookup {

public:

	NatExp_LUT() : Lookup() {}

	NatExp_LUT(double xmin, double xmax, double precision, double a) : Lookup() {

		init(xmin, xmax, precision, a);
	}

	virtual void init(double xmin, double xmax, double precision, double a) {

		std::vector<double> params(1);
		params[0] = a;
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		return exp(params[0] * x);
	}

	~NatExp_LUT() {}

};

class NatLog_LUT : public Lookup {

public:

	NatLog_LUT() : Lookup() {}

	NatLog_LUT(double xmin, double xmax, double precision, double a) : Lookup() {

		init(xmin, xmax, precision, a);
	}

	virtual void init(double xmin, double xmax, double precision, double a) {

		std::vector<double> params(1);
		params[0] = a;
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		return log(params[0] * x);
	}

	~NatLog_LUT() {}

};

class TFac_LUT : public Lookup {

public:

	TFac_LUT() : Lookup() {}

	TFac_LUT(double xmin, double xmax, double a, double b, double precision) : Lookup() {

		init(xmin, xmax, a, b, precision);
	}

	void init(double xmin, double xmax, double a, double b, double precision) {

		std::vector<double> params(2);
		params[0] = a;
		params[1] = b;
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		return exp(params[0] * (x - params[1]));
	}

	~TFac_LUT() {}

};

class TempMod_LUT : public Lookup {

public:

	TempMod_LUT() : Lookup() {}

	TempMod_LUT(double xmin, double xmax, double a, double b, double c, double d, double e, double precision) : Lookup() {

		init(xmin, xmax, a, b, c, d, e, precision);
	}

	void init(double xmin, double xmax, double a, double b, double c, double d, double e, double precision) {

		std::vector<double> params(5);
		params[0] = a;
		params[1] = b;
		params[2] = c;
		params[3] = d;
		params[4] = e;
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		return max(0.0, params[0] + params[1] * pow(x, params[2]) - pow(x / params[3], params[4]));
	}

	~TempMod_LUT() {}

};

class TempFreezeMod_LUT : public Lookup {

public:

	TempFreezeMod_LUT() : Lookup() {}

	TempFreezeMod_LUT(double xmin, double xmax, double a, double b, double c,double precision) : Lookup() {

		init(xmin, xmax, a, b, c, precision);
	}

	void init(double xmin, double xmax, double a, double b, double c, double precision) {

		std::vector<double> params(3);
		params[0] = a;
		params[1] = b;
		params[2] = c;
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		return params[0] * pow(params[1], x / params[2]);
	}

	~TempFreezeMod_LUT() {}

};


class XpowA_LUT : public Lookup {

public:

	XpowA_LUT() : Lookup() {}

	XpowA_LUT(double xmin, double xmax, double precision, double a) : Lookup() {

		init(xmin, xmax, precision, a);
	}

	virtual void init(double xmin, double xmax, double precision, double a) {

		std::vector<double> params(2);
		params[0] = a;
		params[1] = xmax;
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		return std::min(pow(x, params[0]), pow(params[1], params[0]));
	}

	~XpowA_LUT() {}

};

class X01powA_LUT : public XpowA_LUT {

public:

	X01powA_LUT() : XpowA_LUT() {}

	X01powA_LUT(double a, double precision) : XpowA_LUT(0., 1., precision, a) {}

	void init(double a, double precision) {

		std::vector<double> params(2);
		params[0] = a;
		params[1] = 1;
		initialize(0., 1., precision, params);
	}

	~X01powA_LUT() {}
};


class ApowX_LUT : public Lookup {

public:

	ApowX_LUT() : Lookup() {}

	ApowX_LUT(double xmin, double xmax, double precision, double a) : Lookup() {

		init(xmin, xmax, precision, a);
	}

	virtual void init(double xmin, double xmax, double precision, double a) {

		std::vector<double> params(1);
		params[0] = a;
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		return pow(params[0], x);
	}

	~ApowX_LUT() {}

};


class ApowX01_LUT : public ApowX_LUT {

public:

	ApowX01_LUT() : ApowX_LUT() {}

	ApowX01_LUT(double a, double precision) : ApowX_LUT(0., 1., precision, a) {}

	void init(double a, double precision) {

		std::vector<double> params(1);
		params[0] = a;
		initialize(0., 1., precision, params);
	}

	~ApowX01_LUT() {}
};


class AlogX_LUT : public Lookup {

public:

	AlogX_LUT() : Lookup() {}
	
	AlogX_LUT(double xmin, double xmax, double precision, double a) {

		init(xmin, xmax, precision, a);
	}

	virtual void init(double xmin, double xmax, double precision, double a) {

		std::vector<double> params(1);
		params[0] = a;
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		return params[0]*log(x);
	}

	~AlogX_LUT() {}


};


class Kersten_LUT : public Lookup {

public:

	Kersten_LUT() : Lookup() {}
	
	Kersten_LUT(double xmin, double xmax, double precision) {

		init(xmin, xmax, precision);
	}

	void init(double xmin, double xmax, double precision) {

		std::vector<double> params(0);
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		// This is for temp_s >= 0.;
		// with x = s;
		// ke = temp_s >= 0. ? std::max(log(s) + 1.,0.) : s;
		return std::max(log(x) + 1.,0.);
	}

	~Kersten_LUT() {}

};


class Tscal_LUT : public Lookup {

public:
	
	Tscal_LUT() : Lookup() {}

	Tscal_LUT(double temp_min, double temp_low, double temp_high, double temp_max, double precision) : Lookup() {

		init(temp_min, temp_low, temp_high, temp_max, precision);
	}

	void init(double temp_min, double temp_low, double temp_high, double temp_max, double precision) {

		std::vector<double> params(4);
		params[0] = temp_min;
		params[1] = temp_low;
		params[2] = temp_high;
		params[3] = temp_max;
		initialize(temp_min, temp_max, precision, params);
	}

	double fn(double temp, const std::vector<double>& params) {

		double temp_min = params[0];
		double temp_low = params[1];
		double temp_high = params[2];
		double temp_max = params[3];
		double k1 = (temp_min+temp_low) / 2.0;

		return (1. - .01*exp(4.6/(temp_max-temp_high)*(temp-temp_high))) / (1.0+exp((k1-temp)/(k1-temp_min)*4.6));
	}

	~Tscal_LUT();

};


class GreenAmpt_LUT : public Lookup {

public:

	GreenAmpt_LUT() : Lookup() {}
	
	GreenAmpt_LUT(double xmin, double xmax, double precision, double dt_step, double mp_sat, double vol_sat, double k_sat) {

		init(xmin, xmax, precision, dt_step, mp_sat, vol_sat, k_sat);
	}

	void init(double xmin, double xmax, double precision, double dt_step, double mp_sat, double vol_sat, double k_sat) {

		std::vector<double> params(4);
		params[0] = dt_step;
		params[1] = mp_sat;
		params[2] = vol_sat;
		params[3] = k_sat;
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		double EPS = 1.e-4; // Error
		int MAXITER = 6; // Maximum number of iterations

		double dt_step	= params[0];
		double mp_sat = params[1];
		double vol_sat = params[2];
		double k_sat = params[3];
		double vol_wc = x;

		double m = mp_sat*(vol_sat - vol_wc);
		double ks = k_sat*1000.; // mm/s

		double inf = 0.05;
		double dinf = EPS + 1.;
		int niter = 0;
		while (fabs(dinf) > EPS && niter<MAXITER) {
			niter ++;
			double f = inf/ks + m/ks*log(1.-inf/m) - dt_step;
			double df = 1./ks*(1. - 1./(1.-inf/m));
			dinf = -f / df;
			inf += dinf;
		}

		return inf;
	}

	~GreenAmpt_LUT() {}

};


class Test_LUT : public Lookup {

public:


	Test_LUT() : Lookup() {}
	
	Test_LUT(double xmin, double xmax, double precision) {

		init(xmin, xmax, precision);
	}

	void init(double xmin, double xmax, double precision) {

		std::vector<double> params(0);
		initialize(xmin, xmax, precision, params);
	}

	double fn(double x, const std::vector<double>& params) {

		//std::cout << x << "\n";
		return x*x;
	}

	~Test_LUT() {}
};


#endif // LPJ_GUESS_LOOKUP_H
