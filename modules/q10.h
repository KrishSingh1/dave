// Constants required for Q10 lookup tables used by photosynthesis

const double LOOKUPQ10_MINTEMP=-70;
	// minimum temperature ever (deg C)
const double LOOKUPQ10_MAXTEMP=70;
	// maximum temperature ever (deg C)
const double LOOKUPQ10_PRECISION=0.01;
	// rounding precision for temperature in Q10 lookup tables
const int LOOKUPQ10_NDATA=(LOOKUPQ10_MAXTEMP-LOOKUPQ10_MINTEMP+1.0)/
	LOOKUPQ10_PRECISION+0.5;
	// maximum number of values to store in each lookup table
	

// Definition of Q10 lookup table class

class LookupQ10 {

private:
	double* data;

public:
	inline int element(double& temp) {

		// Returns element number corresponding to a particular temperature

		if (temp<LOOKUPQ10_MINTEMP) temp=LOOKUPQ10_MINTEMP;
		else if (temp>LOOKUPQ10_MAXTEMP) temp=LOOKUPQ10_MAXTEMP;

		return (temp-LOOKUPQ10_MINTEMP)/LOOKUPQ10_PRECISION+0.5;
	}

	LookupQ10(double q10,double base25) {
		
		// Constructor (initialises lookup table)
		
		double temp;

		data=new double[LOOKUPQ10_NDATA];
		if (!data) fail("LookupQ10: out of memory creating array");

		for (temp=LOOKUPQ10_MINTEMP;temp<=LOOKUPQ10_MAXTEMP;
			temp+=LOOKUPQ10_PRECISION) {

			data[element(temp)]=base25*pow(q10,(temp-25.0)/10.0);
		}
	}

	double& operator[](double& temp) {
		
		// "Array element" operator (returns temperature-adjusted value
		// based on Q10 and 25-degree base value)
		
		return data[element(temp)];
	}

	// guess2008 - new destructor added
	~LookupQ10() {
		
		delete[] data;
	}
};


// Constants for parameters with Q10 temperature responses used in photosynthesis
// calculations

const double Q10KO=1.2;
	// Q10 for temperature dependency of Michaelis constant for O2 (ko)
const double Q10KC=2.1;
	// Q10 for temperature dependency of Michaelis constant for CO2 (kc)
const double Q10TAU=0.57;
	// Q10 for temperature dependency of CO2/O2 specificity ratio (tau)
const double KO25=3.0E4; // value of ko at 25 deg C (Pa)
const double KC25=30.0; // value of kc at 25 deg C (Pa)
const double TAU25=2600.0; // value of tau at 25 deg C

