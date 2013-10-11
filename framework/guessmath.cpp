////////////////////////////////////////////////////////////////////////////////
/// \file guessmath.cpp
/// \brief Implementation of various mathematical functions
///
/// $Date$
///
////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guessmath.h"

double variation_coefficient(double data[], int n) {
	// 0 and 1 will give division with zero.
	if(n>1){
		  double avg,dev = 0, varcoe = 0, sum = 0;
		  int i;
		  double std = 0;

		  for (i=0; i<n; i++)
			  sum += data[i];
		  avg = fabs(sum / n);
		  for (i=0; i<n; i++)
			  dev += (data[i]-avg) * (data[i] - avg);
		  std = sqrt(fabs(dev / (n-1)));

		  if (std > 0 && avg > 0)	// check that data appear in the array
			  varcoe = std / avg;

		  return varcoe;
	}
	else
		return -1.0;
}
