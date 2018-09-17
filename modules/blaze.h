///////////////////////////////////////////////////////////////////////////////////////
/// \file blaze.cpp
/// \brief BLAZE fire simulation and combustion
///
/// \author Lars Nieradzik
/// $Date: 2017-01-24 16:02:51 +0100 (Tue, 24 Jan 2017) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_BLAZE_H
#define LPJ_GUESS_BLAZE_H

#include "guess.h"
#include "driver.h" //	for randfrac
#include "growth.h" //	for allometry
#include "somdynam.h" // for lignin_to_n_ratio & metabolic_litter_fraction
#include "simfire.h" 
#include "gfed31.h"
#include "plib.h"

//WK explain a bit what this is
const double turnoverfract[12][5] = {
	{ .0 , .0 , .05, .2 , .2 }, //   0 Stems       -> ATM
       	{ .0 , .0 , .15, .2 , .2 }, //   1 Branches    -> ATM
        { .03, .13, .25, .5 , .5 }, //   2 Bark        -> ATM
        { .02, .05, .1 , .6 , .6 }, //   3 Leaves      -> ATM
        { .0 , .0 , .05, .2 , .8 }, //   4 Stems       -> Litter (CWD) !corrected*
        { .0 , .02, .07, .2 , .8 }, //   5 Branches    -> Litter (CWD) !corrected*
        { .03, .13, .25, .5 , .5 }, //   6 Bark        -> Litter (str)
        { .05, .1 , .15, .3 , .4 }, //   7 Leaves      -> Litter (str)
        { .0 , .02, .02, .04, .04}, //   8 FDEAD roots -> ATM
        { .5 , .75, .75, .8 , .8 }, //   9 CWD         -> ATM
	{ .6 , .65, .85, 1. , 1. }, //  10 Bark Litter -> ATM
	{ .6 , .65, .85, 1. , 1. }  //  11 Leaf Litter -> ATM*/
};

//WK explain a bit what this is
const double fbranch = 0.05; // 
const double fbark   = 0.01; // 

const double kg2g    = 1000.;
const double min_fuel= 120.; // available fuel threshold[gC/m2]
	
double pixelsize(double latpos,double longsize,double latsize,int postype);

void blaze_accounting_gridcell(Climate& climate);
	
double available_fuel (Patch& patch, int flix);

int get_fli_index(double fli, bool is_sprouter);

//double get_firelineintensity(Patch patch, Climate climate);
void get_firelineintensity(Patch& patch, Climate climate);

bool burntime();

//WK explain a bit what this is
double surv_prob_boreal(double fli) ;
double surv_prob_temp_nl(double dbh, double fli, double mass_cwd) ;
double surv_prob_temp_bl(double dbh, double fli, bool res) ;
double surv_prob_Savanna(double height, double fli) ;
double surv_prob_OzSavanna(double height, double fli) ;
double surv_prob_tropics(double dbh, double fli) ;

double survival_probability(Patch& patch, Individual& indiv, Climate& climate);

void get_combustion_rates(Patch& patch, int flix);

void combust(Patch& patch, Climate& climate);

//		void Individual::blaze_reduce_biomass(double frac_survive); 
//		!CLN defined in guess.h

void blaze_ignition(Climate& climate);

void blaze(Patch& patch, Climate& climate);

#endif // LPJ_GUESS_BLAZE_H


