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
#include "plib.h"

const double turnoverfract[13][5] = {
	{ .0 , .0 , .05, .2 , .2 }, //   0 Stems       -> ATM
       	{ .0 , .0 , .15, .2 , .2 }, //   1 Branches    -> ATM
        { .03, .13, .25, .5 , .5 }, //   2 Bark        -> ATM
        { .02, .05, .1 , .6 , .6 }, //   3 Leaves      -> ATM
        { .0 , .0 , .05, .2 , .8 }, //   4 Stems       -> Litter (DWD) !corrected*
        { .0 , .02, .07, .2 , .8 }, //   5 Branches    -> Litter (CWD) !corrected*
        { .03, .13, .25, .5 , .5 }, //   6 Bark        -> Litter (str)
        { .05, .1 , .15, .3 , .4 }, //   7 Leaves      -> Litter (str)
        { .0 , .02, .02, .04, .04}, //   8 FDEAD roots -> ATM
	{ .5 , .75, .75, .8 , .8 }, //   9 CWD         -> ATM
//  	{ .12, .15, .18, .2 , .2 }, //   9 CWD         -> ATM
	{ .6 , .65, .85, 1. , 1. }, //  10 Bark Litter -> ATM
	{ .6 , .65, .85, 1. , 1. }, //  11 Leaf Litter -> ATM*
	{ .0 , .0 , .1 , .8 , .8 }, //  12 Deadwood    -> ATM
};

/// tuning faktors for litter ready for combustion

// boreal
const double k_tun_bor_lit = 0.1 ;

// temperate region
const double k_tun_tmp_lit = 0.5 ;

// tropics
const double k_tun_trp_lit = 0.75 ;

// fraction of life woody biomass that is branch
const double fbranch   = 0.05;

// fraction of life woody biomass that is bark
const double fbark     = 0.01;

// conversion kg -> g
const double kg2g      = 1000.;

// min. available fuel to start a fire [gC/m2]
const double min_fuel  = 120.; 
	
/// Subroutines will be described in their respective headers in blaze.cpp

double pixelsize(double latpos,double longsize,double latsize,int postype);

void blaze_accounting_gridcell(Climate& climate);
	
double available_fuel (Patch& patch, int fli_index, double k_tun_litter);

int get_fli_index(double fli, bool is_sprouter);

void get_firelineintensity(Patch& patch, Climate climate);

// survival propabilities used  for different biomes 
double surv_prob_boreal(double fli) ;
double surv_prob_temp_nl(double dbh, double fli, double mass_cwd) ;
double surv_prob_temp_bl(double dbh, double fli, bool res) ;
double surv_prob_Savanna(double height, double fli) ;
double surv_prob_Sprouter_Savanna(double height, double fli) ;
double surv_prob_tropics(double dbh, double fli) ;

double survival_probability(Patch& patch, Individual& indiv, Climate& climate);

void get_combustion_rates(Patch& patch, int fli_index, double k_tun_litter);

void blaze(Patch& patch, Climate& climate);

void blaze_burned_area(Climate& climate);

void blaze_driver(Patch& patch, Climate& climate);

#endif // LPJ_GUESS_BLAZE_H


