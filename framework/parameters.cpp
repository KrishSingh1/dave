///////////////////////////////////////////////////////////////////////////////////////
/// \file parameters.cpp
/// \brief Implementation of the parameters module
///
/// \author Joe Siltberg
/// $Date: $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "parameters.h"
#include "shell.h"

// Definitions of parameters defined globally in parameters.h,
// for documentation, see parameters.h

vegmodetype vegmode;
int npatch;
double patcharea;
bool ifbgestab;
bool ifsme;
bool ifstochestab;
bool ifstochmort;
bool iffire;
bool ifdisturb;
bool ifcalcsla;
bool ifcalccton;
int estinterval;
double distinterval;
bool ifcdebt;

bool ifcentury;
bool ifnlim;
int freenyears;
double nrelocfrac;
double nfix_a;
double nfix_b;

bool ifsmoothgreffmort;
bool ifdroughtlimitedestab;
bool ifrainonwetdaysonly;

bool ifbvoc;

wateruptaketype wateruptake;

bool run_landcover;
bool run[NLANDCOVERTYPES];
bool lcfrac_fixed;
bool all_fracs_const;
bool ifslowharvestpool;
int nyear_spinup;

xtring state_path;
bool restart;
bool save_state;
int state_year;


///////////////////////////////////////////////////////////////////////////////////////
// Implementation of the Paramlist class

Paramlist param;

void Paramlist::addparam(xtring name,xtring value) {
	Paramtype* p = find(name);
	if (p == 0) {
		p = &createobj();
	}
	p->name=name.lower();
	p->str=value;
}

void Paramlist::addparam(xtring name,double value) {
	Paramtype* p = find(name);
	if (p == 0) {
		p = &createobj();
	}
	p->name=name.lower();
	p->num=value;
}

Paramtype& Paramlist::operator[](xtring name) {
	Paramtype* param = find(name);

	if (param == 0) {
		fail("Paramlist::operator[]: parameter \"%s\" not found",(char*)name);
	}

	return *param;
}

Paramtype* Paramlist::find(xtring name) {
	name = name.lower();
	firstobj();
	while (isobj) {
		Paramtype& p=getobj();
		if (p.name==name) return &p;
		nextobj();
	}
	// nothing found
	return 0;
}
