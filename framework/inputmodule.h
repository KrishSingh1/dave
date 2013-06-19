///////////////////////////////////////////////////////////////////////////////////////
/// \file inputmodule.h
/// \brief Base class for input modules
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_INPUT_MODULE_H
#define LPJ_GUESS_INPUT_MODULE_H

class Gridcell;

class InputModule {
public:
	virtual void init() = 0;

	virtual bool getgridcell(Gridcell& gridcell) = 0;

	virtual bool getclimate(Gridcell& gridcell) = 0;

	virtual void getlandcover(Gridcell& gridcell) = 0;
};

#endif // LPJ_GUESS_INPUT_MODULE_H
