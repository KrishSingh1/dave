///////////////////////////////////////////////////////////////////////////////////////
/// \file cruinputmodule.h
/// \brief Input module for the CRU TS 3.0 data set
///
/// \author Joe Siltberg
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_GUESSIO_CRU_H
#define LPJ_GUESS_GUESSIO_CRU_H

#include "inputmodule.h"

class CRUInputModule : public InputModule {
public:

	CRUInputModule();

	~CRUInputModule();

	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

	void getlandcover(Gridcell& gridcell);

};

#endif // LPJ_GUESS_GUESSIO_CRU_H
