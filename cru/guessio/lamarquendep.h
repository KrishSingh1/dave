///////////////////////////////////////////////////////////////////////////////////////
/// \file lamarquendep.h
/// \brief Functionality for reading the Lamarque Nitrogen deposition data set
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_LAMARQUENDEP_H
#define LPJ_GUESS_LAMARQUENDEP_H

namespace Lamarque {

/// number of years of historical nitrogen deposition 
const int NYEAR_HISTNDEP = 16;

/// calender year corresponding to first year nitrogen deposition
const int FIRSTHISTYEARNDEP=1850;

/// Retrieves nitrogen deposition for a particular gridcell
/** The values are either taken from a binary archive file or when it's not
 *  provided default to pre-industrial level of 2 kgN/ha/year.
 *
 *  The binary archive files have nitrogen deposition in gN/m2 on a monthly timestep
 *  for 16 years with 10 year interval starting from 1850 (Lamarque et. al., 2011).
 *
 *  \param  file_ndep   Path to binary archive (empty gives pre-industrial values)
 *  \param  lon         Longitude
 *  \param  lat         Latitude
 *  \param  NHxDryDep   Monthly data on daily dry NHx deposition (kgN/m2/day)
 *  \param  NHxWetDep   Monthly data on daily wet NHx deposition (kgN/m2/day)
 *  \param  NOyDryDep   Monthly data on daily dry NOy deposition (kgN/m2/day)
 *  \param  NOyWetDep   Monthly data on daily wet NOy deposition (kgN/m2/day)
 */
void getndep(const char* file_ndep,
             double lon, double lat, 
             double NHxDryDep[NYEAR_HISTNDEP][12],
             double NHxWetDep[NYEAR_HISTNDEP][12],
             double NOyDryDep[NYEAR_HISTNDEP][12],
             double NOyWetDep[NYEAR_HISTNDEP][12]);

}

#endif // LPJ_GUESS_LAMARQUENDEP_H
