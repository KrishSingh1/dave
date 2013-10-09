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

/// Nitrogen deposition forcing for a single grid cell
class NDepData {
public:

	/// Default constructor
	/** Before getndep is called the data will all be set to
	 *  pre-industrial level (same as calling getndep with an empty filename).
	 */
	NDepData();

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
	 */
	void getndep(const char* file_ndep,
		double lon, double lat);

	/// Returns nitrogen deposition for one year
	/** Given a calendar year, this function chooses values from the correct 
	 *  10 year interval, and sums the different types of wet and dry 
	 *  deposition.
	 *
	 *  If calendar_year is earlier than the first year in data set (1850),
	 *  the values for the first year will be used. If later than the last
	 *  year (2009), fail() is called and the program terminated.
	 *
	 *  \param calendar_year The year for which to get ndep data
	 *  \param mndrydep      Monthly values for dry nitrogen deposition (kgN/m2/day)
	 *  \param mnwetdep      Monthly values for wet nitrogen deposition (kgN/m2/day)
	 */
	void get_one_calendar_year(int calendar_year,
		double mndrydep[12],
		double mnwetdep[12]);

private:

	/// Fills all arrays with pre-industrial level of 2 kgN/ha/year
	void set_to_pre_industrial();

	/// Monthly data on daily dry NHx deposition (kgN/m2/day)
	double NHxDryDep[NYEAR_HISTNDEP][12];

	/// Monthly data on daily wet NHx deposition (kgN/m2/day)
	double NHxWetDep[NYEAR_HISTNDEP][12];

	/// Monthly data on daily dry NOy deposition (kgN/m2/day)
	double NOyDryDep[NYEAR_HISTNDEP][12];

	/// Monthly data on daily wet NOy deposition (kgN/m2/day)
	double NOyWetDep[NYEAR_HISTNDEP][12];
};


}

#endif // LPJ_GUESS_LAMARQUENDEP_H
