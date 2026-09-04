#ifndef LPJ_GUESS_NCINPUT_H
#define LPJ_GUESS_NCINPUT_H

#ifdef HAVE_NETCDF

#include <map>
#include <memory>
#include <string>

#include "driver.h"
#include "guess.h"
#include "gutil.h"
// TODO: ncutil.h and ncvariable.h include order is important, but shouldn't be.
#include "ncutil.h"
#include "ncvariable.h"
#include "variableresolver.h"
#include "ncfile.h"

/*
Configuration options which control the behaviour of the NCInput module.
*/
struct NCInputOptions {
    // Path to NetCDF file containing met forcing data.
    std::string file_met_forcing;
    // Path to NetCDF file containing met forcing data to be used for spinup. If
    // not set, the file specified by file_met_forcing will be used for the
    // spinup as well.
    std::string file_met_spinup;
    // Number of timesteps per day to be used in the simulation.
    int nsteps_day;
    // Cache strategy.
    CacheStrategy cache_strategy;
    // Whether to perform strict checks on the input data.
    bool nc_check_inputs;
    // Whether wind speed data is required.
    bool wind_required;
    // Whether daily diurnal temperature range (Tmax/Tmin) data is required.
    bool dtr_required;
    // Whether to allow constant air pressure.
    bool allow_fixed_pressure;
    // Whether to allow NetCDF CO2 data.
    bool allow_netcdf_co2;
    // Whether to allow NetCDF NDEP data.
    bool allow_netcdf_ndep;
    // Constant air pressure to be used (Pa). This is only used if
    // allow_fixed_pressure is true and the input files don't contain surface
    // air pressure.
    double constant_ps;
    // If false, raw co2 values from NetCDF file (if it contains CO2 data) will
    // be used. If true, the mean annual co2 will be used each day instead.
    bool use_mean_annual_co2;
    // The method to use for estimating saturation vapor pressure, which is
    // required for calculating VPD from RH or SH.
    esat_method_t esat_method;
    /// Whether thome should be calculated from the input dataset and stored in
    /// the climate object.
    bool calculate_thome;
};

class NCInput {
public:
    /*
    Constructor. Initializes member variables, but does not open any files or
    read any data.
    */
    NCInput();

    /*
    Destructor. Closes any open file handles.
    */
    ~NCInput();

    /*
    Initialise the NCInput module with the given configuration options.
    */
    void init(const NCInputOptions &options);

    /*
    Load data for the given gridcell. Data for this gridcell will be used for
    all subsequent calls to read_day() until begin_gridcell() is called again
    with a different gridcell.

    @param lon The longitude of the gridcell.
    @param lat The latitude of the gridcell.
    */
    void begin_gridcell(double lon, double lat);

    /*
    Read the climate data for the current day from the NetCDF file, and populate
    the given Climate object with this data.

    @param year    The calendar year of the day to read.
    @param day     The day of the year to read (0-based).
    @param climate The Climate object to populate with the data read from the
                   NetCDF file.
    @return        True if the data was successfully read and the Climate
                   object was populated, false if there was no more data to
                   read (e.g. end of file reached).

    TODO: refactor this so that it returns some kind of daily climate
    structure. Passing the entire climate object is not ideal.
    */
    bool read_day(int year, int day, Climate &climate);

    /*
    Get the complete first year in the spinup period.
    */
    int get_first_spinup_year() const;

    /*
    Get the last complete spinup year used for spinup climate cycling. If a
    separate spinup file overlaps the main forcing, this excludes the overlap.
    */
    int get_last_spinup_year_before_main() const;

    /*
    Get the first complete year in the "main" simulation timeseries.
    */
    int get_first_main_year() const;

    /*
    Get the total number of years in the main simulation.
    */
    int get_simulation_length() const;

    /*
    Check if the input data contains CO2 information.
    */
    bool has_co2() const;

    /*
    Read one year of precipitation data (in mm) for the current gridcell.

    @param year The calendar year to read.
    @return A vector of 365 * timestep subdaily precipitation values,
    */
    std::vector<double> read_annual_precip(int year);

    /*
    Get the start year of the main simulation.
    */
    int get_simulation_start_year() const;
private:
    /*
    Close all open file handles.
    */
    void close();

    /*
    Open the specified file.

    @param plan The input management plan specifying which variable to open and how to obtain it.
    @param var The variable to open.
    @param limits The valid data limits for this variable.
    @param units The valid units for this variable.
    @param variable The NCVariable object to store the opened variable in.
    @return A shared pointer to the opened NCFile.
    */
    void open_file(const ResolvedPlan &plan,
                   Variable var,
                   bounds limits,
                   const std::vector<xtring> &units);

    /*
    Open the specified file for reading.

    @param file The path to the file to open.
    @return A shared pointer to the opened NCFile.
    */
    std::shared_ptr<NCFile> open_file(const std::string &file);

    /*
    Get the NCVariable object for the specified variable. This will throw if the
    currently-open dataset does not contain the specified variable.

    @param var The variable to get.
    @return A reference to the NCVariable object for the specified variable.
    */
    NCVariable &get(Variable var) const;

    /*
    The different datasets that can be used in the simulation.
    */
    enum class Dataset {
        Spinup,
        Main,
        None
    };

    /*
    Time coverage information for a dataset.
    */
    struct TimeCoverage {
        timepoint literal_start;
        timepoint literal_end;
        int first_complete_year;
        int last_complete_year;
        int nyear() const {
            return max(0, last_complete_year - first_complete_year + 1);
        }
        bool contains_complete_year(int year) const {
            return year >= first_complete_year && year <= last_complete_year;
        }
    };

    /*
    Open all files required by the given input plan, store the file handles in
    the files map, and store the variable objects in the corresponding
    NCVariable members.

    @param plan The input management plan specifying which variables are
                required and how to obtain them.
    @param dataset The dataset encapsulated by the plan.
    */
    void open_files(const ResolvedPlan &plan, Dataset dataset);

    /*
    Verify that the given timestep is valid for the simulation configuration.

    @param timestep The actual number of timesteps per day in input data.
    @throws std::runtime_error if the timestep is invalid.
    */
    void verify_timestep(int timestep) const;

	/*
	Verify that the specified time range covers at least one complete year.

	@param start: Start of time range.
	@param end: End of time range.
	@param timestep: Timestep of time range.
	@param calendar: The CF calendar type used in this file.
	*/
    void verify_data_length(timepoint start, timepoint end, int timestep,
                            GuessNC::CF::CalendarType calendar) const;

	/*
	Get the number of complete years of data which exist in the input file.

	@param start: Start of time range.
	@param end: End of time range.
	@param day_length: Number of timesteps per day in the specified time range.
	@param calendar: The calendar type used in the file.
	*/
	int get_nyear_data(timepoint start, timepoint end, int day_length,
                       GuessNC::CF::CalendarType calendar) const;

	/*
	Set the insoltype in the gridcell to that specified by the insolation
	variable in the currently-open NetCDF file.
	*/
    void set_insoltype(Climate &climate) const;

    /*
    Set the CO2 value in the given Climate object according to the specified
    input plan. This may involve reading from the NetCDF file.

    @param plan The input management plan specifying how to obtain the CO2 data.
    @param year The current calendar year.
    @param day The current day of the year (0-based).
    @param climate The Climate object to set the CO2 value in.
    */
    void set_co2(const ResolvedPlan &plan, int year, int day, Climate &climate);

    /*
    Distribute this year's N deposition data to the daily timestep, based on
    daily precipitation. This will populate the nh4_dep and no3_dep arrays with
    the daily N deposition values for the current year in units of kgN/m2/day.

    @param plan The input management plan specifying how to obtain the N deposition data.
    @param year The current calendar year.
    @param day The current day of the year (0-based).
    @param climate The Climate object to set the N deposition value in.
    */
    void distribute_ndep(const ResolvedPlan &plan, int year, int day);

    /*
    Get the subdaily surface air pressure values for the specified day (Pa).

    @param plan The input management plan specifying how to obtain the surface
                air pressure data.
    @param year The current calendar year.
    @param day The current day of the year (0-based).
    @return A vector of surface air pressure values for the specified day, with one value per subdaily timestep.
    */
    std::vector<double> get_ps(const ResolvedPlan &plan, int year, int day);

    /*
    Get the VPD values for the current day according to the specified input plan.
    This may involve reading from the NetCDF file, or calculating from other
    variables.

    @param plan The input management plan specifying how to obtain the VPD data.
    @param year The current calendar year.
    @param day The current day of the year (0-based).
    @return A vector of VPD values for the current day.
    */
    std::vector<double> get_vpd(const ResolvedPlan &plan, int year, int day);

    /*
    Calculate the thome value for the given gridcell from the currently-loaded
    dataset.

    @param lon The longitude of the gridcell.
    @param lat The latitude of the gridcell.
    */
    double calculate_thome(double lon, double lat) const;

    /*
    Configuration options which control the behaviour of this module.
    */
    NCInputOptions options;

    /*
    Input management plan for the spinup. This will be the same as
    simulation_plan if not using a separate dataset for the spinup.
    */
    ResolvedPlan spinup_plan;

    /*
    Input management plan for the main simulation.
    */
    ResolvedPlan simulation_plan;

    /*
    A map of file names to open files.
    */
    std::map<std::string, std::shared_ptr<NCFile>> files;

    /*
    A map of variables to NCVariable objects. This may not contain an entry for
    every variable, but will contain entries for all variables required by the
    input plan corresponding to the currently-open dataset.
    */
    std::map<Variable, std::unique_ptr<NCVariable>> variables;

    /*
    Longitude of the current gridcell. This is set in begin_gridcell() and used
    in read_day() when reading data from the NetCDF file.
    */
    double lon;

    /*
    Latitude of the current gridcell. This is set in begin_gridcell() and used
    in read_day() when reading data from the NetCDF file.
    */
    double lat;

    /*
    The timestep used in the input data for the main simulation.
    */
    int simulation_timestep;

    /*
    The calendar type used in the input data for the main simulation. This is
    used when calculating the number of years of data available in the input
    file, and for determining leap years when updating the date.
    */
    GuessNC::CF::CalendarType simulation_calendar;

	/*
	The timestep used in the spinup data file, if a separate file is being used
	for spinup.
	*/
	int spinup_timestep;

    /*
    Is the spinup being run with a separate dataset to the main simulation?
    */
    bool separate_spinup;

	/*
	Are we using a separate input file for the spinup, and if so, should we use
	the data for this file as part of the "main" simulation after the spinup?
	This can only happen if the last complete year in the spinup dataset is
	either adjacent to or later than, the first complete year in the "main"
	simulation dataset. In that case, the simulation will be run as:

	1. Spinup (recycling first N years of data from spinup file)
	2. Main simulation (using timeseries from spinup file)
	3. Main simulation continued (using entire timeseries from "main" data file)

	This parameter is not a user input, but is automatically enabled if the
	input files are configured appropriately. In that case, a log message will
	be written. This is provided as an optimisation, to avoid having to
	concatenate large datasets.
	*/
	bool spinup_in_simulation;

	/*
	The calendar used in the spinup data file, if a separate file is being used
	for spinup.
	*/
	GuessNC::CF::CalendarType spinup_calendar;

    /*
    Time coverage of spinup dataset, both effective and literal.
    */
    TimeCoverage spinup_coverage;

    /*
    Time coverage of main simulation dataset, both effective and literal. This
    will be the same as spinup_coverage unless using a separate dataset for the
    spinup.
    */
    TimeCoverage simulation_coverage;

	/*
	Daily NH4 deposition data for the current year.
	*/
	double nh4_dep[Date::MAX_YEAR_LENGTH];

	/*
	Daily NO3 deposition data for the current year.
	*/
	double no3_dep[Date::MAX_YEAR_LENGTH];

    /*
    The currently-open dataset.
    */
    Dataset open_dataset;
};

#endif // HAVE_NETCDF

#endif // LPJ_GUESS_NCINPUT_H
