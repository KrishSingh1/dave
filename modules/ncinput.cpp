#ifdef HAVE_NETCDF

#include "driver.h"
#include "guessnc.h"
#include "pathresolver.h"
#include "util.h"
#include "ncinput.h"
#include "ncinventory.h"
#include "ncutil.h"
#include "ncvariable.h"
#include "variableresolver.h"

ResolveOptions to_resolve_options(const NCInputOptions &options) {
    ResolveOptions resolve_options;
    resolve_options.blaze_enabled = options.blaze_enabled;
    resolve_options.allow_fixed_pressure = options.allow_fixed_pressure;
    resolve_options.allow_netcdf_co2 = options.allow_netcdf_co2;
    resolve_options.allow_netcdf_ndep = options.allow_netcdf_ndep;
    resolve_options.subdaily = options.nsteps_day > 1;
    return resolve_options;
}

ResolvedPlan get_forcing_plan(const std::string &file_path,
                              const NCInputOptions &options) {
    PathResolver path_resolver;
    std::vector<std::string> files = path_resolver.resolve(file_path);
    if (files.empty()) {
        throw std::runtime_error("No files found matching path: " + file_path);
    }

    ResolveOptions resolve_options = to_resolve_options(options);
    NCInventory inventory(files);
    VariableResolver resolver;
    return resolver.resolve(inventory, resolve_options);
}

void open_file(const ResolvedPlan &plan, Variable var,
               CacheStrategy cache_strategy, bool sensibility_checks,
               bounds limits, const std::vector<xtring> &units,
               int &nc_id, std::unique_ptr<NCVariable> &variable) {
    // Get variable info from bindings map.
    auto it = plan.bindings.find(var);
    if (it == plan.bindings.end()) {
        // This should never happen. Such a case should be caught by
        // ResolvedPlan::validate() before we get here.
        throw std::runtime_error("Variable " + std::to_string(static_cast<int>(var)) + " is required by the input plan but is not bound to any variable in the inventory");
    }

    // TODO: this is really hacky.
    // 1. Coordinate variables and dimensions may not have the same name.
    // 2. This discovery is surely not a responsibility of NCInput.
    // 3. We are double-parsing (this parsing has already occurred).
    // Need to refactor this requirement out of the NCVariable contructor.
    const VariableInfo &var_info = it->second;
    nc_id = GuessNC::open_ncdf(var_info.file.c_str());
    variable = std::make_unique<NCVariable>(nc_id, var_info.name.c_str(),
                                            cache_strategy, sensibility_checks,
                                            limits, units, true);
}

NCInput::NCInput() {
    nc_tas = -1;
    nc_rs = -1;
    nc_pr = -1;
    nc_vpd = -1;
    nc_rh = -1;
    nc_sh = -1;
    nc_ps = -1;
    nc_co2 = -1;
    nc_nhx = -1;
    nc_noy = -1;
    nc_tmax = -1;
    nc_tmin = -1;
    nc_ws = -1;
    open_dataset = Dataset::None;
}

void NCInput::init(const NCInputOptions &options) {
    this->options = options;

    std::string file_met_spinup = options.file_met_spinup.empty() ?
        options.file_met_forcing :
        options.file_met_spinup;

    spinup_plan = get_forcing_plan(file_met_spinup, options);
    simulation_plan = get_forcing_plan(options.file_met_forcing, options);

    separate_spinup = file_met_spinup != options.file_met_forcing;
    open_files(spinup_plan, Dataset::Spinup);

    // Ensure that input timestep is an integer multiple of simulation
    // timestep. Note that if using a separate .nc file for spinup and
    // simulation, the timestep may be different between the two files. In
    // this case, we use the input timestep for the spinup, and then use
    // the nsteps_day timestep for the simulation.
    int input_timestep = tas->get_day_length();
    spinup_timestep = input_timestep;
    spinup_calendar = tas->get_calendar();

    simulation_calendar = spinup_calendar;

    spinup_start = first_timestep = tas->get_first_date();
    spinup_end = final_timestep = tas->get_last_date();

    // Ensure that there is at least one year of data for the spinup.
    verify_data_length(first_timestep, final_timestep, input_timestep,
                       spinup_calendar);

    // Record length of input data.
    if (separate_spinup) {
        // Temporarily open a file for the main simulation period.
        int nc_tmp;
        std::unique_ptr<NCVariable> tmp;
        open_file(simulation_plan, Variable::Temp, CacheStrategy::Year, true,
                  { -273.15, 1e4 }, {"degC", "degc", "°C", "℃"}, nc_tmp, tmp);

        // TODO: RAII.
        try {
            // Update simulation timesteps.
            first_timestep = tmp->get_first_date();
            final_timestep = tmp->get_last_date();

            // Set input_timestep to the value from file_met_forcing.
            input_timestep = tmp->get_day_length();
            simulation_calendar = tmp->get_calendar();

            // These three parameters all now come from file_met_forcing; last
            // time we called this function, they came from file_met_spinup.
            // Therefore we should check the main file. Better to fail early.
            verify_data_length(first_timestep, final_timestep, input_timestep, tmp->get_calendar());

            GuessNC::close_ncdf(nc_tmp);
            nc_tmp = -1;
            tmp.reset();
        }
        catch (const std::exception &e) {
            if (nc_tmp != -1) {
                GuessNC::close_ncdf(nc_tmp);
            }
            throw;
        }
    }

    // input_timestep is the timestep from file_met_forcing. Ensure that it
    // may be aggregated to nsteps_day.
    simulation_timestep = input_timestep;
    verify_timestep(input_timestep);

    // Now we must determine the first calendar year of the simulation (not
    // the spinup).
    int first_year = first_timestep.year;

    // If using a separate file for spinup, and this file ends in the year
    // before the "main" file, we can include it as part of the main
    // simulation after running the spinup.
    spinup_in_simulation = false;
    if (separate_spinup && spinup_end.year + 1 == first_timestep.year && is_end_of_year(spinup_end, spinup_timestep, spinup_calendar)) {
        first_year = spinup_start.year;
        spinup_in_simulation = true;
        dprintf("file_met_spinup will be included in main simulation, which will run from %04d-%02d-%02d - %04d-%02d-%02d.\n",
            spinup_start.year, get_month(spinup_start, spinup_calendar) + 1, get_day_of_month(spinup_start, spinup_calendar) + 1,
            final_timestep.year, get_month(final_timestep, simulation_calendar) + 1, get_day_of_month(final_timestep, simulation_calendar) + 1);
    }

    lon = 0;
    lat = 0;

    // TODO: verify consistency of coordinates across all input files.
    // E.g. longitude, latitude and time.
    // At the very least, all files should have all spatial coordinates required
    // by the gridlist, and they should have at least 1 year of overlapping time
    // coverage. We could conditionally be more strict and require identical
    // coordinates in all files.
}

NCInput::~NCInput() {
    close();
}

void NCInput::begin_gridcell(double lon, double lat) {
    this->lon = lon;
    this->lat = lat;

    // Revert to spinup dataset when necessary.
    if (separate_spinup && open_dataset != Dataset::Spinup) {
        open_files(spinup_plan, Dataset::Spinup);
    }
}

void close_if_open(int &nc_id, std::unique_ptr<NCVariable> &variable) {
    if (nc_id != -1) {
        GuessNC::close_ncdf(nc_id);
        nc_id = -1;
    }
    variable.reset();
}

void NCInput::close() {
    close_if_open(nc_tas, tas);
    close_if_open(nc_rs, rs);
    close_if_open(nc_pr, pr);
    close_if_open(nc_ws, ws);
    close_if_open(nc_co2, co2);
    close_if_open(nc_tmax, tmax);
    close_if_open(nc_tmin, tmin);
    close_if_open(nc_vpd, vpd);
    close_if_open(nc_sh, sh);
    close_if_open(nc_rh, rh);
    close_if_open(nc_ps, ps);
    close_if_open(nc_nhx, nhx);
    close_if_open(nc_noy, noy);
    open_dataset = Dataset::None;
}

void NCInput::open_files(const ResolvedPlan &plan, Dataset dataset) {
    // If files aren't already open, this will do nothing.
    close();

    // TODO: automatic unit conversion (where feasible).

    std::vector<xtring> temp_units = {"degC", "degc", "°C", "℃"};
    bounds temp_bounds = {-100, 100};
    open_file(plan, Variable::Temp, options.cache_strategy,
              options.nc_check_inputs, temp_bounds, temp_units, nc_tas, tas);

	// TODO: Should provide units for other insol types as well.
    std::vector<xtring> rs_units = {"W/m2", "W m-2", "w/m^2"};
	bounds bounds_rs = {0, 1e6};
    open_file(plan, Variable::Rs, options.cache_strategy, options.nc_check_inputs,
              bounds_rs, rs_units, nc_rs, rs);
	std::vector<xtring> pr_units = {"mm", "kg m-2", "kg/m2"};
	bounds bounds_pr = {0, 1000};
    open_file(plan, Variable::Pr, options.cache_strategy, options.nc_check_inputs,
              bounds_pr, pr_units, nc_pr, pr);

    if (plan.ws_mode == WsMode::Direct) {
        std::vector<xtring> ws_units = {"m/s", "m s-1", "m s^-1"};
        bounds bounds_ws = {0, 100};
        open_file(plan, Variable::Ws, options.cache_strategy, options.nc_check_inputs,
                  bounds_ws, ws_units, nc_ws, ws);
    }

    if (plan.co2_mode == Co2Mode::Direct) {
        std::vector<xtring> co2_units = {"ppm"};
        bounds bounds_co2 = {0, 900};
        open_file(plan, Variable::CO2, options.cache_strategy, options.nc_check_inputs,
                  bounds_co2, co2_units, nc_co2, co2);
    }

    if (plan.vpd_mode == VpdMode::Direct) {
        std::vector<xtring> vpd_units = {"kPa"};
        // CS updated max from 200 to 100, min from 0 to 0.01 to avoid log(0)
        // issues.
        bounds bounds_vpd = {0.01, 100};
        open_file(plan, Variable::Vpd, options.cache_strategy, options.nc_check_inputs,
                  bounds_vpd, vpd_units, nc_vpd, vpd);
    } else if (plan.vpd_mode == VpdMode::FromSH) {
        std::vector<xtring> sh_units = {"1", "kg/kg"};
        bounds bounds_sh = {0, 1};
        open_file(plan, Variable::SH, options.cache_strategy, options.nc_check_inputs,
                  bounds_sh, sh_units, nc_sh, sh);
    } else if (plan.vpd_mode == VpdMode::FromRH) {
        std::vector<xtring> rh_units = {"1"};
        bounds bounds_rh = {0, 1};
        open_file(plan, Variable::RH, options.cache_strategy, options.nc_check_inputs,
                  bounds_rh, rh_units, nc_rh, rh);
    } else {
        throw std::runtime_error("Invalid VPD mode in input plan: " + std::to_string(static_cast<int>(plan.vpd_mode)));
    }

    if (plan.ndep_mode == NdepMode::Direct) {
        // We will read nitrogen deposition data directly from the input file.
        std::vector<xtring> ndep_units = {"kgN/m2/day"};
        bounds bounds_ndep = {0, 1000};
        open_file(plan, Variable::Nhx, options.cache_strategy, options.nc_check_inputs,
                  bounds_ndep, ndep_units, nc_nhx, nhx);
        open_file(plan, Variable::Noy, options.cache_strategy, options.nc_check_inputs,
                  bounds_ndep, ndep_units, nc_noy, noy);
    }

    if (plan.ps_mode == PsMode::Direct) {
        std::vector<xtring> ps_units = {"Pa"};
        bounds bounds_ps = {50000, 150000};
        open_file(plan, Variable::Ps, options.cache_strategy, options.nc_check_inputs,
                  bounds_ps, ps_units, nc_ps, ps);
    }

    open_dataset = dataset;
}

bool NCInput::read_day(int year, int day, Climate &climate) {
	// First, check if we have enough data to simulate another day.
	if (date.year > nyear_spinup && year > final_timestep.year) {
		return false;
	}

	if (day == 0 && date.get_calendar_year() == final_timestep.year) {
		if (final_timestep.day < Date::MAX_YEAR_LENGTH - 1) {
			// We are in the final year, which is missing some days of data.
			// Because the N deposition disaggregation requires an entire years'
			// worth of precipitation data, we can't simulate part of a year.
            dprintf("Final year (%d) of data is missing some days (contains %d but requires %d for a full year). Because the N deposition disaggregation requires an entire years' worth of precipitation data, we can't simulate this year.\n",
                    final_timestep.year, final_timestep.day + 1,
                    Date::MAX_YEAR_LENGTH);
			return false;
		}

		if (final_timestep.timestep < tas->get_day_length() - 1) {
			// We are in the final year, and the final day is missing some
			// timesteps.
            dprintf("Final year (%d) of data is missing some timesteps (contains %d but requires %d for a full day). Because the N deposition disaggregation requires an entire years' worth of precipitation data, we can't simulate this year.\n",
                    final_timestep.year, final_timestep.timestep + 1,
                    tas->get_day_length());
			return false;
		}
	}

    if (day == 0 && date.year == 0) {
        // TODO: find a way to improve the insoltype situation. This is a mess.
		set_insoltype(climate);
    }

    if (day == 0 && date.get_calendar_year() == first_timestep.year && separate_spinup) {
        // Open the NetCDF files used for the main simulation.
        open_files(simulation_plan, Dataset::Main);

        // Update the cliamte insoltype (as it could have changed).
        set_insoltype(climate);
    }

    // TODO: is this correct? What about spinup_in_simulation?
    ResolvedPlan plan = date.year < nyear_spinup ?
                                spinup_plan :
                                simulation_plan;

	distribute_ndep(plan, year, day);
    set_co2(plan, year, day, climate);

	// Read data from variables (not necessarily from the file - these variable
	// objects will handle caching internally).
	// Note that these arrays will still be on the timestep of the input data,
	// which is not necessarily the same as that of the nsteps_day parameter.
	std::vector<double> temps = tas->read(lon, lat, year, day, 1);
	std::vector<double> precs = pr->read(lon, lat, year, day, 1);
	std::vector<double> insols = rs->read(lon, lat, year, day, 1);
	std::vector<double> pressure = get_ps(plan, year, day);
    std::vector<double> vpds = get_vpd(plan, year, day);

	if ((tas->has_missing_value() && contains(tas->get_missing_value(), temps)) ||
		(pr->has_missing_value() && contains(pr->get_missing_value(), precs)) ||
		(rs->has_missing_value() && contains(rs->get_missing_value(), insols)) ||
		(vpd && vpd->has_missing_value() && contains(vpd->get_missing_value(), vpds)) ||
		(ps && ps->has_missing_value() && contains(ps->get_missing_value(), pressure))) {

		dprintf("WARNING: Encountered missing value in climate data for gridcell (%f, %f) on Year %04d day (of year) %03d. Simulation will now terminate.\n",
                lon, lat, year, day);
		return false;
	}

	// Aggregate subdaily values to daily where necessary, and store them in the
	// climate object. Note that sum inputs are rates and sum are quantities;
	// We can use the mean of the rates and the sum of the quantities. Check the
	// units for more info.
	climate.temp = mean(temps);
	climate.prec = sum(precs);
	climate.insol = mean(insols);

	// Aggregate from one subdaily timescale to another. Note that this will
	// only work if the input timestep is larger than the output timestep, a
	// condition which we have already verified is true in init().
	//
	// nsteps_day is the desired simulation timestep.
	// temp->get_day_length is the input data timestep.

    // Get min/max subdaily temperature before aggregation to target timestep.
	double maxt = max(temps);
	double mint = min(temps);

    // TODO: should be handled by NCVariable.
	if (options.nsteps_day != tas->get_day_length()) {
		temps = aggregate(temps, options.nsteps_day, mean);
		insols = aggregate(insols, options.nsteps_day, mean);
		vpds = aggregate(vpds, options.nsteps_day, mean);
		pressure = aggregate(pressure, options.nsteps_day, mean);
		precs = aggregate(precs, options.nsteps_day, sum);
		// If required, we can aggregate using other functions as well. E.g.
		// precs = aggregate(precs, options.nsteps_day, sum);
	}

	if (plan.ws_mode == WsMode::Direct) {
		std::vector<double> winds = ws->read(lon, lat, year, day, 1);
		climate.u10 = mean(winds);
	}

	for (size_t i = 0; i < vpds.size(); i++) {
		if (vpds[i] < 0.01) {
			if (options.nc_check_inputs) {
				dprintf("Year %d day %d timestep %zu vpd is below 0.01: %f\n", year, day, i, vpds[i]);
			}
			vpds[i] = 0.01;
		}
	}

	climate.temps = temps;
	climate.insols = insols;
	climate.vpds = vpds;
	climate.precs = precs;
	climate.gridcell.dNH4dep = nh4_dep[day];
	climate.gridcell.dNO3dep = no3_dep[day];
	climate.pressure = pressure;

	if (plan.dtr_mode == DtrMode::FromSubdailyTemp) {
		climate.dtr = maxt - mint;
		climate.tmax = maxt;
		climate.tmin = mint;
	} else if (plan.dtr_mode == DtrMode::Direct) {
		// It only really makes sense to provide separate variables for tmin and
		// tmax if input data is daily. Therefore these arrays will usually be
		// of length 1.
		std::vector<double> tmaxs = tmax->read(lon, lat, year, day, 1);
		std::vector<double> tmins = tmin->read(lon, lat, year, day, 1);

		climate.tmax = mean(tmaxs);
		climate.tmin = mean(tmins);
		climate.dtr = climate.tmax - climate.tmin;
	} else if (plan.dtr_mode != DtrMode::NotNeeded) {
        // Should never happen - this is covered by ResolvedPlan::validate().
		throw std::runtime_error("Missing tmax/tmin data, and day length is 1");
	}

    return true;
}

void NCInput::verify_timestep(int timestep) const {
	if (timestep < options.nsteps_day) {
		throw_xtring("Input timestep (%d/day) must be >= simulation timestep (%d/day)"
		, timestep, options.nsteps_day);
	}
	if (timestep % options.nsteps_day != 0) {
		throw_xtring("Input timestep (%d/day) must be an integer multiple of simulation timestep (%d/day)"
		, timestep, options.nsteps_day);
	}
}

int NCInput::get_nyear_data(timepoint start, timepoint end, int day_length,
                            GuessNC::CF::CalendarType calendar) const {
	int nyear = end.year - start.year;

	// If final timestep is the last timestep in a year, then we can add that
	// simulate that entire year (so increment nyear).
	if (is_end_of_year(end, day_length, calendar)) {
		nyear++;
	}

	// If first timestep is not the first timestep in the year, then we don't
	// have a full year of data for the first year. Therefore decrement nyear.
	if (!is_start_of_year(start)) {
		nyear--;
	}

	// Clamp to lower bound of 0.
	return max(nyear, 0);
}

void NCInput::verify_data_length(timepoint start, timepoint end, int timestep,
                                 GuessNC::CF::CalendarType calendar) const {
	// Ensure the simulation file contains at least 1 year of data.
	int nyear_data = get_nyear_data(start, end, timestep, calendar);
	if (nyear_data < 1) {
		throw std::runtime_error("Insufficient met forcing data (<1 year)");
	}
}

void NCInput::set_insoltype(Climate &climate) const {
	// We can get instype from insol's standard_name attribute.
	xtring insol_name = rs->get_standard_name();
	climate.instype = parse_insoltype((const char*)insol_name);
}

void NCInput::set_co2(const ResolvedPlan &plan, int year, int day, Climate &climate) {
    
    if (plan.co2_mode != Co2Mode::Direct) {
        // If CO2 isn't obtained from NC data, it's not out responsibility.
        return;
    }

    if (options.use_mean_annual_co2) {
        // Mean annual co2 - update once per year.
        if (date.day == 0) {
            int year_len = get_year_len(year, co2->get_calendar());
            int co2_year = date.get_calendar_year() < spinup_start.year ?
                                spinup_start.year :
                                date.get_calendar_year();
            std::vector<double> co2_annual = co2->read(lon, lat, co2_year, 0, year_len);
            climate.co2 = mean(co2_annual);
        }
    } else {
        // Mean of subdaily co2.
        climate.co2 = mean(co2->read(lon, lat, year, day, 1));
    }
}

void NCInput::distribute_ndep(const ResolvedPlan &plan, int year, int day) {
    // Nothing to do if NDep isn't coming from NetCDF input.
    if (plan.ndep_mode != NdepMode::Direct) {
        return;
    }

    // TODO: nh4_dep/no3_dep arrays are always size 365, even if they'll never
    // be used. We could reduce memory usage by only allocating these arrays
    // when we'll actually need them.

    if (date.day == 0) {
        // Generate daily N deposition values on first day of year.
        // TODO: leap years.
        int year_len = Date::MAX_YEAR_LENGTH;

        // Read n deposition data from NetCDF input.
        std::vector<double> nhx_subdaily = nhx->read(lon, lat, year, 0, year_len);
        std::vector<double> noy_subdaily = noy->read(lon, lat, year, 0, year_len);

        // Aggregate to daily timestep. Input units are mgN/m2/day, so we take
        // the mean rate over the course of a day, and that will give us the
        // amount of N/m2 deposited in a day. nhx and noy then are vectors with
        // 1 value per day, and units of kgN/m2/day.
        std::vector<double> nhx = aggregate(nhx_subdaily, year_len, mean);
        std::vector<double> noy = aggregate(noy_subdaily, year_len, mean);

        // Copy data into the nh4_dep and no3_dep arrays. nhx = nh4, noy = no3a
        memcpy(nh4_dep, &nhx.front(), year_len * sizeof(double));
        memcpy(no3_dep, &noy.front(), year_len * sizeof(double));
    }
}

std::vector<double> NCInput::get_vpd(const ResolvedPlan &plan, int year, int day) {
    if (plan.vpd_mode == VpdMode::Direct) {
        return vpd->read(lon, lat, year, day, 1);
    }

    // Read subdaily air temperature - required for both SH and RH methods.
    std::vector<double> temp_values = tas->read(lon, lat, year, day, 1);

    // TODO: should look into ways to vectorise this.
    std::vector<double> vpds(temp_values.size());

    if (plan.vpd_mode == VpdMode::FromSH) {
        // Estimate VPD from specific humidity.
        std::vector<double> sh_values = sh->read(lon, lat, year, day, 1);
        std::vector<double> ps_values = get_ps(plan, year, day);
        for (size_t i = 0; i < temp_values.size(); i++) {
            vpds[i] = sh_to_vpd(sh_values[i], temp_values[i], ps_values[i], options.esat_method);
        }
        return vpds;
    }

    if (plan.vpd_mode == VpdMode::FromRH) {
        // Estimate VPD from relative humidity.
        std::vector<double> rh_values = rh->read(lon, lat, year, day, 1);
        for (size_t i = 0; i < temp_values.size(); i++) {
            vpds[i] = rh_to_vpd(rh_values[i], temp_values[i], options.esat_method);
        }
        return vpds;
    }

    throw std::runtime_error("Invalid VPD mode in input plan: " + std::to_string(static_cast<int>(plan.vpd_mode)));
}

std::vector<double> NCInput::get_ps(const ResolvedPlan &plan, int year, int day) {
    if (plan.ps_mode == PsMode::Direct) {
        return ps->read(lon, lat, year, day, 1);
    }

    std::vector<double> pressure(tas->get_day_length());
    for (size_t i = 0; i < pressure.size(); i++) {
        pressure[i] = options.constant_ps;
    }

    return pressure;
}

int NCInput::get_simulation_length() const {
    int nyear = get_nyear_data(first_timestep, final_timestep, simulation_timestep, simulation_calendar);
    if (spinup_in_simulation) {
        nyear += get_nyear_data(spinup_start, spinup_end, spinup_timestep, spinup_calendar);
    }
    return nyear;
}

bool NCInput::has_co2() const {
    // Could check ResolvedPlan, but then we'd need to know year, because it
    // can be different for main simulation vs spinup. This implementation is
    // simpler, although it does assume that input files are already open.
    return co2 != nullptr;
}

int NCInput::get_spinup_length() const {
    return get_nyear_data(spinup_start, spinup_end, spinup_timestep, spinup_calendar);
}

std::vector<double> NCInput::read_annual_precip(int year) {
    // Read subdaily precipitation data for the entire year.
    int year_len = Date::MAX_YEAR_LENGTH;
    return pr->read(lon, lat, year, 0, year_len);
}

timepoint NCInput::get_spinup_start() const {
    return spinup_start;
}

timepoint NCInput::get_spinup_end() const {
    return spinup_end;
}

timepoint NCInput::get_first_timestep() const {
    return first_timestep;
}

int NCInput::get_simulation_start_year() const {
    if (spinup_in_simulation) {
        return spinup_start.year;
    }
    return first_timestep.year;
}

#endif // HAVE_NETCDF
