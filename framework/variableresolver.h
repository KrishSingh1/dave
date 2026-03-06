#ifndef LPJ_GUESS_VARIABLE_RESOLVER_H
#define LPJ_GUESS_VARIABLE_RESOLVER_H

#ifdef HAVE_NETCDF

#include <map>
#include "ncinventory.h"

/*
Encapsulates the different variables which may be read from climate datasets.
*/
enum class Variable {
    // Longitude
    Lon,
    // Latitude
    Lat,
    // Time
    Time,
    // Air temperature
    Temp,
    // Precipitation
    Pr,
    // Insolation (shortwave radiation)
    Rs,
    // Wind speed
    Ws,
    // Vapor pressure deficit
    Vpd,
    // Relative humidity
    RH,
    // Specific humidity
    SH,
    // Surface pressure
    Ps,
    // Carbon dioxide concentration
    CO2,
    // Nitrogen deposition - NHx
    Nhx,
    // Nitrogen deposition - NOy
    Noy,
    // Maximum temperature
    Tmax,
    // Minimum temperature
    Tmin
};

/*
How should VPD be obtained?
*/
enum class VpdMode {
    // Read directly from the input file.
    Direct,
    // Calculate (estimate) from relative humidity.
    FromRH,
    // Calculate (estimate) from specific humidity + temperature + pressure.
    FromSH
};

/*
How should CO2 be obtained?
*/
enum class Co2Mode {
    // Read directly from the input file.
    Direct,
    // Read from an external source.
    External
};

/*
How should nitrogen deposition be obtained?
*/
enum class NdepMode {
    // Read directly from the input file.
    Direct,
    // Read from an external source.
    External
};

enum class PsMode {
    // Read directly from the input file.
    Direct,
    // Use a fixed/constant pressure value.
    Fixed
};

/*
How should diurnal temperature range be obtained?
*/
enum class DtrMode {
    // Read Tmax and Tmin directly from the input file.
    Direct,
    // Calculate from subdaily temperature data.
    FromSubdailyTemp,
    // Not available but fortunately not needed either.
    NotNeeded
};

/*
How should wind speed be obtained?
*/
enum class WsMode {
    // Read directly from the input file.
    Direct,
    // Not available but fortunately not needed either.
    NotNeeded
};

/*
Options for how to resolve variables from the inventory.
*/
struct ResolveOptions {
    // Whether BLAZE fire module is enabled.
    bool blaze_enabled;
    // Whether to allow using a fixed/constant pressure value, if available. If
    // this is false, air pressure is always a hard requirement. Note that this
    // constant air pressure may later be used to estimate VPD.
    bool allow_fixed_pressure;
    // Whether to allow CO2 data to be read from netcdf inputs. 
    bool allow_netcdf_co2;
    // Whether to allow nitrogen deposition to be read from netcdf inputs.
    bool allow_netcdf_ndep;
    // Whether the input data is subdaily.
    bool subdaily;
};

/*
Encapsulates the resolved variable bindings and modes for how to obtain the
required data from the input files.
*/
struct ResolvedPlan {
    // Variable bindings for each required variable. Every variable required by
    // the below data modes will have an entry. Other variables may not have an
    // entry. This will not contain entries for coordinate variables.
    std::map<Variable, VariableInfo> bindings;

    // How to obtain VPD data.
    VpdMode vpd_mode;

    // How to obtain CO2 data.
    Co2Mode co2_mode;

    // How to obtain nitrogen deposition data.
    NdepMode ndep_mode;

    // How to obtain surface pressure data.
    PsMode ps_mode;

    // How to obtain diurnal temperature range data.
    DtrMode dtr_mode;

    // How to obtain wind speed data.
    WsMode ws_mode;

    // Validate the resolved plan, throwing an exception if any required
    // variables are missing, or if there are any incompatible modes.
    void validate() const;
private:
    void assert_binding_present(Variable var, const std::string &name) const;
};

class VariableResolver {
public:
    VariableResolver() { }

    /*
    Generate a plan for how to resolve the required variables from the
    inventory, based on the available variables and the specified options.
    */
    ResolvedPlan resolve(const NCInventory &inventory,
                         const ResolveOptions &options) const;
};

#endif // HAVE_NETCDF

#endif // LPJ_GUESS_VARIABLE_RESOLVER_H
