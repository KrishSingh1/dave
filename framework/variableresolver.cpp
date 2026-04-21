#ifdef HAVE_NETCDF

#include <algorithm>
#include <stdarg.h>
#include <string>
#include <vector>
#include <iterator>
#include <sstream>
#include "shell.h"

#include "variableresolver.h"

// Standard name of the temperature variable, as per CF conventions. Canonical
// units are K.
static const std::string STD_TEMP = "air_temperature";

// Standard names of the insolation variable, as per CF conventions. Canonical
// units are W m-2.
static std::vector<std::string> STDS_INSOL = {
	"surface_downwelling_shortwave_flux_in_air",
	"surface_downwelling_shortwave_flux"
};

// Standard name of the precipitation variable, as per CF conventions. Canonical
// units are kg m-2.
static std::vector<std::string> STDS_PR = {
	"thickness_of_rainfall_amount", // m
	"precipitation_amount", // kg m-2
	// TODO: automatic unit conversions required to handle this case
	"precipitation_flux" // kg m-2 s-1
};

// Optional standard name of the wind speed variable, as per CF conventions.
// Canonical units are m s-1.
static const std::string STD_WIND = "wind_speed";

// Standard name of the CO2 variable, as per CF conventions. Canonical units are
// 1.
static const std::string STD_CO2 = "mole_fraction_of_carbon_dioxide_in_air";

// Standard names of the VPD variable, as per CF conventions. Canonical units
// are Pa.
static std::vector<std::string> STDS_VPD = {
	"water_vapor_saturation_deficit_in_air",
	"water_vapor_saturation_deficit",
};

// Standard name of the relative humidity variable, as per CF conventions.
// Canonical units are Pa.
static const std::string STD_PS = "surface_air_pressure";

// Required standard name of maximum air temperature. This is not part of the
// CF conventions. Expected units are °C.
static const std::string STD_TMAX = "maximum_air_temperature";

// Required standard name of minimum air temperature. This is not part of the
// CF conventions. Expected units are °C.
static const std::string STD_TMIN = "minimum_air_temperature";

// Standard name of the relative humidity variable, as per CF conventions.
// Canonical units are 1.
static const std::string STD_RH = "relative_humidity";

// Standard name of the specific humidity variable, as per CF conventions.
// Canonical units are 1.
static const std::string STD_SH = "specific_humidity";

// Standard name of the NHX nitrogen deposition variable. This is not part of
// the CF conventions. Expected units are kgN m-2 day-1.
static const std::string STD_NDEP_NHX = "ndep_nhx";

// Standard name of the NOY nitrogen deposition variable. This is not part of
// the CF conventions. Expected units are kgN m-2 day-1.
static const std::string STD_NDEP_NOY = "ndep_noy";

enum class LogLevel {
    Debug,
    Diagnostic,
    Info,
    Warning,
    Error
};

std::string pfx(LogLevel level) {
    switch (level) {
        case LogLevel::Debug:
            return "[DBG]";
        case LogLevel::Diagnostic:
            return "[DIA]";
        case LogLevel::Info:
            return "[INF]";
        case LogLevel::Warning:
            return "[WRN]";
        case LogLevel::Error:
            return "[ERR]";
        default:
            throw std::runtime_error("Invalid log level");
    }
}

void vlog(const std::string &message, LogLevel level, va_list args) {
    lprintf("%s ", pfx(level).c_str());
    vlprintf(message.c_str(), args);
    lprintf("\n");
}

// Stub logging function for now. TODO: implement proper logging.
void log(const std::string &message, LogLevel level, ...) {
    va_list args;
    va_start(args, level);
    vlog(message, level, args);
    va_end(args);
}

class PressureWriter {
public:
    explicit PressureWriter(ResolvedPlan &plan) : plan(plan) { }
    void set_direct(const VariableInfo &ps) {
        plan.bindings[Variable::Ps] = ps;
        plan.ps_mode = PsMode::Direct;
    }
    void set_fixed() {
        plan.ps_mode = PsMode::Fixed;
    }
private:
    ResolvedPlan &plan;
};

class WindSpeedWriter {
public:
    explicit WindSpeedWriter(ResolvedPlan &plan) : plan(plan) { }
    void set_direct(const VariableInfo &ws) {
        plan.bindings[Variable::Ws] = ws;
        plan.ws_mode = WsMode::Direct;
    }
    void set_not_needed() {
        plan.ws_mode = WsMode::NotNeeded;
    }
private:
    ResolvedPlan &plan;
};

class CO2Writer {
public:
    explicit CO2Writer(ResolvedPlan &plan) : plan(plan) { }
    void set_direct(const VariableInfo &co2) {
        plan.bindings[Variable::CO2] = co2;
        plan.co2_mode = Co2Mode::Direct;
    }
    void set_external() {
        plan.co2_mode = Co2Mode::External;
    }
private:
    ResolvedPlan &plan;
};

class NDepWriter {
public:
    explicit NDepWriter(ResolvedPlan &plan) : plan(plan) { }
    void set_direct(const VariableInfo &nhx, const VariableInfo &noy) {
        plan.bindings[Variable::Nhx] = nhx;
        plan.bindings[Variable::Noy] = noy;
        plan.ndep_mode = NdepMode::Direct;
    }
    void set_external() {
        plan.ndep_mode = NdepMode::External;
    }
private:
    ResolvedPlan &plan;
};

class DtrWriter {
public:
    explicit DtrWriter(ResolvedPlan &plan) : plan(plan) { }
    void set_direct(const VariableInfo &tmax, const VariableInfo &tmin) {
        plan.bindings[Variable::Tmax] = tmax;
        plan.bindings[Variable::Tmin] = tmin;
        plan.dtr_mode = DtrMode::Direct;
    }
    void set_from_subdaily() {
        plan.dtr_mode = DtrMode::FromSubdailyTemp;
    }
    void set_not_needed() {
        plan.dtr_mode = DtrMode::NotNeeded;
    }
private:
    ResolvedPlan &plan;
};

class VpdWriter {
public:
    explicit VpdWriter(ResolvedPlan &plan) : plan(plan) { }
    void set_direct(const VariableInfo &vpd) {
        plan.bindings[Variable::Vpd] = vpd;
        plan.vpd_mode = VpdMode::Direct;
    }
    void set_from_rh(const VariableInfo &rh) {
        plan.bindings[Variable::RH] = rh;
        plan.vpd_mode = VpdMode::FromRH;
    }
    void set_from_sh(const VariableInfo &sh) {
        plan.bindings[Variable::SH] = sh;
        plan.vpd_mode = VpdMode::FromSH;
    }
private:
    ResolvedPlan &plan;
};

bool try_get_binding(const NCInventory &inventory, const std::vector<std::string> &aliases, VariableInfo &binding) {
    for (const std::string &std_name : aliases) {
        if (inventory.try_get_variable_name(std_name, binding)) {
            return true;
        }
    }
    return false;
}

VariableInfo get_var_info(const NCInventory &inventory, const std::vector<std::string> &aliases) {
    for (const std::string &std_name : aliases) {
        VariableInfo binding;
        if (inventory.try_get_variable_name(std_name, binding)) {
            return binding;
        }
    }
    std::ostringstream msg;
    std::copy(aliases.begin(), aliases.end(), std::ostream_iterator<std::string>(msg, ", "));
    throw std::runtime_error("No variable found with any of the following standard names: " + msg.str());
}

void resolve_pressure(const NCInventory &inventory,
                      const ResolveOptions &options, PressureWriter &writer) {
    VariableInfo ps;
    if (inventory.try_get_variable_name(STD_PS, ps)) {
        log("Found surface pressure variable: %s (file: %s)", LogLevel::Info,
            ps.name.c_str(), ps.file.c_str());
        writer.set_direct(ps);
    } else if (options.allow_fixed_pressure) {
        log("No variable found with standard name '%s', but constant pressure is allowed, and will therefore be used.",
            LogLevel::Info, STD_PS.data());
        writer.set_fixed();
    } else {
        throw std::runtime_error("No variable found with standard name '" + STD_PS + "' and fixed pressure is not allowed");
    }
}

void resolve_wind_speed(const NCInventory &inventory,
                        const ResolveOptions &options, WindSpeedWriter &writer) {
    VariableInfo ws;
    if (inventory.try_get_variable_name(STD_WIND, ws)) {
        log("Found wind speed variable: %s (file: %s)", LogLevel::Info,
            ws.name.c_str(), ws.file.c_str());
        writer.set_direct(ws);
    } else if (options.blaze_enabled) {
        throw std::runtime_error("Blaze fire module is enabled but no variable found with standard name '" + STD_WIND + "'");
    } else {
        log("No variable found with standard name '%s', but blaze fire module is not enabled. Wind speed will not be available, but this is not a problem since the fire module is disabled.", LogLevel::Diagnostic, STD_WIND.data());
        writer.set_not_needed();
    }
}

void resolve_co2(const NCInventory &inventory,
                 const ResolveOptions &options, CO2Writer &writer) {
    VariableInfo co2;
    if (options.allow_netcdf_co2 && inventory.try_get_variable_name(STD_CO2, co2)) {
        log("Found CO2 variable: %s (file: %s)", LogLevel::Info,
            co2.name.c_str(), co2.file.c_str());
        writer.set_direct(co2);
    } else {
        if (options.allow_netcdf_co2) {
            log("No variable found with standard name '%s'. CO2 data will be read from an external source.",
                LogLevel::Warning, STD_CO2.data());
        }
        writer.set_external();
    }
}

void resolve_ndep(const NCInventory &inventory,
                  const ResolveOptions &options, NDepWriter &writer) {
    VariableInfo nhx;
    VariableInfo noy;
    bool has_nhx = options.allow_netcdf_ndep && inventory.try_get_variable_name(STD_NDEP_NHX, nhx);
    bool has_noy = options.allow_netcdf_ndep && inventory.try_get_variable_name(STD_NDEP_NOY, noy);
    if (has_nhx && has_noy) {
        log("Found NHX variable: %s (file: %s)", LogLevel::Info,
            nhx.name.c_str(), nhx.file.c_str());
        log("Found NOY variable: %s (file: %s)", LogLevel::Info,
            noy.name.c_str(), noy.file.c_str());
        writer.set_direct(nhx, noy);
    } else {
        if (options.allow_netcdf_ndep && !has_nhx) {
            log("NHX Variable with standard name '%s' not found. NDep data will therefore need to be obtained from external sources.",
                LogLevel::Warning, STD_NDEP_NHX.data());
        }
        if (options.allow_netcdf_ndep && !has_noy) {
            log("NOY Variable with standard name '%s' not found. NDep data will therefore need to be obtained from external sources.",
                LogLevel::Warning, STD_NDEP_NOY.data());
        }
        writer.set_external();
    }
}

void resolve_dtr(const NCInventory &inventory,
                 const ResolveOptions &options, DtrWriter &writer) {
    if (options.subdaily) {
        log("Input data is subdaily, so Tmax and Tmin will be calculated from subdaily temperature data.", LogLevel::Info);
        writer.set_from_subdaily();
    } else {
        VariableInfo tmax;
        VariableInfo tmin;
        bool has_tmax = inventory.try_get_variable_name(STD_TMAX, tmax);
        bool has_tmin = inventory.try_get_variable_name(STD_TMIN, tmin);
        if (has_tmax && has_tmin) {
            log("Found Tmax variable: %s (file: %s)", LogLevel::Info,
                tmax.name.c_str(), tmax.file.c_str());
            log("Found Tmin variable: %s (file: %s)", LogLevel::Info,
                tmin.name.c_str(), tmin.file.c_str());
            writer.set_direct(tmax, tmin);
        } else if (!options.blaze_enabled) {
            log("Data is daily, and tmin/tmax variables couldn't both be located. However, BLAZE fire module is not enabled, so diurnal temperature range will be treated as optional.", LogLevel::Warning);
            writer.set_not_needed();
        } else {
            if (!has_tmax) {
                log("Variable with standard name '%s' not found. Tmax data will not be available, which is required for the BLAZE fire module.",
                    LogLevel::Error, STD_TMAX.data());
            }
            if (!has_tmin) {
                log("Variable with standard name '%s' not found. Tmin data will not be available, which is required for the BLAZE fire module.",
                    LogLevel::Error, STD_TMIN.data());
            }
            throw std::runtime_error("BLAZE fire module is enabled but Tmax and/or Tmin variables are missing");
        }
    }
}

void resolve_vpd(const NCInventory &inventory,
                 VpdWriter &writer) {
    VariableInfo var;
    if (try_get_binding(inventory, STDS_VPD, var)) {
        // VPD has been provided directly.
        log("Found VPD variable: %s (file: %s)", LogLevel::Info,
            var.name.c_str(), var.file.c_str());
        writer.set_direct(var);
    } else if (inventory.try_get_variable_name(STD_RH, var)) {
        // VPD can be calculated from relative humidity.
        log("No variable found with standard names for VPD, but found variable with standard name '%s'. VPD will be calculated from relative humidity.",
            LogLevel::Info, STD_RH.data());
        writer.set_from_rh(var);
    } else if (inventory.try_get_variable_name(STD_SH, var)) {
        // VPD can be calculated from specific humidity.
        log("No variable found with standard names for VPD, but found variable with standard name '%s', and temperature and pressure data are available. VPD will be calculated from specific humidity, temperature, and pressure.",
            LogLevel::Info, STD_SH.data());
        writer.set_from_sh(var);
    } else {
        throw std::runtime_error("No variables found with standard names for VPD, relative humidity or specific humidity variable");
    }
}

ResolvedPlan VariableResolver::resolve(const NCInventory &inventory,
                                       const ResolveOptions &options) const {
    ResolvedPlan plan;

    plan.bindings[Variable::Pr] = get_var_info(inventory, STDS_PR);
    plan.bindings[Variable::Rs] = get_var_info(inventory, STDS_INSOL);
    plan.bindings[Variable::Temp] = inventory.get_variable_name(STD_TEMP);

    WindSpeedWriter ws_writer(plan);
    resolve_wind_speed(inventory, options, ws_writer);

    // TODO: should we allow the user to *force* the usage of fixed pressure,
    // even if pressure data is available?
    PressureWriter ps_writer(plan);
    resolve_pressure(inventory, options, ps_writer);

    CO2Writer co2_writer(plan);
    resolve_co2(inventory, options, co2_writer);

    NDepWriter ndep_writer(plan);
    resolve_ndep(inventory, options, ndep_writer);

    DtrWriter dtr_writer(plan);
    resolve_dtr(inventory, options, dtr_writer);

    VpdWriter vpd_writer(plan);
    resolve_vpd(inventory, vpd_writer);

    plan.validate();
    return plan;
}

void ResolvedPlan::assert_binding_present(Variable var, const std::string &name) const {
    if (bindings.find(var) == bindings.end()) {
        throw std::runtime_error("Variable '" + name + "' is required but was not found in the inventory");
    }
}

void ResolvedPlan::validate() const {
    assert_binding_present(Variable::Pr, "Precipitation");
    assert_binding_present(Variable::Rs, "SW Radiation");
    assert_binding_present(Variable::Temp, "Temperature");

    // If VPD is to be calculated from specific humidity, temperature and pressure data must be available.
    if (vpd_mode == VpdMode::Direct) {
        assert_binding_present(Variable::Vpd, "VPD");
    } else if (vpd_mode == VpdMode::FromRH) {
        assert_binding_present(Variable::RH, "Relative Humidity");
    } else if (vpd_mode == VpdMode::FromSH) {
        assert_binding_present(Variable::SH, "Specific Humidity");
    }

    if (ps_mode == PsMode::Direct) {
        assert_binding_present(Variable::Ps, "Surface Pressure");
    }
    if (co2_mode == Co2Mode::Direct) {
        assert_binding_present(Variable::CO2, "CO2");
    }
    if (ndep_mode == NdepMode::Direct) {
        assert_binding_present(Variable::Nhx, "NHX nitrogen deposition");
        assert_binding_present(Variable::Noy, "NOY nitrogen deposition");
    }
    if (dtr_mode == DtrMode::Direct) {
        assert_binding_present(Variable::Tmax, "Tmax");
        assert_binding_present(Variable::Tmin, "Tmin");
    }
    if (ws_mode == WsMode::Direct) {
        assert_binding_present(Variable::Ws, "Wind Speed");
    }
}

#endif // HAVE_NETCDF
