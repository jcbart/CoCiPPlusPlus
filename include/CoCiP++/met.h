#ifndef MET_H
#define MET_H

#include <cmath>
#include <string>
#include <CoCiP++/CoCiPLog.h>

// Forward declaration
struct Params;
struct CoCiPTime;

// Virtual struct to hold meteorological inputs to CoCiP and variables derived from them
// Is derived in child types depending on met calculation method
struct IMet {
    const double dz_m; // Difference in altitude used to calculate _lower values and by extension dtheta_dz and wind shear (m)

    // Single values from met input
    double tnsr; // TOA net downward shortwave radiation (SDR - RSR) (W m-2)
    double olr; // TOA outgoing (upward) longwave radiation (W m-2)

    // Values taken from get_local_values
    double air_pressure; // Pressure (Pa)
    double air_temperature; // Air temperature (K)
    double specific_humidity; // Specific humidity of water vapor (kg (kg moist air)-1)
    double u_wind; // Eastward wind (m s-1)
    double v_wind; // Northward wind (m s-1)

    double air_pressure_old; // Pressure at last time step (Pa)
    double air_temperature_old; // Temperature inside contrail at last time step (K)
    double specific_humidity_old; // Specific humidity of water vapor at last time step (kg (kg moist air)-1)

    double air_pressure_lower; // Pressure (Pa) at altitude - dz_m (i.e. grid cell below)
    double air_temperature_lower; // Temperature (K) at grid cell below
    double u_wind_lower; // Eastward wind (m s-1) at grid cell below
    double v_wind_lower; // Northward wind (m s-1) at grid cell below

    double effective_vertical_resolution; // Effective vertical resolution of met data (m)

    // Values calculated
    double dtheta_dz = 0; // Potential temperature gradient (K m-1) (dT_dz in pycontrails)
    double ds_dz = 0; // Wind shear (m s-1 m-1)
    double dsn_dz = 0; // Wind shear normal (m s-1 m-1)
    double rho_air = 0; // Dry air mass (kg m-3)
    double rh_i = 0; // Relative humidity w.r.t. ice as a decimal (not percentage)
    double tau_cirrus = 0; // Optical depth of cirrus above the contrail
    double sdr = 0; // Solar direct radiation (W m-2); from met or calculated in calc_shortwave_radiation
    double rsr = 0; // Reflected solar radiation (W m-2); calculated in calc_shortwave_radiation

    IMet(double dz_m) : dz_m(dz_m) {}

    virtual ~IMet() = default;

    void calc_variables(const CoCiPTime& datetime, const double altitude, const double longitude,
        const double latitude, const double cos_a, const double sin_a, const double depth,
        const double cumul_heat, const Params& params);

    // Virtual method called by parent which uses a sedimented altitude to save the local values
    // at that altitude using type-specific method
    virtual void get_sedimented_values(const double altitude_sed, double& air_pressure_sed,
        double& air_temperature_sed, double& specific_humidity_sed) = 0;

private:
    // Virtual method called by calc_variables to update the local variables using type-specific method
    virtual void get_local_values(const double altitude) = 0;
    
    // Virtual method to calculate cirrus optical depth at contrail level ()
    virtual double calc_tau_cirrus(const double altitude) const = 0;
};

// Array met type: takes 1D (vertical) arrays at the contrail's position (lon, lat)
template <typename arrayType>
struct ArrayMet : public IMet {

    // Vertical arrays from met input
    const int vsize = 0; // Unstaggered length of each vertical array
    arrayType* P = nullptr; // Air pressure (Pa)
    arrayType* T_POT = nullptr; // Potential air temperature (K)
    arrayType* QV = nullptr; // Specific humidity of water vapor (kg (kg moist air)-1)
    arrayType* U = nullptr; // Eastward wind (m s-1)
    arrayType* V = nullptr; // Northward wind (m s-1)
    arrayType* CIWC = nullptr; // Cloud ice water mixing ratio (kg (kg dry air)-1)
    arrayType* Z = nullptr; // Altitude at grid cell centres (m)
    arrayType* Z_AT_W = nullptr; // Altitude at grid cell boundaries (m); array len is 1 larger

    // Initialise Met object with unstaggered vertical length
    ArrayMet(double dz_m, int vsize) : IMet(dz_m), vsize(vsize) {}

    // Checks if the array pointers have valid values
    void check_valid_arrays() const;

    // Find the vertical index k of the grid cell in which the contrail is located
    // Will fail if altitude is not in Z_AT_W range
    int find_k_inside(const double altitude) const {
        for (int k_trial = 0; k_trial < vsize; k_trial++) {
            if ((altitude >= Z_AT_W[k_trial]) && (altitude < Z_AT_W[k_trial+1])) {
                return k_trial;
            }
        }
        std::string msg = std::format(
            "altitude {} m is not in Z_AT_W range (min: {}, max: {})",
            altitude, Z_AT_W[0], Z_AT_W[vsize]
        );
        CoCiP_RaiseError(msg, __FILE__, __LINE__);
        return -1;
    }

    // Find the vertical index k of the grid cell centre below the contrail altitude
    // Will fail if altitude is not in Z range
    int find_k_below(const double altitude) const {
        for (int k_trial = 0; k_trial < vsize-1; k_trial++) {
            if ((altitude >= Z[k_trial]) && (altitude < Z[k_trial+1])) {
                return k_trial;
            }
        }
        std::string msg = std::format(
            "altitude {} m is not in Z range (min: {}, max: {})",
            altitude, Z[0], Z[vsize-1]
        );
        CoCiP_RaiseError(msg, __FILE__, __LINE__);
        return -1;
    }

    // Calculate the height fraction of altitude between Z[k_below] and Z[k_below + 1]
    constexpr double calc_interp_fraction(const double altitude, const int k_below) {
        return (altitude - Z[k_below]) / (Z[k_below + 1] - Z[k_below]);
    }

    // Interpolate P array logarithmically given grid cell index below altitude (k_below) and
    // height fraction of altitude between grid cell centres (interp_fraction)
    constexpr double interp_P(const int k_below, const double interp_fraction) {
        return (
            std::pow(P[k_below], 1 - interp_fraction)
            * std::pow(P[k_below + 1], interp_fraction)
        );
    }

    // Interpolate T_POT array linearly given grid cell index below altitude (k_below) and
    // height fraction of altitude between grid cell centres (interp_fraction)
    constexpr double interp_T_POT(const int k_below, const double interp_fraction) {
        return (T_POT[k_below] + interp_fraction * (T_POT[k_below + 1] - T_POT[k_below]));
    }

    // Interpolate QV array linearly given grid cell index below altitude (k_below) and
    // height fraction of altitude between grid cell centres (interp_fraction)
    constexpr double interp_QV(const int k_below, const double interp_fraction) {
        return (QV[k_below] + interp_fraction * (QV[k_below + 1] - QV[k_below]));
    }

    // Interpolate U array linearly given grid cell index below altitude (k_below) and
    // height fraction of altitude between grid cell centres (interp_fraction)
    constexpr double interp_U(const int k_below, const double interp_fraction) {
        return (U[k_below] + interp_fraction * (U[k_below + 1] - U[k_below]));
    }

    // Interpolate V array linearly given grid cell index below altitude (k_below) and
    // height fraction of altitude between grid cell centres (interp_fraction)
    constexpr double interp_V(const int k_below, const double interp_fraction) {
        return (V[k_below] + interp_fraction * (V[k_below + 1] - V[k_below]));
    }

    void get_sedimented_values(const double altitude_sed, double& air_pressure_sed,
        double& air_temperature_sed, double& specific_humidity_sed) override;

private:
    void get_local_values(const double altitude) override;

    // Equivalent to pycontrails' tau_cirrus
    double calc_tau_cirrus(const double altitude) const override;
};

// Array met type: uses simple formulae to calculate local variables
struct SimpleMet : public IMet {

    double z0;
    double P0;
    double T0;
    double lapse_rate;
    double rh_i1;
    double rh_i0;
    double D1;
    double DT;
    double ds_dz_cross_track;

    SimpleMet(double dz_m, double z0, double P0, double T0, double lapse_rate, double rh_i1,
        double rh_i0, double D1, double DT, double ds_dz_cross_track)
        : IMet(dz_m), z0(z0), P0(P0), T0(T0), lapse_rate(lapse_rate), rh_i1(rh_i1), rh_i0(rh_i0),
          D1(D1), DT(DT), ds_dz_cross_track(ds_dz_cross_track) {
    
        effective_vertical_resolution = 0;
    }

    void get_sedimented_values(const double altitude_sed, double& air_pressure_sed,
        double& air_temperature_sed, double& specific_humidity_sed) override;

private:
    void get_local_values(const double altitude) override;

    double calc_tau_cirrus(const double /*altitude (unused)*/) const override {
        return 0;
    }
};

#endif