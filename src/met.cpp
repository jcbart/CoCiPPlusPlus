#include <cmath>
#include <string>
#include <format>
#include <CoCiP++/met.h>
#include <CoCiP++/params.h>
#include <CoCiP++/CoCiPTime.h>
#include "wind_shear.h"
#include "tau_cirrus.h"
#include "thermo.h"
#include "geo.h"
#include "constants.h"

// IMet

void IMet::calc_variables(const CoCiPTime& datetime, const double altitude, const double longitude,
    const double latitude, const double cos_a, const double sin_a, const double depth,
    const double cumul_heat, const Params& params) {
    
    // Save required values from last time step
    air_pressure_old = air_pressure;
    air_temperature_old = air_temperature;
    specific_humidity_old = specific_humidity;

    // Get met type-derived values
    get_local_values(altitude);

    air_temperature += cumul_heat;

    // Potential temperature gradient
    dtheta_dz = thermo::T_potential_gradient(air_temperature, air_pressure,
        air_temperature_lower, air_pressure_lower, dz_m);
    
    // Wind shear
    ds_dz = wind_shear::wind_shear(u_wind, u_wind_lower, v_wind, v_wind_lower, dz_m);

    // Wind shear normal to contrail heading
    dsn_dz = wind_shear::wind_shear_normal(u_wind, u_wind_lower, v_wind, v_wind_lower, cos_a,
        sin_a, dz_m);
    
    // Wind shear enhancement (in calc_contrail_properties in pycontrails)
    double shear_enhancement = wind_shear::wind_shear_enhancement_factor(depth,
        effective_vertical_resolution, params.wind_shear_enhancement_exponent);
    ds_dz *= shear_enhancement;
    dsn_dz *= shear_enhancement;

    rho_air = thermo::rho_d(air_temperature, air_pressure);

    rh_i = thermo::rh_i(specific_humidity, air_temperature, air_pressure);

    tau_cirrus = calc_tau_cirrus(altitude);

    sdr = geo::solar_direct_radiation(longitude, latitude, datetime, 0.01);

    rsr = std::max(sdr - tnsr, 0.);

    // Calculate values at sedimented altitude

}

// ArrayMet

// Types to compile
template struct ArrayMet<float>;
template struct ArrayMet<double>;

template <typename arrayType>
void ArrayMet<arrayType>::check_valid_arrays() const {
    std::string msg;

    // Check for null pointers
    if (P == nullptr) { CoCiP_RaiseError("P array pointer is null", __FILE__, __LINE__); }
    if (T_POT == nullptr) { CoCiP_RaiseError("T_POT array pointer is null", __FILE__, __LINE__); }
    if (QV == nullptr) { CoCiP_RaiseError("QV array pointer is null", __FILE__, __LINE__); }
    if (U == nullptr) { CoCiP_RaiseError("U array pointer is null", __FILE__, __LINE__); }
    if (V == nullptr) { CoCiP_RaiseError("V array pointer is null", __FILE__, __LINE__); }
    if (CIWC == nullptr) { CoCiP_RaiseError("CIWC array pointer is null", __FILE__, __LINE__); }
    if (Z == nullptr) { CoCiP_RaiseError("Z array pointer is null", __FILE__, __LINE__); }
    if (Z_AT_W == nullptr) { CoCiP_RaiseError("Z_AT_W array pointer is null", __FILE__, __LINE__); }

    // P
    for (int i = 0; i < vsize-1; i++) {
        if (P[i] <= 0) {
            msg = std::format("P array invalid - P[{}] = {}", i, P[i]);
            CoCiP_RaiseError(msg, __FILE__, __LINE__);
        }
        if (P[i] <= P[i+1]) {
            msg = std::format("P array invalid - found P[{}] = {} and P[{}] = {}",
                i, P[i], i+1, P[i+1]);
            CoCiP_RaiseError(msg, __FILE__, __LINE__);
        }
    }
    // T_POT
    for (int i = 0; i < vsize-1; i++) {
        if (T_POT[i] <= 0) {
            msg = std::format("T_POT array invalid - T_POT[{}] = {}", i, T_POT[i]);
            CoCiP_RaiseError(msg, __FILE__, __LINE__);
        }
    }
    // QV
    for (int i = 0; i < vsize-1; i++) {
        if (QV[i] < 0) {
            msg = std::format("QV array invalid - QV[{}] = {}", i, QV[i]);
            CoCiP_RaiseError(msg, __FILE__, __LINE__);
        }
    }
    // CIWC
    for (int i = 0; i < vsize-1; i++) {
        if (CIWC[i] < 0) {
            msg = std::format("CIWC array invalid - CIWC[{}] = {}", i, CIWC[i]);
            CoCiP_RaiseError(msg, __FILE__, __LINE__);
        }
    }
    // Z
    for (int i = 0; i < vsize-1; i++) {
        if (Z[i] >= Z[i+1]) {
            msg = std::format("Z array invalid - found Z[{}] = {} and Z[{}] = {}",
                i, Z[i], i+1, Z[i+1]);
            CoCiP_RaiseError(msg, __FILE__, __LINE__);
        }
    }
    // Z_AT_W
    for (int i = 0; i < vsize; i++) {
        if (Z_AT_W[i] >= Z_AT_W[i+1]) {
            msg = std::format("Z_AT_W array invalid - found Z_AT_W[{}] = {} and Z_AT_W[{}] = {}",
                i, Z_AT_W[i], i+1, Z_AT_W[i+1]);
            CoCiP_RaiseError(msg, __FILE__, __LINE__);
        }
    }
    // Z with Z_AT_W
    for (int i = 0; i < vsize; i++) {
        if (Z[i] < Z_AT_W[i] || Z[i] >= Z_AT_W[i+1]) {
            msg = std::format("Z or Z_AT_W array invalid - found Z[{}] = {}, "
                "but Z_AT_W[{}] = {} and Z_AT_W[{}] = {}",
                i, Z[i], i, Z_AT_W[i], i+1, Z_AT_W[i+1]);
            CoCiP_RaiseError(msg, __FILE__, __LINE__);
        }
    }
}

template <typename arrayType>
void ArrayMet<arrayType>::get_sedimented_values(const double altitude_sed,
    double& air_pressure_sed, double& air_temperature_sed, double& specific_humidity_sed) {
    
    int k_below_sed = find_k_below(altitude_sed);
    double interp_fraction_sed = calc_interp_fraction(altitude_sed, k_below_sed);
    air_pressure_sed = interp_P(k_below_sed, interp_fraction_sed);
    air_temperature_sed = thermo::T_from_T_potential(
        interp_T_POT(k_below_sed, interp_fraction_sed), air_pressure_sed
    );
    specific_humidity_sed = interp_QV(k_below_sed, interp_fraction_sed);
}

template <typename arrayType>
void ArrayMet<arrayType>::get_local_values(const double altitude) {
    
    //check_valid_arrays();

    // Find index of grid cell centre below altitude
    int k_below = find_k_below(altitude);

    // Find height fraction of altitude between grid cell centres
    double interp_fraction = calc_interp_fraction(altitude, k_below);

    // Interpolate values
    air_pressure = interp_P(k_below, interp_fraction);
    air_temperature = thermo::T_from_T_potential(
        interp_T_POT(k_below, interp_fraction), air_pressure
    );
    specific_humidity = interp_QV(k_below, interp_fraction);
    u_wind = interp_U(k_below, interp_fraction);
    v_wind = interp_V(k_below, interp_fraction);

    // Find index of grid cell centre dz_m below altitude
    int k_below_lower = find_k_below(altitude - dz_m);

    // Find height fraction of altitude between grid cell centres
    double interp_fraction_lower = calc_interp_fraction(altitude - dz_m, k_below_lower);

    // Interpolate lower values
    air_pressure_lower = interp_P(k_below_lower, interp_fraction_lower);
    air_temperature_lower = thermo::T_from_T_potential(
        interp_T_POT(k_below_lower, interp_fraction_lower), air_pressure_lower
    );
    u_wind_lower = interp_U(k_below_lower, interp_fraction_lower);
    v_wind_lower = interp_V(k_below_lower, interp_fraction_lower);

    effective_vertical_resolution = Z[k_below + 1] - Z[k_below];
}

template <typename arrayType>
double ArrayMet<arrayType>::calc_tau_cirrus(const double altitude) const {
    // Find grid cell index of contrail
    int k = find_k_inside(altitude);
    // Integrate from contrail altitude to TOA
    double cumsum = 0;
    // Cirrus optical depth is defined as cirrus optical depth above the contrail
    for (int i = k; i < vsize; i++) {
        double beta_e = tau_cirrus::cirrus_effective_extinction_coef(
            CIWC[i],
            thermo::T_from_T_potential(T_POT[i], P[i]),
            P[i]
        );
        double dz = tau_cirrus::height_to_geopt_height(Z_AT_W[i+1])
                    - tau_cirrus::height_to_geopt_height(Z_AT_W[i]);
        cumsum += beta_e * dz;
    }
    return cumsum;
}

// SimpleMet

void SimpleMet::get_sedimented_values(const double altitude_sed, double& air_pressure_sed,
    double& air_temperature_sed, double& specific_humidity_sed) {
    
    air_temperature_sed = T0 - lapse_rate * (altitude_sed - z0);

    air_pressure_sed = P0 * std::pow(
        air_temperature_sed / T0,
        constants::M_d * constants::GRAVITY / (constants::R * lapse_rate)
    );

    double rh_i_temp;
    if (std::abs(altitude_sed - z0) < D1) {
        rh_i_temp = rh_i1;
    }
    else if (std::abs(altitude_sed - z0) < D1 + DT) {
        rh_i_temp = rh_i1 + (std::abs(altitude_sed - z0) - D1) * (rh_i0 - rh_i1) / DT;
    }
    else {
        rh_i_temp = rh_i0;
    }

    specific_humidity_sed = rh_i_temp * thermo::q_sat_ice(air_temperature_sed, air_pressure_sed);
}

void SimpleMet::get_local_values(const double altitude) {
    // Remove cmath and constants and move to thermo

    air_temperature = T0 - lapse_rate * (altitude - z0);

    air_temperature_lower = T0 - lapse_rate * (altitude - dz_m - z0);

    air_pressure = P0 * std::pow(
        air_temperature / T0,
        constants::M_d * constants::GRAVITY / (constants::R * lapse_rate)
    );

    air_pressure_lower = P0 * std::pow(
        air_temperature_lower / T0,
        constants::M_d * constants::GRAVITY / (constants::R * lapse_rate)
    );

    // This is calculated in reverse from specific_humidity immediately afterwards in calc_variables
    // except using air_temperature + cumul_heat
    double rh_i_temp;
    if (std::abs(altitude - z0) < D1) {
        rh_i_temp = rh_i1;
    }
    else if (std::abs(altitude - z0) < D1 + DT) {
        rh_i_temp = rh_i1 + (std::abs(altitude - z0) - D1) * (rh_i0 - rh_i1) / DT;
    }
    else {
        rh_i_temp = rh_i0;
    }

    specific_humidity = rh_i_temp * thermo::q_sat_ice(air_temperature, air_pressure);

    u_wind = 0;
    v_wind = 0;
    u_wind_lower = 0;
    v_wind_lower = -ds_dz_cross_track * dz_m;
}