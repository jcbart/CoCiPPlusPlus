#include <cmath>
#include <string>
#include <CoCiP++/met.h>
#include <CoCiP++/params.h>
#include <CoCiP++/CoCiPTime.h>
#include "wind_shear.h"
#include "tau_cirrus.h"
#include "thermo.h"
#include "geo.h"
#include "constants.h"
#include "CoCiPLog.h"

// IMet

void IMet::calc_variables(double altitude, double cumul_heat, double depth, double cos_a, double sin_a,
    double longitude, double latitude, CoCiPTime& datetime, const Params* params) {
    
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
        effective_vertical_resolution, params->wind_shear_enhancement_exponent);
    ds_dz *= shear_enhancement;
    dsn_dz *= shear_enhancement;

    rho_air = thermo::rho_d(air_temperature, air_pressure);

    rh_i = thermo::rh_i(specific_humidity, air_temperature, air_pressure);

    tau_cirrus = calc_tau_cirrus();

    sdr = geo::solar_direct_radiation(longitude, latitude, datetime, 0.01);

    rsr = std::max(sdr - tnsr, 0.);
}

// ArrayMet

// Types to compile
template struct ArrayMet<float>;
template struct ArrayMet<double>;

template <typename arrayType>
void ArrayMet<arrayType>::get_local_values(double altitude) {
    k = find_k_index(altitude);

    k_lower = std::max(k-1, 0);

    air_pressure = P[k];
    air_temperature = thermo::T_from_T_potential(T_POT[k], air_pressure);
    specific_humidity = QV[k];
    u_wind = U[k];
    v_wind = V[k];
    air_pressure_lower = P[k_lower];
    air_temperature_lower = thermo::T_from_T_potential(T_POT[k_lower], air_pressure_lower);
    u_wind_lower = U[k_lower];
    v_wind_lower = V[k_lower];
    ciwc = CIWC[k];
    effective_vertical_resolution = Z[k] - Z[k_lower];
    dz_m = Z[k] - Z[k_lower];
}

template <typename arrayType>
int ArrayMet<arrayType>::find_k_index(double altitude) const {
    for (int k_trial = 0; k_trial < vsize; k_trial++) {
        if ((altitude >= Z_AT_W[k_trial]) && (altitude < Z_AT_W[k_trial+1])) {
            return k_trial;
        }
    }
    std::string msg = "altitude " + std::to_string(altitude)
        + " m is not in Z_AT_W range (min: " + std::to_string(Z_AT_W[0])
        + ", max: " + std::to_string(Z_AT_W[vsize]) + ")";
    CoCiP_RaiseError(msg, __FILE__, __LINE__);
    return -1;
}

template <typename arrayType>
double ArrayMet<arrayType>::calc_tau_cirrus() const {
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

SimpleMet::SimpleMet(double z0, double P0, double T0, double lapse_rate, double rh_i1,
    double rh_i0, double D1, double DT, double ds_dz_cross_track) {
    
    this->z0 = z0;
    this->P0 = P0;
    this->T0 = T0;
    this->lapse_rate = lapse_rate;
    this->rh_i1 = rh_i1;
    this->rh_i0 = rh_i0;
    this->D1 = D1;
    this->DT = DT;
    this->ds_dz_cross_track = ds_dz_cross_track;
    effective_vertical_resolution = 0;
    dz_m = 100; // Value for calculating _lower variables
}

void SimpleMet::get_local_values(double altitude) {
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
        // Check this
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
    ciwc = 0;
}

double SimpleMet::calc_tau_cirrus() const {
    return 0;
}