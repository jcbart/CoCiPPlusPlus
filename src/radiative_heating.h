#ifndef RADIATIVE_HEATING_H
#define RADIATIVE_HEATING_H

#include <cmath>
#include <algorithm>
#include "constants.h"

namespace radiative_heating {

// Calculate the convective velocity scale, i.e., vertical mixing rate (m s-1)
constexpr double convective_velocity_scale(double depth_eff,
    double eff_heat_rate, double air_temperature) {
    
    return std::pow(
        (constants::GRAVITY * depth_eff*depth_eff * std::max(-eff_heat_rate, 0.)) / air_temperature,
        1./3.
    );
}

// Calculate effective heating rate (K s-1)
constexpr double effective_heating_rate(double d_heat_rate, double cumul_rad_heat,
    double dtheta_dz, double depth) {

    if (cumul_rad_heat <= 0) {
        return 0;
    }
    double heat_denom = std::max(0.5 * dtheta_dz * depth, 0.);
    double heat_ratio = cumul_rad_heat / (cumul_rad_heat + heat_denom);
    return (d_heat_rate * heat_ratio);
}

// Calculate the differential heating rate affecting the contrail plume (K s-1)
double differential_heating_rate(double air_temperature, double rh_i,
    double rho_air, double r_ice_vol, double depth_eff, double tau_contrail, double tau_cirrus,
    double sd0, double sdr, double rsr, double olr);

// Calculate shortwave differential heating rate (K s-1)
double differential_heating_rate_shortwave(double cp_contrail,
    double r_ice_vol_um, double tau_contrail, double tau_cirrus, double sd0, double sdr,
    double rsr);

// Calculate longwave differential heating rate (K s-1)
double differential_heating_rate_longwave(double air_temperature, double rh_i,
    double cp_contrail, double r_ice_vol_um, double tau_contrail, double tau_cirrus, double olr);

// Calculate the heating rate affecting the contrail plume (K s-1)
double heating_rate(double air_temperature, double rh_i, double rho_air,
    double r_ice_vol, double depth_eff, double tau_contrail, double tau_cirrus, double sd0,
    double sdr, double rsr, double olr);

// Calculate shortwave heating rate (K s-1)
double heating_rate_shortwave(double cp_contrail, double r_ice_vol_um,
    double tau_contrail, double tau_cirrus, double sd0, double sdr, double rsr);

// Calculate longwave heating rate (K s-1)
double heating_rate_longwave(double air_temperature, double rh_i,
    double cp_contrail, double r_ice_vol_um, double tau_contrail, double tau_cirrus, double olr);

// Calculate contrail heat capacity per unit length and width (J K-1 m-2)
constexpr double contrail_heat_capacity(double rho_air, double depth_eff) {
    return (depth_eff * rho_air * constants::c_pd);
}

}

#endif