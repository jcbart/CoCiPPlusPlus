#ifndef WAKE_VORTEX_H
#define WAKE_VORTEX_H

#include "constants.h"

namespace wake_vortex {

// Calculate the maximum contrail downward displacement after the wake vortex phase (m)
double max_downward_displacement(double wingspan, double true_airspeed,
    double aircraft_mass, double air_temperature, double dtheta_dz, double ds_dz,
    double air_pressure, double effective_vertical_resolution,
    double wind_shear_enhancement_exponent, double turbulent_vertical_velocity_scale);

// Calculate the effective time scale of the wake vortex (s)
constexpr double effective_time_scale(double wingspan, double true_airspeed,
    double aircraft_mass, double rho_air) {
    
    return (
        std::pow(constants::PI, 4) / 32 * std::pow(wingspan, 3) * rho_air * true_airspeed
        / (aircraft_mass * constants::GRAVITY)
    );
}

// Calculate the maximum contrail downward displacement under strongly stratified conditions (m)
constexpr double downward_displacement_strongly_stratified(double wingspan,
    double true_airspeed, double aircraft_mass, double rho_air, double n_bv) {
    
    return (
        ((1.49 * 16) / (2 * constants::PI*constants::PI*constants::PI)
        * aircraft_mass * constants::GRAVITY)
        / (wingspan*wingspan * rho_air * true_airspeed * n_bv)
    );
}

// Calculate the maximum contrail downward displacement under weakly/stably stratified conditions (m)
double downward_displacement_weakly_stratified(double wingspan, double true_airspeed,
    double aircraft_mass, double rho_air, double n_bv, double dz_max_strong, double ds_dz,
    double t_0, double effective_vertical_resolution, double wind_shear_enhancement_exponent,
    double turbulent_vertical_velocity_scale);

// Calculate the wake vortex separation (m)
constexpr double wake_vortex_separation(double wingspan) {
    return (constants::PI * wingspan / 4);
}

// Calculate the turbulent kinetic energy dissipation rate (epsilon; m2 s-3)
// The shear enhancement factor is used to account for any sub-grid scale turbulence
constexpr double turbulent_kinetic_energy_dissipation_rate(double ds_dz,
    double shear_enhancement_factor, double turbulent_vertical_velocity_scale) {
    
    return (0.5 * turbulent_vertical_velocity_scale*turbulent_vertical_velocity_scale
            * (ds_dz * shear_enhancement_factor*shear_enhancement_factor));
}

// Calculate the normalized dissipation rate of the sinking wake vortex (units?)
constexpr double normalized_dissipation_rate(double epsilon, double wingspan,
    double true_airspeed, double aircraft_mass, double rho_air) {

    double c = std::pow(constants::PI/4, 1./3.) * constants::PI*constants::PI*constants::PI / 8;
    double numer = c * std::pow(epsilon * wingspan, 1./3.) * wingspan*wingspan * rho_air * true_airspeed;
    double epsn_st = numer / (constants::GRAVITY * aircraft_mass);
    return std::min(epsn_st, 0.36);
}

// Calculate the initial contrail width (m)
constexpr double initial_contrail_width(double wingspan) {
    return (constants::PI/4 * wingspan);
}

// Calculate the initial contrail depth (m)
constexpr double initial_contrail_depth(double dz_max, double initial_wake_vortex_depth) {
    return (dz_max * initial_wake_vortex_depth);
}

}

#endif