#include <cmath>
#include <algorithm>
#include "wake_vortex.h"
#include "wind_shear.h"
#include "thermo.h"
#include "constants.h"

double wake_vortex::max_downward_displacement(double wingspan, double true_airspeed,
    double aircraft_mass, double air_temperature, double dtheta_dz, double ds_dz,
    double air_pressure, double effective_vertical_resolution,
    double wind_shear_enhancement_exponent) {
    
    double rho_air = thermo::rho_d(air_temperature, air_pressure);
    double n_bv = thermo::brunt_vaisala_frequency(air_pressure, air_temperature, dtheta_dz);
    double t_0 = effective_time_scale(wingspan, true_airspeed, aircraft_mass, rho_air);
    bool is_weakly_stratified = (n_bv * t_0 < 0.8);
    double dz_max_strong = downward_displacement_strongly_stratified(wingspan, true_airspeed, aircraft_mass, rho_air, n_bv);
    double dz_max;
    if (is_weakly_stratified) {
        dz_max = downward_displacement_weakly_stratified(wingspan, true_airspeed,
            aircraft_mass, rho_air, n_bv, dz_max_strong, ds_dz, t_0, effective_vertical_resolution,
            wind_shear_enhancement_exponent);
    }
    else {
        dz_max = dz_max_strong;
    }
    return dz_max;    
}

double wake_vortex::downward_displacement_weakly_stratified(double wingspan, double true_airspeed,
    double aircraft_mass, double rho_air, double n_bv, double dz_max_strong, double ds_dz,
    double t_0, double effective_vertical_resolution, double wind_shear_enhancement_exponent) {
    
    double b_0 = wake_vortex_separation(wingspan);
    double dz_max = std::max(dz_max_strong, 10.);
    double shear_enhancement_factor = wind_shear::wind_shear_enhancement_factor(dz_max,
        effective_vertical_resolution, wind_shear_enhancement_exponent);
    double epsn = turbulent_kinetic_energy_dissipation_rate(ds_dz, shear_enhancement_factor);
    double epsn_st = normalized_dissipation_rate(epsn, wingspan, true_airspeed, aircraft_mass, rho_air);
    return (b_0 * (7.68 * (1 - 4.07 * epsn_st + 5.67 * epsn_st*epsn_st) * (0.79 - n_bv * t_0) + 1.88));
}