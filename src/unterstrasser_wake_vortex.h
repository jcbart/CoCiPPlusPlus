#ifndef UNTERSTRASSER_WAKE_VORTEX_H
#define UNTERSTRASSER_WAKE_VORTEX_H

#include <cmath>
#include <algorithm>
#include "constants.h"

namespace unterstrasser_wake_vortex {

// Calculate area of the wake-vortex plume (m2)
constexpr double plume_area(double wingspan) {
    double r_plume = 1.5 + 0.314 * wingspan;
    return (2 * constants::PI * r_plume*r_plume);
}

// Calculate the total length-scale effect of the wake vortex downwash (m)
constexpr double z_total_length_scale(double z_atm, double z_emit, double z_desc,
    double true_airspeed, double fuel_flow, double aei_n, double wingspan) {
    
    double fuel_dist = fuel_flow / true_airspeed;
    double n_ice_dist = fuel_dist * aei_n;
    double n_ice_per_vol = n_ice_dist / plume_area(wingspan);
    double n_ice_per_vol_ref = 3.38e12 / plume_area(60.3);
    double psi = std::pow(n_ice_per_vol_ref / n_ice_per_vol, 0.16);
    return (psi * (1.27 * z_atm + 0.42 * z_emit) - 0.49 * z_desc);
}

// Calculate the length-scale effect of ambient supersaturation on the ice crystal mass budget (m)
constexpr double z_atm_length_scale_analytical(double air_temperature, double rh_i) {
    
    // Only perform operation when the ambient condition is supersaturated w.r.t. ice
    // Otherwise, z_atm = 0
    double s_i_clamped = std::max(rh_i - 1, 0.);
    double z_atm = 607.46 * std::pow(s_i_clamped, 0.897) * std::pow(air_temperature/205, 2.225);
    return z_atm;
}

// Calculate aircraft-emitted water vapour concentration in the plume (kg m-3)
constexpr double emitted_water_vapour_concentration(double ei_h2o, double wingspan,
    double true_airspeed, double fuel_flow) {
    
    double h2o_per_dist = (ei_h2o * fuel_flow) / true_airspeed;
    double area_p = plume_area(wingspan);
    return (h2o_per_dist / area_p);
}

// Calculate the length-scale effect of water vapour emissions on the ice crystal mass budget (m)
constexpr double z_emit_length_scale_analytical(double rho_emit, double air_temperature) {
    
    double t_205 = air_temperature - 205;
    double z_emit = 1106.6 * std::pow(rho_emit * 1e5, 0.678 + 0.0116 * t_205)
                   * std::exp(-(0.0807 + 0.000428 * t_205) * t_205);
    return z_emit;
}

// Calculate fraction of ice particle number surviving the wake vortex phase
constexpr double survival_fraction_from_length_scale(double z_total) {
    double f_surv = 0.42 + (1.31 / constants::PI) * std::atan(-1 + z_total/100);
    f_surv = std::max(0., std::min(1., f_surv));
    return f_surv;
}

// Calculate fraction of ice particle number surviving the wake vortex phase and required inputs
// Based on Unterstrasser et al. (2016)
double ice_particle_number_survival_fraction(double air_temperature, double rh_i, double ei_h2o,
    double wingspan, double true_airspeed, double fuel_flow, double aei_n, double z_desc) {
    
    double rho_emit = emitted_water_vapour_concentration(ei_h2o, wingspan, true_airspeed, fuel_flow);

    // Analytical
    double z_atm = z_atm_length_scale_analytical(air_temperature, rh_i);
    double z_emit = z_emit_length_scale_analytical(rho_emit, air_temperature);

    double z_total = z_total_length_scale(z_atm, z_emit, z_desc, true_airspeed, fuel_flow, aei_n,
        wingspan);
    return survival_fraction_from_length_scale(z_total);
}

}

#endif