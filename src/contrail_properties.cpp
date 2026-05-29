#include <cmath>
#include <algorithm>
#include <CoCiP++/params.h>
#include "contrail_properties.h"
#include "radiative_heating.h"
#include "constants.h"
#include "thermo.h"

double contrail_properties::iwc_adiabatic_heating(double air_temperature_pre_vortex,
    double air_pressure_pre_vortex, double air_pressure_post_vortex) {
    
    double e_sat_ice_pre_vortex = thermo::e_sat_ice(air_temperature_pre_vortex);
    double air_temperature_post_vortex = temperature_adiabatic_heating(air_temperature_pre_vortex,
        air_pressure_pre_vortex, air_pressure_post_vortex);
    
    double e_sat_ice_post_vortex = thermo::e_sat_ice(air_temperature_post_vortex);

    double delta_q = (constants::EPSILON) * (
        (e_sat_ice_post_vortex / air_pressure_post_vortex) 
        - (e_sat_ice_pre_vortex / air_pressure_pre_vortex)
    );
    return std::max(delta_q, 0.);
}

bool contrail_properties::contrail_persistent(double tau_contrail, double n_ice_per_m3,
    const Params& params) {
    
    if (tau_contrail < params.min_tau) { return false; }
    else if (tau_contrail > params.max_tau) { return false; }
    else if (n_ice_per_m3 < params.min_n_ice_per_m3) { return false; }
    else if (n_ice_per_m3 > params.max_n_ice_per_m3) { return false; }
    return true;
}

double contrail_properties::ice_particle_terminal_fall_speed(double air_pressure,
    double air_temperature, double r_ice_vol) {
    
    double ipm = ice_particle_mass(r_ice_vol);
    double alpha;
    if (ipm >= 4.264e-8) {
        alpha = 8.80 * std::pow(ipm, 0.096);
    }
    else if (ipm >= 2.166e-9) {
        alpha = 329.8 * std::pow(ipm, 0.31);
    }
    else if (ipm >= 2.146e-13) {
        alpha = 63292.4 * std::pow(ipm, 0.57);
    }
    else {
        alpha = 735.4 * std::pow(ipm, 0.42);
    }
    return (alpha * std::pow(3e4 / air_pressure, 0.178) * std::pow(233 / air_temperature, 0.394));
}

double contrail_properties::vertical_diffusivity(double air_pressure, double air_temperature,
    double dtheta_dz, double depth_eff, double terminal_fall_speed,
    double turbulent_vertical_velocity_scale, double sedimentation_impact_factor,
    std::optional<double> eff_heat_rate, double max_vertical_diffusivity) {

    double n_bv = thermo::brunt_vaisala_frequency(air_pressure, air_temperature, dtheta_dz);
    n_bv = std::max(0.001, n_bv);

    double w_prime;
    if (eff_heat_rate.has_value()) {
        w_prime = radiative_heating::convective_velocity_scale(depth_eff, *eff_heat_rate,
            air_temperature);
        w_prime = std::max(turbulent_vertical_velocity_scale, w_prime);
    }
    else {
        w_prime = turbulent_vertical_velocity_scale;
    }
    
    double d_v = w_prime*w_prime / n_bv
        + sedimentation_impact_factor * terminal_fall_speed * depth_eff;
    return std::min(d_v, max_vertical_diffusivity);
}

double contrail_properties::particle_losses_turbulence(double width, double depth, double depth_eff,
    double diffuse_h, double diffuse_v, double turb_efficiency) {
    
    double inner_term = (diffuse_h / std::pow(std::max(width, depth), 2))
                        + (diffuse_v / (depth_eff*depth_eff));
    return (turb_efficiency * std::abs(inner_term));
}

double contrail_properties::contrail_optical_depth(double r_ice_vol, double n_ice_per_m,
    double width) {
    
    double q_ext = scattering_extinction_efficiency(r_ice_vol);
    double tau_contrail = constants::c_r * constants::PI * r_ice_vol*r_ice_vol
                          * (n_ice_per_m / width) * q_ext;
    if (r_ice_vol <= 1e-9) {
        tau_contrail = 0;
    }
    return std::max(0., tau_contrail);
}

double contrail_properties::new_ice_particle_number(double n_ice_per_m, double dn_dt_agg,
    double dn_dt_turb, double length_ratio, double dt_s) {
    
    n_ice_per_m = std::max(n_ice_per_m, 0.);

    double n_ice_per_m_new;

    // Small loss
    if (dn_dt_turb * dt_s < 1e-5) {
        n_ice_per_m_new = n_ice_per_m / (1 + (dn_dt_agg * dt_s * n_ice_per_m));
    }
    else {
        double exp_term = (dn_dt_turb * dt_s < 80) ? std::exp(-dn_dt_turb * dt_s) : 0;
        double numerator = dn_dt_turb * n_ice_per_m * exp_term;
        double denominator = dn_dt_turb + (dn_dt_agg * n_ice_per_m * (1 - exp_term));
        n_ice_per_m_new = (numerator / denominator) * length_ratio;
    }
    return n_ice_per_m_new;
}