#include <cmath>
#include <algorithm>
#include <CoCiP++/params.h>
#include "contrail_properties.h"
#include "radiative_heating.h"
#include "constants.h"
#include "thermo.h"

double contrail_properties::initial_iwc(double air_temperature, double specific_humidity,
    double air_pressure, double fuel_dist, double width, double depth, double ei_h2o) {
    
    double q_sat = thermo::q_sat_ice(air_temperature, air_pressure);
    double q_exh = q_exhaust(air_temperature, air_pressure, fuel_dist, width, depth, ei_h2o);
    return std::max(q_exh + specific_humidity - q_sat, 0.);
}

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
    const Params* params) {
    
    if (tau_contrail < params->min_tau) { return false; }
    else if (tau_contrail > params->max_tau) { return false; }
    else if (n_ice_per_m3 < params->min_n_ice_per_m3) { return false; }
    else if (n_ice_per_m3 > params->max_n_ice_per_m3) { return false; }
    return true;
}

double contrail_properties::plume_effective_cross_sectional_area(double width, double depth,
    double sigma_yz) {
    
    double sigma_yy = 0.125 * width*width;
    double sigma_zz = 0.125 * depth*depth;
    return new_effective_area_from_sigma(sigma_yy, sigma_zz, sigma_yz);
}

double contrail_properties::ice_particle_volume_mean_radius(double iwc, double n_ice_per_kg_air) {
    // Force out negative iwc instead of masking
    double total_ice_volume = std::max(0., iwc) / constants::RHO_ICE; // m3 per kg air
    double r_ice_vol = std::pow(
        3. / (4. * constants::PI) * total_ice_volume / n_ice_per_kg_air,
        1./3.
    );
    return std::max(1e-10, r_ice_vol);
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
    double sedimentation_impact_factor, double eff_heat_rate, double max_vertical_diffusivity) {

    double n_bv = thermo::brunt_vaisala_frequency(air_pressure, air_temperature, dtheta_dz);
    n_bv = std::max(0.01, n_bv);
    double cvs = radiative_heating::convective_velocity_scale(depth_eff, eff_heat_rate,
        air_temperature);
    cvs = std::max(0.01, cvs);
    
    double d_v = cvs / n_bv + sedimentation_impact_factor * terminal_fall_speed * depth_eff;
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

double contrail_properties::new_ice_water_content(double iwc_old, double q_old, double q_new,
    double q_sat_old, double q_sat_new, double plume_mass_per_m_old, double plume_mass_per_m_new) {
    
    double q_mean = 0.5 * (q_old + q_new);
    // Total mass of H2O (ice + vapor) per m
    double mass_h2o_old = plume_mass_per_m_old * (iwc_old + q_sat_old);
    double mass_h2o_new = mass_h2o_old + (plume_mass_per_m_new - plume_mass_per_m_old) * q_mean;
    // IWC is assumed to be total H2O specific humidity - saturation vapor specific humdity
    double iwc_new = std::max(0., (mass_h2o_new / plume_mass_per_m_new) - q_sat_new);
    return iwc_new;
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