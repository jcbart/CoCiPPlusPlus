#ifndef CONTRAIL_PROPERTIES_H
#define CONTRAIL_PROPERTIES_H

#include <cmath>
#include <algorithm>
#include <optional>
#include "thermo.h"
#include "constants.h"

// Forward declaration
struct Params;

namespace contrail_properties {

// Calculate the specific humidity released by water vapor from aircraft emissions
constexpr double q_exhaust(double air_temperature, double air_pressure,
    double fuel_dist, double width, double depth, double ei_h2o) {
    
    return (
        (ei_h2o * fuel_dist)
        / ((constants::PI / 4) * width * depth * thermo::rho_d(air_temperature, air_pressure))
    );
}

// Estimate the initial contrail ice water content (iwc; kg (kg air)-1) before the wake vortex
// phase
constexpr double initial_iwc(double air_temperature, double specific_humidity,
    double air_pressure, double fuel_dist, double width, double depth, double ei_h2o) {

    double q_sat = thermo::q_sat_ice(air_temperature, air_pressure);
    double q_exh = q_exhaust(air_temperature, air_pressure, fuel_dist, width, depth, ei_h2o);
    return std::max(q_exh + specific_humidity - q_sat, 0.);
}

// Calculate the change in ice water content (kg (kg air)-1) due to adiabatic heating from the wake
// vortex phase
double iwc_adiabatic_heating(double air_temperature_pre_vortex,
    double air_pressure_pre_vortex, double air_pressure_post_vortex);

// Calculate the ambient air temperature (K) after the wake vortex phase
constexpr double temperature_adiabatic_heating(double air_temperature_pre_vortex,
    double air_pressure_pre_vortex, double air_pressure_post_vortex) {
    
    return air_temperature_pre_vortex * std::pow(
        (air_pressure_post_vortex / air_pressure_pre_vortex),
        (constants::GAMMA - 1) / constants::GAMMA
    );
}

// Calculate the ice water content after the wake vortex phase
constexpr double iwc_post_wake_vortex(double iwc, double iwc_ad) {
    return std::max(iwc - iwc_ad, 0.);
}

// Calculate the initial number of ice particles per distance after the wake vortex phase
constexpr double initial_ice_particle_number(double aei, double fuel_dist,
    double min_aei) {
    
    return fuel_dist * std::max(aei, min_aei);
}

// Calculate the activation rate of black carbon particles to contrail ice crystals
// No vPM contribution
constexpr double ice_particle_activation_rate(double air_temperature,
    double T_crit_sac) {
    
    // Ignore rounding line
    return (-0.661 * std::exp(std::min(air_temperature - T_crit_sac, 0.)) + 1);
}

// Determine if contrail is initially "persistent" (not the same as rh_i > 100%)
constexpr bool initial_persistent(double iwc_post_vortex, double rh_i) {
    bool persistent = ((iwc_post_vortex > 1e-12) && (iwc_post_vortex < 1e10)
                       && (rh_i > 0) && (rh_i < 1e10));
    return persistent;
}

// Determine if contrail passes survival threshold criteria
// Does not include age, altitude, or latitude conditions
bool contrail_persistent(double tau_contrail, double n_ice_per_m3,
    const Params& params);

// Calculate effective cross-sectional area of contrail plume from sigma parameters (m2)
// This function calculates the same output as plume_effective_cross_sectional_area, but
// calculated with different input parameters
constexpr double new_effective_area_from_sigma(double sigma_yy,
    double sigma_zz, double sigma_yz) {
    
    return (2 * constants::PI * std::sqrt(sigma_yy * sigma_zz - sigma_yz*sigma_yz));
}

// Calculate the effective cross-sectional area of the contrail plume (m2)
constexpr double plume_effective_cross_sectional_area(double width, double depth,
    double sigma_yz) {

    double sigma_yy = 0.125 * width*width;
    double sigma_zz = 0.125 * depth*depth;
    return new_effective_area_from_sigma(sigma_yy, sigma_zz, sigma_yz);
}

// Calculate the effective depth of the contrail plume (m)
constexpr double plume_effective_depth(double width, double area_eff) {
    return (area_eff / width);
}

// Calculate the contrail plume mass per unit length (kg m-1)
// rho_air should be the total air density, but the input is always dry air density?
constexpr double plume_mass_per_distance(double area_eff, double rho_air) {
    return (area_eff * rho_air);
}

// Calculate the number of contrail ice particles per volume of plume (# m-3)
constexpr double ice_particle_number_per_volume_of_plume(double n_ice_per_m,
    double area_eff) {

    return (n_ice_per_m / area_eff);
}

// Calculate the number of contrail ice particles per mass of air (# kg-1)
constexpr double ice_particle_number_per_mass_of_air(double n_ice_per_vol,
    double rho_air) {
        
    return (n_ice_per_vol / rho_air);
}

// Calculate the ice particle volume mean radius (m)
constexpr double ice_particle_volume_mean_radius(double iwc, double n_ice_per_kg_air) {
    // Force out negative iwc instead of masking
    double total_ice_volume = std::max(0., iwc) / constants::RHO_ICE; // m3 per kg air
    double r_ice_vol = std::pow(
        3. / (4. * constants::PI) * total_ice_volume / n_ice_per_kg_air,
        1./3.
    );
    return std::max(1e-10, r_ice_vol);
}

// Calculate the terminal fall speed of contrail ice particles (m s-1)
double ice_particle_terminal_fall_speed(double air_pressure,
    double air_temperature, double r_ice_vol);

// Calculate the mass of a contrail ice particle (kg)
constexpr double ice_particle_mass(double r_ice_vol) {
    return (4./3. * constants::PI * r_ice_vol*r_ice_vol*r_ice_vol * constants::RHO_ICE);
}

// Calculate contrail horizontal diffusivity (m2 s-1)
constexpr double horizontal_diffusivity(double ds_dz, double depth,
    double max_horizontal_diffusivity) {
    
    return std::min(0.1 * ds_dz * depth*depth, max_horizontal_diffusivity);
}

// Calculate the contrail phase relaxation rate (s-1) - the inverse of the time scale over which
// specific humidity inside a contrail relaxes toward saturation due to sublimation or deposition
constexpr double phase_relaxation_rate(double r_ice_vol, double n_ice_per_vol,
    double diffusivity_water_vapor) {
    
    return 4 * constants::PI * r_ice_vol * n_ice_per_vol * diffusivity_water_vapor;
}

// Calculate contrail vertical diffusivity (m2 s-1)
double vertical_diffusivity(double air_pressure, double air_temperature,
    double dtheta_dz, double depth_eff, double terminal_fall_speed,
    double turbulent_vertical_velocity_scale, double sedimentation_impact_factor,
    std::optional<double> eff_heat_rate, double max_vertical_diffusivity);

// Calculate the rate of contrail ice particle losses due to sedimentation-induced aggregation (# s-1)
constexpr double particle_losses_aggregation(double r_ice_vol,
    double terminal_fall_speed, double area_eff, double agg_efficiency = 1) {

    return ((8 * agg_efficiency * constants::PI * r_ice_vol*r_ice_vol * terminal_fall_speed)
            / area_eff);
}

// Calculate the rate of contrail ice particle losses due to plume-internal turbulence (# s-1)
double particle_losses_turbulence(double width, double depth, double depth_eff,
    double diffuse_h, double diffuse_v, double turb_efficiency = 0.1);

// Calculate the contrail optical depth ()
double contrail_optical_depth(double r_ice_vol, double n_ice_per_m,
    double width);

// Calculate the phase delay of the light wave passing through the contrail ice particle (rad)
constexpr double light_wave_phase_delay(double r_ice_vol) {
    double phase_delay = (4 * constants::PI * (constants::MU_ICE - 1) / constants::LAMBDA_LIGHT)
                         * r_ice_vol;
    return std::min(100., phase_delay);
}

// Calculate the scattering extinction efficiency based on Mie-theory ()
constexpr double scattering_extinction_efficiency(double r_ice_vol) {
    double phase_delay = light_wave_phase_delay(r_ice_vol);
    return (
        2 - (4 / phase_delay)
            * (std::sin(phase_delay) - ((1 - std::cos(phase_delay)) / phase_delay))
    );
}

// Calculate water mass per unit length (kg m-1) taken from the atmosphere due to plume dilution
// In pycontrails, this is included within new_ice_water_content
constexpr double delta_mass_h2o_dilution(double q_old, double q_new, double plume_mass_per_m_old,
    double plume_mass_per_m_new) {

    double q_mean = 0.5 * (q_old + q_new);
    return (plume_mass_per_m_new - plume_mass_per_m_old) * q_mean;
}

// Calculate the new contrail ice water content (kg (kg air)-1) after the time integration step
constexpr double new_ice_water_content(double iwc_old, double q_sat_old, double q_sat_new,
    double plume_mass_per_m_old, double plume_mass_per_m_new, double delta_mass_h2o_dil) {

    // Total mass of H2O (ice + vapor) per m
    double mass_h2o_old = plume_mass_per_m_old * (iwc_old + q_sat_old);
    double mass_h2o_new = mass_h2o_old + delta_mass_h2o_dil;
    // IWC is assumed to be total H2O specific humidity - saturation vapor specific humdity
    double iwc_new = std::max(0., (mass_h2o_new / plume_mass_per_m_new) - q_sat_new);
    return iwc_new;
}

// Calculate water mass per unit length (kg m-1) taken from the atmosphere due to sedimentation
// In pycontrails, this is included within new_ice_water_content_revised
constexpr double delta_mass_h2o_sedimentation_revised(double q_old, double q_sed, double q_sat_old,
    double q_sat_sed, double plume_mass_per_m_old, double plume_mass_per_m_sed,
    double depth_eff, double terminal_fall_speed, double phase_relax_rate, double dt_s) {
    
    double qa = 0.5 * (q_old + q_sed);
    double qs = 0.5 * (q_sat_old + q_sat_sed);
    double m = 0.5 * (plume_mass_per_m_old + plume_mass_per_m_sed);
    double delta_mass_phase_relax = (
        m
        * (qs - qa)
        * std::expm1(-depth_eff * phase_relax_rate / terminal_fall_speed)
        * terminal_fall_speed
        * dt_s
        / depth_eff
    );

    return (plume_mass_per_m_sed * q_sat_sed - plume_mass_per_m_old * q_sat_old)
        + delta_mass_phase_relax;
}

// Calculate water mass per unit length (kg m-1) taken from the atmosphere due to plume dilution
// In pycontrails, this is included within new_ice_water_content_revised
constexpr double delta_mass_h2o_dilution_revised(double q_sed, double q_new, double plume_mass_per_m_sed,
    double plume_mass_per_m_new) {

    double qa = 0.5 * (q_sed + q_new);
    return (plume_mass_per_m_new - plume_mass_per_m_sed) * qa;
}

// Calculate the new contrail ice water content (kg (kg air)-1) after the time integration step
constexpr double new_ice_water_content_revised(double iwc_old, double q_sat_old, double q_sat_new,
    double plume_mass_per_m_old, double plume_mass_per_m_new, double delta_mass_h2o_sed,
    double delta_mass_h2o_dil) {

    double mass_h2o_old = plume_mass_per_m_old * (iwc_old + q_sat_old);
    double mass_h2o_new = mass_h2o_old + delta_mass_h2o_dil + delta_mass_h2o_sed;
    double iwc_new = std::max(0., (mass_h2o_new / plume_mass_per_m_new) - q_sat_new);
    return iwc_new;
}

// Calculate the number of ice particles per distance at the end of the time step (# m-1)
double new_ice_particle_number(double n_ice_per_m, double dn_dt_agg,
    double dn_dt_turb, double length_ratio, double dt_s);

// Calculate contrail altitude after sedimenting for time step dt (s)
constexpr double altitude_after_sedimentation(double altitude,
    double terminal_fall_speed, double dt_s) {
    
    return (altitude - terminal_fall_speed * dt_s);
}

}

#endif