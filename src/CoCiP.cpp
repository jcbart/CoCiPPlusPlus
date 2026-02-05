#include <cmath>
#include <algorithm>
#include <vector>
#include <CoCiP++/CoCiP.h>
#include <CoCiP++/met.h>
#include <CoCiP++/params.h>
#include "constants.h"
#include "sac.h"
#include "wake_vortex.h"
#include "unterstrasser_wake_vortex.h"
#include "contrail_properties.h"
#include "radiative_heating.h"
#include "radiative_forcing.h"
#include "thermo.h"
#include "geo.h"
#include "CoCiPLog.h"

void CoCiP::formation() {
    // Calculate local met variables
    update_met_calculations();
    double G = sac::slope_mixing_line(met->specific_humidity, met->air_pressure, engine_efficiency,
        ei_h2o, q_fuel);
    double T_sat_liq = sac::T_sat_liquid(G);
    double rh_crit_sac = sac::rh_critical_sac(met->air_temperature, T_sat_liq, G);
    double rh = thermo::rh_l(met->specific_humidity, met->air_temperature, met->air_pressure);
    if (std::isnan(rh) || std::isnan(rh_crit_sac)) {
        sac = false;
    }
    else {
        sac = (rh > rh_crit_sac);
    }
    T_crit_sac = sac::T_critical_sac(T_sat_liq, rh, G);
}

void CoCiP::simulate_wake_vortex_downwash() {
    update_met_calculations();

    // Initial contrail width, depth, and downward displacement
    double dz_max = wake_vortex::max_downward_displacement(wingspan, true_airspeed, aircraft_mass,
        met->air_temperature, met->dtheta_dz, met->ds_dz, met->air_pressure,
        met->effective_vertical_resolution, params->wind_shear_enhancement_exponent);

    width = wake_vortex::initial_contrail_width(wingspan);
    depth = wake_vortex::initial_contrail_depth(dz_max, params->initial_wake_vortex_depth);
}

void CoCiP::initial_properties() {
    double fuel_dist = fuel_flow / true_airspeed; // (kg m-1)
    altitude -= 0.5 * depth;

    // Save pre-vortex values required below since pre- and post- vortex are both used
    double air_pressure_pre_vortex = met->air_pressure;
    double air_temperature_pre_vortex = met->air_temperature;
    double specific_humidity_pre_vortex = met->specific_humidity;
    // Update local meteorology after change in altitude
    update_met_calculations();

    // Ignoring humidity scaling

    //double q_sat = thermo::q_sat_ice(met->air_temperature, met->air_pressure);
    //double rho_air = thermo::rho_d(met->air_temperature, met->air_pressure);

    // Uses pre-vortex values
    double iwc_pre_vortex = contrail_properties::initial_iwc(air_temperature_pre_vortex,
        specific_humidity_pre_vortex, air_pressure_pre_vortex, fuel_dist, width, depth, ei_h2o);

    double iwc_ad = contrail_properties::iwc_adiabatic_heating(air_temperature_pre_vortex,
        air_pressure_pre_vortex, met->air_pressure);

    iwc = contrail_properties::iwc_post_wake_vortex(iwc_pre_vortex, iwc_ad);

    // Either Extended K15 model
    //double aei = droplet_apparent_emission_index(specific_humidity, air_temperature, T_exhaust,
    //    air_pressure, nvpm_ei_n, G);
    
    // Or...
    // Uses air temperature pre-vortex
    double f_activation = contrail_properties::ice_particle_activation_rate(
        air_temperature_pre_vortex, T_crit_sac);
    double aei = nvpm_ei_n * f_activation;

    // Joint
    double n_ice_per_m_pre_vortex = contrail_properties::initial_ice_particle_number(aei,
        fuel_dist, params->min_ice_particle_number_nvpm_ei_n);

    // Unterstrasser version
    double rh_i_pre_vortex = thermo::rh_i(specific_humidity_pre_vortex, air_temperature_pre_vortex,
        air_pressure_pre_vortex);
    // Why pass nvpm_ei_n instead of aei? This seems like a mistake
    double f_surv = unterstrasser_wake_vortex::ice_particle_number_survival_fraction(
        air_temperature_pre_vortex, rh_i_pre_vortex, ei_h2o, wingspan, true_airspeed, fuel_flow,
        nvpm_ei_n, 0.5 * depth);

    n_ice_per_m = n_ice_per_m_pre_vortex * f_surv;
    
    persistent = contrail_properties::initial_persistent(iwc, met->rh_i);

    // Persistent => is downwash flight
}

void CoCiP::set_heading(const double a) {
    double a_rad = constants::RAD_PER_DEG * a;
    sin_a = std::sin(a_rad);
    cos_a = std::cos(a_rad);
}

void CoCiP::update_met_calculations() {
    met->calc_variables(altitude, cumul_heat, depth, cos_a, sin_a, longitude, latitude, dayOfYear,
        params);
}

void CoCiP::calc_contrail_properties() {
    // Shear enhancement moved to met->calc_variables

    area_eff = contrail_properties::plume_effective_cross_sectional_area(width, depth, sigma_yz);
    double depth_eff = contrail_properties::plume_effective_depth(width, area_eff);

    n_ice_per_vol = contrail_properties::ice_particle_number_per_volume_of_plume(n_ice_per_m,
        area_eff);
    double n_ice_per_kg_air = contrail_properties::ice_particle_number_per_mass_of_air(
        n_ice_per_vol, met->rho_air);
    plume_mass_per_m = contrail_properties::plume_mass_per_distance(area_eff, met->rho_air);
    r_ice_vol = contrail_properties::ice_particle_volume_mean_radius(iwc, n_ice_per_kg_air);
    tau_contrail = contrail_properties::contrail_optical_depth(r_ice_vol, n_ice_per_m, width);
    terminal_fall_speed =  contrail_properties::ice_particle_terminal_fall_speed(
        met->air_pressure, met->air_temperature, r_ice_vol);
    diffuse_h = contrail_properties::horizontal_diffusivity(met->ds_dz, depth,
        params->max_horizontal_diffusivity);

    // Radiative heating
    double theta_rad = geo::orbital_position(dayOfYear);
    double sd0 = geo::solar_constant(theta_rad);
    heat_rate = radiative_heating::heating_rate(met->air_temperature, met->rh_i, met->rho_air,
        r_ice_vol, depth_eff,tau_contrail, met->tau_cirrus, sd0, met->sdr, met->rsr, met->olr);
    d_heat_rate = radiative_heating::differential_heating_rate(met->air_temperature, met->rh_i,
        met->rho_air, r_ice_vol, depth_eff, tau_contrail, met->tau_cirrus, sd0, met->sdr, met->rsr,
        met->olr);
    double eff_heat_rate = radiative_heating::effective_heating_rate(d_heat_rate,
        cumul_differential_heat, met->dtheta_dz, depth);

    diffuse_v = contrail_properties::vertical_diffusivity(met->air_pressure, met->air_temperature,
        met->dtheta_dz, depth_eff, terminal_fall_speed, params->sedimentation_impact_factor,
        eff_heat_rate, params->max_vertical_diffusivity);

    dn_dt_agg = contrail_properties::particle_losses_aggregation(r_ice_vol, terminal_fall_speed,
        area_eff);
    dn_dt_turb = contrail_properties::particle_losses_turbulence(width, depth, depth_eff,
        diffuse_h, diffuse_v);
}

void CoCiP::calc_radiative_properties() {
    double theta_rad = geo::orbital_position(dayOfYear);
    double sd0 = geo::solar_constant(theta_rad);

    double r_vol_um = r_ice_vol * 1e6;
    std::vector<double> habit_weights = radiative_forcing::calc_habit_weights(r_vol_um,
        params->habit_distributions, params->radius_threshold_um);
    
    rf_lw = radiative_forcing::longwave_radiative_forcing(r_vol_um, met->olr, met->air_temperature,
        tau_contrail, met->tau_cirrus, habit_weights);
    rf_sw = radiative_forcing::shortwave_radiative_forcing(r_vol_um, met->sdr, met->rsr, sd0,
        tau_contrail, met->tau_cirrus, habit_weights);

    double rf_lw_scaled = rf_lw * params->rf_lw_enhancement_factor;
    double rf_sw_scaled = rf_sw * params->rf_sw_enhancement_factor;
    rf_net = radiative_forcing::net_radiative_forcing(rf_lw_scaled, rf_sw_scaled);
}

void CoCiP::process_downwash_flight() {
    update_met_calculations();
    calc_contrail_properties();
    calc_radiative_properties();
}

void CoCiP::plume_temporal_evolution(const double length_ratio, const double dt_s) {
    double sigma_yy = 0.125 * width*width;
    double sigma_zz = 0.125 * depth*depth;
    double max_sigma_zz = 0.125 * params->max_depth*params->max_depth;
    double max_diffuse_v = (max_sigma_zz - sigma_zz) / (2 * dt_s);
    diffuse_v = std::min(diffuse_v, max_diffuse_v);

    // To avoid redundancy
    double dsn_dz_2 = met->dsn_dz*met->dsn_dz;
    double dt_s_2 = dt_s*dt_s;
    double dt_s_3 = dt_s_2 * dt_s;

    double sigma_yy_new = (
        ((2. / 3.) * dsn_dz_2 * diffuse_v * dt_s_3)
        + (dsn_dz_2 * sigma_zz * dt_s_2)
        + (2 * (diffuse_h + met->dsn_dz * sigma_yz) * dt_s)
        + sigma_yy
    ) * (length_ratio*length_ratio);

    double sigma_zz_new = (2 * diffuse_v * dt_s) + sigma_zz;

    double sigma_yz_new = (
        (met->dsn_dz * diffuse_v * dt_s_2) + (met->dsn_dz * sigma_zz * dt_s) + sigma_yz
    ) * length_ratio;

    sigma_yy = sigma_yy_new;
    sigma_zz = sigma_zz_new;
    sigma_yz = sigma_yz_new;

    width = std::sqrt(8 * sigma_yy);
    depth = std::sqrt(8 * sigma_zz);

    area_eff = contrail_properties::new_effective_area_from_sigma(sigma_yy, sigma_zz, sigma_yz);
}

void CoCiP::calc_timestep_contrail_evolution(const double length_ratio, const double dt_s) {

    // Radiative heating
    // Like pycontrails, uses heat rates from previous time step
    cumul_heat = std::min(1.5, cumul_heat + heat_rate * dt_s);
    cumul_differential_heat -= d_heat_rate * dt_s;

    // Update plume parameters including width, depth, and area_eff
    // Like pycontrails, uses the values of diffuse_h/_v from previous time step
    plume_temporal_evolution(length_ratio, dt_s);

    // Calculate new plume mass per distance using new area_eff and rho_air
    double plume_mass_per_m_new = contrail_properties::plume_mass_per_distance(area_eff,
        met->rho_air);
    
    double q_sat = thermo::q_sat_ice(met->air_temperature, met->air_pressure);
    double q_sat_old = thermo::q_sat_ice(met->air_temperature_old, met->air_pressure_old);
    // The q_sat calculations should use temperature inside the contrail, as they do,
    // but the iwc calculation should use ambient specific humidity, i.e. with ambient temperature
    double iwc_new = contrail_properties::new_ice_water_content(iwc, met->specific_humidity_old,
        met->specific_humidity, q_sat_old, q_sat, plume_mass_per_m, plume_mass_per_m_new);
    double n_ice_per_m_new = contrail_properties::new_ice_particle_number(n_ice_per_m, dn_dt_agg,
        dn_dt_turb, length_ratio, dt_s);
    
    // Update saved values
    plume_mass_per_m = plume_mass_per_m_new; // overwritten in calc_contrail_properties anyway
    iwc = iwc_new;
    n_ice_per_m = n_ice_per_m_new;

    calc_contrail_properties();
    calc_radiative_properties();
    
    persistent = contrail_properties::contrail_persistent(tau_contrail, n_ice_per_vol, params);
    
    // No energy forcing because requires length, can be calculated externally if desired
}

void CoCiP::evolve(const double length_ratio, const double a, const double dt_s) {
    // calc_continuous not required (do not need to know segment chain)
    
    // Set heading angle; required for calculating wind shear normal used in
    // plume_temporal_evolution 
    set_heading(a);

    // Altitude is updated here before update_met_calculations rather than in
    // calc_timestep_contrail_evolution to better suit the order of calculations
    altitude = contrail_properties::altitude_after_settling(altitude, terminal_fall_speed, dt_s);

    update_met_calculations();

    calc_timestep_contrail_evolution(length_ratio, dt_s);

    // Ignore contrail contrail overlapping at least for now
}