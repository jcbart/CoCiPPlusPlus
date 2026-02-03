#ifndef PARAMS_H
#define PARAMS_H

#include <vector>

struct Params {
    bool isInitialised = false;

    // Initial wake vortex depth scaling factor
    // This factor scales max contrail downward displacement after the wake vortex phase
    // to set the initial contrail depth
    double initial_wake_vortex_depth;

    // Sedimentation impact factor
    double sedimentation_impact_factor;

    // n in Schumann (2012) eq. 39
    double wind_shear_enhancement_exponent;

    // Lower bound for nvpm_ei_n to account for ambient aerosol particles for newer engines (kg-1)
    double min_ice_particle_number_nvpm_ei_n;

    // Upper bound for contrail plume depth (m)
    double max_depth;

    // Upper bound for contrail horizontal plume diffusivities (m2 s-1)
    double max_horizontal_diffusivity;

    // Upper bound for contrail vertical plume diffusivities (m2 s-1)
    // Set to nullptr to for no max
    double max_vertical_diffusivity;

    std::vector<double> radius_threshold_um;

    std::vector<std::vector<double>> habit_distributions;

    double rf_sw_enhancement_factor;

    double rf_lw_enhancement_factor;

    // Minimum contrail optical depth
    double min_tau;

    // Maximum contrail optical depth
    double max_tau;

    // Minimum contrail ice particle number per volume of air
    double min_n_ice_per_m3;

    // Maximum contrail ice particle number per volume of air
    double max_n_ice_per_m3;

    void readYAML();
};

#endif