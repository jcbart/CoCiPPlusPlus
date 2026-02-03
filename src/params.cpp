#include <limits>
#include <vector>
#include <string>
#include <yaml-cpp/yaml.h>
#include <CoCiP++/params.h>
#include "constants.h"
#include "CoCiPLog.h"

void Params::readYAML() {
    YAML::Node node = YAML::LoadFile("CoCiP-params.yaml");

    initial_wake_vortex_depth = node["initial_wake_vortex_depth"].as<double>();

    sedimentation_impact_factor = node["sedimentation_impact_factor"].as<double>();

    wind_shear_enhancement_exponent = node["wind_shear_enhancement_exponent"].as<double>();

    min_ice_particle_number_nvpm_ei_n = (
        node["min_ice_particle_number_nvpm_ei_n"]
        && !node["min_ice_particle_number_nvpm_ei_n"].IsNull()
        ) ?
        node["min_ice_particle_number_nvpm_ei_n"].as<double>()
        : std::numeric_limits<double>::infinity();

    max_depth = node["max_depth"] ?
        node["max_depth"].as<double>() : std::numeric_limits<double>::infinity();

    max_horizontal_diffusivity = (
        node["max_horizontal_diffusivity"]
        && !node["max_horizontal_diffusivity"].IsNull()
        ) ?
        node["max_horizontal_diffusivity"].as<double>() : std::numeric_limits<double>::infinity();

    max_vertical_diffusivity = (
        node["max_vertical_diffusivity"]
        && !node["max_vertical_diffusivity"].IsNull()
        ) ?
        node["max_vertical_diffusivity"].as<double>() : std::numeric_limits<double>::infinity();

    radius_threshold_um = node["radius_threshold_um"].as<std::vector<double>>();

    habit_distributions = node["habit_distributions"].as<std::vector<std::vector<double>>>();

    // Ensure habit_distributions and radius_threshold_um match in size
    if (habit_distributions.size() != 1 + radius_threshold_um.size()) {
        std::string msg = "habit_distributions size must be 1 + radius_threshold_um size\n";
        msg    += "habit_distributions: " + std::to_string(habit_distributions.size())
                  + ", radius_threshold_um: " + std::to_string(radius_threshold_um.size());
        CoCiP_RaiseError(msg, __FILE__, __LINE__);
    }

    // Ensure each row of habit_distributions is of size rf_const::num_habits and sums to 1
    for (const std::vector<double>& row : habit_distributions) {
        if (row.size() != rf_const::num_habits) {
            std::string msg = "habit_distributions row size = " + std::to_string(row.size())
                + " does not equal num_habits = " + std::to_string(rf_const::num_habits);
            CoCiP_RaiseError(msg, __FILE__, __LINE__);
        }
        double sum = 0;
        for (size_t i = 0; i < row.size(); i++) {
            sum += row[i];
        }
        if (std::abs(sum - 1) > 1e-3) {
            std::string msg = "habit_distributions row must sum to 1, but sums to "
                + std::to_string(sum);
            CoCiP_RaiseError(msg, __FILE__, __LINE__);
        }
    }

    rf_sw_enhancement_factor = node["rf_sw_enhancement_factor"].as<double>();

    rf_lw_enhancement_factor = node["rf_lw_enhancement_factor"].as<double>();

    min_tau = node["min_tau"].as<double>();

    max_tau = node["max_tau"].as<double>();

    min_n_ice_per_m3 = node["min_n_ice_per_m3"].as<double>();

    max_n_ice_per_m3 = node["max_n_ice_per_m3"].as<double>();

    isInitialised = true;
}