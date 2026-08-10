#include <limits>
#include <vector>
#include <string>
#include <yaml-cpp/yaml.h>
#include <CoCiP++/CoCiPLog.h>
#include <CoCiP++/params.h>
#include "constants.h"

void Params::readYAML() {
    // Lambda function to parse a node which contains infinity by default
    auto parseDefaultInfNode = []<typename T>(const YAML::Node& node) -> T {
        if (node && !node.IsNull()) {
            return node.as<T>();
        }
        return std::numeric_limits<T>::infinity();
    };
    // Lambda function to parse a node and return std::nullopt if not found
    auto parseOptionalNode = []<typename T>(const YAML::Node& node) -> std::optional<T> {
        if (node && !node.IsNull()) {
            return node.as<T>();
        }
        return std::nullopt;
    };

    YAML::Node config = YAML::LoadFile("CoCiP-params.yaml");

    interp_with_pressure = config["interp_with_pressure"].as<bool>();

    dz_m = config["dz_m"].as<double>();

    effective_vertical_resolution = parseOptionalNode.operator()<double>(
        config["effective_vertical_resolution"]
    );

    initial_wake_vortex_depth = config["initial_wake_vortex_depth"].as<double>();

    turbulent_vertical_velocity_scale = config["turbulent_vertical_velocity_scale"].as<double>();

    sedimentation_impact_factor = config["sedimentation_impact_factor"].as<double>();

    wind_shear_enhancement_exponent = config["wind_shear_enhancement_exponent"].as<double>();

    min_ice_particle_number_nvpm_ei_n = parseDefaultInfNode.operator()<double>(
        config["min_ice_particle_number_nvpm_ei_n"]
    );

    max_depth = parseDefaultInfNode.operator()<double>(config["max_depth"]);

    max_horizontal_diffusivity = parseDefaultInfNode.operator()<double>(
        config["max_horizontal_diffusivity"]
    );

    max_vertical_diffusivity = parseDefaultInfNode.operator()<double>(
        config["max_vertical_diffusivity"]
    );

    radiative_heating_effects = config["radiative_heating_effects"].as<bool>();

    radius_threshold_um = config["radius_threshold_um"].as<std::vector<double>>();

    habit_distributions = config["habit_distributions"].as<std::vector<std::vector<double>>>();

    // Ensure habit_distributions and radius_threshold_um match in size
    if (habit_distributions.size() != 1 + radius_threshold_um.size()) {
        std::string msg = "habit_distributions size must be 1 + radius_threshold_um size\n";
        msg += "habit_distributions: " + std::to_string(habit_distributions.size())
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

    revised_contrail_ice_budget = config["revised_contrail_ice_budget"].as<bool>();

    rf_sw_enhancement_factor = config["rf_sw_enhancement_factor"].as<double>();

    rf_lw_enhancement_factor = config["rf_lw_enhancement_factor"].as<double>();

    min_tau = config["min_tau"].as<double>();

    max_tau = config["max_tau"].as<double>();

    min_n_ice_per_m3 = config["min_n_ice_per_m3"].as<double>();

    max_n_ice_per_m3 = config["max_n_ice_per_m3"].as<double>();

    isInitialised = true;
}