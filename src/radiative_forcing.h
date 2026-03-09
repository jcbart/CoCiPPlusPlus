#ifndef RADIATIVE_FORCING_H
#define RADIATIVE_FORCING_H

#include <cmath>
#include <algorithm>
#include <vector>

namespace radiative_forcing {

// Assigns weights to each ice particle habit
// Equivalent of pycontrails' habit_weights
std::vector<double> calc_habit_weights(double r_vol_um,
    const std::vector<std::vector<double>>& habit_distributions,
    const std::vector<double>& radius_threshold_um);

// Determine regime of ice particle habits based on contrail ice particle volume mean radius
int habit_weight_regime_idx(double r_vol_um,
    const std::vector<double>& radius_threshold_um);

// Calculates a vector of effective radii, one for each ice crystal habit (um)
// Uses habit_weights to ignore where weight is zero
std::vector<double> effective_radius_by_habit(double r_vol_um,
    const std::vector<double>& habit_weights);

// Effective radius of an ice crystal with a sphere habit (um)
constexpr double effective_radius_sphere(double r_vol_um) {
    return std::min(r_vol_um, 25.);
}

// Effective radius of an ice crystal with a solid column habit (um)
constexpr double effective_radius_solid_column(double r_vol_um) {
    double r_eff_um = (
        0.2588 * std::exp(-(6.912e-3 * r_vol_um)) + 0.6372 * std::exp(-(3.142e-4 * r_vol_um))
    ) * r_vol_um;
    if (r_vol_um <= 42.2) {
        r_eff_um = 0.824 * r_vol_um;
    }
    return std::min(r_eff_um, 45.);
}

// Effective radius of an ice crystal with a hollow column habit (um)
constexpr double effective_radius_hollow_column(double r_vol_um) {
    double r_eff_um = (
        0.2281 * std::exp(-(7.359e-3 * r_vol_um)) + 0.5651 * std::exp(-(3.350e-4 * r_vol_um))
    ) * r_vol_um;
    if (r_vol_um <= 39.7) {
        r_eff_um = 0.729 * r_vol_um;
    }
    return std::min(r_eff_um, 45.);
}

// Effective radius of an ice crystal with a rough aggregate habit (um)
constexpr double effective_radius_rough_aggregate(double r_vol_um) {
    return std::min(0.574 * r_vol_um, 45.);
}

// Effective radius of an ice crystal with a rosette habit (um)
constexpr double effective_radius_rosette(double r_vol_um) {
    double r_eff_um = r_vol_um * (
        0.1770 * std::exp(-(2.144e-2 * r_vol_um)) + 0.4267 * std::exp(-(3.562e-4 * r_vol_um))
    );
    return std::min(r_eff_um, 45.);
}

// Effective radius of an ice crystal with a plate habit (um)
constexpr double effective_radius_plate(double r_vol_um) {
    double r_eff_um = r_vol_um * (
        0.1663 + 0.3713 * std::exp(-(0.0336 * r_vol_um)) + 0.3309 * std::exp(-(0.0035 * r_vol_um))
    );
    return std::min(r_eff_um, 45.);
}

// Effective radius of an ice crystal with a droxtal habit (um)
constexpr double effective_radius_droxtal(double r_vol_um) {
    return std::min(0.94f * r_vol_um, 45.);
}

// Effective radius of an ice crystal with a myhre habit (um)
constexpr double effective_radius_myhre(double r_vol_um) {
    return std::min(r_vol_um, 45.);
}

// Calculate the local contrail longwave radiative forcing (positive) (W m-2)
// Ignoring r_eff_um parameter
double longwave_radiative_forcing(double r_vol_um, double olr,
    double air_temperature, double tau_contrail, double tau_cirrus,
    const std::vector<double>& habit_weights);

// Calculate the local contrail shortwave radiative forcing (negative) (W m-2)
double shortwave_radiative_forcing(double r_vol_um, double sdr, double rsr,
    double sd0, double tau_contrail, double tau_cirrus, const std::vector<double>& habit_weights);

// Calculate the local contrail net radiative forcing (rf_lw + rf_sw) (W m-2)
constexpr double net_radiative_forcing(double rf_lw, double rf_sw) {
    return (rf_lw + rf_sw);
}

// Calculate relative reduction in outgoing longwave radiation (OLR) due to the presence of natural
// cirrus for each ice particle habit
std::vector<double> olr_reduction_natural_cirrus(double tau_cirrus);

// Calculate the effective emissivity of the contrail for each ice particle habit
std::vector<double> contrail_effective_emissivity(
    std::vector<double>& r_eff_um
);

// Calculate total albedo [0 - 1]
// Equivalent of pycontrails' albedo
constexpr double calc_albedo(double sdr, double rsr) {
    return std::max(0., std::min(1., rsr / sdr));
}

// Calculate contrail albedo for each ice particle habit
std::vector<double> contrail_albedo(double tau_contrail, double mue,
    std::vector<double>& r_eff_um);

// Calculate the effective optical depth of natural cirrus above the contrail for each ice particle
// habit
std::vector<double> effective_tau_cirrus(double tau_cirrus, double mue);

}

#endif