#include <cmath>
#include <algorithm>
#include <vector>
#include "radiative_forcing.h"
#include "constants.h"

std::vector<double> radiative_forcing::calc_habit_weights(double r_vol_um,
    const std::vector<std::vector<double>>& habit_distributions,
    const std::vector<double>& radius_threshold_um) {

    int idx = habit_weight_regime_idx(r_vol_um, radius_threshold_um);
    return habit_distributions[idx];
}

std::vector<double> radiative_forcing::effective_radius_by_habit(double r_vol_um,
    const std::vector<double>& habit_weights) {
    
    std::vector<double> r_eff_um(rf_const::num_habits);
    r_eff_um[0] = (habit_weights[0] > 0) ? effective_radius_sphere(r_vol_um) : 0;
    r_eff_um[1] = (habit_weights[1] > 0) ? effective_radius_solid_column(r_vol_um) : 0;
    r_eff_um[2] = (habit_weights[2] > 0) ? effective_radius_hollow_column(r_vol_um) : 0;
    r_eff_um[3] = (habit_weights[3] > 0) ? effective_radius_rough_aggregate(r_vol_um) : 0;
    r_eff_um[4] = (habit_weights[4] > 0) ? effective_radius_rosette(r_vol_um) : 0;
    r_eff_um[5] = (habit_weights[5] > 0) ? effective_radius_plate(r_vol_um) : 0;
    r_eff_um[6] = (habit_weights[6] > 0) ? effective_radius_droxtal(r_vol_um) : 0;
    r_eff_um[7] = (habit_weights[7] > 0) ? effective_radius_myhre(r_vol_um) : 0;
    return r_eff_um;
}

double radiative_forcing::longwave_radiative_forcing(double r_vol_um, double olr,
    double air_temperature, double tau_contrail, double tau_cirrus,
    const std::vector<double>& habit_weights) {
    
    std::vector<double> r_eff_um = effective_radius_by_habit(r_vol_um, habit_weights);

    std::vector<double> e_lw = olr_reduction_natural_cirrus(tau_cirrus);
    std::vector<double> f_lw = contrail_effective_emissivity(r_eff_um);
    
    double rf_lw_sum = 0;
    
    // Find rf_lw due to each habit, weight it, and add to rf_lw_sum
    for (int i = 0; i < rf_const::num_habits; i++) {
        double rf_lw_habit = (
            (olr - rf_const::k_t[i] * (air_temperature - rf_const::T_0[i]))
            * e_lw[i]
            * (1 - std::exp(-rf_const::delta_t[i] * f_lw[i] * tau_contrail))
        );

        rf_lw_sum += std::max(0., rf_lw_habit) * habit_weights[i];
    }

    return rf_lw_sum;
}

double radiative_forcing::shortwave_radiative_forcing(double r_vol_um, double sdr, double rsr,
    double sd0, double tau_contrail, double tau_cirrus, const std::vector<double>& habit_weights) {
    
    // If not daytime, return zero
    if (!(sdr > 0)) {
        return 0;
    }

    double albedo = calc_albedo(sdr, rsr);
    double mue = std::min(sdr / sd0, 1.);

    // Calculated again??
    std::vector<double> r_eff_um = effective_radius_by_habit(r_vol_um, habit_weights);

    std::vector<double> alpha_c = contrail_albedo(tau_contrail, mue, r_eff_um);

    std::vector<double> e_sw = effective_tau_cirrus(tau_cirrus, mue);

    double rf_sw_sum = 0;
    
    // Find rf_sw due to each habit, weight it, and add to rf_sw_sum
    for (int i = 0; i < rf_const::num_habits; i++) {
        double rf_sw_habit = -sdr * std::pow(rf_const::t_a[i] - albedo, 2) * alpha_c[i] * e_sw[i];
        rf_sw_sum += std::min(0., rf_sw_habit) * habit_weights[i];
    }

    return rf_sw_sum;
}

std::vector<double> radiative_forcing::olr_reduction_natural_cirrus(double tau_cirrus) {
    std::vector<double> result(rf_const::num_habits);
    for (int i = 0; i < rf_const::num_habits; i++) {
        result[i] = std::exp(-rf_const::delta_lc[i] * tau_cirrus);
    }
    return result;
}

std::vector<double> radiative_forcing::contrail_effective_emissivity(
    std::vector<double>& r_eff_um
) {
    std::vector<double> result(rf_const::num_habits);
    for (int i = 0; i < rf_const::num_habits; i++) {
        result[i] = 1 - std::exp(-rf_const::delta_lr[i] * r_eff_um[i]);
    }
    return result;
}

std::vector<double> radiative_forcing::contrail_albedo(double tau_contrail, double mue,
    std::vector<double>& r_eff_um) {

    std::vector<double> alpha_c(rf_const::num_habits);

    for (int i = 0; i < rf_const::num_habits; i++) {
        double tau_aps =
            tau_contrail * (1 -
                            rf_const::F_r[i]
                            * (1 - std::exp(-rf_const::delta_sr[i] * r_eff_um[i]))
                           );
        double tau_eff = tau_aps / (mue + 1e-6);
        double r_c = 1 - std::exp(-rf_const::gamma_upper[i] * tau_eff);
        double r_c_aps = std::exp(-rf_const::gamma_lower[i] * tau_eff);
        double f_mu = std::pow(2 * (1 - mue), rf_const::B_mu[i]) - 1;
        alpha_c[i] = r_c * (rf_const::C_mu[i] + (rf_const::A_mu[i] * r_c_aps * f_mu));
    }
    return alpha_c;
}

std::vector<double> radiative_forcing::effective_tau_cirrus(double tau_cirrus, double mue) {
    double tau_cirrus_eff = tau_cirrus / (mue + 1e-6);

    std::vector<double> e_sw(rf_const::num_habits);

    for (int i = 0; i < rf_const::num_habits; i++) {
        e_sw[i] = std::exp(tau_cirrus * rf_const::delta_sc_aps[i]
                           - tau_cirrus_eff * rf_const::delta_sc[i]);
    }
    return e_sw;
}