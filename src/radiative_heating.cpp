#include <cmath>
#include <algorithm>
#include "radiative_heating.h"
#include "constants.h"

double radiative_heating::effective_heating_rate(double d_heat_rate, double cumul_rad_heat,
    double dtheta_dz, double depth) {

    if (cumul_rad_heat <= 0) {
        return 0;
    }
    double heat_denom = std::max(0.5 * dtheta_dz * depth, 0.);
    double heat_ratio = cumul_rad_heat / (cumul_rad_heat + heat_denom);
    return (d_heat_rate * heat_ratio);
}

double radiative_heating::differential_heating_rate(double air_temperature, double rh_i,
    double rho_air, double r_ice_vol, double depth_eff, double tau_contrail, double tau_cirrus,
    double sd0, double sdr, double rsr, double olr) {
    
    double r_ice_vol_um = r_ice_vol * 1e6;
    double cp_contrail = contrail_heat_capacity(rho_air, depth_eff);
    double d_heat_rate_sw = differential_heating_rate_shortwave(cp_contrail, r_ice_vol_um,
        tau_contrail, tau_cirrus, sd0, sdr, rsr);
    double d_heat_rate_lw = differential_heating_rate_longwave(air_temperature, rh_i, cp_contrail,
        r_ice_vol_um, tau_contrail, tau_cirrus, olr);
    return std::min(d_heat_rate_sw + d_heat_rate_lw, 0.);
}

double radiative_heating::differential_heating_rate_shortwave(double cp_contrail,
    double r_ice_vol_um, double tau_contrail, double tau_cirrus, double sd0, double sdr,
    double rsr) {
    
    using namespace rad_heat;

    if (sdr <= 0) {
        return 0;
    }

    double mue = std::min(sdr / sd0, 1.);
    double tau_eff = tau_contrail / (mue + 1e-6);

    return (
        (1 - std::exp(-dqsw * r_ice_vol_um))
        * (dtt * sdr - ddcth * rsr)
        * tau_contrail
        * (-d_gamma_s * tau_contrail + tau_eff)
        * (1 / cp_contrail)
        * (1 - dacth * mue + d_gamma * mue*mue + (dacth3 - 1) * mue*mue*mue)
        * std::exp(dbcth * tau_cirrus - dccth * tau_cirrus / (mue + 1e-6))
        * std::exp(-dgalbs * tau_contrail * std::pow(1 - mue, dexalb))
        * std::pow(mue, draddsw)
        * (1 + dfrsw * std::pow(1 - mue, 2))
    );
}

double radiative_heating::differential_heating_rate_longwave(double air_temperature, double rh_i,
    double cp_contrail, double r_ice_vol_um, double tau_contrail, double tau_cirrus, double olr) {
    
    using namespace rad_heat;

    double cool = dsigma * std::pow(air_temperature, dak);
    double epsc = 1 - std::exp(-ddelta * (tau_contrail + tau_cirrus));
    return (
        -dfrlw
        * (1 / cp_contrail)
        * (olr - cool)
        * (epsc / ddelta)
        * tau_contrail
        * std::exp(-dqlw * tau_cirrus)
        * std::max(1 - draddlw * 10 / (r_ice_vol_um + 30), 0.)
        * (1 - std::exp(-dqrlw * r_ice_vol_um))
        * std::exp(-(rh_i - 0.9) * dcrhi)
    );
}

double radiative_heating::heating_rate(double air_temperature, double rh_i, double rho_air,
    double r_ice_vol, double depth_eff, double tau_contrail, double tau_cirrus, double sd0,
    double sdr, double rsr, double olr) {
    
    double r_ice_vol_um = r_ice_vol * 1e6;
    double cp_contrail = contrail_heat_capacity(rho_air, depth_eff);
    double heat_rate_sw = heating_rate_shortwave(cp_contrail, r_ice_vol_um, tau_contrail,
        tau_cirrus, sd0, sdr, rsr);
    double heat_rate_lw = heating_rate_longwave(air_temperature, rh_i, cp_contrail, r_ice_vol_um,
        tau_contrail, tau_cirrus, olr);
    return (heat_rate_sw + heat_rate_lw);
}

double radiative_heating::heating_rate_shortwave(double cp_contrail, double r_ice_vol_um,
    double tau_contrail, double tau_cirrus, double sd0, double sdr, double rsr) {
    
    using namespace rad_heat;

    double mue = std::min(sdr / sd0, 1.);
    double tau_eff = tau_contrail / (mue + 1e-6);
    double heat_rate_sw = (
        (1 - std::exp(-q_sw * r_ice_vol_um))
        * (ttt * sdr + dcth * rsr)
        * tau_eff
        * (1 / cp_contrail)
        * (1 - acth * mue + gamma_r * mue*mue)
        * std::exp(bcth * tau_cirrus - ccth * tau_cirrus / (mue + 1e-6))
        * std::exp(-tau_contrail * std::pow(1 - mue, 2.))
        * mue
    );
    return std::max(heat_rate_sw, 0.);
}

double radiative_heating::heating_rate_longwave(double air_temperature, double rh_i,
    double cp_contrail, double r_ice_vol_um, double tau_contrail, double tau_cirrus, double olr) {
    
    using namespace rad_heat;

    double fzlw = std::exp(-(rh_i - 0.9) * czlw);
    double cool = sigma * std::pow(air_temperature, ak);
    double epsc = 1 - std::exp(-delta * (tau_contrail + tau_cirrus));
    double heat_rate_lw = (
        std::exp(-(rh_i - 0.9) * crhi)
        * fr_lw
        * (1 / cp_contrail)
        * (olr / fzlw - cool * fzlw)
        * (epsc / delta)
        * std::exp(-q_lw * tau_cirrus)
        * (1 + radd_lw / (r_ice_vol_um + 30))
    );
    return std::max(heat_rate_lw, 0.);
}