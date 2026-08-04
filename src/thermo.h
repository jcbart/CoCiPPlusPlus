#ifndef THERMO_H
#define THERMO_H

#include <cmath>
#include "constants.h"

namespace thermo {

// Calculate air density (kg m-3) for (T, P) assuming dry air
constexpr double rho_d(double T, double P) {
    return (P / (constants::R_d * T));
}

// Calculate isobaric heat capacity of moist air (J kg-1 K-1) given specific humidity (kg (kg air)-1)
constexpr double c_pm(double q) {
    return (constants::c_pd * (1 - q) + constants::c_pv * q);
}

// Calculate molecular diffusivity of water vapor (m2 s-1)
constexpr double diffusivity_water_vapor(double T, double P) {
    // Clip between -40 and 40 C
    T = std::max(-40 - constants::ABSOLUTE_ZERO, std::min(40 - constants::ABSOLUTE_ZERO, T));
    constexpr double T0 = 273.15;
    constexpr double P0 = 101325;
    return 0.0000211 * std::pow(T / T0, 1.94) * (P0 / P);
}

// Calculate saturation pressure of water vapor over ice (Pa) given T (K) according to Sonntag (1994)
constexpr double e_sat_ice(double T) {
    return (100.0 * std::exp(
        (-6024.5282 / T)
        + 24.7219
        + (0.010613868 * T)
        - (1.3198825e-5 * (T*T))
        - 0.49382577 * std::log(T)
    ));
}

// Calculate saturation pressure of water vapor over liquid water (Pa) given T (K) using Murphy and Koop (2005)
constexpr double mk05_e_sat_liquid(double T) {
    return std::exp(
        54.842763
        - 6763.22 / T
        - 4.21 * std::log(T)
        + 0.000367 * T
        + std::tanh(0.0415 * (T - 218.8))
        * (53.878 - 1331.22 / T - 9.44523 * std::log(T) + 0.014025 * T)
    );
}

// Calculate saturation water vapor pressure with respect to liquid water (Pa) given T (K)
// Configured to return mk05_e_sat_liquid(T)
constexpr double e_sat_liquid(double T) {
    return mk05_e_sat_liquid(T);
}

// Calculate the derivate of mk05_e_sat_liquid (Pa K-1)
constexpr double mk05_e_sat_liquid_prime(double T) {
    double tanh_term = std::tanh(0.0415 * (T - 218.8));
    return mk05_e_sat_liquid(T) * (
        6763.22 / (T*T)
        - 4.21 / T
        + 0.000367
        + 0.0415 * (1 - tanh_term*tanh_term) * (
            53.878 - 1331.22 / T - 9.44523 * std::log(T) + 0.014025 * T
        )
        + tanh_term * (1331.22 / (T*T) - 9.44523 / T + 0.014025)
    );
}

// Calculate the derivative of saturation water vapor pressure with respect to liquid water
// (Pa K-1) given T (K)
// Configured to return mk05_e_sat_liquid_prime(T)
constexpr double e_sat_liquid_prime(double T) {
    return mk05_e_sat_liquid_prime(T);
}

// Calculate saturation specific humidity over ice (kg (kg air)-1)
constexpr double q_sat_ice(double T, double P) {
    return (constants::EPSILON * e_sat_ice(T) / P);
}

// Calculate the relative humidity with respect to liquid water ()
// given q (kg (kg air)-1), T (K), and P (Pa)
// Equivalent to pycontrails' rh
constexpr double rh_l(double q, double T, double P) {
    return ((q * P) / (constants::EPSILON * e_sat_liquid(T)));
}

// Calculate the relative humidity with respect to ice ()
// given q (kg (kg air)-1), T (K), and P (Pa)
// Equivalent to pycontrails' rhi
constexpr double rh_i(double q, double T, double P) {
    return ((q * P) / (constants::EPSILON * e_sat_ice(T)));
}

// Calculate potential temperature (K) given temperature (K) and air pressure (Pa)
constexpr double T_potential(double T, double P) {
    return T * std::pow(constants::P_SURFACE / P, constants::R_d / constants::c_pd);
}

// Calculate temperature (K) given potential temperature (K) and air pressure (Pa)
constexpr double T_from_T_potential(double T_pot, double P) {
     return T_pot * std::pow(P / constants::P_REF, constants::R_d / constants::c_pd);
}

// Calculate the potential temperature gradient (K m-1) between two altitudes
constexpr double T_potential_gradient(double T_top, double P_top, double T_btm,
    double P_btm, double dz) {
    
    return (T_potential(T_top, P_top) - T_potential(T_btm, P_btm)) / dz;
}

// Calculate the Brunt-Vaisala frequency (s-1) where dtheta_dz is potential temperature gradient
constexpr double brunt_vaisala_frequency(double P, double T, double dtheta_dz) {
    return std::sqrt(std::max(1e-6, dtheta_dz) * constants::GRAVITY / T_potential(T, P));
}

}

#endif