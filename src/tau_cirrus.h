#ifndef TAU_CIRRUS_H
#define TAU_CIRRUS_H

#include <cmath>
#include <algorithm>
#include "thermo.h"
#include "constants.h"

namespace tau_cirrus {

// calc_tau_cirrus is in met.h/.cpp since it depends on met type

// Approximate conversion of height (m) to geopotential height (m)
constexpr double height_to_geopt_height(double z) {
    return (z * constants::RADIUS_EARTH / (constants::RADIUS_EARTH + z));
}

// Calculate the effective extinction coefficient for spectral range 0.2-0.69 um (m-1)
// for layer with cloud ice water content ciwc (kg (kg dry air)-1), temperature T (K),
// and pressure P (Pa)
inline double cirrus_effective_extinction_coef(double ciwc, double T, double P) {
    ciwc = std::max(0., ciwc);
    const double a_0_beta = -1.30817e-4;
    const double a_1_beta = 2.52883e0;
    double rho_air = thermo::rho_d(T, P);
    double riwc = ciwc * rho_air * 1000;
    double tiwc = T + constants::ABSOLUTE_ZERO + 190;
    double d_eff = 45.8966 * std::pow(riwc, 0.2214) + 0.7957 * tiwc * std::pow(riwc, 0.2535);
    d_eff = std::max(d_eff, 10.);
    return (riwc * (a_0_beta + a_1_beta / d_eff));
}

}

#endif