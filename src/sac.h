#ifndef SAC_H
#define SAC_H

#include <cmath>
#include <algorithm>
#include <limits>
#include "thermo.h"
#include "constants.h"

namespace sac {

// Calculate the slope of the mixing line in a temperature-humidity diagram (Pa K-1)
constexpr double slope_mixing_line(double specific_humidity, double air_pressure,
    double engine_efficiency, double ei_h2o, double q_fuel) {
    
    double c_pm = thermo::c_pm(specific_humidity);
    double G = (ei_h2o * c_pm * air_pressure)
                / (constants::EPSILON * q_fuel * (1. - engine_efficiency));
    return G;
}

// Calculate temperature at which liquid saturation curve has slope G (K)
constexpr double T_sat_liquid(double G) {
    double log_ = std::log(G - 0.053);
    double T_sat_liquid_ = -46.46 - constants::ABSOLUTE_ZERO + 9.43*log_ + 0.72*log_*log_;
    return T_sat_liquid_;
}

// Calculate critical relative humidity threshold of contrail formation ()
constexpr double rh_critical_sac(double air_temperature, double T_sat_liq, double G) {
    double e_sat_T_sat_liq = thermo::e_sat_liquid(T_sat_liq);
    double e_sat_T = thermo::e_sat_liquid(air_temperature);
    double rh_crit = (G * (air_temperature - T_sat_liq) + e_sat_T_sat_liq) / e_sat_T;
    rh_crit = std::max(0., std::min(1., rh_crit));

    if (air_temperature > T_sat_liq) {
        rh_crit = std::numeric_limits<double>::infinity();
    }
    return rh_crit;
}

// Estimate temperature threshold for persistent contrail formation (K)
double T_critical_sac(double T_LM, double rh, double G);

}

#endif