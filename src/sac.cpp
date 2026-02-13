#include <cmath>
#include <algorithm>
#include <functional>
#include <limits>
#include "sac.h"
#include "thermo.h"
#include "constants.h"

double sac::slope_mixing_line(double specific_humidity, double air_pressure, double engine_efficiency,
    double ei_h2o, double q_fuel) {
    
    double c_pm = thermo::c_pm(specific_humidity);
    double G = (ei_h2o * c_pm * air_pressure)
              / (constants::EPSILON * q_fuel * (1. - engine_efficiency));
    return G;
}

double sac::T_sat_liquid(double G) {
    double log_ = std::log(G - 0.053);
    double T_sat_liquid_ = -46.46 - constants::ABSOLUTE_ZERO + 9.43*log_ + 0.72*log_*log_;
    return T_sat_liquid_;
}

double sac::rh_critical_sac(double air_temperature, double T_sat_liq, double G) {
    double e_sat_T_sat_liq = thermo::e_sat_liquid(T_sat_liq);
    double e_sat_T = thermo::e_sat_liquid(air_temperature);
    double rh_crit = (G * (air_temperature - T_sat_liq) + e_sat_T_sat_liq) / e_sat_T;
    rh_crit = std::max(0., std::min(1., rh_crit));

    if (air_temperature > T_sat_liq) {
        rh_crit = std::numeric_limits<double>::infinity();
    }
    return rh_crit;
}

// A Newton-Raphson solver for finding the root of function f given its derivative fprime
double newtonRaphsonSolver(std::function<double(double)> f,
                           std::function<double(double)> fprime,
                           double x0, double tol = 1e-7, int maxIter = 100) {
    double x = x0;
    for (int i = 0; i < maxIter; i++) {
        double f_x = f(x);
        double fprime_x = fprime(x);
        if (std::abs(fprime_x) < 1e-12) {
            // Avoid division by 0
            break;
        }
        double x_next = x - f_x / fprime_x;
        if (std::abs(x_next - x) < tol) {
            return x_next;
        }
        x = x_next;
    }
    return x;
}

double sac::T_critical_sac(double T_LM, double rh, double G) {
    if (!(rh < 0.999 && std::isfinite(T_LM))) {
        return T_LM;
    }
    double e_L_of_T_LM = thermo::e_sat_liquid(T_LM);

    auto Schumann_eq11 = [T_LM, e_L_of_T_LM, rh, G](double T) -> double {
        return (T - T_LM + (e_L_of_T_LM - rh * thermo::e_sat_liquid(T)) / G);
    };

    auto Schumann_eq11_prime = [rh, G](double T) -> double {
        return (1. - rh * thermo::e_sat_liquid_prime(T) / G);
    };

    double init_guess = T_LM - 1;
    double tol = 1e-7;
    int maxIter = 10;
    double T_crit_sac = newtonRaphsonSolver(
        Schumann_eq11, Schumann_eq11_prime, init_guess, tol, maxIter
    );
    return T_crit_sac;
}