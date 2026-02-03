#ifndef SAC_H
#define SAC_H

namespace sac {

// Calculate the slope of the mixing line in a temperature-humidity diagram (Pa K-1)
double slope_mixing_line(double specific_humidity, double air_pressure, double engine_efficiency,
    double ei_h2o, double q_fuel);

// Calculate temperature at which liquid saturation curve has slope G (K)
double T_sat_liquid(double G);

// Calculate critical relative humidity threshold of contrail formation ()
double rh_critical_sac(double air_temperature, double T_sat_liq, double G);

// Estimate temperature threshold for persistent contrail formation (K)
double T_critical_sac(double T_LM, double rh, double G);

}

#endif