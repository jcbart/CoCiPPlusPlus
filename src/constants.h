#ifndef CONSTANTS_H
#define CONSTANTS_H

// General constants
namespace constants {

constexpr double PI = 3.141592653589793238462643383279502884197169399;

// Radians per degree (pi / 180)
constexpr double RAD_PER_DEG = PI / 180;

// Degrees per radian (180 / pi)
constexpr double DEG_PER_RAD = 1 / RAD_PER_DEG;

// Absolute zero in Celcius
constexpr double ABSOLUTE_ZERO = -273.15;

// Acceleration due to gravity (m s-2)
constexpr double GRAVITY = 9.80665;

// Radius of Earth (m)
constexpr double RADIUS_EARTH = 6371229;

// Standard surface pressure (Pa)
constexpr double P_SURFACE = 101325;

// Reference pressure for potential temperature (Pa)
constexpr double P_REF = 1e5;

// Isobaric heat capacity of dry air (J kg-1 K-1)
constexpr double c_pd = 1004;

// Isobaric heat capacity of water vapor (J kg-1 K-1)
constexpr double c_pv = 1870;

// Molar mass of dry air (kg mol-1)
constexpr double M_d = 28.9647e-3;

// Ratio of specific heat capacity at constant pressure to that at constant volume for a diatomic
// ideal gas 
// Used in adiabatic heating
constexpr double GAMMA = 1.4;

// Molar gas constant (J mol-1 K-1)
constexpr double R = 8.314462618;

// Gas constant of dry air (J kg-1 K-1)
constexpr double R_d = 287.05;

// Gas constant of water vapor (J kg-1 K-1)
constexpr double R_v = 461.51;

// Ratio of gas constant for dry air / gas constant for water vapor (R_d / R_v)
constexpr double EPSILON = R_d / R_v;

// Density of ice (kg m-3)
constexpr double RHO_ICE = 917;

// Annual average incident solar radiation (W m-2)
constexpr double SOLAR_CONSTANT = 1361;

// Real refractive index of ice
constexpr double MU_ICE = 1.31;

// (Approximate) wavelength of visible light
constexpr double LAMBDA_LIGHT = 550e-9;

// Ratio between the volume mean radius and the effective radius (uncertainty +/- 0.3)
constexpr double c_r = 0.9;

}

// Radiative heating constants
namespace rad_heat {

// Coefficients for shortwave differential heating rate

constexpr double dacth = 0.205747e01;
constexpr double dacth3 = 0.898366e00;
constexpr double dbcth = 0.791045e00;
constexpr double dccth = 0.612725e00;
constexpr double ddcth = 0.517342e-02;
constexpr double dexalb = 0.267568e01;
constexpr double dfrsw = 0.139286e01;
constexpr double dgalbs = 0.178497e01;
constexpr double d_gamma = 0.142104e01;
constexpr double d_gamma_s = 0.882497e00;
constexpr double dqsw = 0.631427e-01;
constexpr double draddsw = 0.261780e00;
constexpr double dtt = 0.339171e-01;

// Coefficients for longwave differential heating rate

constexpr double dak = 0.357181e01;
constexpr double dcrhi = 0.623019e-01;
constexpr double ddelta = 0.198000e01;
constexpr double dfrlw = 0.609262e00;
constexpr double dqlw = 0.100000e-05;
constexpr double dqrlw = 0.160286e00;
constexpr double draddlw = 0.898529e-05;
constexpr double dsigma = 0.159884e-06;

// Coefficients for shortwave heating rate

constexpr double acth = 0.156899e01;
constexpr double bcth = 0.875130e00;
constexpr double ccth = 0.112445e01;
constexpr double dcth = 0.236688e-01;
constexpr double exal_b = 0.410705e00;
constexpr double fr_sw = 0.537577e01;
constexpr double gamma_r = 0.762254e00;
constexpr double q_sw = 0.454176e-01;
constexpr double radd_sw = 0.991554e00;
constexpr double ttt = 0.985031e-01;

// Coefficients for longwave heating rate

constexpr double ak= 0.294930e01;
constexpr double crhi = 0.174422e00;
constexpr double czlw = 0.393884e-01;
constexpr double delta = 0.860746e00;
constexpr double fr_lw = 0.760423e00;
constexpr double q_lw = 0.152075e02;
constexpr double radd_lw = 0.308486e-02;
constexpr double sigma = 0.253499e-04;

}

// Radiative forcing constants
namespace rf_const {

constexpr int num_habits = 8;

constexpr double k_t[num_habits] = {1.93466, 1.95456, 1.95994, 1.95906, 1.94397, 1.95123, 2.30363, 1.94611};

constexpr double T_0[num_habits] = {152.237, 152.724, 152.923, 152.360, 151.879, 152.318, 165.692, 153.073};

constexpr double delta_t[num_habits] = {0.940846, 0.808397, 0.736222, 0.675591, 0.748757, 0.708515, 0.927592, 0.795527};

constexpr double delta_lr[num_habits] = {0.211276, 0.341194, 0.325496, 0.255921, 0.170265, 1.65441, 0.201949, 0};

constexpr double delta_lc[num_habits] = {0.159942, 0.0958129, 0.0924850, 0.0462023, 0.132925, 0.0870067, 0.0626339, 0.0665289};

constexpr double t_a[num_habits] = {0.879119, 0.901701, 0.881812, 0.899144, 0.879896, 0.883212, 0.899096, 1.00744};

constexpr double A_mu[num_habits] = {0.361226, 0.294072, 0.343894, 0.317866, 0.337227, 0.310978, 0.342593, 0.269179};

constexpr double C_mu[num_habits] = {0.709300, 0.678016, 0.687546, 0.675315, 0.712041, 0.713317, 0.660267, 0.545716};

constexpr double delta_sr[num_habits] = {0.149851, 0.0254270, 0.0238836, 0.0463724, 0.0478892, 0.0700234, 0.0517942, 0};

constexpr double F_r[num_habits] = {0.511852, 0.576911, 0.597351, 0.225750, 0.550734, 0.817858, 0.249004, 0};

constexpr double gamma_lower[num_habits] = {0.323166, 0.392598, 0.356189, 0.345040, 0.407515, 0.523604, 0.310853, 0.274741};

constexpr double gamma_upper[num_habits] = {0.241507, 0.347023, 0.288452, 0.296813, 0.327857, 0.437560, 0.274710, 0.208154};

constexpr double B_mu[num_habits] = {1.67592, 1.55687, 1.71065, 1.55843, 1.70782, 1.71789, 1.56399, 1.59015};

constexpr double delta_sc[num_habits] = {0.157017, 0.143274, 0.167995, 0.148547, 0.173036, 0.162442, 0.171855, 0.213488};

constexpr double delta_sc_aps[num_habits] = {0.229574, 0.197611, 0.245036, 0.204875, 0.248328, 0.254029, 0.244051, 0.302246};

}

#endif