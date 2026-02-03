#ifndef WIND_SHEAR_H
#define WIND_SHEAR_H

#include <cmath>

namespace wind_shear {

// Calculate the magnitude of the wind shear (s-1) between two altitudes 
inline double wind_shear(double u_wind_top, double u_wind_btm, double v_wind_top,
    double v_wind_btm, double dz) {
    
    double du_dz = (u_wind_top - u_wind_btm) / dz;
    double dv_dz = (v_wind_top - v_wind_btm) / dz;
    return std::sqrt(du_dz*du_dz + dv_dz*dv_dz);
}

// Calculate the wind shear normal to the contrail heading (s-1) between two altitudes 
inline double wind_shear_normal(double u_wind_top, double u_wind_btm, double v_wind_top,
    double v_wind_btm, double cos_a, double sin_a, double dz) {
    
    double du_dz = (u_wind_top - u_wind_btm) / dz;
    double dv_dz = (v_wind_top - v_wind_btm) / dz;
    return (dv_dz * cos_a - du_dz * sin_a);
}

// Calculate the multiplication factor to enhance the wind shear based on contrail depth
// If effective_vertical_resolution (m) or wind_shear_enhancement_exponent is zero,
// there is no enhancement
constexpr double wind_shear_enhancement_factor(double contrail_depth,
    double effective_vertical_resolution, double wind_shear_enhancement_exponent) {

    return (
        (contrail_depth > 0) ?
        0.5 * (1 + std::pow(effective_vertical_resolution / contrail_depth,
                            wind_shear_enhancement_exponent))
        : 1
    );
}

}

#endif