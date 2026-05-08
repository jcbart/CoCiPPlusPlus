#ifndef GEO_H
#define GEO_H

#include <CoCiP++/CoCiPTime.h>
#include "constants.h"

namespace geo {

// Calculate the orbital position of Earth relative to the start of the year (radians)
constexpr double orbital_position(double day_of_year) {
    return (constants::RAD_PER_DEG * (360. * (day_of_year / 365.25)));
}

// Calculate the solar electromagnetic radiation per unit area from orbital position (W m-2)
constexpr double solar_constant(double theta_rad) {
    return constants::SOLAR_CONSTANT * (
        1.00011
        + (0.034221 * std::cos(theta_rad))
        + (0.001280 * std::sin(theta_rad))
        + (0.000719 * std::cos(theta_rad * 2))
        + (0.000077 * std::sin(theta_rad * 2))
    );
}

// Calculate the solar declination angle (radians) from the orbital position in radians (theta_rad)
constexpr double solar_declination_angle(double theta_rad) {
    return (
        0.396372
        - (22.91327 * std::cos(theta_rad))
        + (4.02543 * std::sin(theta_rad))
        - (0.387205 * std::cos(2 * theta_rad))
        + (0.051967 * std::sin(2 * theta_rad))
        - (0.154527 * std::cos(3 * theta_rad))
        + (0.084798 * std::sin(3 * theta_rad))
    );
}

// Calculate correction to the solar hour angle due to Earth's orbital location (degrees)
constexpr double orbital_correction_for_solar_hour_angle(double theta_rad) {
    return (
        0.004297
        + (0.107029 * std::cos(theta_rad))
        - (1.837877 * std::sin(theta_rad))
        - (0.837378 * std::cos(2 * theta_rad))
        - (2.340475 * std::sin(2 * theta_rad))
    );
}

// Calculate the Sun's East to West angular displacement around the polar axis (degrees)
constexpr double solar_hour_angle(double longitude, double hour, double theta_rad) {
    return (((hour - 12) * 15) + longitude + orbital_correction_for_solar_hour_angle(theta_rad));
}

// Calculate the cosine of the solar zenith angle (radians)
constexpr double cosine_solar_zenith_angle(double longitude, double latitude,
    const CoCiPTime& datetime, double theta_rad) {
    
    double lat_rad = constants::RAD_PER_DEG * latitude;
    double sdec_rad = constants::RAD_PER_DEG * solar_declination_angle(theta_rad);
    double sha_rad = constants::RAD_PER_DEG * solar_hour_angle(longitude, datetime.hour_of_day(), theta_rad);
    return (
        std::sin(lat_rad) * std::sin(sdec_rad) + (
        std::cos(lat_rad) * std::cos(sdec_rad) * std::cos(sha_rad)
        )
    );
}

// Calculate the instantaneous theoretical solar direct radiation (SDR; W m-2)
constexpr double solar_direct_radiation(double longitude, double latitude, const CoCiPTime& datetime,
    double threshold_cos_sza = 0) {
    
    double theta_rad = orbital_position(datetime.day_of_year());
    double solar_const = solar_constant(theta_rad);
    double cos_sza = cosine_solar_zenith_angle(longitude, latitude, datetime, theta_rad);
    return (cos_sza < threshold_cos_sza) ? 0 : (cos_sza * solar_const);
}

}

#endif