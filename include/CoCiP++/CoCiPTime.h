#ifndef COCIPTIME_H
#define COCIPTIME_H

#include <chrono>

// Struct for storing a timepoint
// Wrapper for a std::chrono::sys_time
struct CoCiPTime {
    // Time point - uses milliseconds as precision
    std::chrono::sys_time<std::chrono::milliseconds> timepoint;

    // Empty constructor
    CoCiPTime() {}

    // Constructor from timepoint
    template<typename Duration>
    explicit CoCiPTime(std::chrono::sys_time<Duration> timepoint)
        : timepoint(std::chrono::time_point_cast<std::chrono::milliseconds>(timepoint)) {}

    // Constructor from values
    CoCiPTime(int year, int month, int day, int hours, int minutes, double seconds) {
        set(year, month, day, hours, minutes, seconds);
    }

    // Set from timepoint
    template<typename Duration>
    void set(std::chrono::sys_time<Duration> timepoint) {
        this->timepoint = std::chrono::time_point_cast<std::chrono::milliseconds>(timepoint);
    }

    // Set from values
    void set(int year, int month, int day, int hours, int minutes, float seconds) {
        std::chrono::year_month_day ymd = std::chrono::year{year} / month / day;
        std::chrono::sys_days date_part{ymd};
        auto time_part = std::chrono::hours{hours}
            + std::chrono::minutes(minutes)
            + std::chrono::duration<double>{seconds};
        timepoint = std::chrono::time_point_cast<std::chrono::milliseconds>(date_part + time_part);
    }

    // Returns number of days (including partial) passed since the start of the year
    constexpr double day_of_year() const {
        // Time point at start of day
        std::chrono::sys_time<std::chrono::days> datepoint
            = std::chrono::floor<std::chrono::days>(timepoint);
        
        // Start of day as year_month_day
        std::chrono::year_month_day ymd{datepoint};

        // 1st January of same year as a timepoint
        auto jan1 = std::chrono::sys_days{ymd.year() / std::chrono::January / 1};

        // Duration between start of day and timepoint
        std::chrono::milliseconds time_of_day = timepoint - datepoint;

        // Add number of whole days passed and fraction of day
        // (using ratio instead of duration_cast avoids truncation)
        return (
            (datepoint - jan1).count()
            + std::chrono::duration<double, std::ratio<86400>>(time_of_day).count()
        );
    }

    // Returns hours (including partial) since the start of the day
    constexpr double hour_of_day() const {
        // Time point at start of day
        std::chrono::sys_time<std::chrono::days> datepoint
            = std::chrono::floor<std::chrono::days>(timepoint);
        
        // Duration between start of day and timepoint
        std::chrono::milliseconds time_of_day = timepoint - datepoint;

        // Convert to hours as a double
        // (using ratio instead of duration_cast avoids truncation)
        return std::chrono::duration<double, std::ratio<3600>>(time_of_day).count();
    }
};

#endif