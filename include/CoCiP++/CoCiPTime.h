#ifndef COCIPTIME_H
#define COCIPTIME_H

// Struct for storing a datetime: yy-mm-dd h:m:s
// Does not check if datetime is valid
struct CoCiPTime {
    int yy; // year
    int mm; // month
    int dd; // day
    int h; // hour
    int m; // minute
    float s; // seconds

    // Sets the internal variables
    void set(int yy, int mm, int dd, int h, int m, float s) {
        this->yy = yy;
        this->mm = mm;
        this->dd = dd;
        this->h = h;
        this->m = m;
        this->s = s;
    }

    // Returns true if internal year is a leap year
    constexpr bool is_leap() const {
        return ((yy % 4 == 0 && yy % 100 != 0) || (yy % 400 == 0));
    }

    // Returns days (including partial) since the start of the year
    constexpr double day_of_year() const {
        constexpr int DAYS_IN_MONTH[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        constexpr int DAYS_IN_MONTH_LEAP[12] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        
        double doy = 0;
        if (is_leap()) {
            for (int i = 0; i < mm-1; i++) {
                doy += DAYS_IN_MONTH_LEAP[i];
            }
        }
        else {
            for (int i = 0; i < mm-1; i++) {
                doy += DAYS_IN_MONTH[i];
            }
        }
        doy += dd-1; // 1st day of month means 0 full days have passed in month
        doy += (h / 24.) + (m / 1440.) + (s / 86400.);
        return doy;
    }

    // Returns hours (including partial) since the start of the day
    constexpr double hour_of_day() const {
        return h + (m / 60.) + (s / 3600.);
    }
};

#endif