#ifndef MET_H
#define MET_H

// Forward declaration
struct Params;
struct CoCiPTime;

// Virtual struct to hold meteorological inputs to CoCiP and variables derived from them
// Is derived in child types depending on met calculation method
struct IMet {

    // Single values from met input
    double tnsr; // TOA net downward shortwave radiation (SDR - RSR) (W m-2)
    double olr; // TOA outgoing (upward) longwave radiation (W m-2)

    // Values taken from get_local_values
    double air_pressure; // Pressure (Pa)
    double air_temperature; // Air temperature (K)
    double specific_humidity; // Specific humidity of water vapor (kg (kg moist air)-1)
    double u_wind; // Eastward wind (m s-1)
    double v_wind; // Northward wind (m s-1)
    double ciwc; // Cloud ice water mixing ratio (kg (kg dry air)-1)

    double air_pressure_old; // Pressure at last time step (Pa)
    double air_temperature_old; // Temperature inside contrail at last time step (K)
    double specific_humidity_old; // Specific humidity of water vapor at last time step (kg (kg moist air)-1)

    double air_pressure_lower; // Pressure (Pa) at altitude - dz_m (i.e. grid cell below)
    double air_temperature_lower; // Temperature (K) at grid cell below
    double u_wind_lower; // Eastward wind (m s-1) at grid cell below
    double v_wind_lower; // Northward wind (m s-1) at grid cell below

    double effective_vertical_resolution; // Effective vertical resolution of met data (m)
    double dz_m; // Difference in altitude between centre of contrail grid cell and grid cell below (m)

    // Values calculated
    double dtheta_dz = 0; // Potential temperature gradient (K m-1) (dT_dz in pycontrails)
    double ds_dz = 0; // Wind shear (m s-1 m-1)
    double dsn_dz = 0; // Wind shear normal (m s-1 m-1)
    double rho_air = 0; // Dry air mass (kg m-3)
    double rh_i = 0; // Relative humidity w.r.t. ice as a decimal (not percentage)
    double tau_cirrus = 0; // Optical depth of cirrus above the contrail
    double sdr = 0; // Solar direct radiation (W m-2); from met or calculated in calc_shortwave_radiation
    double rsr = 0; // Reflected solar radiation (W m-2); calculated in calc_shortwave_radiation

    virtual ~IMet() = default;

    void calc_variables(double altitude, double cumul_heat, double depth, double cos_a, double sin_a,
        double longitude, double latitude, CoCiPTime& datetime, const Params* params);

private:
    // Virtual method called by calc_variables to update the local variables using type-specific
    // method
    virtual void get_local_values(double altitude) = 0;
    
    // Virtual method to calculate cirrus optical depth at contrail level ()
    virtual double calc_tau_cirrus() const = 0;
};

// Array met type: takes 1D (vertical) arrays at the contrail's position (lon, lat)
template <typename arrayType>
struct ArrayMet : public IMet {

    // Vertical arrays from met input
    int vsize = 0; // Unstaggered length of each vertical array
    arrayType* T_POT = nullptr; // Potential air temperature (K)
    arrayType* P = nullptr; // Air pressure (Pa)
    arrayType* QV = nullptr; // Specific humidity of water vapor (kg (kg moist air)-1)
    arrayType* U = nullptr; // Eastward wind (m s-1)
    arrayType* V = nullptr; // Northward wind (m s-1)
    arrayType* CIWC = nullptr; // Cloud ice water mixing ratio (kg (kg dry air)-1)
    arrayType* Z = nullptr; // Altitude at grid cell centres (m)
    arrayType* Z_AT_W = nullptr; // Altitude at grid cell boundaries (m); array len is 1 larger

    int k;
    int k_lower;

    // Initialise Met object with unstaggered vertical length
    ArrayMet(int vsize) : vsize(vsize) { }

private:
    void get_local_values(double altitude) override;

    // Equivalent to pycontrails' tau_cirrus
    double calc_tau_cirrus() const override;

    // Find the vertical index k of the grid cell in which the contrail is located
    // len(Z_AT_W) = klen + 1
    // Will fail if altitude is not in Z_AT_W range
    int find_k_index(double altitude) const;
};

// Array met type: uses simple formulae to calculate local variables
struct SimpleMet : public IMet {

    double z0;
    double P0;
    double T0;
    double lapse_rate;
    double rh_i1;
    double rh_i0;
    double D1;
    double DT;
    double ds_dz_cross_track;

    SimpleMet(double z0, double P0, double T0, double lapse_rate, double rh_i1, double rh_i0,
        double D1, double DT, double ds_dz_cross_track);

private:
    void get_local_values(double altitude) override;

    double calc_tau_cirrus() const override;
};

#endif