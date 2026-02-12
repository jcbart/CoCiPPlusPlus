#ifndef COCIP_H
#define COCIP_H

#include <CoCiP++/CoCiPTime.h>

// Forward declarations

struct Params;

// Class for a single CoCiP instance
template <typename MetType>
class CoCiP {
private:
    // Find met variables derived from met inputs
    // Equivalent to pycontrails' calc_timestep_meteorology plus calc_shortwave_radiation and
    // calc_outgoing_longwave_radiation
    // Must only call once per timestep so _old values are correct
    // Ambient and contrail air temperatures (and so RH etc) are conflated
    void update_met_calculations();

    // Calculate contrail properties
    // Only called at end of time step in line with pycontrails meaning variables used in time step
    // are from previous time step even though the meteorology has changed
    void calc_contrail_properties();

    // Calculate radiative properties for contrail
    void calc_radiative_properties();

    // Updates sin_a and cos_a internally where a is angle in degrees between segment and
    // longitudinal axis
    // See pycontrails physics.geo.segment_angle
    void set_heading(const double a);

    // Calculate the contrail evolution across time step dt (s)
    // Advection through the time step should have already occurred giving
    // length_ratio = old_length / new_length
    void calc_timestep_contrail_evolution(const double length_ratio, const double dt_s);

    // Calculate the temporal evolution of the contrail plume parameters
    // length_ratio = old_length/new_length
    // Equivalent to pycontrails' plume_temporal_evolution plus new_contrail_dimensions
    void plume_temporal_evolution(const double length_ratio, const double dt_s);

public:
    // Pointer to meteorology object (derived IMet struct)
    MetType* met;

    // Pointer to params
    Params* params;

    double longitude; // (degrees); only for calculating solar direction radiation and must be driven externally
    double latitude; // (degrees); only for calculating solar direction radiation and must be driven externally
    double altitude; // (m); updated by initial_properties and altitude_after_settling, but should also be driven externally by advection
    CoCiPTime datetime;

    double sin_a; // Sin of angle between contrail and longitude axis; set in evolve
    double cos_a; // Cos of angle between contrail and longitude axis; set in evolve

    // Flight inputs
    double engine_efficiency;
    double ei_h2o;
    double q_fuel;
    double aircraft_mass;
    double wingspan;
    double true_airspeed;
    double fuel_flow;
    double T_exhaust;
    double nvpm_ei_n;

    // Formation
    bool sac = false;
    double T_crit_sac = 0;

    // Contrail properties
    double width = 0; // (m); first set in simulate_wake_vortex_downwash, updated in plume_temporal_evolution
    double depth = 0; // (m); first set in simulate_wake_vortex_downwash, updated in plume_temporal_evolution
    double area_eff = 0; // (m2); first set in calc_contrail_properties
    double sigma_yz = 0; // What is it? Initially zero, updated in plume_temporal_evolution
    double n_ice_per_m = 0; // Ice crystal number per unit length (# m-1); first set in initial_properties, updated in calc_timestep_contrail_evolution
    double n_ice_per_vol = 0; // Ice crystal number per unit volume (# m-3); calculated in calc_contrail_properties
    double iwc = 0; // Ice water content (kg (kg air)-1); first set in initial_properties, updated in calc_timestep_contrail_evolution
    double plume_mass_per_m = 0; // Plume mass per metre (kg m-1); updated in calc_timestep_contrail_evolution
    double r_ice_vol = 0; // Ice particle volume mean radius (m); calculated in calc_contrail_properties
    double tau_contrail = 0; // Contrail optical depth (); calculated in calc_contrail_properties
    double terminal_fall_speed = 0; // Ice particle terminal fall speed (m s-1); calculated in calc_contrail_properties
    double diffuse_h = 0; // Horizontal diffusivity (m2 s-1); calculated in calc_contrail_properties
    double diffuse_v = 0; // Vertical diffusivity (m2 s-1); calculated in calc_contrail_properties
    double dn_dt_agg = 0; // Rate of contrail ice particle losses due to sedimentation-induced aggregation (# s-1); calculated in calc_contrail_properties
    double dn_dt_turb = 0; // Rate of contrail ice particle losses due to plume-internal turbulence (# s-1); calculated in calc_contrail_properties
    double heat_rate = 0; // Heating rate (K s-1); calculated in calc_contrail_properties
    double d_heat_rate = 0; // Differential heating rate (K s-1); calculated in calc_contrail_properties
    double rf_sw = 0; // Shortwave radiative forcing (negative value) (W m-2); calculated in calc_radiative_properties
    double rf_lw = 0; // Longwave radiative forcing (positive value) (W m-2); calculated in calc_radiative_properties
    double rf_net = 0; // Net radiative forcing (rf_lw + rf_sw) (W m-2); calculated in calc_radiative_properties
    double cumul_heat = 0; // Cumulative heat of contrail relative to ambient air due to absorbing radiation (K); updated in calc_timestep_contrail_evolution
    double cumul_differential_heat = 0; // Cumulative differential heat of contrail (K); updated in calc_timestep_contrail_evolution
    // No energy forcing because requires length, can be calculated externally if desired
    bool persistent = false; // True if contrail survives time step; first set in initial_proprties, updated in calc_timestep_contrail_evolution

    // Determine if contrail forms; CoCiP::sac is true if so
    // Must be called externally once
    void formation();

    // Apply wake vortex model to calculate initial contrail width and depth
    // Must be called externally once
    void simulate_wake_vortex_downwash();

    // Calculate the initial contrail properties at the end of the wake vortex phase
    // Equivalent to pycontrails' _find_initial_persistent_contrails
    // Must be called externally once
    void initial_properties();

    // Initialise and calculate properties of a contrail which survives the downwash vortex
    // a is angle in degrees between segment and longitudinal axis (to update heading)
    // Similar to pycontrails' _process_downwash_flight
    // Must be called externally once
    void process_downwash_flight(const double a);

    // Evolve contrail by time step dt (s)
    // Advection through the time step should have already occurred giving
    // length_ratio = old_length / new_length
    // a is angle in degrees between segment and longitudinal axis (to update heading)
    // Similar to the contents of the loop in pycontrails' _simulate_contrail_evolution
    // Must be driven externally each time step
    void evolve(const double length_ratio, const double a, const double dt_s);
};

#endif