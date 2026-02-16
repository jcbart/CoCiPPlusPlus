#include <fstream>
#include <yaml-cpp/yaml.h>
#include <CoCiP++/CoCiP.h>
#include <CoCiP++/met.h>
#include <CoCiP++/params.h>
#include "CoCiPLog.h"

struct Config {
    double duration_m;
    double dt_s;
    double save_interval_m;

    // Time and location
    int year;
    int month;
    int day;
    int hour;
    int minute;
    float seconds;
    double longitude;
    double latitude;

    // Inputs to met
    double z0;
    double P0;
    double T0;
    double lapse_rate;
    double rh_i1;
    double rh_i0;
    double D1;
    double DT;
    double ds_dz;
    double tnsr;
    double olr;

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

    void readYAML(std::string& configPath) {
        YAML::Node node = YAML::LoadFile(configPath);

        duration_m = node["duration_m"].as<double>();
        dt_s = node["dt_s"].as<double>();
        save_interval_m = node["save_interval_m"].as<double>();
        year = node["year"].as<int>();
        month = node["month"].as<int>();
        day = node["day"].as<int>();
        hour = node["hour"].as<int>();
        minute = node["minute"].as<int>();
        seconds = node["seconds"].as<float>();
        longitude = node["longitude"].as<double>();
        latitude = node["latitude"].as<double>();
        z0 = node["z0"].as<double>();
        P0 = node["P0"].as<double>();
        T0 = node["T0"].as<double>();
        lapse_rate = node["lapse_rate"].as<double>();
        rh_i1 = node["rh_i1"].as<double>();
        rh_i0 = node["rh_i0"].as<double>();
        D1 = node["D1"].as<double>();
        DT = node["DT"].as<double>();
        ds_dz = node["ds_dz"].as<double>();
        tnsr = node["tnsr"].as<double>();
        olr = node["olr"].as<double>();
        engine_efficiency = node["engine_efficiency"].as<double>();
        ei_h2o = node["ei_h2o"].as<double>();
        q_fuel = node["q_fuel"].as<double>();
        aircraft_mass = node["aircraft_mass"].as<double>();
        wingspan = node["wingspan"].as<double>();
        true_airspeed = node["true_airspeed"].as<double>();
        fuel_flow = node["fuel_flow"].as<double>();
        T_exhaust = node["T_exhaust"].as<double>();
        nvpm_ei_n = node["nvpm_ei_n"].as<double>();
    }
};

template <typename MetType>
void write_to_csv(CoCiP<MetType>& cocip, double time_elapsed_m, bool first_write) {
    std::ofstream file;

    if (first_write) {
        file.open("cocip.out");
        file << "Time (m), Altitude (m), Width (m), Depth (m), sigma_yz (m2), "
             << "n_ice_per_m3 (# m-3), n_ice_per_m (# m-1), IWC (kg kg-1), IWC (kg m-1), "
             << "r_ice_vol (m), Optical depth (), Cumulative heat (K), Air temperature (K), "
             << "RHi (), RF SW (W m-2), RF LW (W m-2), RF net (W m-2)\n";
    }
    else {
        file.open("cocip.out", std::ios::app);
    }
    file << time_elapsed_m << ", " << cocip.altitude << ", " << cocip.width << ", " << cocip.depth
         << ", " << cocip.sigma_yz << ", " << cocip.n_ice_per_vol << ", " << cocip.n_ice_per_m
         << ", " << cocip.iwc << ", " << cocip.iwc * cocip.plume_mass_per_m << ", "
         << cocip.r_ice_vol << ", " << cocip.tau_contrail << ", " << cocip.cumul_heat << ", "
         << cocip.met->air_temperature << ", " << cocip.met->rh_i << ", " << cocip.rf_sw << ", "
         << cocip.rf_lw << ", " << cocip.rf_net
         << "\n";
    
    file.close();
}

// CoCiP standalone program
int main(int argc, char* argv[]) {

    // Read config file path
    std::string configPath = "CoCiP-config.yaml";
    if (argc > 1) {
        configPath = argv[1];
    }

    // Read config containing simulation parameters and flight and met inputs
    Config config;
    config.readYAML(configPath);

    CoCiP<SimpleMet> cocip;

    // Give pointer to Params object and read YAML
    cocip.params = new Params;
    cocip.params->readYAML();

    // Initialise cocip.met
    cocip.met = new SimpleMet(config.z0, config.P0, config.T0, config.lapse_rate, config.rh_i1,
        config.rh_i0, config.D1, config.DT, config.ds_dz);
    
    cocip.met->tnsr = config.tnsr;
    cocip.met->olr = config.olr;

    // Give datetime and location
    cocip.datetime.set(
        config.year,
        config.month,
        config.day,
        config.hour,
        config.minute,
        config.seconds
    );
    cocip.altitude = config.z0;
    cocip.longitude = config.longitude;
    cocip.latitude = config.latitude;

    // Give flight inputs
    cocip.engine_efficiency = config.engine_efficiency;
    cocip.ei_h2o = config.ei_h2o;
    cocip.q_fuel = config.q_fuel;
    cocip.aircraft_mass = config.aircraft_mass;
    cocip.wingspan = config.wingspan;
    cocip.true_airspeed = config.true_airspeed;
    cocip.fuel_flow = config.fuel_flow;
    cocip.T_exhaust = config.T_exhaust;
    cocip.nvpm_ei_n = config.nvpm_ei_n;
    
    cocip.formation();

    if (!cocip.sac) {
        CoCiP_LogWrite("Finished: No contrail formed");
        exit(EXIT_SUCCESS);
    }
    
    cocip.simulate_wake_vortex_downwash();
    
    cocip.initial_properties();

    if (!cocip.persistent) {
        CoCiP_LogWrite("Finished: Not initially persistent");
        exit(EXIT_SUCCESS);
    }

    cocip.process_downwash_flight(0);
    
    write_to_csv(cocip, 0, true);

    double time_elapsed_m = 0;
    double next_save_m = config.save_interval_m;
    
    while (time_elapsed_m < config.duration_m) {
        time_elapsed_m += (config.dt_s / 60);

        cocip.evolve(1, 0, config.dt_s);

        if (time_elapsed_m >= next_save_m) {
            write_to_csv(cocip, time_elapsed_m, false);
            next_save_m += config.save_interval_m;
        }

        if (!cocip.persistent) {
            CoCiP_LogWrite("Finished: Not persistent at " + std::to_string(time_elapsed_m) + " m");
            exit(EXIT_SUCCESS);
        }
    }

    return 0;
}