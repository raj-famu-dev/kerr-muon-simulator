#pragma once

namespace monograph {

    // Geometric units: G = c = M = 1. One unit of time is GM/c^3.
    constexpr double kSolarMassSeconds = 4.925490947e-6;  // GM_sun / c^3, in seconds
    constexpr double kSolarMassMeters  = 1476.625;        // GM_sun / c^2, in meters
    constexpr double kMuonLifetimeSec  = 2.2e-6;          // muon proper lifetime

    // Seconds in one geometric time unit for a black hole of the given mass.
    inline double seconds_per_M(double mass_in_solar_masses) {
        return kSolarMassSeconds * mass_in_solar_masses;
    }

    inline double seconds_to_M(double seconds, double mass_in_solar_masses) {
        return seconds / seconds_per_M(mass_in_solar_masses);
    }

    inline double M_to_seconds(double time_in_M, double mass_in_solar_masses) {
        return time_in_M * seconds_per_M(mass_in_solar_masses);
    }

    inline double M_to_meters(double length_in_M, double mass_in_solar_masses) {
        return length_in_M * kSolarMassMeters * mass_in_solar_masses;
    }

    // The muon's 2.2 microsecond lifetime expressed in geometric units.
    inline double muon_lifetime_M(double mass_in_solar_masses) {
        return seconds_to_M(kMuonLifetimeSec, mass_in_solar_masses);
    }

}  // namespace monograph
