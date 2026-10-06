#pragma once
#include "monograph/state.hpp"

namespace monograph {

    // Direction of an orbit relative to the black hole's spin.
    enum class Direction { Prograde = +1, Retrograde = -1 };

    // A circular orbit exists at radius r only if
    //   r^(3/2) - 3 r^(1/2) + 2 s a > 0     (s = +1 prograde, -1 retrograde)
    // Below the photon orbit no circular orbit of a massive particle exists.
    bool circular_orbit_exists(double r, double a, Direction dir);

    // Specific energy E and angular momentum L of the circular orbit at radius r.
    // Retrograde orbits have negative L. Only valid if circular_orbit_exists().
    Constants circular_orbit(double r, double a, Direction dir);

    // ISCO radius from Bardeen's closed-form expression (the ground truth).
    double isco_radius_analytic(double a, Direction dir);

    // Numerical stability test using the RK4 engine.
    // Launches the particle on the circular orbit at radius r, gives it a tiny
    // inward nudge, and integrates for tau_max. Returns true if it stays near r
    // (stable), false if it drifts away or no circular orbit exists (unstable).
    bool is_orbit_stable(double r, double a, Direction dir,
                         double tau_max = 20000.0, double dtau = 0.05);

    // Bisection on is_orbit_stable() to locate the ISCO numerically.
    double isco_radius_numerical(double a, Direction dir,
                                 double tolerance = 1e-3,
                                 double tau_max = 20000.0, double dtau = 0.05);

}  // namespace monograph