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

// Everything needed to start a plunging particle.
struct Launch {
    Constants c;  // a, E, L (L already reduced)
    State s0;     // t = 0, r = r0, phi = 0, r_dot < 0 (moving inward)
};

// Start on the circular orbit at r0, cut the angular momentum by l_fraction
// (0.3 means |L| is reduced by 30%), and give the particle the inward radial
// speed that the conserved E and L demand at r0. A smaller L means a weaker
// centrifugal barrier, so the orbit decays. Returns false if no circular orbit
// exists at r0 or the resulting state is not allowed.
bool make_plunge_launch(double r0, double a, Direction dir, double l_fraction, Launch& out);

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