#pragma once
#include "monograph/state.hpp"

namespace monograph {

    // Outer event horizon: r+ = 1 + sqrt(1 - a^2)
    double horizon_radius(double a);

    // Delta = r^2 - 2r + a^2 (zero at the horizon)
    double delta(double r, double a);

    // Geodesic derivatives with respect to proper time (equatorial plane, M = 1).
    // They depend only on r and the constants E, L, a.
    double dt_dtau(double r, const Constants& c);     // dt/dtau
    double dphi_dtau(double r, const Constants& c);   // dphi/dtau
    double d2r_dtau2(double r, const Constants& c);   // d^2 r / dtau^2

    // The value of (dr/dtau)^2 that E and L demand at radius r.
    double radial_potential(double r, const Constants& c);

    // r_dot^2 - radial_potential(r). Zero for an exact solution, so its drift
    // during a run measures the integrator's numerical error.
    double constraint_error(const State& s, const Constants& c);

}  // namespace monograph
