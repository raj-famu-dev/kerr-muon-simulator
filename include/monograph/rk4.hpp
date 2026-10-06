#pragma once
#include <functional>

#include "monograph/state.hpp"

namespace monograph {

// --- Vector-space operations on State so RK4 can combine slopes -------------
// A State doubles as a "slope": its fields hold dt/dtau, dr/dtau, dphi/dtau,
// d(r_dot)/dtau.
inline State operator+(const State& x, const State& y) {
    return State{x.t + y.t, x.r + y.r, x.phi + y.phi, x.r_dot + y.r_dot};
}

inline State operator*(double k, const State& x) {
    return State{k * x.t, k * x.r, k * x.phi, k * x.r_dot};
}

// --- Generic RK4 step --------------------------------------------------------
// Works for any type S that supports (S + S) and (double * S), for example the
// State struct or a tiny harmonic-oscillator struct in a unit test.
// f(y) must return dy/dx. The system is autonomous (no explicit dependence on
// the independent variable), which is true of the geodesic equations in tau.
template <class S, class F>
S rk4_step(const S& y, double h, F&& f) {
    const S k1 = f(y);
    const S k2 = f(y + (0.5 * h) * k1);
    const S k3 = f(y + (0.5 * h) * k2);
    const S k4 = f(y + h * k3);
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// --- Engine-specific wrappers (defined in rk4.cpp) ---------------------------

// The four rates of change of the state, from the Kerr geodesic equations.
State geodesic_derivative(const State& s, const Constants& c);

// Advance the particle by one tick dtau of its own proper time.
State advance(const State& s, const Constants& c, double dtau);

// --- Adaptive integration loop -----------------------------------------------

struct IntegratorOptions {
    double dtau_max = 2e-3;      // largest step in proper time (units of M)
    double dtau_min = 1e-12;     // smallest allowed step
    double eta = 0.01;           // near the horizon: step <= eta * (r - r+) / |dr/dtau|
    double eps = 1e-6;           // stop when r <= r+ + eps
    double tau_max = 5000.0;     // give up after this much proper time
    double r_escape = 1e3;       // stop if the particle flies out this far
    double land_on_tau = -1.0;   // if > 0, one step ends exactly at this tau
};

struct IntegrationResult {
    State final_state;
    double tau;                  // proper time elapsed
    long steps;
    bool reached_horizon;        // true if r <= r+ + eps was reached
    double max_constraint_error; // largest |r_dot^2 - R(r)| seen
};

// Called with the initial state and again after every step.
using Observer = std::function<void(double tau, const State& s)>;

// The step shrinks in proportion to the remaining distance to the horizon, so
// the integrator can follow dt/dtau -> infinity without losing accuracy.
double adaptive_dtau(const State& s, double a, const IntegratorOptions& o);

IntegrationResult integrate(const State& s0, const Constants& c,
                            const IntegratorOptions& o, const Observer& observer = nullptr);

}  // namespace monograph