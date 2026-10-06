#include <cmath>
#include <cstdio>

#include "monograph/kerr.hpp"
#include "monograph/state.hpp"
#include "monograph/units.hpp"

using namespace monograph;

static int failures = 0;

static void check(bool ok, const char* name) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++failures;
}

static bool close(double x, double y, double tol = 1e-9) {
    return std::fabs(x - y) <= tol * (1.0 + std::fabs(y));
}

// Bardeen circular-orbit constants. sign = +1 prograde, -1 retrograde.
static Constants circular_orbit(double r, double a, int sign) {
    const double s = static_cast<double>(sign);
    const double sq = std::sqrt(r);
    const double denom = std::pow(r, 0.75) *
                         std::sqrt(r * sq - 3.0 * sq + 2.0 * s * a);
    Constants c;
    c.a = a;
    c.E = (r * sq - 2.0 * sq + s * a) / denom;
    c.L = s * (r * r - 2.0 * s * a * sq + a * a) / denom;
    return c;
}

int main() {
    // 1. Horizon radius
    check(close(horizon_radius(0.0), 2.0), "horizon a=0 is 2");
    check(close(horizon_radius(1.0), 1.0), "horizon a=1 is 1");
    check(close(horizon_radius(0.9), 1.0 + std::sqrt(0.19)), "horizon a=0.9");

    // 2. Schwarzschild (a = 0) closed forms at an arbitrary point
    {
        const double r = 10.0;
        const Constants c{0.0, 0.95, 4.0};
        check(close(dt_dtau(r, c), c.E / (1.0 - 2.0 / r)), "a=0: dt/dtau = E/(1-2/r)");
        check(close(dphi_dtau(r, c), c.L / (r * r)), "a=0: dphi/dtau = L/r^2");
        const double expected = -1.0 / (r * r) + c.L * c.L / (r * r * r)
                              - 3.0 * c.L * c.L / (r * r * r * r);
        check(close(d2r_dtau2(r, c), expected), "a=0: d2r/dtau2 matches Schwarzschild");
    }

    // 3. Schwarzschild circular orbit at r = 10
    {
        const double r = 10.0;
        Constants c;
        c.a = 0.0;
        c.L = std::sqrt(r * r / (r - 3.0));
        c.E = (1.0 - 2.0 / r) / std::sqrt(1.0 - 3.0 / r);
        check(std::fabs(radial_potential(r, c)) < 1e-12, "a=0 circular: (dr/dtau)^2 = 0");
        check(std::fabs(d2r_dtau2(r, c)) < 1e-12, "a=0 circular: acceleration = 0");
    }

    // 4. Kerr circular orbits, prograde and retrograde (independent formulas)
    for (int sign = 1; sign >= -1; sign -= 2) {
        const double a = 0.9, r = 5.0;
        const Constants c = circular_orbit(r, a, sign);
        const char* tag = sign > 0 ? "prograde" : "retrograde";
        char name[96];

        std::snprintf(name, sizeof name, "Kerr %s circular: (dr/dtau)^2 = 0", tag);
        check(std::fabs(radial_potential(r, c)) < 1e-12, name);

        std::snprintf(name, sizeof name, "Kerr %s circular: acceleration = 0", tag);
        check(std::fabs(d2r_dtau2(r, c)) < 1e-12, name);

        // Orbital angular velocity: Omega = sign / (r^1.5 + sign*a)
        const double omega = dphi_dtau(r, c) / dt_dtau(r, c);
        const double expected = sign / (std::pow(r, 1.5) + sign * a);
        std::snprintf(name, sizeof name, "Kerr %s circular: Omega = dphi/dt", tag);
        check(close(omega, expected), name);
    }

    // 5. Units
    check(close(seconds_per_M(10.0), 4.925490947e-5), "10 solar masses: 49.25 microseconds per M");
    check(std::fabs(muon_lifetime_M(10.0) - 0.04466) < 1e-4, "muon lifetime ~0.0447 M at 10 Msun");

    std::printf("\n%s\n", failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED");
    return failures == 0 ? 0 : 1;
}
