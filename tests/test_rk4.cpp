#include <algorithm>
#include <cmath>
#include <cstdio>

#include "monograph/kerr.hpp"
#include "monograph/rk4.hpp"
#include "monograph/state.hpp"

using namespace monograph;

static const double kPi = 3.14159265358979323846;

static int failures = 0;

static void check(bool ok, const char* name) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++failures;
}

// ---- Test 1: harmonic oscillator x'' = -x, exact solution cos(t) ----------
struct Osc {
    double x, v;
};
static Osc operator+(const Osc& a, const Osc& b) { return {a.x + b.x, a.v + b.v}; }
static Osc operator*(double k, const Osc& a) { return {k * a.x, k * a.v}; }

static double oscillator_error(double h) {
    Osc y{1.0, 0.0};
    const double T = 10.0;
    const int n = static_cast<int>(std::lround(T / h));
    for (int i = 0; i < n; ++i) {
        y = rk4_step(y, h, [](const Osc& s) { return Osc{s.v, -s.x}; });
    }
    return std::fabs(y.x - std::cos(T));
}

// ---- Helpers for the geodesic tests ----------------------------------------
static Constants circular_orbit(double r, double a, int sign) {
    const double s = static_cast<double>(sign);
    const double sq = std::sqrt(r);
    const double denom = std::pow(r, 0.75) * std::sqrt(r * sq - 3.0 * sq + 2.0 * s * a);
    Constants c;
    c.a = a;
    c.E = (r * sq - 2.0 * sq + s * a) / denom;
    c.L = s * (r * r - 2.0 * s * a * sq + a * a) / denom;
    return c;
}

// Plunge from r0 with the given constants, stopping just outside the horizon.
// Returns the largest constraint error seen. Also reports final t, tau, r.
struct PlungeResult {
    double max_constraint;
    double t, tau, r;
    bool reached_horizon;
};

static PlungeResult plunge(const Constants& c, double r0, double h) {
    const double r_plus = horizon_radius(c.a);
    State s{0.0, r0, 0.0, -std::sqrt(radial_potential(r0, c))};
    double tau = 0.0, worst = 0.0;
    while (s.r > r_plus + 0.01 && tau < 2000.0) {
        s = advance(s, c, h);
        tau += h;
        worst = std::max(worst, std::fabs(constraint_error(s, c)));
    }
    return {worst, s.t, tau, s.r, s.r <= r_plus + 0.01};
}

int main() {
    // 1. RK4 is 4th order: halving h cuts the error by about 2^4 = 16
    {
        const double e1 = oscillator_error(0.1);
        const double e2 = oscillator_error(0.05);
        const double e3 = oscillator_error(0.025);
        const double ratio1 = e1 / e2;
        const double ratio2 = e2 / e3;
        std::printf("  oscillator errors: %.3e  %.3e  %.3e\n", e1, e2, e3);
        std::printf("  error ratios:      %.2f  %.2f  (expect ~16)\n", ratio1, ratio2);
        check(ratio1 > 14.0 && ratio1 < 18.0, "RK4 order: ratio h -> h/2 is ~16");
        check(ratio2 > 14.0 && ratio2 < 18.0, "RK4 order: ratio h/2 -> h/4 is ~16");
    }

    // 2. A circular Schwarzschild orbit stays circular for a full revolution
    {
        const double r = 10.0;
        const Constants c = circular_orbit(r, 0.0, +1);
        State s{0.0, r, 0.0, 0.0};
        const double dphi_rate = dphi_dtau(r, c);
        const double period = 2.0 * kPi / dphi_rate;
        const double h = 0.01;
        const int n = static_cast<int>(period / h);
        double max_dr = 0.0;
        for (int i = 0; i < n; ++i) {
            s = advance(s, c, h);
            max_dr = std::max(max_dr, std::fabs(s.r - r));
        }
        std::printf("  max |r - r0| over one orbit: %.3e, phi = %.6f\n", max_dr, s.phi);
        check(max_dr < 1e-6, "a=0 circular orbit: radius stays at 10");
        check(std::fabs(s.phi - 2.0 * kPi) < 1e-3, "a=0 circular orbit: phi advances 2 pi");
    }

    // 3. Kerr circular orbits (prograde and retrograde) hold their radius
    for (int sign = 1; sign >= -1; sign -= 2) {
        const double a = 0.9, r = 5.0;
        const Constants c = circular_orbit(r, a, sign);
        State s{0.0, r, 0.0, 0.0};
        double max_dr = 0.0, tau = 0.0;
        const double h = 0.01;
        for (int i = 0; i < 20000; ++i) {
            s = advance(s, c, h);
            tau += h;
            max_dr = std::max(max_dr, std::fabs(s.r - r));
        }
        std::printf("  Kerr %s: max |r - r0| = %.3e after tau = %.0f\n",
                    sign > 0 ? "prograde" : "retrograde", max_dr, tau);
        check(max_dr < 1e-6, sign > 0 ? "Kerr prograde circular orbit holds radius"
                                      : "Kerr retrograde circular orbit holds radius");
    }

    // 4. Plunges: constraint drift stays tiny all the way to the horizon
    {
        const Constants schw{0.0, 0.98, 3.0};
        const PlungeResult p = plunge(schw, 15.0, 0.005);
        std::printf("  a=0 plunge: drift %.3e, tau %.2f, t %.2f, r %.4f\n",
                    p.max_constraint, p.tau, p.t, p.r);
        check(p.reached_horizon, "a=0 plunge reaches the horizon");
        check(p.max_constraint < 1e-8, "a=0 plunge: constraint drift < 1e-8");

        const Constants kerr{0.9, 0.95, 2.0};
        const PlungeResult q = plunge(kerr, 15.0, 0.005);
        std::printf("  a=0.9 plunge: drift %.3e, tau %.2f, t %.2f, r %.4f\n",
                    q.max_constraint, q.tau, q.t, q.r);
        check(q.reached_horizon, "a=0.9 plunge reaches the horizon");
        check(q.max_constraint < 1e-8, "a=0.9 plunge: constraint drift < 1e-8");
        check(q.t > q.tau, "a=0.9 plunge: coordinate time exceeds proper time");
    }

    std::printf("\n%s\n", failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED");
    return failures == 0 ? 0 : 1;
}