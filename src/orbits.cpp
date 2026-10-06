#include "monograph/orbits.hpp"

#include <cmath>

#include "monograph/kerr.hpp"
#include "monograph/rk4.hpp"

namespace monograph {

static double sgn(Direction d) { return d == Direction::Prograde ? 1.0 : -1.0; }

bool circular_orbit_exists(double r, double a, Direction dir) {
    const double sq = std::sqrt(r);
    return r > horizon_radius(a) && (r * sq - 3.0 * sq + 2.0 * sgn(dir) * a) > 0.0;
}

Constants circular_orbit(double r, double a, Direction dir) {
    const double s = sgn(dir);
    const double sq = std::sqrt(r);
    const double denom = std::pow(r, 0.75) * std::sqrt(r * sq - 3.0 * sq + 2.0 * s * a);
    Constants c;
    c.a = a;
    c.E = (r * sq - 2.0 * sq + s * a) / denom;
    c.L = s * (r * r - 2.0 * s * a * sq + a * a) / denom;
    return c;
}

double isco_radius_analytic(double a, Direction dir) {
    const double z1 = 1.0 + std::cbrt(1.0 - a * a) * (std::cbrt(1.0 + a) + std::cbrt(1.0 - a));
    const double z2 = std::sqrt(3.0 * a * a + z1 * z1);
    const double root = std::sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2));
    return dir == Direction::Prograde ? 3.0 + z2 - root : 3.0 + z2 + root;
}

bool is_orbit_stable(double r, double a, Direction dir, double tau_max, double dtau) {
    if (!circular_orbit_exists(r, a, dir)) return false;

    const Constants c = circular_orbit(r, a, dir);
    // Tiny inward nudge. If the orbit is unstable the deviation grows
    // exponentially; if stable it only oscillates with amplitude ~ nudge/omega.
    State s{0.0, r, 0.0, -1e-6};

    const double r_plus = horizon_radius(a);
    const double lower = r - 0.02;   // drifted inward: plunging
    const double upper = r + 0.5;    // drifted outward: escaping
    const long steps = static_cast<long>(tau_max / dtau);

    for (long i = 0; i < steps; ++i) {
        s = advance(s, c, dtau);
        if (s.r < lower || s.r > upper || s.r < r_plus + 0.01) return false;
    }
    return true;
}

double isco_radius_numerical(double a, Direction dir, double tolerance,
                             double tau_max, double dtau) {
    double lo = horizon_radius(a) + 1e-3;  // unstable (or no circular orbit)
    double hi = 12.0;                      // comfortably stable
    while (hi - lo > tolerance) {
        const double mid = 0.5 * (lo + hi);
        if (is_orbit_stable(mid, a, dir, tau_max, dtau)) hi = mid; else lo = mid;
    }
    return 0.5 * (lo + hi);
}

}  // namespace monograph
