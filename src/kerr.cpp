
#include "monograph/kerr.hpp"

#include <cmath>

namespace monograph {

    double horizon_radius(double a) {
        return 1.0 + std::sqrt(1.0 - a * a);
    }

    double delta(double r, double a) {
        return r * r - 2.0 * r + a * a;
    }

    double dt_dtau(double r, const Constants& c) {
        const double a = c.a;
        const double numerator = (r * r + a * a + 2.0 * a * a / r) * c.E
                               - (2.0 * a / r) * c.L;
        return numerator / delta(r, a);
    }

    double dphi_dtau(double r, const Constants& c) {
        const double numerator = (1.0 - 2.0 / r) * c.L + (2.0 * c.a / r) * c.E;
        return numerator / delta(r, c.a);
    }

    double d2r_dtau2(double r, const Constants& c) {
        const double K = c.a * c.a * (c.E * c.E - 1.0) - c.L * c.L;
        const double J = (c.L - c.a * c.E) * (c.L - c.a * c.E);
        return -1.0 / (r * r) - K / (r * r * r) - 3.0 * J / (r * r * r * r);
    }

    double radial_potential(double r, const Constants& c) {
        const double K = c.a * c.a * (c.E * c.E - 1.0) - c.L * c.L;
        const double J = (c.L - c.a * c.E) * (c.L - c.a * c.E);
        return c.E * c.E - 1.0 + 2.0 / r + K / (r * r) + 2.0 * J / (r * r * r);
    }

    double constraint_error(const State& s, const Constants& c) {
        return s.r_dot * s.r_dot - radial_potential(s.r, c);
    }

}  // namespace monograph