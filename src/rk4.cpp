#include "monograph/rk4.hpp"

#include <algorithm>
#include <cmath>

#include "monograph/kerr.hpp"

namespace monograph {

    State geodesic_derivative(const State& s, const Constants& c) {
        return State{
            dt_dtau(s.r, c),    // dt/dtau
            s.r_dot,            // dr/dtau
            dphi_dtau(s.r, c),  // dphi/dtau
            d2r_dtau2(s.r, c)   // d(r_dot)/dtau
        };
    }

    State advance(const State& s, const Constants& c, double dtau) {
        return rk4_step(s, dtau, [&c](const State& y) { return geodesic_derivative(y, c); });
    }

    double adaptive_dtau(const State& s, double a, const IntegratorOptions& o) {
        const double gap = s.r - horizon_radius(a);
        const double speed = std::fabs(s.r_dot);
        double h = o.dtau_max;
        if (speed > 1e-12) h = std::min(h, o.eta * gap / speed);
        return std::max(h, o.dtau_min);
    }

    IntegrationResult integrate(const State& s0, const Constants& c,
                                const IntegratorOptions& o, const Observer& observer) {
        State s = s0;
        double tau = 0.0;
        long steps = 0;
        double worst = std::fabs(constraint_error(s, c));
        const double r_stop = horizon_radius(c.a) + o.eps;
        bool reached = false;

        if (observer) observer(tau, s);

        while (tau < o.tau_max) {
            if (s.r <= r_stop) { reached = true; break; }
            if (s.r > o.r_escape) break;

            double h = adaptive_dtau(s, c.a, o);
            if (o.land_on_tau > tau && tau + h > o.land_on_tau) h = o.land_on_tau - tau;

            s = advance(s, c, h);
            tau = (o.land_on_tau > 0.0 && std::fabs(tau + h - o.land_on_tau) < 1e-9)
                      ? o.land_on_tau : tau + h;
            ++steps;
            worst = std::max(worst, std::fabs(constraint_error(s, c)));
            if (observer) observer(tau, s);
        }
        if (s.r <= r_stop) reached = true;
        return IntegrationResult{s, tau, steps, reached, worst};
    }

}  // namespace monograph