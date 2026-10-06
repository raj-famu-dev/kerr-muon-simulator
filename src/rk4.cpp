#include "monograph/rk4.hpp"

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

}  // namespace monograph
