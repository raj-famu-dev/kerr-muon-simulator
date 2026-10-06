
#pragma once

namespace monograph {

    // Dynamic state of the particle. Proper time tau is NOT stored here:
    // it is the independent variable (the integration clock).
    struct State {
        double t;      // coordinate time (distant observer's clock), units of M
        double r;      // Boyer-Lindquist radius, units of M
        double phi;    // azimuthal angle, radians
        double r_dot;  // dr/dtau
    };

    // Black hole spin and the two conserved quantities of the orbit.
    // Computed once at launch, then never changed.
    struct Constants {
        double a;  // spin parameter, 0 <= a < 1
        double E;  // specific energy
        double L;  // specific angular momentum (negative = retrograde)
    };

}  // namespace monograph