#include <cmath>
#include <cstdio>
#include <initializer_list>

#include "monograph/kerr.hpp"
#include "monograph/orbits.hpp"

using namespace monograph;

static int failures = 0;

static void check(bool ok, const char* name) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++failures;
}

static bool near(double x, double y, double tol) { return std::fabs(x - y) <= tol; }

int main() {
    const Direction pro = Direction::Prograde, ret = Direction::Retrograde;

    // 1. Bardeen's closed form reproduces the known targets
    check(near(isco_radius_analytic(0.0, pro), 6.0, 1e-12), "analytic ISCO a=0 is 6 M");
    check(near(isco_radius_analytic(0.9, pro), 2.3209, 1e-4), "analytic prograde ISCO a=0.9 is 2.32 M");
    check(near(isco_radius_analytic(0.9, ret), 8.7174, 1e-4), "analytic retrograde ISCO a=0.9 is 8.72 M");
    check(near(isco_radius_analytic(0.998, pro), 1.2370, 1e-4), "analytic prograde ISCO a=0.998 is 1.24 M");
    check(near(isco_radius_analytic(0.998, ret), 8.9944, 1e-4), "analytic retrograde ISCO a=0.998 is 8.99 M");
    check(near(isco_radius_analytic(1.0, pro), 1.0, 1e-9), "analytic prograde ISCO a=1 is 1 M");
    check(near(isco_radius_analytic(1.0, ret), 9.0, 1e-9), "analytic retrograde ISCO a=1 is 9 M");

    // 2. Circular-orbit constants at the Schwarzschild ISCO: E = sqrt(8/9), L = sqrt(12)
    {
        const Constants c = circular_orbit(6.0, 0.0, pro);
        check(near(c.E, std::sqrt(8.0 / 9.0), 1e-12), "a=0 ISCO: E = sqrt(8/9)");
        check(near(c.L, std::sqrt(12.0), 1e-12), "a=0 ISCO: L = sqrt(12)");
    }

    // 3. Existence of circular orbits (photon orbit at r = 3 for a = 0)
    check(!circular_orbit_exists(2.9, 0.0, pro), "no circular orbit inside r = 3 (a=0)");
    check(circular_orbit_exists(3.1, 0.0, pro), "circular orbit exists at r = 3.1 (a=0)");

    // 4. Clear-cut stability: well outside versus well inside the ISCO
    for (double a : {0.0, 0.9}) {
        for (Direction d : {pro, ret}) {
            const double isco = isco_radius_analytic(a, d);
            char n1[96], n2[96];
            std::snprintf(n1, sizeof n1, "a=%.1f %s: stable at ISCO + 0.5", a, d == pro ? "prograde" : "retrograde");
            std::snprintf(n2, sizeof n2, "a=%.1f %s: unstable at ISCO - 0.5", a, d == pro ? "prograde" : "retrograde");
            check(is_orbit_stable(isco + 0.5, a, d), n1);
            if (isco - 0.5 > horizon_radius(a)) check(!is_orbit_stable(isco - 0.5, a, d), n2);
        }
    }

    // 5. The engine finds the ISCO numerically, matching Bardeen
    struct Case { double a; Direction d; const char* label; };
    const Case cases[] = {
        {0.0, pro, "numerical ISCO a=0"},
        {0.9, pro, "numerical prograde ISCO a=0.9"},
        {0.9, ret, "numerical retrograde ISCO a=0.9"},
        {0.998, pro, "numerical prograde ISCO a=0.998"},
        {0.998, ret, "numerical retrograde ISCO a=0.998"},
    };
    for (const Case& c : cases) {
        const double an = isco_radius_analytic(c.a, c.d);
        const double nu = isco_radius_numerical(c.a, c.d);
        std::printf("  %-34s analytic %.4f  numerical %.4f\n", c.label, an, nu);
        check(near(nu, an, 5e-3), c.label);
    }

    // 6. The headline result: prograde survives closer than retrograde
    {
        const double p = isco_radius_numerical(0.9, pro);
        const double r = isco_radius_numerical(0.9, ret);
        check(p < r - 5.0, "a=0.9: prograde ISCO is more than 5 M inside the retrograde one");
    }

    std::printf("\n%s\n", failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED");
    return failures == 0 ? 0 : 1;
}
