// Phase 1 executable: locate the Innermost Stable Circular Orbit (ISCO) using
// the RK4 engine and compare it with Bardeen's closed-form result.
//
// Usage: monograph_isco [--spin 0.9] [--out isco_scan.csv] [--step 0.05]

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <string>

#include "monograph/kerr.hpp"
#include "monograph/orbits.hpp"

using namespace monograph;

static const double kTauMax = 20000.0;
static const double kDtau = 0.05;

static const char* name(Direction d) {
    return d == Direction::Prograde ? "prograde" : "retrograde";
}

int main(int argc, char** argv) {
    double sweep_spin = 0.9;
    double step = 0.05;
    std::string out_path = "isco_scan.csv";

    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--spin") && i + 1 < argc) sweep_spin = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--out") && i + 1 < argc) out_path = argv[++i];
        else if (!std::strcmp(argv[i], "--step") && i + 1 < argc) step = std::atof(argv[++i]);
    }

    // ---- Part 1: numerical versus analytic ISCO for several spins -----------
    std::printf("ISCO search with the RK4 engine (tau_max = %.0f M, dtau = %.2f M)\n\n", kTauMax, kDtau);
    std::printf("  spin a   direction    analytic   numerical      error\n");
    std::printf("  ------   ----------   --------   ---------   ---------\n");
    const double spins[] = {0.0, 0.5, 0.9, 0.998};
    for (double a : spins) {
        for (Direction d : {Direction::Prograde, Direction::Retrograde}) {
            if (a == 0.0 && d == Direction::Retrograde) continue;  // identical at a = 0
            const double an = isco_radius_analytic(a, d);
            const double nu = isco_radius_numerical(a, d, 1e-3, kTauMax, kDtau);
            std::printf("  %6.3f   %-10s   %8.4f   %9.4f   %+9.5f\n", a, name(d), an, nu, nu - an);
        }
    }

    // ---- Part 2: radius sweep, which particles survive? ---------------------
    std::printf("\nRadius sweep for a = %.3f (r_horizon = %.4f M)\n", sweep_spin, horizon_radius(sweep_spin));

    std::ofstream csv(out_path);
    csv << "r,prograde_stable,retrograde_stable\n";

    std::string row_pro, row_ret, ruler;
    const double r_start = horizon_radius(sweep_spin) + 0.02;
    int k = 0;
    for (double r = r_start; r <= 12.0 + 1e-9; r += step, ++k) {
        const bool sp = is_orbit_stable(r, sweep_spin, Direction::Prograde, kTauMax, kDtau);
        const bool sr = is_orbit_stable(r, sweep_spin, Direction::Retrograde, kTauMax, kDtau);
        csv << r << ',' << (sp ? 1 : 0) << ',' << (sr ? 1 : 0) << '\n';
        // Console map: one character per ~0.2 M
        if (k % static_cast<int>(0.2 / step + 0.5) == 0) {
            row_pro += sp ? '#' : '.';
            row_ret += sr ? '#' : '.';
        }
    }
    std::printf("\n  '#' = particle survives, '.' = particle plunges   (left = near the horizon, right = r = 12 M)\n\n");
    std::printf("  prograde   %s\n", row_pro.c_str());
    std::printf("  retrograde %s\n", row_ret.c_str());

    const double pro = isco_radius_numerical(sweep_spin, Direction::Prograde, 1e-3, kTauMax, kDtau);
    const double ret = isco_radius_numerical(sweep_spin, Direction::Retrograde, 1e-3, kTauMax, kDtau);
    std::printf("\n  Prograde ISCO   ~ %.3f M\n  Retrograde ISCO ~ %.3f M\n", pro, ret);
    if (sweep_spin > 0.0)
        std::printf("  A prograde particle survives %.2f M closer to the black hole than a retrograde one.\n", ret - pro);

    std::printf("\nSweep data written to: %s\n", std::filesystem::absolute(out_path).string().c_str());
    return 0;
}
