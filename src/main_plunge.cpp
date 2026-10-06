// Phase 2 executable: launch a muon on a decaying orbit around a Kerr black
// hole and record every frame of its plunge to a CSV file.
//
// Usage: monograph_plunge [--spin 0.9] [--direction prograde|retrograde]
//                         [--r0 6] [--lfrac 0.3] [--mass 0.1]
//                         [--dtau 0.002] [--eps 1e-6] [--taumax 5000]
//                         [--out plunge.csv]

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <string>

#include "monograph/csv_writer.hpp"
#include "monograph/kerr.hpp"
#include "monograph/orbits.hpp"
#include "monograph/rk4.hpp"
#include "monograph/units.hpp"

using namespace monograph;

static const double kPi = 3.14159265358979323846;

int main(int argc, char** argv) {
    double a = 0.9, r0 = 6.0, lfrac = 0.3, mass = 0.1;
    Direction dir = Direction::Prograde;
    IntegratorOptions opt;
    std::string out_path = "plunge.csv";

    for (int i = 1; i < argc; ++i) {
        const bool more = i + 1 < argc;
        if (!std::strcmp(argv[i], "--spin") && more) a = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--r0") && more) r0 = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--lfrac") && more) lfrac = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--mass") && more) mass = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--dtau") && more) opt.dtau_max = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--eps") && more) opt.eps = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--taumax") && more) opt.tau_max = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--out") && more) out_path = argv[++i];
        else if (!std::strcmp(argv[i], "--direction") && more) {
            const std::string d = argv[++i];
            if (d == "prograde") dir = Direction::Prograde;
            else if (d == "retrograde") dir = Direction::Retrograde;
            else { std::fprintf(stderr, "--direction must be prograde or retrograde\n"); return 1; }
        } else { std::fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]); return 1; }
    }

    if (a < 0.0 || a >= 1.0 || mass <= 0.0) {
        std::fprintf(stderr, "Need 0 <= spin < 1 and mass > 0.\n");
        return 1;
    }

    Launch launch;
    if (!make_plunge_launch(r0, a, dir, lfrac, launch)) {
        std::fprintf(stderr,
            "Cannot launch from r0 = %.3f M with lfrac = %.2f: no circular orbit exists there, "
            "or the reduced angular momentum gives a forbidden state. Try a larger --r0.\n", r0, lfrac);
        return 1;
    }
    const Constants& c = launch.c;
    const double r_plus = horizon_radius(a);
    const double tau_decay = muon_lifetime_M(mass);
    const double us_per_M = seconds_per_M(mass) * 1e6;
    opt.land_on_tau = tau_decay;  // so one frame lands exactly on the decay time

    CsvWriter csv(out_path, {"tau", "t", "r", "phi", "r_dot", "x", "y", "dt_dtau",
                             "constraint_error", "muon_alive", "tau_us", "t_us"});
    if (!csv.is_open()) {
        std::fprintf(stderr, "Cannot open %s for writing.\n", out_path.c_str());
        return 1;
    }

    // Record the first frame at or after the decay time.
    bool decayed = false;
    State decay_state{};
    double decay_tau = 0.0;

    auto observer = [&](double tau, const State& s) {
        const bool alive = tau < tau_decay - 1e-12;
        if (!alive && !decayed) { decayed = true; decay_state = s; decay_tau = tau; }
        csv.write_row({tau, s.t, s.r, s.phi, s.r_dot, s.r * std::cos(s.phi), s.r * std::sin(s.phi),
                       dt_dtau(s.r, c), constraint_error(s, c), alive ? 1.0 : 0.0,
                       tau * us_per_M, s.t * us_per_M});
    };

    const IntegrationResult res = integrate(launch.s0, c, opt, observer);

    // ---------------------------------------------------------------- report
    std::printf("Launch\n");
    std::printf("  spin a = %.3f, %s, r0 = %.3f M, angular momentum reduced by %.0f%%\n",
                a, dir == Direction::Prograde ? "prograde" : "retrograde", r0, lfrac * 100.0);
    std::printf("  E = %.5f   L = %.5f   initial dr/dtau = %.5f\n", c.E, c.L, launch.s0.r_dot);
    std::printf("  horizon r+ = %.5f M      black hole mass = %g solar masses (1 M = %.4f us)\n\n",
                r_plus, mass, us_per_M);

    std::printf("Journey\n");
    if (res.reached_horizon) {
        std::printf("  Reached the horizon (r = r+ + %.0e) after %ld steps.\n", opt.eps, res.steps);
        std::printf("  Proper time   tau = %9.4f M   (%.4f us)\n", res.tau, res.tau * us_per_M);
        std::printf("  Coordinate t      = %9.4f M   (%.4f us)   t/tau = %.2f\n",
                    res.final_state.t, res.final_state.t * us_per_M, res.final_state.t / res.tau);
        std::printf("  Angle swept       = %.2f turns\n", res.final_state.phi / (2.0 * kPi));
        std::printf("  Note: t grows without limit as r -> r+, so this value depends on --eps.\n");
    } else {
        std::printf("  The particle did NOT reach the horizon within tau = %.0f M\n", res.tau);
        std::printf("  (it is on a bound orbit). Try a larger --lfrac.\n");
    }
    std::printf("  Frames written: %ld   Largest constraint error: %.2e\n\n",
                csv.rows_written(), res.max_constraint_error);

    std::printf("Muon decay (2.2 us of its own time = %.4f M for this mass)\n", tau_decay);
    if (decayed) {
        std::printf("  The muon decays at tau = %.4f M:\n", decay_tau);
        std::printf("    radius r = %.4f M   angle = %.1f deg   coordinate time t = %.4f M (%.3f us)\n",
                    decay_state.r, decay_state.phi * 180.0 / kPi, decay_state.t,
                    decay_state.t * us_per_M);
        if (res.reached_horizon)
            std::printf("    that is %.1f%% of the way (in proper time) to the horizon.\n",
                        100.0 * decay_tau / res.tau);
    } else if (res.reached_horizon) {
        std::printf("  The muon reaches the horizon alive.\n");
    }

    if (res.reached_horizon && tau_decay < 0.02 * res.tau) {
        const double mid = kMuonLifetimeSec / (kSolarMassSeconds * 0.5 * res.tau);
        std::printf("\n  NOTE: at this mass the muon decays almost at launch. For the decay to happen\n"
                    "  mid-plunge, try a much smaller black hole, e.g. --mass %.3g\n", mid);
    }

    std::printf("\nFrame data written to: %s\n", std::filesystem::absolute(out_path).string().c_str());
    return 0;
}
