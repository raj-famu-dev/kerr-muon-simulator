#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "monograph/csv_writer.hpp"
#include "monograph/kerr.hpp"
#include "monograph/orbits.hpp"
#include "monograph/rk4.hpp"

using namespace monograph;

static int failures = 0;

static void check(bool ok, const char* name) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++failures;
}

static bool near(double x, double y, double tol) { return std::fabs(x - y) <= tol; }

int main() {
    const Direction pro = Direction::Prograde, ret = Direction::Retrograde;

    // 1. Launcher: the launch state satisfies the constraint exactly
    {
        Launch l;
        check(make_plunge_launch(6.0, 0.9, pro, 0.3, l), "launcher succeeds at r0 = 6");
        check(std::fabs(constraint_error(l.s0, l.c)) < 1e-12, "launch state satisfies the constraint");
        check(l.s0.r_dot < 0.0, "launch velocity points inward");
        check(!make_plunge_launch(2.5, 0.0, pro, 0.3, l), "launcher refuses r0 inside the photon orbit");
    }

    // 2. A full prograde Kerr plunge: constraint holds to the horizon, t diverges
    {
        Launch l;
        make_plunge_launch(6.0, 0.9, pro, 0.3, l);
        IntegratorOptions opt;
        const double r_plus = horizon_radius(0.9);

        // Record t when the gap to the horizon first drops below 1e-2, 1e-3, ...
        std::vector<double> t_at(6, -1.0);
        auto obs = [&](double, const State& s) {
            const double gap = s.r - r_plus;
            for (int k = 2; k <= 5; ++k)
                if (t_at[k] < 0.0 && gap <= std::pow(10.0, -k)) t_at[k] = s.t;
        };
        const IntegrationResult res = integrate(l.s0, l.c, opt, obs);

        std::printf("  tau = %.4f M, t = %.4f M, steps = %ld, drift = %.2e\n",
                    res.tau, res.final_state.t, res.steps, res.max_constraint_error);
        check(res.reached_horizon, "a=0.9 prograde plunge reaches the horizon");
        check(res.max_constraint_error < 1e-8, "constraint drift < 1e-8 all the way to the horizon");
        check(res.final_state.t > 3.0 * res.tau, "coordinate time far exceeds proper time");

        // Near the horizon t grows like -(r+^2 + a^2)/(r+ - r-) * ln(gap)
        // (for a = 0 this is the familiar t = -2M ln(r - 2M)), so every decade
        // of shrinking gap adds ln(10) * (r+^2 + a^2)/(r+ - r-) to t.
        const double r_minus = 1.0 - std::sqrt(1.0 - 0.9 * 0.9);
        const double expected = std::log(10.0) * (r_plus * r_plus + 0.9 * 0.9) / (r_plus - r_minus);
        const double d1 = t_at[4] - t_at[3], d2 = t_at[5] - t_at[4];
        std::printf("  t gained per decade of gap: %.3f, %.3f   (theory = %.3f)\n", d1, d2, expected);
        check(std::fabs(d2 / expected - 1.0) < 0.02, "t diverges logarithmically at the predicted rate");
    }

    // 2b. NFR-1 across spins and directions with the default step controls
    {
        double worst = 0.0;
        int count = 0;
        for (double a : {0.0, 0.5, 0.9, 0.998}) {
            for (Direction d : {pro, ret}) {
                if (a == 0.0 && d == ret) continue;
                for (double r0 : {6.0, 10.0}) {
                    Launch l;
                    if (!make_plunge_launch(r0, a, d, 0.5, l)) continue;
                    IntegratorOptions opt;
                    opt.tau_max = 500.0;
                    const IntegrationResult res = integrate(l.s0, l.c, opt);
                    if (res.reached_horizon) {
                        worst = std::max(worst, res.max_constraint_error);
                        ++count;
                    }
                }
            }
        }
        std::printf("  %d plunges across spins and directions, worst drift = %.2e\n", count, worst);
        check(count >= 10, "enough plunges were checked");
        check(worst < 1e-8, "constraint drift < 1e-8 for every spin and direction");
    }

    // 3. Landing exactly on the decay time
    {
        Launch l;
        make_plunge_launch(6.0, 0.9, pro, 0.3, l);
        IntegratorOptions opt;
        opt.land_on_tau = 0.7123456;
        bool landed = false;
        auto obs = [&](double tau, const State&) { if (std::fabs(tau - 0.7123456) < 1e-12) landed = true; };
        integrate(l.s0, l.c, opt, obs);
        check(landed, "a frame lands exactly on the requested proper time");
    }

    // 4. At a = 0 prograde and retrograde launches are mirror images
    {
        Launch lp, lr;
        make_plunge_launch(6.0, 0.0, pro, 0.3, lp);
        make_plunge_launch(6.0, 0.0, ret, 0.3, lr);
        IntegratorOptions opt;
        const IntegrationResult rp = integrate(lp.s0, lp.c, opt);
        const IntegrationResult rr = integrate(lr.s0, lr.c, opt);
        std::printf("  a=0: prograde tau = %.6f, retrograde tau = %.6f\n", rp.tau, rr.tau);
        check(near(rp.tau, rr.tau, 1e-9), "a=0: same proper time to the horizon either way");
        check(near(rp.final_state.phi, -rr.final_state.phi, 1e-9), "a=0: angle is mirrored");
    }

    // 5. Frame dragging: same launch, different spin direction, different fate
    {
        Launch lp, lr;
        make_plunge_launch(6.0, 0.9, pro, 0.1, lp);
        make_plunge_launch(6.0, 0.9, ret, 0.1, lr);
        IntegratorOptions opt;
        opt.tau_max = 500.0;
        opt.dtau_max = 0.02;
        const IntegrationResult rp = integrate(lp.s0, lp.c, opt);
        const IntegrationResult rr = integrate(lr.s0, lr.c, opt);
        check(!rp.reached_horizon, "a=0.9, 10% kick: prograde particle stays in orbit");
        check(rr.reached_horizon, "a=0.9, 10% kick: retrograde particle plunges");
    }

    // 6. CSV round trip
    {
        const std::string path = "test_csv_roundtrip.csv";
        {
            CsvWriter w(path, {"a", "b", "c"});
            w.write_row({1.0, 2.5, -3.0});
            w.write_row({4.0, 5.5, 6.0});
            check(w.rows_written() == 2, "CsvWriter counts rows");
        }
        std::ifstream in(path);
        std::string header, row1;
        std::getline(in, header);
        std::getline(in, row1);
        check(header == "a,b,c", "CSV header is written");
        check(row1 == "1,2.5,-3", "CSV row is written");
        std::remove(path.c_str());
    }

    std::printf("\n%s\n", failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED");
    return failures == 0 ? 0 : 1;
}
