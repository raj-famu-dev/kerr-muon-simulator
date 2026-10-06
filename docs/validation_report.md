# Validation Report

Every number below was produced by this repository's own code: the C++ engine
(`test_kerr`, `test_rk4`, `test_isco`, `test_constraint`, `monograph_isco`,
`monograph_plunge`) and the Python scripts in `python/`. Units are geometric
(G = c = M = 1) unless a time is given in microseconds.

## 1. Kerr equations

`test_kerr` (16 checks) confirms that the geodesic equations

* reduce to the Schwarzschild closed forms when `a = 0`,
* give zero radial velocity and zero radial acceleration on circular orbits
  at `a = 0.9`, both prograde and retrograde (Bardeen's independent formulas),
* satisfy `dφ/dt = ±1 / (r^{3/2} ± a)` on those orbits, which ties the
  `dt/dτ` and `dφ/dτ` equations to each other.

## 2. RK4 convergence

Integrating the harmonic oscillator `x'' = -x` to `t = 10`, halving the step
should cut the global error by 2⁴ = 16.

| step h | error | ratio to previous |
|---|---|---|
| 0.1 | 3.935 × 10⁻⁶ | |
| 0.05 | 2.649 × 10⁻⁷ | 14.86 |
| 0.025 | 1.714 × 10⁻⁸ | 15.46 |

The ratios approach 16, the signature of a 4th-order method.

## 3. Innermost stable circular orbit

The engine finds the ISCO by launching a particle on the circular orbit at each
radius, nudging it inward, integrating for 20 000 M, and bisecting on
survive / plunge. Ground truth is Bardeen's closed form.

| Spin a | Direction | Analytic | Engine | Error |
|---|---|---|---|---|
| 0 | either | 6.0000 | 6.0005 | +0.0005 |
| 0.5 | prograde | 4.2330 | 4.2336 | +0.0006 |
| 0.5 | retrograde | 7.5546 | 7.5554 | +0.0008 |
| 0.9 | prograde | 2.3209 | 2.3211 | +0.0002 |
| 0.9 | retrograde | 8.7174 | 8.7187 | +0.0013 |
| 0.998 | prograde | 1.2370 | 1.2368 | −0.0002 |
| 0.998 | retrograde | 8.9944 | 8.9954 | +0.0010 |

All errors are within the bisection tolerance (10⁻³ M) plus the nudge-growth
resolution. For `a = 0.9` a prograde particle survives **6.40 M** closer to the
hole than a retrograde one.

![ISCO survival map](../assets/isco_survival.png)

## 4. Numerical error along a plunge

The radial first integral `(dr/dτ)² = E² − 1 + 2/r + K/r² + 2J/r³` is not used
by the integrator, so its violation measures the integration error directly.

| Run (r₀ = 6 M, L reduced 30%) | Largest `|ṙ² − R(r)|` |
|---|---|
| a = 0 | 2.0 × 10⁻¹³ |
| a = 0.9 prograde | 5.6 × 10⁻¹⁵ |
| a = 0.9 retrograde | 1.0 × 10⁻¹⁰ |

Across 14 plunges covering a ∈ {0, 0.5, 0.9, 0.998}, both directions, and
r₀ ∈ {6, 10} M, the worst drift is 6.6 × 10⁻¹⁰, below the 10⁻⁸ requirement (NFR-1).

## 5. Coordinate time diverges logarithmically

Near the horizon `t ≈ −[(r₊² + a²) / (r₊ − r₋)] · ln(r − r₊)`. For `a = 0.9`
every factor-of-ten decrease of `r − r₊` should therefore add 7.585 M to `t`.
The engine measures:

| Decade | Measured Δt | Theory |
|---|---|---|
| 10⁻³ → 10⁻⁴ | 7.588 M | 7.585 M |
| 10⁻⁴ → 10⁻⁵ | 7.582 M | 7.585 M |

The proper time at the same moment stays finite, so `dt/dτ → ∞` while `τ` does not.

![Time dilation](../assets/time_dilation.png)

## 6. Where the muon dies

The muon lives 2.2 μs of its own time. For a black hole of 0.1 M☉
(1 M = 0.4925 μs, so 2.2 μs = 4.47 M) launched from r₀ = 6 M with its angular
momentum reduced by 30%:

| | a = 0 | a = 0.9 prograde | a = 0.9 retrograde |
|---|---|---|---|
| Muon clock at the horizon | 4.22 μs | 5.99 μs | 3.32 μs |
| Distant clock at the horizon (r − r₊ = 10⁻⁶) | 21.15 μs | 33.02 μs | 29.05 μs |
| Angle swept | 0.23 turns | 2.67 turns | 1.97 turns |
| Radius where the muon decays | 4.32 M | 4.55 M | 3.81 M |
| Share of the plunge completed at decay | 52% | 37% | 66% |

Frame dragging does not change the muon's lifetime. It changes **where** along
the path that lifetime runs out, and how much coordinate time a distant
observer records.

![Comparison](../assets/compare_spins.png)

## 7. Limits and caveats

* **Coordinate time at the horizon depends on the cutoff.** `t` diverges, so the
  value reported at `r = r₊ + ε` grows by about 7.6 M (for a = 0.9) every time
  `ε` shrinks tenfold. Proper time and angle are insensitive to `ε` in
  comparison.
* **Validated range.** The 10⁻⁸ drift requirement is checked for launches from
  r₀ ≥ 6 M. Launches very close to the photon orbit of a near-extremal hole
  (for example a = 0.998 retrograde from r₀ = 4 M) have energies above 16 and
  need a much smaller `--dtau` to reach the same accuracy.
* **Mass scale.** For a 10 M☉ hole the muon's 2.2 μs is only 0.045 M, so it
  decays almost at launch. The default of 0.1 M☉ puts the decay mid-plunge.
* **Coordinates.** Boyer–Lindquist `x = r cos φ`, `y = r sin φ` are plotting
  coordinates; BL radius is not a Euclidean distance close to the hole.
