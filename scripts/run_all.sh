#!/usr/bin/env bash
# Build the engine, run all three phases, and regenerate every figure.
# Usage (from the repository root):  bash scripts/run_all.sh
set -e

bash scripts/build.sh
mkdir -p data assets

echo "== Phase 1: ISCO =="
./build/monograph_isco --spin 0.9 --out data/isco_scan.csv

echo "== Phase 2: muon plunges =="
./build/monograph_plunge --spin 0   --out data/schwarzschild.csv
./build/monograph_plunge --spin 0.9 --direction prograde   --out data/prograde.csv
./build/monograph_plunge --spin 0.9 --direction retrograde --out data/retrograde.csv

echo "== Phase 3: figures =="
python python/plot_isco.py data/isco_scan.csv --spin 0.9 --out assets/isco_survival.png
python python/plot_trajectory.py data/prograde.csv --out assets/trajectory.png
python python/plot_time_dilation.py data/prograde.csv --spin 0.9 --out assets/time_dilation.png
python python/compare_spins.py data/schwarzschild.csv data/prograde.csv data/retrograde.csv \
    --labels "a = 0" "a = 0.9 prograde" "a = 0.9 retrograde" --out assets/compare_spins.png
python python/animate_plunge.py data/prograde.csv --out assets/plunge.gif

echo "Done. Figures are in assets/."