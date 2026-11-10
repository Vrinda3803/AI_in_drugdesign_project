#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RECEPTOR_DIR="$PROJECT_ROOT/data/external/ftsZ/prepared"
RUN_ROOT="$PROJECT_ROOT/results/docking/cabsdock/runs"
command -v CABSdock >/dev/null 2>&1 || { echo "CABSdock is not installed or not on PATH" >&2; exit 127; }
mkdir -p "$RUN_ROOT"
echo "CABSdock path: $(command -v CABSdock)"
echo "Run plan: $PROJECT_ROOT/results/docking/cabsdock/run_matrix.csv"

if [[ -e "$RUN_ROOT/B2RP_parent__22HR_chainA_apo__seed104729" ]]; then echo "Output already exists: B2RP_parent__22HR_chainA_apo__seed104729" >&2; exit 2; fi
echo "Starting B2RP_parent__22HR_chainA_apo__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWDTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_parent__22HR_chainA_apo__seed104729"

if [[ -e "$RUN_ROOT/B2RP_parent__22HR_chainA_apo__seed130363" ]]; then echo "Output already exists: B2RP_parent__22HR_chainA_apo__seed130363" >&2; exit 2; fi
echo "Starting B2RP_parent__22HR_chainA_apo__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWDTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_parent__22HR_chainA_apo__seed130363"

if [[ -e "$RUN_ROOT/B2RP_parent__22HR_chainA_apo__seed155921" ]]; then echo "Output already exists: B2RP_parent__22HR_chainA_apo__seed155921" >&2; exit 2; fi
echo "Starting B2RP_parent__22HR_chainA_apo__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWDTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_parent__22HR_chainA_apo__seed155921"

if [[ -e "$RUN_ROOT/B2RP_parent__22HQ_chainA_GDPbound__seed104729" ]]; then echo "Output already exists: B2RP_parent__22HQ_chainA_GDPbound__seed104729" >&2; exit 2; fi
echo "Starting B2RP_parent__22HQ_chainA_GDPbound__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWDTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_parent__22HQ_chainA_GDPbound__seed104729"

if [[ -e "$RUN_ROOT/B2RP_parent__22HQ_chainA_GDPbound__seed130363" ]]; then echo "Output already exists: B2RP_parent__22HQ_chainA_GDPbound__seed130363" >&2; exit 2; fi
echo "Starting B2RP_parent__22HQ_chainA_GDPbound__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWDTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_parent__22HQ_chainA_GDPbound__seed130363"

if [[ -e "$RUN_ROOT/B2RP_parent__22HQ_chainA_GDPbound__seed155921" ]]; then echo "Output already exists: B2RP_parent__22HQ_chainA_GDPbound__seed155921" >&2; exit 2; fi
echo "Starting B2RP_parent__22HQ_chainA_GDPbound__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWDTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_parent__22HQ_chainA_GDPbound__seed155921"

if [[ -e "$RUN_ROOT/D4K__22HR_chainA_apo__seed104729" ]]; then echo "Output already exists: D4K__22HR_chainA_apo__seed104729" >&2; exit 2; fi
echo "Starting D4K__22HR_chainA_apo__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K__22HR_chainA_apo__seed104729"

if [[ -e "$RUN_ROOT/D4K__22HR_chainA_apo__seed130363" ]]; then echo "Output already exists: D4K__22HR_chainA_apo__seed130363" >&2; exit 2; fi
echo "Starting D4K__22HR_chainA_apo__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K__22HR_chainA_apo__seed130363"

if [[ -e "$RUN_ROOT/D4K__22HR_chainA_apo__seed155921" ]]; then echo "Output already exists: D4K__22HR_chainA_apo__seed155921" >&2; exit 2; fi
echo "Starting D4K__22HR_chainA_apo__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K__22HR_chainA_apo__seed155921"

if [[ -e "$RUN_ROOT/D4K__22HQ_chainA_GDPbound__seed104729" ]]; then echo "Output already exists: D4K__22HQ_chainA_GDPbound__seed104729" >&2; exit 2; fi
echo "Starting D4K__22HQ_chainA_GDPbound__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K__22HQ_chainA_GDPbound__seed104729"

if [[ -e "$RUN_ROOT/D4K__22HQ_chainA_GDPbound__seed130363" ]]; then echo "Output already exists: D4K__22HQ_chainA_GDPbound__seed130363" >&2; exit 2; fi
echo "Starting D4K__22HQ_chainA_GDPbound__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K__22HQ_chainA_GDPbound__seed130363"

if [[ -e "$RUN_ROOT/D4K__22HQ_chainA_GDPbound__seed155921" ]]; then echo "Output already exists: D4K__22HQ_chainA_GDPbound__seed155921" >&2; exit 2; fi
echo "Starting D4K__22HQ_chainA_GDPbound__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K__22HQ_chainA_GDPbound__seed155921"

if [[ -e "$RUN_ROOT/L18K__22HR_chainA_apo__seed104729" ]]; then echo "Output already exists: L18K__22HR_chainA_apo__seed104729" >&2; exit 2; fi
echo "Starting L18K__22HR_chainA_apo__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWDTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18K__22HR_chainA_apo__seed104729"

if [[ -e "$RUN_ROOT/L18K__22HR_chainA_apo__seed130363" ]]; then echo "Output already exists: L18K__22HR_chainA_apo__seed130363" >&2; exit 2; fi
echo "Starting L18K__22HR_chainA_apo__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWDTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18K__22HR_chainA_apo__seed130363"

if [[ -e "$RUN_ROOT/L18K__22HR_chainA_apo__seed155921" ]]; then echo "Output already exists: L18K__22HR_chainA_apo__seed155921" >&2; exit 2; fi
echo "Starting L18K__22HR_chainA_apo__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWDTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18K__22HR_chainA_apo__seed155921"

if [[ -e "$RUN_ROOT/L18K__22HQ_chainA_GDPbound__seed104729" ]]; then echo "Output already exists: L18K__22HQ_chainA_GDPbound__seed104729" >&2; exit 2; fi
echo "Starting L18K__22HQ_chainA_GDPbound__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWDTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18K__22HQ_chainA_GDPbound__seed104729"

if [[ -e "$RUN_ROOT/L18K__22HQ_chainA_GDPbound__seed130363" ]]; then echo "Output already exists: L18K__22HQ_chainA_GDPbound__seed130363" >&2; exit 2; fi
echo "Starting L18K__22HQ_chainA_GDPbound__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWDTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18K__22HQ_chainA_GDPbound__seed130363"

if [[ -e "$RUN_ROOT/L18K__22HQ_chainA_GDPbound__seed155921" ]]; then echo "Output already exists: L18K__22HQ_chainA_GDPbound__seed155921" >&2; exit 2; fi
echo "Starting L18K__22HQ_chainA_GDPbound__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWDTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18K__22HQ_chainA_GDPbound__seed155921"

if [[ -e "$RUN_ROOT/D4K_L18K__22HR_chainA_apo__seed104729" ]]; then echo "Output already exists: D4K_L18K__22HR_chainA_apo__seed104729" >&2; exit 2; fi
echo "Starting D4K_L18K__22HR_chainA_apo__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18K__22HR_chainA_apo__seed104729"

if [[ -e "$RUN_ROOT/D4K_L18K__22HR_chainA_apo__seed130363" ]]; then echo "Output already exists: D4K_L18K__22HR_chainA_apo__seed130363" >&2; exit 2; fi
echo "Starting D4K_L18K__22HR_chainA_apo__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18K__22HR_chainA_apo__seed130363"

if [[ -e "$RUN_ROOT/D4K_L18K__22HR_chainA_apo__seed155921" ]]; then echo "Output already exists: D4K_L18K__22HR_chainA_apo__seed155921" >&2; exit 2; fi
echo "Starting D4K_L18K__22HR_chainA_apo__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18K__22HR_chainA_apo__seed155921"

if [[ -e "$RUN_ROOT/D4K_L18K__22HQ_chainA_GDPbound__seed104729" ]]; then echo "Output already exists: D4K_L18K__22HQ_chainA_GDPbound__seed104729" >&2; exit 2; fi
echo "Starting D4K_L18K__22HQ_chainA_GDPbound__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18K__22HQ_chainA_GDPbound__seed104729"

if [[ -e "$RUN_ROOT/D4K_L18K__22HQ_chainA_GDPbound__seed130363" ]]; then echo "Output already exists: D4K_L18K__22HQ_chainA_GDPbound__seed130363" >&2; exit 2; fi
echo "Starting D4K_L18K__22HQ_chainA_GDPbound__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18K__22HQ_chainA_GDPbound__seed130363"

if [[ -e "$RUN_ROOT/D4K_L18K__22HQ_chainA_GDPbound__seed155921" ]]; then echo "Output already exists: D4K_L18K__22HQ_chainA_GDPbound__seed155921" >&2; exit 2; fi
echo "Starting D4K_L18K__22HQ_chainA_GDPbound__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18K__22HQ_chainA_GDPbound__seed155921"

if [[ -e "$RUN_ROOT/D4R_hypothesis__22HR_chainA_apo__seed104729" ]]; then echo "Output already exists: D4R_hypothesis__22HR_chainA_apo__seed104729" >&2; exit 2; fi
echo "Starting D4R_hypothesis__22HR_chainA_apo__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWRTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_hypothesis__22HR_chainA_apo__seed104729"

if [[ -e "$RUN_ROOT/D4R_hypothesis__22HR_chainA_apo__seed130363" ]]; then echo "Output already exists: D4R_hypothesis__22HR_chainA_apo__seed130363" >&2; exit 2; fi
echo "Starting D4R_hypothesis__22HR_chainA_apo__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWRTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_hypothesis__22HR_chainA_apo__seed130363"

if [[ -e "$RUN_ROOT/D4R_hypothesis__22HR_chainA_apo__seed155921" ]]; then echo "Output already exists: D4R_hypothesis__22HR_chainA_apo__seed155921" >&2; exit 2; fi
echo "Starting D4R_hypothesis__22HR_chainA_apo__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWRTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_hypothesis__22HR_chainA_apo__seed155921"

if [[ -e "$RUN_ROOT/D4R_hypothesis__22HQ_chainA_GDPbound__seed104729" ]]; then echo "Output already exists: D4R_hypothesis__22HQ_chainA_GDPbound__seed104729" >&2; exit 2; fi
echo "Starting D4R_hypothesis__22HQ_chainA_GDPbound__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWRTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_hypothesis__22HQ_chainA_GDPbound__seed104729"

if [[ -e "$RUN_ROOT/D4R_hypothesis__22HQ_chainA_GDPbound__seed130363" ]]; then echo "Output already exists: D4R_hypothesis__22HQ_chainA_GDPbound__seed130363" >&2; exit 2; fi
echo "Starting D4R_hypothesis__22HQ_chainA_GDPbound__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWRTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_hypothesis__22HQ_chainA_GDPbound__seed130363"

if [[ -e "$RUN_ROOT/D4R_hypothesis__22HQ_chainA_GDPbound__seed155921" ]]; then echo "Output already exists: D4R_hypothesis__22HQ_chainA_GDPbound__seed155921" >&2; exit 2; fi
echo "Starting D4R_hypothesis__22HQ_chainA_GDPbound__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWRTIKSMGKVFAGKILQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_hypothesis__22HQ_chainA_GDPbound__seed155921"

if [[ -e "$RUN_ROOT/L18R_hypothesis__22HR_chainA_apo__seed104729" ]]; then echo "Output already exists: L18R_hypothesis__22HR_chainA_apo__seed104729" >&2; exit 2; fi
echo "Starting L18R_hypothesis__22HR_chainA_apo__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWDTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18R_hypothesis__22HR_chainA_apo__seed104729"

if [[ -e "$RUN_ROOT/L18R_hypothesis__22HR_chainA_apo__seed130363" ]]; then echo "Output already exists: L18R_hypothesis__22HR_chainA_apo__seed130363" >&2; exit 2; fi
echo "Starting L18R_hypothesis__22HR_chainA_apo__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWDTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18R_hypothesis__22HR_chainA_apo__seed130363"

if [[ -e "$RUN_ROOT/L18R_hypothesis__22HR_chainA_apo__seed155921" ]]; then echo "Output already exists: L18R_hypothesis__22HR_chainA_apo__seed155921" >&2; exit 2; fi
echo "Starting L18R_hypothesis__22HR_chainA_apo__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWDTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18R_hypothesis__22HR_chainA_apo__seed155921"

if [[ -e "$RUN_ROOT/L18R_hypothesis__22HQ_chainA_GDPbound__seed104729" ]]; then echo "Output already exists: L18R_hypothesis__22HQ_chainA_GDPbound__seed104729" >&2; exit 2; fi
echo "Starting L18R_hypothesis__22HQ_chainA_GDPbound__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWDTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18R_hypothesis__22HQ_chainA_GDPbound__seed104729"

if [[ -e "$RUN_ROOT/L18R_hypothesis__22HQ_chainA_GDPbound__seed130363" ]]; then echo "Output already exists: L18R_hypothesis__22HQ_chainA_GDPbound__seed130363" >&2; exit 2; fi
echo "Starting L18R_hypothesis__22HQ_chainA_GDPbound__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWDTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18R_hypothesis__22HQ_chainA_GDPbound__seed130363"

if [[ -e "$RUN_ROOT/L18R_hypothesis__22HQ_chainA_GDPbound__seed155921" ]]; then echo "Output already exists: L18R_hypothesis__22HQ_chainA_GDPbound__seed155921" >&2; exit 2; fi
echo "Starting L18R_hypothesis__22HQ_chainA_GDPbound__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWDTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/L18R_hypothesis__22HQ_chainA_GDPbound__seed155921"

if [[ -e "$RUN_ROOT/D4K_L18R_hypothesis__22HR_chainA_apo__seed104729" ]]; then echo "Output already exists: D4K_L18R_hypothesis__22HR_chainA_apo__seed104729" >&2; exit 2; fi
echo "Starting D4K_L18R_hypothesis__22HR_chainA_apo__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18R_hypothesis__22HR_chainA_apo__seed104729"

if [[ -e "$RUN_ROOT/D4K_L18R_hypothesis__22HR_chainA_apo__seed130363" ]]; then echo "Output already exists: D4K_L18R_hypothesis__22HR_chainA_apo__seed130363" >&2; exit 2; fi
echo "Starting D4K_L18R_hypothesis__22HR_chainA_apo__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18R_hypothesis__22HR_chainA_apo__seed130363"

if [[ -e "$RUN_ROOT/D4K_L18R_hypothesis__22HR_chainA_apo__seed155921" ]]; then echo "Output already exists: D4K_L18R_hypothesis__22HR_chainA_apo__seed155921" >&2; exit 2; fi
echo "Starting D4K_L18R_hypothesis__22HR_chainA_apo__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18R_hypothesis__22HR_chainA_apo__seed155921"

if [[ -e "$RUN_ROOT/D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed104729" ]]; then echo "Output already exists: D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed104729" >&2; exit 2; fi
echo "Starting D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed104729"

if [[ -e "$RUN_ROOT/D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed130363" ]]; then echo "Output already exists: D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed130363" >&2; exit 2; fi
echo "Starting D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed130363"

if [[ -e "$RUN_ROOT/D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed155921" ]]; then echo "Output already exists: D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed155921" >&2; exit 2; fi
echo "Starting D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4K_L18R_hypothesis__22HQ_chainA_GDPbound__seed155921"

if [[ -e "$RUN_ROOT/D4R_L18K_hypothesis__22HR_chainA_apo__seed104729" ]]; then echo "Output already exists: D4R_L18K_hypothesis__22HR_chainA_apo__seed104729" >&2; exit 2; fi
echo "Starting D4R_L18K_hypothesis__22HR_chainA_apo__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWRTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18K_hypothesis__22HR_chainA_apo__seed104729"

if [[ -e "$RUN_ROOT/D4R_L18K_hypothesis__22HR_chainA_apo__seed130363" ]]; then echo "Output already exists: D4R_L18K_hypothesis__22HR_chainA_apo__seed130363" >&2; exit 2; fi
echo "Starting D4R_L18K_hypothesis__22HR_chainA_apo__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWRTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18K_hypothesis__22HR_chainA_apo__seed130363"

if [[ -e "$RUN_ROOT/D4R_L18K_hypothesis__22HR_chainA_apo__seed155921" ]]; then echo "Output already exists: D4R_L18K_hypothesis__22HR_chainA_apo__seed155921" >&2; exit 2; fi
echo "Starting D4R_L18K_hypothesis__22HR_chainA_apo__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWRTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18K_hypothesis__22HR_chainA_apo__seed155921"

if [[ -e "$RUN_ROOT/D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed104729" ]]; then echo "Output already exists: D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed104729" >&2; exit 2; fi
echo "Starting D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWRTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed104729"

if [[ -e "$RUN_ROOT/D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed130363" ]]; then echo "Output already exists: D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed130363" >&2; exit 2; fi
echo "Starting D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWRTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed130363"

if [[ -e "$RUN_ROOT/D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed155921" ]]; then echo "Output already exists: D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed155921" >&2; exit 2; fi
echo "Starting D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWRTIKSMGKVFAGKIKQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18K_hypothesis__22HQ_chainA_GDPbound__seed155921"

if [[ -e "$RUN_ROOT/D4R_L18R_hypothesis__22HR_chainA_apo__seed104729" ]]; then echo "Output already exists: D4R_L18R_hypothesis__22HR_chainA_apo__seed104729" >&2; exit 2; fi
echo "Starting D4R_L18R_hypothesis__22HR_chainA_apo__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWRTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18R_hypothesis__22HR_chainA_apo__seed104729"

if [[ -e "$RUN_ROOT/D4R_L18R_hypothesis__22HR_chainA_apo__seed130363" ]]; then echo "Output already exists: D4R_L18R_hypothesis__22HR_chainA_apo__seed130363" >&2; exit 2; fi
echo "Starting D4R_L18R_hypothesis__22HR_chainA_apo__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWRTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18R_hypothesis__22HR_chainA_apo__seed130363"

if [[ -e "$RUN_ROOT/D4R_L18R_hypothesis__22HR_chainA_apo__seed155921" ]]; then echo "Output already exists: D4R_L18R_hypothesis__22HR_chainA_apo__seed155921" >&2; exit 2; fi
echo "Starting D4R_L18R_hypothesis__22HR_chainA_apo__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWRTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18R_hypothesis__22HR_chainA_apo__seed155921"

if [[ -e "$RUN_ROOT/D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed104729" ]]; then echo "Output already exists: D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed104729" >&2; exit 2; fi
echo "Starting D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWRTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed104729"

if [[ -e "$RUN_ROOT/D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed130363" ]]; then echo "Output already exists: D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed130363" >&2; exit 2; fi
echo "Starting D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWRTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed130363"

if [[ -e "$RUN_ROOT/D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed155921" ]]; then echo "Output already exists: D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed155921" >&2; exit 2; fi
echo "Starting D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWRTIKSMGKVFAGKIRQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/D4R_L18R_hypothesis__22HQ_chainA_GDPbound__seed155921"

if [[ -e "$RUN_ROOT/B2RP_top_published_lead_control__22HR_chainA_apo__seed104729" ]]; then echo "Output already exists: B2RP_top_published_lead_control__22HR_chainA_apo__seed104729" >&2; exit 2; fi
echo "Starting B2RP_top_published_lead_control__22HR_chainA_apo__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGAIKQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_top_published_lead_control__22HR_chainA_apo__seed104729"

if [[ -e "$RUN_ROOT/B2RP_top_published_lead_control__22HR_chainA_apo__seed130363" ]]; then echo "Output already exists: B2RP_top_published_lead_control__22HR_chainA_apo__seed130363" >&2; exit 2; fi
echo "Starting B2RP_top_published_lead_control__22HR_chainA_apo__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGAIKQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_top_published_lead_control__22HR_chainA_apo__seed130363"

if [[ -e "$RUN_ROOT/B2RP_top_published_lead_control__22HR_chainA_apo__seed155921" ]]; then echo "Output already exists: B2RP_top_published_lead_control__22HR_chainA_apo__seed155921" >&2; exit 2; fi
echo "Starting B2RP_top_published_lead_control__22HR_chainA_apo__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HR_chainA_apo.pdb" -p "GIWKTIKSMGKVFAGAIKQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_top_published_lead_control__22HR_chainA_apo__seed155921"

if [[ -e "$RUN_ROOT/B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed104729" ]]; then echo "Output already exists: B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed104729" >&2; exit 2; fi
echo "Starting B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed104729"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGAIKQNL" -a 20 -y 50 -s 200 -z 104729 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed104729"

if [[ -e "$RUN_ROOT/B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed130363" ]]; then echo "Output already exists: B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed130363" >&2; exit 2; fi
echo "Starting B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed130363"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGAIKQNL" -a 20 -y 50 -s 200 -z 130363 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed130363"

if [[ -e "$RUN_ROOT/B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed155921" ]]; then echo "Output already exists: B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed155921" >&2; exit 2; fi
echo "Starting B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed155921"
CABSdock -i "$RECEPTOR_DIR/22HQ_chainA_GDPbound.pdb" -p "GIWKTIKSMGKVFAGAIKQNL" -a 20 -y 50 -s 200 -z 155921 -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w "$RUN_ROOT/B2RP_top_published_lead_control__22HQ_chainA_GDPbound__seed155921"

echo "All CABSdock runs finished. Inspect outputs and logs before importing any score."
