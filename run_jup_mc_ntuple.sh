#!/bin/bash
# Run JUP MC Ntuple analysis with xrootd input on T2_CN_Beijing
# Mirrors run_analysis.sh paths and run_jjp_mc_ntuple.sh interface

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Default parameters
MODE="DPS_1"          # SPS, DPS_1, DPS_2, DPS_3, TPS
MAX_EVENTS=-1
JPSI_MUON_ID="soft"
UPS_MUON_ID="tight"
INPUT_DIR=""
OUTPUT=""
JOBS=1

# Default MC base on T2_CN_Beijing (IHEP) via xrootd
XROOTD_BASE="root://cceos.ihep.ac.cn//eos/ihep/cms/store/user/xcheng/MC_Production/output"

print_help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -m, --mode MODE       SPS | DPS_1 | DPS_2 | DPS_3 | TPS (default: DPS_1)"
    echo "  -i, --input-dir DIR   Override input directory (xrootd/EOS/local)"
    echo "  -n, --max-events N    Maximum events to process (-1=all)"
    echo "  --jpsi-muon-id TYPE   Muon ID for J/psi (soft/medium/tight/loose/none)"
    echo "  --ups-muon-id TYPE    Muon ID for Upsilon (tight/medium/loose/soft/none)"
    echo "  -o, --output FILE     Output ROOT file (default: output/jup_mc_<MODE>_correlations.root)"
    echo "  -j, --jobs N          Parallel processes (default: 1)"
    echo "  -h, --help            Show this help"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--mode)
            MODE="$2"; shift 2;;
        -n|--max-events)
            MAX_EVENTS="$2"; shift 2;;
        --jpsi-muon-id)
            JPSI_MUON_ID="$2"; shift 2;;
        --ups-muon-id)
            UPS_MUON_ID="$2"; shift 2;;
        -i|--input-dir)
            INPUT_DIR="$2"; shift 2;;
        -o|--output)
            OUTPUT="$2"; shift 2;;
        -j|--jobs)
            JOBS="$2"; shift 2;;
        -h|--help)
            print_help; exit 0;;
        *)
            echo "Unknown option: $1"; print_help; exit 1;;
    esac
done

# Normalize mode to uppercase
MODE_UP=$(echo "$MODE" | tr '[:lower:]' '[:upper:]')

# Choose input directory if not provided
if [[ -z "$INPUT_DIR" ]]; then
    case "$MODE_UP" in
        SPS)
            INPUT_DIR="${XROOTD_BASE}/JUP_SPS";;
        DPS_1|DPS1)
            INPUT_DIR="${XROOTD_BASE}/JUP_DPS1";;
        DPS_2|DPS2)
            INPUT_DIR="${XROOTD_BASE}/JUP_DPS2";;
        DPS_3|DPS3)
            INPUT_DIR="${XROOTD_BASE}/JUP_DPS3";;
        TPS)
            INPUT_DIR="${XROOTD_BASE}/JUP_TPS";;
        *)
            echo "Invalid mode: $MODE (choose SPS, DPS_1, DPS_2, DPS_3, TPS)"; exit 1;;
    esac
fi

# Derive default output if not set
if [[ -z "$OUTPUT" ]]; then
    mkdir -p output
    OUTPUT="output/jup_mc_${MODE_UP}_correlations.root"
fi

PLOT_DIR="$(dirname "$OUTPUT")/plots_JUP_MC_${MODE_UP}"

cat <<EOF
==========================================
JUP MC Ntuple Correlation Analysis
==========================================
Mode: ${MODE_UP}
Input dir: ${INPUT_DIR}
Max events: ${MAX_EVENTS}
J/psi muon ID: ${JPSI_MUON_ID}
Upsilon muon ID: ${UPS_MUON_ID}
Output: ${OUTPUT}
Plots: ${PLOT_DIR}
Jobs: ${JOBS}
==========================================
EOF

echo "Running analysis..."
CMD=(python3 analyze_ntuple_JUP.py -n "$MAX_EVENTS" --jpsi-muon-id "$JPSI_MUON_ID" --ups-muon-id "$UPS_MUON_ID" -i "$INPUT_DIR" -o "$OUTPUT" -j "$JOBS")
"${CMD[@]}"
STATUS=$?

if [[ $STATUS -ne 0 ]]; then
    echo "Analysis failed (exit $STATUS)"; exit $STATUS
fi

echo "Creating plots..."
mkdir -p "$PLOT_DIR"
python3 plot_ntuple_results.py -i "$OUTPUT" -o "$PLOT_DIR" -p JUP
EXIT_PLOT=$?

if [[ $EXIT_PLOT -ne 0 ]]; then
    echo "Plotting failed (exit $EXIT_PLOT)"; exit $EXIT_PLOT
fi

echo "=========================================="
echo "MC Ntuple analysis complete"
echo "Output histograms: $OUTPUT"
echo "Output plots: $PLOT_DIR/"
echo "=========================================="
