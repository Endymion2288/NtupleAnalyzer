#!/bin/bash
# Run JJP MC Ntuple analysis with local or default IHEP paths

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Default parameters
MODE="DPS"          # DPS or TPS
MAX_EVENTS=-1
MUON_ID="soft"
INPUT_DIR=""
OUTPUT=""
JOBS=1

# Default MC directories (IHEP)
MC_DPS_DIR="/eos/user/x/xcheng/learn_MC/JJP_DPS_MC_output/Ntuple/"
MC_TPS_DIR="/eos/user/x/xcheng/learn_MC/JJP_TPS_MC_output/Ntuple/"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--mode)
            MODE="$2"
            shift 2
            ;;
        -n|--max-events)
            MAX_EVENTS="$2"
            shift 2
            ;;
        --muon-id)
            MUON_ID="$2"
            shift 2
            ;;
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -j|--jobs)
            JOBS="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  -m, --mode MODE       MC sample mode: DPS or TPS (default: DPS)"
            echo "  -i, --input-dir DIR   Override input directory (local/EOS)"
            echo "  -n, --max-events N    Maximum events to process (-1=all)"
            echo "  --muon-id TYPE        Muon ID requirement (soft/medium/tight/loose/none)"
            echo "  -o, --output FILE     Output ROOT file"
            echo "  -j, --jobs N          Parallel processes"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Choose input directory if not provided
if [ -z "$INPUT_DIR" ]; then
    case ${MODE^^} in
        DPS)
            INPUT_DIR="$MC_DPS_DIR"
            ;;
        TPS)
            INPUT_DIR="$MC_TPS_DIR"
            ;;
        *)
            echo "Invalid mode: $MODE (choose DPS or TPS)"
            exit 1
            ;;
    esac
fi

# Derive default output if not set
if [ -z "$OUTPUT" ]; then
    mkdir -p output
    OUTPUT="output/jjp_mc_${MODE^^}_correlations.root"
fi

PLOT_DIR="$(dirname "$OUTPUT")/plots_JJP_MC_${MODE^^}"

echo "=========================================="
echo "JJP MC Ntuple Correlation Analysis"
echo "=========================================="
echo "Mode: ${MODE^^}"
echo "Input dir: $INPUT_DIR"
echo "Max events: $MAX_EVENTS"
echo "Muon ID: $MUON_ID"
echo "Output: $OUTPUT"
echo "Plots: $PLOT_DIR"
echo "=========================================="

echo "Running analysis..."
CMD=(python3 analyze_ntuple_JJP.py -n "$MAX_EVENTS" --muon-id "$MUON_ID" -i "$INPUT_DIR" -o "$OUTPUT" -j "$JOBS")
"${CMD[@]}"
STATUS=$?

if [ $STATUS -ne 0 ]; then
    echo "Analysis failed (exit $STATUS)"
    exit $STATUS
fi

echo "Creating plots..."
mkdir -p "$PLOT_DIR"
python3 plot_ntuple_results.py -i "$OUTPUT" -o "$PLOT_DIR" -p JJP

EXIT_PLOT=$?
if [ $EXIT_PLOT -ne 0 ]; then
    echo "Plotting failed (exit $EXIT_PLOT)"
    exit $EXIT_PLOT
fi

echo "=========================================="
echo "MC Ntuple analysis complete"
echo "Output histograms: $OUTPUT"
echo "Output plots: $PLOT_DIR/"
echo "=========================================="
