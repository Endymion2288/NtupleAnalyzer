#!/bin/bash
# Run JJP Ntuple analysis

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Default parameters
MAX_EVENTS=-1
MUON_ID="soft"
OUTPUT=""
JOBS=16

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--max-events)
            MAX_EVENTS="$2"
            shift 2
            ;;
        --muon-id)
            MUON_ID="$2"
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
            echo "  -n, --max-events N   Maximum events to process (-1=all)"
            echo "  --muon-id TYPE       Muon ID requirement (soft/medium/tight/loose/none)"
            echo "  -o, --output FILE    Output ROOT file"
            echo "  -j, --jobs N         Parallel processes"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "=========================================="
echo "JJP Ntuple Correlation Analysis"
echo "=========================================="
echo "Max events: $MAX_EVENTS"
echo "Muon ID: $MUON_ID"

# Build command
CMD="python3 analyze_ntuple_JJP.py -n $MAX_EVENTS --muon-id $MUON_ID -j $JOBS"
if [ -n "$OUTPUT" ]; then
    CMD="$CMD -o $OUTPUT"
fi

echo "Running: $CMD"
echo "=========================================="

$CMD

# If analysis succeeded, run plotting
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Creating Plots"
    echo "=========================================="
    
    if [ -n "$OUTPUT" ]; then
        INPUT_FILE="$OUTPUT"
        PLOT_DIR="$(dirname $OUTPUT)/plots_JJP"
    else
        INPUT_FILE="output/jjp_ntuple_correlations.root"
        PLOT_DIR="output/plots_JJP"
    fi
    
    mkdir -p "$PLOT_DIR"
    python3 plot_ntuple_results.py -i "$INPUT_FILE" -o "$PLOT_DIR" -p JJP
    
    echo ""
    echo "=========================================="
    echo "Analysis complete!"
    echo "=========================================="
    echo "Output histograms: $INPUT_FILE"
    echo "Output plots: $PLOT_DIR/"
fi
