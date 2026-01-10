#!/bin/bash
# Run both JJP and JUP analyses

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "Running JJP and JUP Ntuple Analysis"
echo "=========================================="

# Run JJP analysis
echo ""
echo ">>> Running JJP analysis..."
./run_jjp_analysis.sh "$@"

echo ""
echo "=========================================="

# Run JUP analysis
echo ""
echo ">>> Running JUP analysis..."
./run_jup_analysis.sh "$@"

echo ""
echo "=========================================="
echo "All analyses complete!"
echo "=========================================="
