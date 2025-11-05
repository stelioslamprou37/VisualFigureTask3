#!/bin/bash

# ============================================================================
# Forest Plot Analysis - Data Download Script
# Dataset: GSE84437 (Gastric Cancer)
# ============================================================================

echo "============================================================================"
echo "Forest Plot Analysis - GSE84437 Dataset"
echo "============================================================================"
echo ""

# Create directory
mkdir -p input_data

echo "Dataset Information:"
echo "  Name: GSE84437"
echo "  Samples: 483 gastric cancer"
echo "  Platform: Illumina HumanWG-6 v3.0 expression beadchip"
echo "  Paper: Discovery Oncology 2025"
echo ""

echo "Option 1: Automatic Download (via GEOquery)"
echo "  - Run the R script directly"
echo "  - It will download automatically"
echo "  - Command: Rscript expert_solution_forest.R"
echo ""

echo "Option 2: Manual Download from GEO"
echo "  Step 1: Go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84437"
echo "  Step 2: Find 'Supplementary file' section"
echo "  Step 3: Download: GSE84437_series_matrix.txt.gz (~125 MB)"
echo "  Step 4: Save to: ./input_data/"
echo ""

echo "Option 3: Docker (Automated)"
echo "  Command 1: docker build -f Dockerfile_forest -t forest-analysis ."
echo "  Command 2: docker run -v \$(pwd):/workspace forest-analysis"
echo ""

echo "============================================================================"
echo "For more information, see INSTALLATION_GUIDE.md"
echo "============================================================================"
