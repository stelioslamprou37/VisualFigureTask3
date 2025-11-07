#!/bin/bash

# ============================================================================
# Forest Plot Analysis - Automated Data Download Script
# Dataset: GSE84437 (Gastric Cancer)
# Purpose: Download GSE84437 data if not already present
# ============================================================================

set -e  # Exit on error

echo "============================================================================"
echo "Forest Plot Analysis - GSE84437 Dataset Manager"
echo "============================================================================"
echo ""

# Configuration
DATA_DIR="input_data"
FILENAME="GSE84437_series_matrix.txt"
FILEPATH="${DATA_DIR}/${FILENAME}"
GZ_FILENAME="${FILENAME}.gz"
GZ_FILEPATH="${DATA_DIR}/${GZ_FILENAME}"

# Create directory if it doesn't exist
if [ ! -d "${DATA_DIR}" ]; then
    echo "Creating ${DATA_DIR} directory..."
    mkdir -p "${DATA_DIR}"
fi

# ============================================================================
# Check if file already exists
# ============================================================================

if [ -f "${FILEPATH}" ]; then
    echo "✓ File found: ${FILEPATH}"
    FILE_SIZE=$(du -h "${FILEPATH}" | cut -f1)
    echo "  File size: ${FILE_SIZE}"
    echo ""
    echo "Dataset is ready for analysis!"
    echo "Run: Rscript forest_plot_analysis.R"
    exit 0
fi

echo "Dataset file not found: ${FILEPATH}"
echo ""
echo "Attempting automatic download via wget..."
echo ""

# ============================================================================
# Download dataset from GEO
# ============================================================================

GEO_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84437&format=file&file=GSE84437_series_matrix.txt.gz"

# Check if wget is available
if ! command -v wget &> /dev/null; then
    echo "ERROR: wget not found. Please install wget or download manually."
    echo ""
    echo "Manual download option:"
    echo "  1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84437"
    echo "  2. Download: GSE84437_series_matrix.txt.gz (~125 MB)"
    echo "  3. Save to: ./${DATA_DIR}/"
    exit 1
fi

# Download the file
echo "Downloading GSE84437_series_matrix.txt.gz..."
echo "Source: ${GEO_URL}"
echo ""

if wget -q --show-progress -O "${GZ_FILEPATH}" "${GEO_URL}"; then
    echo ""
    echo "✓ Download complete!"
else
    echo "ERROR: Download failed. Please check your internet connection."
    exit 1
fi

# ============================================================================
# Decompress the file
# ============================================================================

echo ""
echo "Decompressing ${GZ_FILENAME}..."

if command -v gunzip &> /dev/null; then
    gunzip -f "${GZ_FILEPATH}"
    echo "✓ Decompression complete!"
elif command -v gzip &> /dev/null; then
    gzip -d -f "${GZ_FILEPATH}"
    echo "✓ Decompression complete!"
else
    echo "ERROR: gunzip or gzip not found. Cannot decompress file."
    exit 1
fi

# ============================================================================
# Verify the extracted file
# ============================================================================

if [ -f "${FILEPATH}" ]; then
    FILE_SIZE=$(du -h "${FILEPATH}" | cut -f1)
    LINE_COUNT=$(wc -l < "${FILEPATH}")
    
    echo ""
    echo "============================================================================"
    echo "✓ Dataset Ready!"
    echo "============================================================================"
    echo "File: ${FILEPATH}"
    echo "Size: ${FILE_SIZE}"
    echo "Lines: ${LINE_COUNT}"
    echo ""
    echo "Next steps:"
    echo "  1. Run analysis: Expert_Solution_RScript.R"
    echo "  2. Output will be saved in output/ directory"
    echo "============================================================================"
else
    echo "ERROR: File extraction failed. Please download manually."
    echo "Manual download: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84437"
    exit 1
fi
