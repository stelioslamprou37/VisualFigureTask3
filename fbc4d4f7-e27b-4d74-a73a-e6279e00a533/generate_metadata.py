#!/usr/bin/env python3

import hashlib
import json
from pathlib import Path

def compute_md5(filepath):
    """Compute MD5 checksum of a file."""
    md5 = hashlib.md5()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            md5.update(chunk)
    return md5.hexdigest()

# Process input_data directory
input_dir = Path('input_data')

if not input_dir.exists():
    print("Error: input_data/ directory not found!")
    print("Please run download.sh first to create input_data/ and download files.")
    exit(1)

print("Generating metadata.jsonl...")
print(f"Scanning {input_dir}/ for files...")

with open('metadata.jsonl', 'w') as out:
    files_found = 0
    for file in sorted(input_dir.iterdir()):
        if file.is_file() and not file.name.startswith('.'):
            print(f"  Processing: {file.name}")
            md5sum = compute_md5(file)
            entry = {"name": file.name, "md5sum": md5sum}
            out.write(json.dumps(entry) + '\n')
            files_found += 1

print(f"\n✓ Created metadata.jsonl with {files_found} files")
print("\nExpected files in input_data/:")
print("  - paper.pdf")
print("  - GSE84437_series_matrix.txt.gz (or .txt after extraction)")
