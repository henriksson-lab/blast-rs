#!/bin/bash
# Download a real NCBI BLAST database for testing.
# Uses 16S_ribosomal_RNA which is small (~50MB) but exercises real-world conditions.
#
# Usage: bash tests/download_test_db.sh /path/to/store/db
#
# After downloading, test with:
#   echo ">primer1
#   CCTACGGGNGGCWGCAG" > /tmp/primer.fa
#   ./target/release/blast-cli blastn --task blastn-short --evalue 5 \
#     --max_target_seqs 10000 --max_hsps 1 --num_threads 8 \
#     --db /path/to/store/db/16S_ribosomal_RNA \
#     --query /tmp/primer.fa --outfmt '6 qseqid saccver bitscore evalue'

set -euo pipefail

DB_DIR="${1:-/husky/henriksson/for_claude/blast}"
mkdir -p "$DB_DIR"

echo "Downloading 16S_ribosomal_RNA database to $DB_DIR..."
if command -v update_blastdb.pl &>/dev/null; then
    cd "$DB_DIR"
    update_blastdb.pl --decompress 16S_ribosomal_RNA
elif command -v wget &>/dev/null; then
    cd "$DB_DIR"
    wget -q "https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz"
    tar xzf 16S_ribosomal_RNA.tar.gz
    rm -f 16S_ribosomal_RNA.tar.gz
else
    echo "Error: need wget or update_blastdb.pl"
    exit 1
fi

echo "Database downloaded to $DB_DIR/16S_ribosomal_RNA"
echo "Files:"
ls -lh "$DB_DIR"/16S_ribosomal_RNA.* 2>/dev/null || echo "(no files found)"
