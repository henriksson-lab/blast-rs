#!/bin/bash
# Download C. elegans genome and build BLAST databases for testing
set -e
cd "$(dirname "$0")"

echo "Downloading C. elegans genome..."
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz" -O celegans.fna.gz
gunzip -f celegans.fna.gz

echo "Building v4 BLAST database..."
MAKEDB="../../ncbi-blast-2.17.0+-src/c++/ReleaseMT/bin/makeblastdb"
if [ ! -f "$MAKEDB" ]; then
    echo "makeblastdb not found at $MAKEDB"
    exit 1
fi
$MAKEDB -in celegans.fna -dbtype nucl -out celegans4 -title "C. elegans WBcel235" -blastdb_version 4

echo "Extracting test queries..."
head -1000 celegans.fna | tail -8 | tr -d '\n' | cut -c1-500 > /tmp/seq.txt
echo ">celegans_query_500bp" > query_500.fa
cat /tmp/seq.txt >> query_500.fa
echo "" >> query_500.fa

head -2000 celegans.fna | tail -28 | tr -d '\n' | cut -c1-2000 > /tmp/seq.txt
echo ">celegans_query_2000bp" > query_2000.fa
cat /tmp/seq.txt >> query_2000.fa
echo "" >> query_2000.fa

echo "Done! Files:"
ls -lh celegans4.* query_*.fa
