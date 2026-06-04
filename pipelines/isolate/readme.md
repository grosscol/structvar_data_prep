# Isolate structural variants to avoid overlap
Generate crams with sequences for multiple structural variant.
Where the variants are non-overlapping by at least 200 bp.

## Pipeline
1. isolate.sh
    Make series of structural variant vcf-like files where each structural variant is isolated by at least 200 bp.
2. extract.sh
    Extract the sequences relevant to each structural variant in the isolated files to produce corresponding crams.
3. index.sh
    Process index made by isolate.sh into simpler format suitable for direct ingest into mongodb.
