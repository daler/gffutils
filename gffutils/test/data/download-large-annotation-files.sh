# Download large annotation files neede for testing
# to gffutils/test/data/ directory.

cd $(dirname $0)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.annotation.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-83/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.83.gff3.gz
gzip -d Saccharomyces_cerevisiae.R64-1-1.83.gff3.gz
