# Download large annotation files neede for testing
# to gffutils/test/data/ directory.

cd $(dirname $0)
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M8/gencode.vM8.annotation.gff3.gz
gzip -d gencode.vM8.annotation.gff3.gz
wget ftp://ftp.ensembl.org/pub/release-83/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.83.gff3.gz
gzip -d Saccharomyces_cerevisiae.R64-1-1.83.gff3.gz
