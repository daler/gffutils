GTF=gencode.v19.annotation.gtf

if [ ! -e "${GTF}" ]; then
    wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
    gunzip ${GTF}.gz
else
    echo "${GTF} up to date"
fi

for n in 10 50 100; do
    head -n ${n}000 ${GTF} > ${GTF}.${n}k
done

