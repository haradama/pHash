ACCESSIONs=(NZ_CP007510 NZ_CP008958)

mkdir tmp
for ACCESSION in ${ACCESSIONs[@]}; do
    curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$ACCESSION&rettype=fasta&retmode=text" > ./tmp/$ACCESSION.fna
done

cat ./tmp/*.fna > testData.fna
