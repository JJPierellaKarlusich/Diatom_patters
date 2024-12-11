#!/bin/bash

# Download the files
mkdir MATOU-v1.5
wget https://www.genoscope.cns.fr/tara/localdata/data/Geneset-v1.5/MATOU-v1.5.taxonomy.tsv.gz MATOU-v1.5/
wget https://www.genoscope.cns.fr/tara/localdata/data/Geneset-v1.5/MATOU-v1.5.pfam.gz MATOU-v1.5/
wget https://www.genoscope.cns.fr/tara/localdata/data/Geneset-v1.5/MATOU-v1.5.metaG.occurrences.gz MATOU-v1.5/
wget https://www.genoscope.cns.fr/tara/localdata/data/Geneset-v1.5/MATOU-v1.5.metaT.occurrences.gz MATOU-v1.5/

mkdir Pfam_per_taxon
mkdir intermediate_files

# List of photosynthetic taxonomy groups
mytaxonlist=("Bacillariophyta" "Pelagophyceae" "Dictyochophyceae" "Chlorophyta" "Haptophyta" "Dinophyceae" "Chlorarachniophyceae" "Euglenida" "Cryptophyta")


for mytaxon in $mytaxonlist
do

echo
echo "Retrieving taxonomy for" $mytaxon
zcat MATOU-v1.5/MATOU-v1.5.taxonomy.tsv.gz | head -n 1 > $mytaxon.taxonomy.tsv
zcat MATOU-v1.5/MATOU-v1.5.taxonomy.tsv.gz | egrep $mytaxon -i -w >> $mytaxon.taxonomy.tsv
echo "Done!"
echo

cut -f1 $mytaxon.taxonomy.tsv | sed '/^geneID$/d' > $mytaxon.unigenelist


echo
echo "Retrieving Pfam assignation for" $mytaxon
zcat MATOU-v1.5/MATOU-v1.5.pfam.gz | head -n 1 > $mytaxon.MATOU-v1.5.Pfam
awk -F "\t" 'FNR==NR {hash[$1]; next} ($2 in hash)' $mytaxon.unigenelist  <(gzip -dc MATOU-v1.5/MATOU-v1.5.pfam.gz) >> $mytaxon.MATOU-v1.5.Pfam
cut -f2-13 $mytaxon.MATOU-v1.5.Pfam > $mytaxon.MATOU-v1.5.Pfam2
mv $mytaxon.MATOU-v1.5.Pfam2 $mytaxon.MATOU-v1.5.Pfam
echo "Done!"
echo


echo
echo "retrieving metaG counts for" $mytaxon
head <(gzip -dc  MATOU-v1.5/MATOU-v1.5.metaG.occurrences.gz) -n 1 > $mytaxon.MATOU-v1.5.metaG.occurrences.tsv
awk -F "\t" 'FNR==NR {hash[$1]; next} $1 in hash' $mytaxon.unigenelist <(gzip -dc  MATOU-v1.5/MATOU-v1.5.metaG.occurrences.gz)  >> $mytaxon.MATOU-v1.5.metaG.occurrences.tsv
echo "Done!"
echo

echo
echo "retrieving metaT counts for" $mytaxon
head <(gzip -dc  MATOU-v1.5/MATOU-v1.5.metaT.occurrences.gz) -n 1 > $mytaxon.MATOU-v1.5.metaT.occurrences.tsv
awk -F "\t" 'FNR==NR {hash[$1]; next} $1 in hash' $mytaxon.unigenelist <(gzip -dc MATOU-v1.5/MATOU-v1.5.metaT.occurrences.gz)  >> $mytaxon.MATOU-v1.5.metaT.occurrences.tsv
echo "Done!"
echo

echo
echo "Rscript"
Rscript data_organizer_MATOUv2_PfamPerTaxon_inputLong.R $mytaxon
echo "Done!"
echo

mv $mytaxon.MATOU-v1.5.Pfam.metaG.tsv Pfam_per_taxon/
mv $mytaxon.MATOU-v1.5.Pfam.metaT.tsv Pfam_per_taxon/

mv $mytaxon.* intermediate_files/

done


