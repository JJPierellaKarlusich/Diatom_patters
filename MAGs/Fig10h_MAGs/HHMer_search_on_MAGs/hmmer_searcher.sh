
wget https://www.genoscope.cns.fr/tara/localdata/data/SMAGs-v1/SMAGs_v1_concat.faa .

my_pfams=$(ls HHMp/ | sed 's/\.hmm//')

for mypfam in $my_pfams

do

echo $mypfam

hmmsearch --domt $mypfam.vs.SMAGs_v1_concat.domt -E 1 HHMp/$mypfam.hmm ../SMAGs/SMAGs_v1_concat.faa > /dev/null

cat $mypfam.vs.SMAGs_v1_concat.domt | sed 's/ .*//; /^#$/d; /^#/d' > $mypfam.vs.SMAGs_v1_concat.protlist

# proteins
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $mypfam.vs.SMAGs_v1_concat.protlist ../SMAGs/SMAGs_v1_concat.faa  > $mypfam.vs.SMAGs_v1_concat.faa

# unigenes
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $mypfam.vs.SMAGs_v1_concat.protlist ../SMAGs/SMAGs_v1_concat.fna  > $mypfam.vs.SMAGs_v1_concat.fna

rm $mypfam.vs.SMAGs_v1_concat.protlist

##### Retriving Pfam domain regions from the translated sequences ##########################

cat $mypfam.vs.SMAGs_v1_concat.domt | sed 's/ \+/\t/g' | cut -f1,20-21 | sed '/#/d; s/\t/-/2'  > $mypfam.translatedUnigene_HMMcoordonates

makeblastdb -in $mypfam.vs.SMAGs_v1_concat.faa -dbtype prot -parse_seqids

cat $mypfam.translatedUnigene_HMMcoordonates | while IFS=$'\t' read -r -a myArray
do
blastdbcmd -entry ${myArray[0]} -range ${myArray[1]} -db $mypfam.vs.SMAGs_v1_concat.faa >> $mypfam.vs.SMAGs_v1_concat.pHMMregion.faa
done

sed 's/>lcl|/>/; s/ unnamed protein product//; s/:/__/' $mypfam.vs.SMAGs_v1_concat.pHMMregion.faa -i

rm $mypfam.vs.SMAGs_v1_concat.faa.p* $mypfam.translatedUnigene_HMMcoordonates


done
