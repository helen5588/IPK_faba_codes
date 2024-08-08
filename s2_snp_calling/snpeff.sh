snpeff='java -jar /filer-dg/agruppen/dg3/zhangh/21.snpeff/snpEff/snpEff.jar'
Config='/filer-dg/agruppen/dg3/zhangh/21.snpeff/snpEff/snpEff.config'

cd /filer-dg/agruppen/dg3/zhangh/21.snpeff/snpEff//data

mkdir faba_newv2
cd faba_newv2

ln -s Hedin2.genome.fa sequences.fa
ln -s Hedin2.gene.gff genes.gff

#Edit snpEff.config and insert your specific database information:


echo "faba_newv2.genome : faba_newv2" >> snpEff.config

#check by
tail snpEff.config

#Build the database
java -jar $snpeff build -gff3 -v faba_newv2

#run the snpeff
for i in *.vcf.gz;do
 aa=`echo $i|awk -F '/' '{print $NF}'|awk -F '.' '{print $1}'`;
 echo "module load jdk;eval $snpeff eff -c $Config faba_newv2 $i | cut -f -8 | gzip > $aa.SNP_snpeff.vcf.gz" >>ww.sh
done
