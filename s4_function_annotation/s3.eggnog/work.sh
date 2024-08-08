pep=Hedin2.v2.transcript.pep.fa
module load eggnog-mapper/2.1.9
diamond blastp --db /filer-dg/agruppen/dg3/zhangh/software/database/emapperdb-5.0.2/eggnog_proteins.dmnd --query $pep --out $pep.eggnog.tab --outfmt 6 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 30 --block-size 20.0 --index-chunks 1 --threads 50
parsing_blast_result.pl -type tab --max-hit-num 1 --db-query $pep --db-subject /filer-dg/agruppen/dg3/zhangh/software/database/emapperdb-5.0.2/eggnog_proteins.dmnd --HSP-num 1 $pep.eggnog.tab | cut -f 1,2,11,12 > $pep.eggnog.tab.seed_orthologs
emapper.py --annotate_hits_table $pep.eggnog.tab.seed_orthologs -o eggNOG --data_dir /filer-dg/agruppen/dg3/zhangh/software/database/emapperdb-5.0.2/
