pep=Hedin2.v2.transcript.pep.fa
###nr database
module load diamond
diamond  blastp -e 1e-5 --threads 8 --outfmt 6  -d /filer-dg/agruppen/dg6/zhangh/software/annotation/nr//Plants.fa.dmnd --query $pep -o $pep.blast.nr # It's better to split the pep to several parts and cat them togther.
perl blast_parser.pl  -tophit 1 $pep.blast.nr >$pep.blast.nr.xls

###nt database
/opt/Bio/blast/2.2.26/bin/blastall  -p tblastn -e 1e-5 -a 8 -F F -d /filer-dg/agruppen/dg6/zhangh/software/annotation/nt//Plants.fa -i $pep -o $pep.blast.nt
perl blast_parser.pl  -tophit 1 $pep.blast.nt >$pep.blast.nt.xls

