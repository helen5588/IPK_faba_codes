pep=Hedin2.v2.transcript.pep.fa
/opt/Bio/blast/2.2.26/bin/blastall  -p blastp -e 1e-5 -a 5 -m 8 -F F -d /filer-dg/agruppen/dg6/zhangh/software/annotation/DNA_annotation/database/uniprot/release-2015_12/swissport/uniprot_sprot.Eukaryota.fasta.simple -i $pep -o $pep.blast.swissprot
perl get_annot_info.pl -tophit 1 -topmatch 1 -id /filer-dg/agruppen/dg6/zhangh/software/annotation/DNA_annotation/database/uniprot/release-2015_12/swissport/uniprot_sprot.Eukaryota.id.annot.xls -input $pep -out $pep.xls
