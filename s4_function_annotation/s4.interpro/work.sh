pep=Hedin2.v2.transcript.pep.fa
export PATH="/opt/Bio/jdk/1.6.0_45/bin/:$PATH"
module load interproscan
sh /opt/Bio/interproscan/5.57-90.0/interproscan.sh -goterms -f tsv  --appl ProDom --appl PRINTS --appl Pfam --appl SMART --appl PANTHER --appl ProSiteProfiles --appl ProSitePatterns  -T ./temp.file -i $pep -o $pep.iprscan
