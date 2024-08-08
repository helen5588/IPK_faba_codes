dir=$1

{
 echo -n -e "#sample\tall_reads\ttrimmed_reads\ttrimmed_base\treads_in_bam\tmapped_reads"
  echo -n -e "\tnon_dup\tq10\tq20\tq30\tq40\tq50\tq60\n"

   find -L $dir -mindepth 1 \
      | grep -E 'bam$|cram$' | grep -v merge | sed 's/.bam$//'|sed 's/.cram$//' | while read base; do
     echo -n -e ${base##*/}"\t"
       cat ${base}_cutadapt.err | grep -m 1 "^Total read pairs processed:" \
   | awk -F: '{printf $2"\t"}' | tr -d ', '|awk -F: '{printf $1*2"\t"}'
	  cat ${base}_cutadapt.err | grep -m1 'Pairs written (passing filters):' \
     | awk -F: '{print $2}' | tr -d ', ' | cut -d '(' -f 1 | tr '\n' '\t'|awk -F: '{printf $1*2"\t"}'
	       cat ${base}_cutadapt.err | grep -m1 'Total written (filtered):' \
    | awk -F: '{print $2}' | tr -d ', ' |sed 's/bp//g'| cut -d '(' -f 1 | tr '\n' '\t'
   cut -f 2 ${base}_mapping_stats.tsv | xargs | tr '  ' '\t'
    done | sort
   } | awk '!p || !/^#/; /^#/ {p=1}' > 240116_$dir\_mapping_stats.tsv

