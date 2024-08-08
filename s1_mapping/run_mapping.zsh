#!/bin/zsh
adapter=$1
shift
minlen=$1
shift
novosort_threads=$1
shift
threads=$1
shift
mem=$1
shift
tmp=$1
shift
ref=$1
shift
mmi=$1
shift
i=$1

minimap='/opt/Bio/minimap2/2.11/bin/minimap2'
samtools='/opt/Bio/samtools/1.8/bin/samtools'
novosort="/opt/Bio/novocraft/V3.06.05/bin/novosort"
cutadapt='/opt/Bio/cutadapt/1.15/bin/cutadapt'

indexsize='50G'
batchmem='5G'

base=$i/$i:t
#bam=${base}.bam
cram=${base}.cram

minimaperr=${base}_minimap.err
samtoolserr=${base}_samtools.err
sorterr=${base}_novosort.err
cutadapterr=${base}_cutadapt.err
cramCoverr=${base}_cramconversion.err

b=$base:t
rgentry="@RG\tID:$b\tPL:ILLUMINA\tPU:$b\tSM:$b"

{
  $cutadapt -f fastq --interleaved -a $adapter -A $adapter -m $minlen \
   <(find $i | egrep '1.fq.gz$|_R1' | grep 'f.*q.gz$' | sort | xargs gzip -cd) \
    <(find $i | egrep '2.fq.gz$|_R2' | grep 'f.*q.gz$' | sort | xargs gzip -cd) \
     2> $cutadapterr \
      | $minimap -ax sr -R $rgentry -t $threads -2 -I $indexsize -K$batchmem -a $mmi /dev/stdin 2> $minimaperr \
       | $samtools view -Su /dev/stdin 2> $samtoolserr \
        | $novosort -c $novosort_threads -t $tmp -m $mem --keepTags --md /dev/stdin 2> $sorterr \
	 | $samtools view -@ $threads -CT $ref /dev/stdin -o $cram && $samtools index -c -@ $threads $cram 2> $cramCoverr

	  echo $pipestatus > ${base}_pipestatus.txt 
}
