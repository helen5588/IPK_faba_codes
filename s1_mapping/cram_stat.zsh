#!/bin/zsh


awkcmd='BEGIN{
  OFS = "\t"
}

!and($2, 2304){
  all++
   if(!and($2, 4)){
      mapped++
        if(!and($2, 1024)){
    nondup++
    if(!and($2, 12)){
       if($5 >= 10){
    q10++
       }
          if($5 >= 20){
       q20++
          }
     if($5 >= 30){
          q30++
     }
        if($5 >= 40){
     q40++
        }
   if($5 >= 50){
        q50++
   }
      if($5 >= 60){
           q60++
      }
         if(and($2, 2)){
      proper++
         }
   }
 }
    }
}

END {
  print "all",  0+all
   print "mapped", 0+mapped
   print "nondup",0+nondup 
   print "q10", 0+q10
     print "q20", 0+q20
      print "q30", 0+q30
       print "q40", 0+q40
        print "q50", 0+q50
 print "q60", 0+q60
}'
samtools='/opt/Bio/samtools/1.16.1/bin/samtools'
i=$1
base=$i/$i:t
bam=${base}.cram

stats=${base}_mapping_stats.tsv
statserr=${base}_mapping_stats.err
$samtools view $bam | awk "$awkcmd" > $stats 2> $statserr
