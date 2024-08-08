#!/bin/awk -f

BEGIN{
 OFS=FS="\t"
 if(header){
  print "##fileformat=VCFv4.2"
  print "##FILTER=<ID=PASS,Description=\"All filters passed\">"
  print "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Average mapping quality\">"
  print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
  print "##FORMAT=<ID=DP,Number=.,Type=Integer,Description=\"Read depth\">"
  print "##FORMAT=<ID=DV,Number=.,Type=Integer,Description=\"Read depth of the alternative allele\">"
 } else {
  h["0/0"]=0
  h["0/1"]=1
  h["1/1"]=2
  if(!dphommax)
   dphommax=100
  if(!dphom)
   dphom=1
  if(!dphet)
   dphet=1
  if(!tol)
   tol=0.2499
  while(getline < bed){
   chr[$1]=$4
   offset[$1]=$5
  }
 }
}

header && /^#(#reference=|CHROM)/

header && /^[^#]/ {
 exit
}

/^#/ || $6 < minqual || $4 == "N" || 0+gensub(".*DP=([0-9]+).*", "\\1", "g", $8) < mindp {
 next
}

{
 o=""
 n=A=B=H=0

 for(i = 10; i <= NF; i++){
  split($i, a, ":")
  dp = 0+a[3]
  dv = 0+a[4]
  if(dp > 0)
   r=dv/dp
  else if($5 == "."){
   o = o"\t0/0:"dp":"dv
   A++
  }
  if(dp >= dphom && dp <=dphommax && r <= tol ){
   o = o"\t0/0:"dp":"dv
   A++
  }
  else if(dp >= dphom && dp <=dphommax && r >= 1-tol){
   o = o"\t1/1:"dp":"dv
   B++
  }
  else if(dp >= dphet && dp <=dphommax && r >= 0.5-tol && r <= 0.5+tol){
   o = o"\t0/1:"dp":"dv
   H++
  }
  else{
   o = o"\t./.:"dp":"dv
   n++
  }
 }
 present = A + B + H
}

!present || A < minhomn || B < minhomn || present < minpresent * (present + n) || A + B < minhomp * present {
 next
}

{
 if(B > A)
  m = (2*A + H) / 2 / present
 else
  m = (2*B + H) / 2 / present
}

m >= minmaf {
 mq = gensub(".*MQ=([0-9]+).*", "\\1", "g", $8)
 print chr[$1], offset[$1]+$2, $3, $4, $5, $6, $7, "MQ="mq, "GT:DP:DV"o
}
