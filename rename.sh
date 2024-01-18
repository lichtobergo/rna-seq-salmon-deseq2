cd data/fastq
pwd

for i in *.fastq.gz
  do
  mv $i $(echo $i | awk '{split($1,a,/_/); print a[1]"_"a[5]"_"a[6]}')
  done

# remove leading project number from fastq filename
rename 's/p1431s//g' *.fastq.gz
#rename 's/Mblood|Mctx|Rblood|Rctx//g' *farename 's/Mblood|Mctx|Rblood|Rctx//g' *fastq.gzstq.gz
rename 's/_001//g' *.fastq.gz
