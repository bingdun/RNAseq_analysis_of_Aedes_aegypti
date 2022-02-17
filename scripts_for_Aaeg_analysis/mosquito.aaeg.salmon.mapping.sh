
cd /media/eightteratwo/backup_fastq/bretta/infection_tissue_20190424/mosquito/aaeg
#LET S TRIM THE SEQUENCES
#go to the folder of fastq.gz files. thses files are downloaded from biotech sequencing facility
cd ~/backup_fastq/

#trim fastq.gz to fastq
for sample in *.fastq.gz; do zcat ${sample} | /home/yoni/bbmap/bbduk.sh in=stdin.fq out=${sample}_trimmed_clean ref=/home/yoni/bbmap/resources/polyA.fa.gz,/home/yoni/bbmap/resources/truseq.fa.gz k=13 ktrim=r forcetrimleft=11 useshortkmers=t mink=5 qtrim=rl trimq=10 minlength=20; done  


#move the trimmed fastq to a new folder
mkdir clean
mv *_clean clean/
cd clean

#rename files. keep only the sample ID
ls *clean |awk -F "_" '{print "mv "$0" "$4"."$7".fq"}' | head
ls *clean |awk -F "_" '{print "mv "$0" "$4"."$7".fq"}' | bash 


#do mapping and count reads with salmon
mkdir salmon_counts

for f in *.fq; 
do echo  ${f}  ${f%%.*}.fq;
salmon quant -i /media/eightteratwo/xiaoli/reference/AaegL5.0/vector_aaegl/salmon_index/aaeg_vector -l A -r  ${f} -o ${f%%.*}.salmon_vector;
#extract salmon column of "NumReads" (the 5th column), put them into separate files
awk '{print $5}'  ${f%%.*}.salmon_vector/quant.sf  > salmon_counts/${f%%.*}.quant.simple.txt ; 
done

#extract salmon column of "Name", put them into a separate txt. this would be the rownames for all samples
awk '{print $1}' C4.salmon_vector/quant.sf > salmon_counts/00.quant.simple.txt ;

#we combine them all together, with "paste" commend. 
#make sure the soft don't mess up (disorder or miss) the gene order
cd salmon_counts
#show the name of samples
#you get two files: the .txt with the counts and the list.txt for sample list.
ls *.simple.txt > list.txt
paste *.quant.simple.txt > mosquito.Aaeg.salmon_vector.txt

