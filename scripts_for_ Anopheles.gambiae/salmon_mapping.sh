export PATH=$HOME/soft/mafft-7.222-without-extensions/core:$PATH
export PATH=$HOME/soft/fastp/:$PATH
export PATH=$HOME/soft/:$PATH
export PATH=$HOME/FastQC/:$PATH
export PATH=$HOME/software/bbmap/:$PATH
export PATH=$HOME/soft/Salmon-latest_linux_x86_64/bin/:$PATH


cd ~/backup_fastq/

mkdir mosquito
cd mosquito/aaeg
ln -s  ../../10663_10644_95220_HCJNKBGXB_Bretta_BF_*.gz  .
rename 's/10663_10644_95220_HCJNKBGXB_Bretta_//'  *.fastq.gz

cd ../agam
ln -s  ../../10663_10644_95220_HCJNKBGXB_Bretta_Ant_*.gz  .
ln -s  ../../10663_10644_95220_HCJNKBGXB_Bretta_Post_*.gz  .
ln -s  ../../10663_10644_95220_HCJNKBGXB_Bretta_Prov_*.gz  .
ln -s  ../../10663_10644_95220_HCJNKBGXB_Bretta_WG_*.gz  .

rename 's/10663_10644_95220_HCJNKBGXB_Bretta_//'  *.fastq.gz


###download referecne from vector base
cd /media/eightteratwo/xiaoli/reference/mosquitos/Agam

wget -c https://www.vectorbase.org/download/anopheles-gambiae-pestchromosomesagamp4fagz
wget -c https://www.vectorbase.org/download/anopheles-gambiae-pestpeptidesagamp412fagz
wget -c https://www.vectorbase.org/download/anopheles-gambiae-pesttranscriptsagamp412fagz
wget -c https://www.vectorbase.org/download/anopheles-gambiae-pestbasefeaturesagamp412gff3gz
wget -c https://www.vectorbase.org/download/anopheles-gambiae-pestbasefeaturesagamp412gtfgz

rename 's/pest/pest\./'  anopheles* 
rename 's/agam/\.agam/'  anopheles* 
rename 's/412/412\./'  anopheles* 
rename 's/gz/\.gz/'  anopheles* 

gunzip -d  anopheles-gambiae-pest.transcripts.agamp412.fa.gz 


#build index of transcripts with salmon
mkdir salmon_index
cd salmon_index
#build a lightweight-alignment (FMD-based) index    use:  "--type fmd "
#quasi-mapping is the default index type in Salmon  ,  "--type quasi"
salmon index -t ../anopheles-gambiae-pest.transcripts.agamp412.fa  -i  agamp412_vector --keepDuplicates
# There were 361 transcripts that would need to be removed to avoid duplicates.
cd ..


#LET S TRIM THE SEQUENCES
#go to the folder of fastq.gz files. thses files are downloaded from biotech sequencing facility
cd ~/backup_fastq/mosquito/agam/

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
salmon quant -i /media/eightteratwo/xiaoli/reference/mosquitos/Agam/salmon_index/agamp412_vector -l A -r  ${f} -o ${f%%.*}.salmon_vector;
#extract salmon column of "NumReads" (the 5th column), put them into separate files
awk '{print $5}'  ${f%%.*}.salmon_vector/quant.sf  > salmon_counts/${f%%.*}.quant.simple.txt ; 
done

#extract salmon column of "Name", put them into a separate txt. this would be the rownames for all samples
awk '{print $1}' A5.salmon_vector/quant.sf > salmon_counts/00.quant.simple.txt ;

#we combine them all together, with "paste" commend. 
#make sure the soft don't mess up (disorder or miss) the gene order
cd salmon_counts
#show the name of samples
#you get two files: the .txt with the counts and the list.txt for sample list.
ls *.simple.txt > list.txt
paste *.quant.simple.txt > mosquito.Agam.salmon_vector.txt

#/media/eightteratwo/backup_fastq/bretta/infection_tissue_20190424/mosquito/agam/clean/salmon_counts


