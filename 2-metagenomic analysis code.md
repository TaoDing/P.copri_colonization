

[TOC]


## 1 Quality control and removal of human host contamination using kneaddata software



~~~ bash
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

gzip -d GRCh38_latest_genomic.fna.gz 



bowtie2-build GRCh38_latest_genomic.fna human_genome 



parallel -j 3 --xapply 'echo {1} {2}'    ::: clean.data/*1.clean.fq.gz ::: clean.data/*2.clean.fq.gz


parallel -j 3 --xapply 'kneaddata -i {1} -i {2} -o kneaddata_out -v \
-db /public1/chenrz/lihui.tiv/metagenomic_108/human_genome  \
--trimmomatic /public/apps/anaconda2/bin/Trimmomatic-0.38/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
-t 9 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output' \
 ::: clean.data/*1.clean.fq.gz ::: clean.data/*2.clean.fq.gz
 
kneaddata_read_count_table --input kneaddata_out --output kneaddata_read_counts.txt


mkdir kneaddata_out/contam_seq
mv kneaddata_out/*_contam*fastq kneaddata_out/contam_seq


mkdir cat_reads
for i in `tail -n+2 design-10.txt | cut -f 1`;do \
cat kneaddata_out/${i}_clean_kneaddata_paired_1.fastq kneaddata_out/${i}_clean_kneaddata_paired_2.fastq kneaddata_out/${i}_clean_kneaddata_unmatched_1.fastq kneaddata_out/${i}_clean_kneaddata_unmatched_2.fastq | awk '{if(NR%4==1) print "@"NR; else print $0}' > cat_reads/${i}.fastq; done

~~~



## 2 Read-based species and functional annotation using humann2 software



~~~ bash

parallel -j 4 'humann2 --threads 32 --nucleotide-database /public/home/chenrz/References/HUMAnN2_database/chocophlan/ --protein-database /public/home/chenrz/References/HUMAnN2_database/uniref --input {} --output humann2_out/{/.}' ::: cat_reads/*fastq




mkdir humann2_final_out

humann2_join_tables -s --input humann2_out/ --file_name pathabundance --output humann2_final_out/humann2_pathabundance.tsv
humann2_join_tables -s --input humann2_out/ --file_name pathcoverage --output humann2_final_out/humann2_pathcoverage.tsv
humann2_join_tables -s --input humann2_out/ --file_name genefamilies --output humann2_final_out/humann2_genefamilies.tsv


humann2_renorm_table --input humann2_final_out/humann2_pathabundance.tsv --units relab --output humann2_final_out/humann2_pathabundance_relab.tsv
 

humann2_renorm_table --input humann2_final_out/humann2_genefamilies.tsv --units relab --output humann2_final_out/humann2_genefamilies_relab.tsv



humann2_split_stratified_table --input humann2_final_out/humann2_pathabundance_relab.tsv --output humann2_final_out
humann2_split_stratified_table --input humann2_final_out/humann2_genefamilies_relab.tsv --output humann2_final_out
humann2_split_stratified_table --input humann2_final_out/humann2_pathcoverage.tsv --output humann2_final_out




hclust2.py -i humann2_final_out/humann2_genefamilies_relab_unstratified.tsv -o humann2_final_out/humann2_genefamilies_relab_unstratified.png \
--ftop 25 --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 \
--slabel_size 6 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300 -c Blues




humann2_renorm_table --input humann2_final_out/humann2_pathabundance.tsv --units cpm --output humann2_final_out/humann2_pathabundance_cpm.tsv
 

humann2_renorm_table --input humann2_final_out/humann2_genefamilies.tsv --units cpm --output humann2_final_out/humann2_genefamilies_cpm.tsv


humann2_split_stratified_table --input humann2_final_out/humann2_pathabundance_cpm.tsv --output humann2_final_out

humann2_split_stratified_table --input humann2_final_out/humann2_genefamilies_cpm.tsv --output humann2_final_out


# Regrouping genes to other functional categories
cd humann2_final_out
mkdir -p KO
cd ../

# Regrouping genes to ko
humann2_regroup_table --input humann2_final_out/humann2_genefamilies_cpm.tsv --output humann2_final_out/KO/humann2_genefamilies_cpm_ko.tsv -c /public/home/chenrz/References/HUMAnN2_database/humann-master/humann2/data/misc/map_ko_uniref90.txt


# Rename ko to ko_name
/public/home/chenrz/References/HUMAnN2_database/humann-master/humann2/tools/rename_table.py \
--input humann2_final_out/KO/humann2_genefamilies_cpm_ko.tsv \
-n kegg-orthology \
--output humann2_final_out/KO/humann2_genefamilies_cpm_ko_name.tsv


humann2_split_stratified_table --input humann2_final_out/KO/humann2_genefamilies_cpm_ko_name.tsv \
--output humann2_final_out/KO

# Regrouping genes to EC
cd humann2_final_out
mkdir EC
cd ../
humann2_regroup_table --input humann2_final_out/humann2_genefamilies_cpm.tsv --output humann2_final_out/EC/humann2_cpm_ec.tsv -c /public/home/chenrz/References/HUMAnN2_database/humann-master/humann2/data/misc/map_level4ec_uniref90.txt.gz


# Rename ko to EC_name
/public/home/chenrz/References/HUMAnN2_database/humann-master/humann2/tools/rename_table.py \
--input humann2_final_out/EC/humann2_cpm_ec.tsv \
-n ec \
--output humann2_final_out/EC/humann2_kegg_ec_name.tsv


humann2_split_stratified_table --input humann2_final_out/EC/humann2_kegg_ec_name.tsv --output humann2_final_out/EC


~~~

## 3 Generating species abundance table using metaphlan2



~~~bash
mkdir metaphlan2_out
cp humann2_out/*/*/*metaphlan_bugs_list.tsv metaphlan2_out/

merge_metaphlan_tables.py metaphlan2_out/*metaphlan_bugs_list.tsv > metaphlan2_out/metaphlan2_merged.txt
sed -i 's/_metaphlan_bugs_list//g' metaphlan2_out/metaphlan2_merged.txt


grep -E "(s__)|(^ID)" metaphlan2_out/metaphlan2_merged.txt | grep -v "t__" | sed 's/^.*s__//g' > metaphlan2_out/metaphlan2_merged_species.txt




grep -E "(g__)|(^ID)" metaphlan2_out/metaphlan2_merged.txt | grep -v "s__" | sed 's/^.*s__//g' > metaphlan2_out/metaphlan2_merged_genus.txt



grep -E "(k__)|(^ID)" metaphlan2_out/metaphlan2_merged.txt | grep -v "p__"  > metaphlan2_out/metaphlan2_merged_kingdom.txt


~~~

## 4 Generating species abundance table using kraken2

~~~bash


kraken2-build  --threads 56 --db ./ --download-taxonomy


kraken2-build  --threads 56 --db ./  --download-library bacteria
kraken2-build  --threads 56 --db ./  --download-library viral
kraken2-build  --threads 56 --db ./  --download-library archaea
kraken2-build  --threads 56 --db ./  --download-library plasmid
kraken2-build  --threads 56 --db ./  --download-library human
kraken2-build  --threads 56 --db ./  --download-library fungi


bracken-build -d /public/home/chenrz/References/kraken2_database/db -t 8 -k 35 -l 150

cd /public1/chenrz/lihui.tiv/metagenomic_107

mkdir kraken2_out

awk '{print $1}' ./design.txt >seq_ID 
sed -i '1d' seq_ID


for i in `tail -n+2 design.txt | cut -f 1`; 
do
   /public/home/chenrz/app/kraken2\
    --db /public/home/chenrz/References/kraken2_database/db/ \
    --threads 32 \
    --report kraken2_out/${i}.report \
	--output kraken2_out/${i}.out  \
 --gzip-compressed   --paired \
clean.data/${i}_R1.fq.gz clean.data/${i}_R2.fq.gz
done




for i in `tail -n+2 design.txt | cut -f 1`; 
do
/public/home/chenrz/app/Bracken-2.5/bracken  -d /public/home/chenrz/References/kraken2_database/db/ -i kraken2_out/${i}.report -o kraken2_out/${i}_bracken.output -r 150 -l S -t32
done


/public/home/chenrz/app/Bracken-2.5/bracken  -d /public/home/chenrz/References/kraken2_database/db/ -i kraken2_out/SeqCRZa_108.report -o kraken2_out/SeqCRZa_108_bracken.output -r 150 -l S -t32



for i in `tail -n+2 design.txt | cut -f 1`;
do
/public/home/chenrz/References/kraken2_database/kraken2-2.0.8-beta/KrakenTools-master/kreport2mpa.py -r kraken2_out/${i}_bracken.report   -o kraken2_out/${i}_bracken.final.report --display-header
done


for i in `tail -n+2 design.txt | cut -f 1`;
do
/public/home/chenrz/References/kraken2_database/kraken2-2.0.8-beta/KrakenTools-master/kreport2mpa.py -r kraken2_out/${i}_bracken.report   -o kraken2_out/${i}_bracken.final.abundance --display-header --percentages
done





/public/home/chenrz/References/kraken2_database/kraken2-2.0.8-beta/KrakenTools-master/combine_mpa.py -i kraken2_out/*_bracken.final.report -o kraken2_out/combined_mpa_report


/public/home/chenrz/References/kraken2_database/kraken2-2.0.8-beta/KrakenTools-master/combine_mpa.py -i kraken2_out/*_bracken.final.abundance -o kraken2_out/combined_mpa_abundance



grep -E "(k__)|(^ID)" kraken2_out/combined_mpa_report | grep -v "p__" | grep -v "f__" | grep -v "g__" | grep -v "o__" | grep -v "c__" | grep -v "s__" | grep -v "t__"> kraken2_out/kraken2_combined_kingdom.txt




grep -E "(s__)|(^ID)" kraken2_out/combined_mpa_report > kraken2_out/kraken2_combined_stratified_species.txt



grep -E "(s__)|(^ID)" kraken2_out/combined_mpa_report  | grep -v "t__" | sed 's/^.*s__//g' > kraken2_out/kraken2_combined_unstratified_species.txt 



grep -E "(g__)|(^ID)" kraken2_out/combined_mpa_report | grep -v "s__" | sed 's/^.*g__//g' > kraken2_out/kraken2_combined_stratified_genus.txt


~~~



## 5 Denovo assembly using megahit software

   

```bash

cd /public1/chenrz/lihui.tiv/metagenomic_107/
mkdir megahit_out


for i in `tail -n+2 design.txt | cut -f 1`; 
do
megahit -1 kneaddata_out/${i}_R1_kneaddata_paired_1.fastq.gz \
-2 kneaddata_out/${i}_R1_kneaddata_paired_2.fastq.gz \
-o megahit_out/${i}  \
--out-prefix ${i}.megahit \
--min-contig-len 1000 \
-m 0.9 -t 32 
done

# Assessing assembly efficiency using quast
quast.py megahit_out/*/*megahit.contigs.fa -o megahit_out/megahit_report -t 32

cd /public1/chenrz/lihui.tiv/metagenomic_107/

mkdir megahit_1000_fasta



for i in `tail -n+2 design.txt | cut -f 1`; 
do
awk '/^>/ {printf("\n%s\t",$0);next;} {printf("%s",$0);} END {printf("\n");}' < megahit_out/${i}/${i}.megahit.contigs.fa | egrep -v '^$' | tr "\t" "\n" | awk '!/^>/ { next } { getline seq } length(seq) > 1000 { print $0 "\n" seq }' > megahit_out/megahit_1000_fasta/${i}.megahit.1000.contigs.fa
done


```



## 6 Gene function annotation using Prokka

~~~bash


for i in `tail -n+2 design.txt | cut -f 1`; 
do
prokka megahit_out/${i}/${i}.megahit.contigs.fa \
--outdir prokka_out/${i}/  \
--prefix ${i} \
--mincontiglen 1000 \
--cpus 46
done
~~~

 

##  7 Non-redundant gene set construction using CD-HIT

~~~bash

mkdir -p prokka_out/CDHIT_out/
cat prokka_out/*/*.faa >prokka_out/CDHIT_out/all.faa

conda activate cd-hit
cd-hit -i prokka_out/CDHIT_out/all.faa -o prokka_out/CDHIT_out/all.cdhit.faa -T 32 -M 0


mkdir -p prokka_out/CDHIT_out/
cat prokka_out/*/*.ffn >prokka_out/CDHIT_out/all.ffn

conda activate cd-hit
cd-hit -i prokka_out/CDHIT_out/all.ffn -o prokka_out/CDHIT_out/all.cdhit.ffn -T 32 -M 0

cd-hit-est -i prokka_out/CDHIT_out/all.ffn -o prokka_out/CDHIT_out/all.cdhit-est.ffn -T 32 -M 0 

~~~



## 8 Rapid estimation of gene abundance without alignment --Salmon

~~~bash

cd /public/home/chenrz/app/
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.8.2/Salmon-0.8.2_linux_x86_64.tar.gz
tar xvfz Salmon-0.8.2_linux_x86_64.tar.gz

vim  ~/.bashrc
export PATH=/public/home/chenrz/app/Salmon-0.8.2_linux_x86_64/bin:$PATH
source ~/.bashrc


cd /public1/chenrz/lihui.tiv/metagenomic_107
mkdir salmon_out
cd salmon_out



ln -fs /public1/chenrz/lihui.tiv/metagenomic_107/prokka_out/CDHIT_out/all.cdhit-est.ffn .
ln -fs /public1/chenrz/lihui.tiv/metagenomic_107/kneaddata_out/*paired_1.fastq.gz .
ln -fs /public1/chenrz/lihui.tiv/metagenomic_107/kneaddata_out/*paired_2.fastq.gz .

rename _R1_kneaddata_paired_1. _R1. * 
rename _R1_kneaddata_paired_2. _R2. *
rename fastq.gz fq.gz *
rename _R1.fq.gz .R1.fq.gz *
rename _R2.fq.gz .R2.fq.gz *

salmon index -t all.cdhit-est.ffn -i transcript_index --type quasi -k 31 -p 32




mkdir quant_out

for i in `tail -n+2 design.txt | cut -f 1`; 
do
salmon quant -i transcript_index --libType IU \
-1 $i.R1.fq.gz -2 $i.R2.fq.gz -o quant_out/$i.quant -p 24;
done


curl -L -O https://raw.githubusercontent.com/ngs-docs/2016-metagenomics-sio/master/gather-counts.py

chmod +x gather-counts.py

./gather-counts.py  



for file in *counts
do
  name=${file%%.*}
  sed -e "s/count/$name/g" $file > tmp
  mv tmp $file
done


for i in `tail -n+2 design.txt | cut -f 1`; 
do
awk '{print $2}' $i.quant.counts > $i.gene_counts
done

awk '{print $1}' SeqCRZa_1.quant.counts >gene.names

paste gene.names *gene_counts >combined.gene_counts

~~~



## 9 Strain level analysis using panphlan3.0

~~~bash
#1 install panphlan3.0
conda create --name panphlan
conda activate panphlan
conda install -c bioconda panphlan 


#2 Download the reference genomes

cd /public/home/chenrz/References/panphlan3.0_database
conda activate panphlan
panphlan_download_pangenome.py -i Prevotella_copri 

tar -xvjf Prevotella_copri.tar.bz2


#3  Map samples against pangenome

cd /public1/chenrz/lihui.tiv/metagenomic_107
conda activate panphlan
mkdir -p panphlan_out/map_results


for i in `tail -n+2 design.txt | cut -f 1`; 
do
panphlan_map.py -i cat_reads/$i.fastq --indexes /public/home/chenrz/References/panphlan3.0_database/Prevotella_copri/Prevotella_copri -p /public/home/chenrz/References/panphlan3.0_database/Prevotella_copri/Prevotella_copri_pangenome.tsv -o panphlan_out/map_results/$i.copri.tsv --nproc 12
done


#4  Profiling strains


panphlan_profiling.py -i panphlan_out/map_results/  --o_matrix panphlan_out/result_profile_copri.tsv --min_coverage 1 --left_max 1.70 --right_min 0.30 -p /public/home/chenrz/References/panphlan3.0_database/Prevotella_copri/Prevotella_copri_pangenome.tsv --add_ref  --func_annot /public/home/chenrz/References/panphlan3.0_database/Prevotella_copri/panphlan_Prevotella_copri_annot.tsv --field 4  --verbose

~~~



## 10 Assembling metagenomic data using MetaWRAP

~~~bash

cd /public1/chenrz/lihui.tiv/metagenomic_107

for i in `head -n 29  design.txt | tail -n +2 | cut -f 1`; 
do
cp kneaddata_out/${i}_R1_kneaddata*.fastq.gz /public5/chenrz/metagenomic_107_feces/metawrap_out/tibetan_2Y/
done

for i in `head -n 60  design.txt | tail -n +30 | cut -f 1`; 
do
cp kneaddata_out/${i}_R1_kneaddata*.fastq.gz /public5/chenrz/metagenomic_107_feces/metawrap_out/tibetan_1Y/
done


for i in `head -n 82  design.txt | tail -n +61 | cut -f 1`; 
do
cp kneaddata_out/${i}_R1_kneaddata*.fastq.gz  /public5/chenrz/metagenomic_107_feces/metawrap_out/tibetan_1M/
done


for i in `head -n 108  design.txt | tail -n +83 | cut -f 1`;
do
cp kneaddata_out/${i}_R1_kneaddata*.fastq.gz  /public5/chenrz/metagenomic_107_feces/metawrap_out/han/
done


cd /public5/chenrz/metagenomic_107_feces/metawrap_out/

cat tibetan_1M/*_1.fastq.gz >tibetan_1M.R1.fastq.gz
cat tibetan_1Y/*_1.fastq.gz >tibetan_1Y.R1.fastq.gz
cat tibetan_2Y/*_1.fastq.gz >tibetan_2Y.R1.fastq.gz
cat han/*_1.fastq.gz >han.R1.fastq.gz


cat tibetan_1M/*_2.fastq.gz >tibetan_1M.R2.fastq.gz
cat tibetan_1Y/*_2.fastq.gz >tibetan_1Y.R2.fastq.gz
cat tibetan_2Y/*_2.fastq.gz >tibetan_2Y.R2.fastq.gz
cat han/*_2.fastq.gz >han.R2.fastq.gz



source activate metawrap-env
cd /public5/chenrz/metagenomic_107_feces/metawrap_out/


metawrap assembly -1 tibetan_1M.R1.fastq.gz -2 tibetan_1M.R2.fastq.gz -m 320 -t 32 -o ./tibetan_1M_metawrap_assembly/ -l 1000


cd /public5/chenrz/metagenomic_107_feces/metawrap_out/fastq
mkdir han_pair
for i in `tail -n+1 han/han.txt | cut -f 1`; 
do
seqkit sort -n han/${i}.fastq.gz > han_pair/${i}.fastq.gz
done


mkdir tibetan_1M_pair
for i in `tail -n+1 tibetan_1M/tibetan_1M.txt | cut -f 1`; 
do
seqkit sort -n tibetan_1M/${i}.fastq.gz > tibetan_1M_pair/${i}.fastq.gz
done

mkdir tibetan_1Y_pair
for i in `tail -n+1 tibetan_1Y/tibetan_1Y.txt | cut -f 1`; 
do
seqkit sort -n tibetan_1Y/${i}.fastq.gz > tibetan_1Y_pair/${i}.fastq.gz
done

mkdir tibetan_2Y_pair
for i in `tail -n+1 tibetan_2Y/tibetan_2Y.txt | cut -f 1`; 
do
seqkit sort -n tibetan_2Y/${i}.fastq.gz > tibetan_2Y_pair/${i}.fastq.gz
done
~~~

###  10.1 Use concoct, metabat2 and maxbin2 for binning

~~~


cd /public5/chenrz/metagenomic_107_feces/metawrap_out/metawrap_assembly_out/
source activate metawrap-env

metawrap binning -a ./tibetan_1M_metawrap_assembly/final_assembly.fasta -t 16 --metabat2 --maxbin2 --concoct -o ./tibetan_1M_metawrap_binning /public5/chenrz/metagenomic_107_feces/metawrap_out/fastq/tibetan_1M_pair/*fastq.gz --run-checkm


metawrap binning -a ./tibetan_1Y_metawrap_assembly/final_assembly.fasta -t 16 --metabat2 --maxbin2 --concoct -o ./tibetan_1Y_metawrap_binning /public5/chenrz/metagenomic_107_feces/metawrap_out/fastq/tibetan_1Y_pair/*fastq.gz --run-checkm

metawrap binning -a ./tibetan_2Y_metawrap_assembly/final_assembly.fasta -t 16 --metabat2 --maxbin2 --concoct -o ./tibetan_2Y_metawrap_binning /public5/chenrz/metagenomic_107_feces/metawrap_out/fastq/tibetan_2Y_pair/*fastq.gz --run-checkm

metawrap binning -a ./tibetan_han_metawrap_assembly/final_assembly.fasta -t 16 --metabat2 --maxbin2 --concoct -o ./tibetan_han_metawrap_binning /public5/chenrz/metagenomic_107_feces/metawrap_out/fastq/tibetan_han_pair/*fastq.gz --run-checkm
~~~

### 10 .2 Bin refinement

~~~bash
cd /public5/chenrz/metagenomic_107_feces/metawrap_out/metawrap_assembly_out/
source activate metawrap-env

metawrap bin_refinement -o ./tibetan_1M_metawrap_binning/metawrap_70_5_bins -t 36 -m 320 -A ./tibetan_1M_metawrap_binning/metabat2_bins/ -B ./tibetan_1M_metawrap_binning/maxbin2_bins/ -C ./tibetan_1M_metawrap_binning/concoct_bins/ -c 70 -x 5

metawrap bin_refinement -o ./tibetan_1Y_metawrap_binning/metawrap_70_5_bins -t 36 -m 320 -A ./tibetan_1Y_metawrap_binning/metabat2_bins/ -B ./tibetan_1Y_metawrap_binning/maxbin2_bins/ -C ./tibetan_1Y_metawrap_binning/concoct_bins/ -c 70 -x 5


metawrap bin_refinement -o ./tibetan_2Y_metawrap_binning/metawrap_70_5_bins -t 36 -m 320 -A ./tibetan_2Y_metawrap_binning/metabat2_bins/ -B ./tibetan_2Y_metawrap_binning/maxbin2_bins/ -C ./tibetan_2Y_metawrap_binning/concoct_bins/ -c 70 -x 5


metawrap bin_refinement -o ./han_metawrap_binning/metawrap_70_5_bins -t 36 -m 320 -A ./han_metawrap_binning/metabat2_bins/ -B ./han_metawrap_binning/maxbin2_bins/ -C ./han_metawrap_binning/concoct_bins/ -c 70 -x 5
~~~



### 10.3 Bins species annotation (based on GTDB database)

~~~bash

wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar zxvf gtdbtk_data.tar.gz


cp /public/home/shuling/.conda/envs/gtdb/gtdbtk_r95_data.tar.gz /public/home/chenrz/.conda/envs/gtdbtk/gtdbtk_r95_data.tar.gz

cd /public/home/chenrz/.conda/envs/gtdbtk/
tar zxvf gtdbtk_r95_data.tar.gz


cd /public/home/chenrz/.conda/envs/gtdbtk/etc/conda/activate.d

vi gtdbtk.sh
/public/home/chenrz/.conda/envs/gtdbtk/gtdb_r95/release95


cd /public5/chenrz/metagenomic_107_feces/metawrap_out/metawrap_assembly_out/


source activate gtdbtk



gtdbtk classify_wf --genome_dir ./han_metawrap_binning/metawrap_70_5_bins/metawrap_70_5_bins --extension fa --out_dir ./han_metawrap_gtdb --cpus 32


gtdbtk classify_wf --genome_dir ./tibetan_1M_metawrap_binning/metawrap_70_5_bins/metawrap_70_5_bins --extension fa --out_dir ./tibetan_1M_metawrap_gtdb --cpus 32

gtdbtk classify_wf --genome_dir ./tibetan_1Y_metawrap_binning/metawrap_70_5_bins/metawrap_70_5_bins --extension fa --out_dir ./tibetan_1Y_metawrap_gtdb --cpus 32

gtdbtk classify_wf --genome_dir ./tibetan_2Y_metawrap_binning/metawrap_70_5_bins/metawrap_70_5_bins --extension fa --out_dir ./tibetan_2Y_metawrap_gtdb --cpus 32



mkdir kraken2_out
/public/home/chenrz/app/kraken2\
    --db /public/home/chenrz/References/kraken2_database/db/ \
    --threads 32 \
    --report kraken2_out/bin.1_han.report \
	--output kraken2_out/bin.1_han.out  \
	all.bins/bin.1_han.fa

~~~

### 10.4 Bins functional gene annotation-prokka



~~~bash

cd /public5/chenrz/metagenomic_107_feces/metawrap_out/metawrap_assembly_out/
mkdir all.bins

cp tibetan_2Y_metawrap_binning/metawrap_70_5_bins/metawrap_70_5_bins/* all.bins
cp tibetan_1Y_metawrap_binning/metawrap_70_5_bins/metawrap_70_5_bins/* all.bins
cp tibetan_1M_metawrap_binning/metawrap_70_5_bins/metawrap_70_5_bins/* all.bins
cp han_metawrap_binning/metawrap_70_5_bins/metawrap_70_5_bins/* all.bins



source activate metawrap_root 

metaWRAP annotate_bins -o ./all.bins_annot -t 36 -b all.bins/



source activate gtdbtk
gtdbtk classify_wf --genome_dir ./all.bins/ --extension fa --out_dir ./all.bins_gtdb --cpus 50


conda activate phylophlan3

phylophlan -i all.bins_annot/bin_translated_genes/ \
 -d phylophlan \
 --databases_folder /public/home/chenrz/.conda/envs/phylophlan3 \
 --diversity low \
 -t a --nproc 60 \
 -o ./all.bins.phylophlan  \
 -f /public/home/chenrz/.conda/envs/phylophlan3/lib/python3.8/site-packages/phylophlan/phylophlan_configs/supermatrix_aa.cfg

~~~

### 10.5 Build an evolutionary tree for all Prevotella 

~~~bash


grep 'Prevotella' gtdbtk.bac120.summary.tsv >bins.prevotella
awk '{print $1}' bins.prevotella >bins.prevotella_ID


cd /public5/chenrz/metagenomic_107_feces/metawrap_out/metawrap_assembly_out/Prevotella
mkdir -p kraken2_out
for i in `tail -n+1 bins.prevotella_ID | cut -f 1`; 
do
   /public/home/chenrz/app/kraken2\
    --db /public/home/chenrz/References/kraken2_database/db/ \
    --threads 32 \
    --report kraken2_out/${i}.report \
	--output kraken2_out/${i}.out  \
  ../all.bins/${i}.fa
done



cd /public5/chenrz/metagenomic_107_feces/metawrap_out/metawrap_assembly_out/Prevotella
for i in `tail -n+1 bins.prevotella_ID | cut -f 1`; 
do
cp -r ../all.bins_annot/bin_translated_genes/${i}.faa faa/
done


cd /public/home/chenrz/References/copri_genome/MAGS_clades_40_test


source activate metawrap_root 
metaWRAP annotate_bins -o ./clades.bins_annot -t 36 -b clades.bins/


cp /public/home/chenrz/References/copri_genome/MAGS_clades_40_test/clades.bins_annot/bin_translated_genes/* \
/public5/chenrz/metagenomic_107_feces/metawrap_out/metawrap_assembly_out/Prevotella/faa/

cd /public/home/chenrz/References/copri_genome/MAGS_clades
source activate phylophlan3

phylophlan -i clades.bins_annot/clades.faa/ \
 -d phylophlan \
 --databases_folder /public/home/chenrz/.conda/envs/phylophlan3 \
 --diversity low \
 -t a --nproc 60 \
 -o ./clades.bins.phylophlan  \
 -f /public/home/chenrz/.conda/envs/phylophlan3/lib/python3.8/site-packages/phylophlan/phylophlan_configs/supermatrix_aa.cfg

~~~



### 10.6 Abundance quantification of Prevotella copri clades in samples

~~~bash
cd /public5/chenrz/metagenomic_107_feces/metawrap_out/metawrap_assembly_out/Prevotella
source activate metawrap_root

metawrap quant_bins -t 54 \
-o prevotella.clades.quant/ \
-b /public/home/chenrz/References/copri_genome/MAGS_clades/clades_bins/  \
/public5/chenrz/metagenomic_107_feces/metawrap_out/fastq/all_feces_pair/*.fastq

~~~



###  10.7 Detection of single nucleotide polymorphisms (SNPs) in Prevotella copri clades



~~~bash

conda create -n snippy
source activate snippy
conda install -c conda-forge -c bioconda -c defaults snippy


# Calling SNPs
# Input Requirements
# a reference genome in FASTA or GENBANK format (can be in multiple contigs)
# sequence read file(s) in FASTQ or FASTA format (can be .gz compressed) format
# a folder to put the results in

cd /public5/chenrz/metagenomic_107_feces/metawrap_out/metawrap_assembly_out/

snippy --cpus 16 --outdir mysnps --ref GCF_000157935.1_ASM15793v1_p.copri18205_genomic.fna --ctgs input.faa/*.faa

~~~

































