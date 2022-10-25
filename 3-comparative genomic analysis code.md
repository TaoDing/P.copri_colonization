# Comparative genomics analysis



## 1 Bacteroides 

It is necessary to download the genomes of all Bacteroides from NCBI Refseq database

~~~bash
cd /public5/chenrz/PC_BO/Bacteroides.all


# 1 Get genome metadata

datasets summary genome taxon 'Bacteroides' >Bacteroides.genome

# 2 Download a genome data package

datasets download genome taxon 'Bacteroides' --filename Bacteroides_dataset.zip

cd Bacteroides.all.genome
ls > ../Bacteroides.id
cd .. 
vim Bacteroides.id
:%s/.fna//g


for i in `tail -n+1 Bacteroides.id | cut -f 1`; 
do
prokka Bacteroides.all.genome/${i}.fna \
--outdir Bacteroides.all.genome_prokka/${i}/  \
--prefix ${i} \
--mincontiglen 1000 \
--cpus 46
done




conda activate /public/home/chenrz/.conda/envs/roary

cd /public5/chenrz/PC_BO/Bacteroides.all

mkdir Bacteroides.all.genome.gff/
cp Bacteroides.all.genome_prokka/*/*.gff  Bacteroides.all.genome.gff/


roary -e --mafft \
-p 50 \
Bacteroides.all.genome.gff/*.gff \
-f Bacteroides.all.genome_roary_out

~~~

 

## 2 Bifidobacterium

It is necessary to download the genomes of all Bifidobacterium from NCBI Refseq database or Download a genome data package by datasets download

~~~bash
cd /public5/chenrz/PC_BO/
mkdir Bifidobacterium.all
cd Bifidobacterium.all

# 1 Get genome metadata

datasets summary genome taxon 'Bifidobacterium' >Bifidobacterium.genome

# 2 Download a genome data package

datasets download genome taxon 'Bifidobacterium' --filename Bifidobacterium_dataset.zip


cd /public5/chenrz/PC_BO/Bifidobacterium.all/Bifidobacterium_dataset/genome



ls > ../Bifidobacterium.id
cd .. 
vim Bifidobacterium.id
:%s/.fna//g


cd /public5/chenrz/PC_BO/Bifidobacterium.all/

for i in `tail -n+1 Bifidobacterium.id | cut -f 1`; 
do
prokka Bifidobacterium_dataset/genome/${i}.fna \
--outdir Bifidobacterium.genome_prokka/${i}/  \
--prefix ${i} \
--mincontiglen 1000 \
--cpus 26 
done



conda activate /public/home/chenrz/.conda/envs/roary

cd /public5/chenrz/PC_BO/Bifidobacterium.all/


mkdir Bifidobacterium.all.genome.gff/
cp  Bifidobacterium.genome_prokka/*/*.gff  Bifidobacterium.all.genome.gff/


roary -e --mafft \
-p 50 \
Bifidobacterium.all.genome.gff/*.gff \
-f Bifidobacterium.all.genome_roary_out

~~~



## 3 Bacillus

~~~bash
cd /public5/chenrz/PC_BO/
mkdir Bacillus.all
cd Bacillus.all

# 1 Get genome metadata

datasets summary genome taxon 'Bacillus' >Bacillus.genome.json

# 2 Download a genome data package

datasets download genome taxon 'Bacillus' --filename Bacillus_dataset.zip

cd /public5/chenrz/PC_BO/Bacillus.all/Bacillus_dataset/genome



ls > ../Bacillus.id
cd .. 
vim Bacillus.id
:%s/.fna//g


cd /public5/chenrz/PC_BO/Bacillus.all/

for i in `tail -n+1 Bacillus.id | cut -f 1`; 
do
prokka Bacillus_dataset/genome/${i}.fna \
--outdir Bacillus.genome_prokka/${i}/  \
--prefix ${i} \
--mincontiglen 1000 \
--cpus 26 
done




conda activate /public/home/chenrz/.conda/envs/roary

cd /public5/chenrz/PC_BO/Bacillus.all/


mkdir Bacillus.all.genome.gff/
cp  Bacillus.genome_prokka/*/*.gff  Bacillus.all.genome.gff/


roary -e --mafft \
-p 50 \
Bacillus.all.genome.gff/*.gff \
-f Bacillus.all.genome_roary_out

~~~

## 4 Prevotella



~~~bash
cd /public5/chenrz/PC_BO/
mkdir Prevotella.all
cd Prevotella.all

# 1 Get genome metadata

datasets summary genome taxon 'Prevotella' >Prevotella.genome.json

# 2 Download a genome data package

datasets download genome taxon 'Prevotella' --filename Prevotella_dataset.zip


cd /public5/chenrz/PC_BO/Bacillus.all/Prevotella_dataset/genome



ls > ../Prevotella.id
cd .. 
vim Prevotella.id
:%s/.fna//g


cd /public5/chenrz/PC_BO/Prevotella.all/



for i in `tail -n+1 Prevotella.id | cut -f 1`; 
do
prokka Prevotella_dataset/genome/${i}.fna \
--outdir Prevotella.genome_prokka/${i}/  \
--prefix ${i} \
--mincontiglen 1000 \
--cpus 26 
done


conda activate /public/home/chenrz/.conda/envs/roary

cd /public5/chenrz/PC_BO/Prevotella.all/



mkdir Prevotella.all.genome.gff/
cp  Prevotella.genome_prokka/*/*.gff  Prevotella.all.genome.gff/



roary -e --mafft \
-p 50 \
Prevotella.all.genome.gff/*.gff \
-f Prevotella.all.genome_roary_out

~~~





## 5 Ruminococcus



~~~bash
cd /public5/chenrz/PC_BO/
mkdir Ruminococcus.all
cd Ruminococcus.all

# 1 Get genome metadata

datasets summary genome taxon 'Ruminococcus' >Ruminococcus.genome.json


# 2 Download a genome data package
datasets download genome taxon 'Ruminococcus' --filename Ruminococcus_dataset.zip


mkdir -p Ruminococcus_dataset/genome

cd Ruminococcus_dataset/genome
ls > ../Ruminococcus.id
cd .. 
vim Ruminococcus.id
:%s/.fna//g


cd /public5/chenrz/PC_BO/Ruminococcus.all/

conda activate /public/apps/anaconda2/envs/prokka

for i in `tail -n+1 Ruminococcus.id | cut -f 1`; 
do
prokka Ruminococcus_dataset/genome/${i}.fna \
--outdir Ruminococcus.genome_prokka/${i}/  \
--prefix ${i} \
--mincontiglen 1000 \
--cpus 26 
done




conda activate /public/home/chenrz/.conda/envs/roary

cd /public5/chenrz/PC_BO/Ruminococcus.all/


mkdir Ruminococcus.all.genome.gff/
cp  Ruminococcus.genome_prokka/*/*.gff  Ruminococcus.all.genome.gff/



roary -e --mafft \
-p 50 \
Ruminococcus.all.genome.gff/*.gff \
-f Ruminococcus.all.genome_roary_out

~~~

## 6 Faecalibacterium 



~~~bash
cd /public5/chenrz/PC_BO/
mkdir Faecalibacterium.all
cd Faecalibacterium.all

# 1 Get genome metadata

datasets summary genome taxon 'Faecalibacterium' >Faecalibacterium.genome.json

# 2 Download a genome data package
datasets download genome taxon 'Faecalibacterium' --filename Faecalibacterium_dataset.zip

mkdir -p Faecalibacterium_dataset/genome

cd Faecalibacterium_dataset/genome
ls > ../Faecalibacterium.id
cd .. 
vim Faecalibacterium.id
:%s/.fna//g


cd /public5/chenrz/PC_BO/Faecalibacterium.all/

conda activate /public/apps/anaconda2/envs/prokka

for i in `tail -n+1 Faecalibacterium.id | cut -f 1`; 
do
prokka Faecalibacterium_dataset/genome/${i}.fna \
--outdir Faecalibacterium.genome_prokka/${i}/  \
--prefix ${i} \
--mincontiglen 1000 \
--cpus 26 
done




conda activate /public/home/chenrz/.conda/envs/roary

cd /public5/chenrz/PC_BO/Faecalibacterium.all/

mkdir Faecalibacterium.all.genome.gff/
cp  Faecalibacterium.genome_prokka/*/*.gff  Faecalibacterium.all.genome.gff/

roary -e --mafft \
-p 50 \
Faecalibacterium.all.genome.gff/*.gff \
-f Faecalibacterium.all.genome_roary_out

~~~

