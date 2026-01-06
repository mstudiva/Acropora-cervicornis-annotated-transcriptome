# Acropora cervicornis Transcriptome Annotation, version December 20, 2023
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com) for use on FAU's HPC (KoKo)


#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster
# terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions -
# please make sure to read them before copy-pasting.


#------------------------------
# setup

# To install Bioperl as a conda environment
conda create -y -n bioperl perl-bioperl

# getting scripts
cd ~/bin
git clone https://github.com/z0on/annotatingTranscriptomes.git
mv annotatingTranscriptomes/* .
rm -rf annotatingTranscriptomes
rm launcher_creator.py

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git
mv emapper_to_GOMWU_KOGMWU/* .
rm -rf emapper_to_GOMWU_KOGMWU

git clone https://github.com/mstudiva/Acropora-cervicornis-annotated-transcriptome.git
mv Acropora-cervicornis-annotated-transcriptome/* .
rm -rf Acropora-cervicornis-annotated-transcriptome

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate


#------------------------------
# getting transcriptomes

# Locatelli et al. (2024): https://zenodo.org/doi/10.5281/zenodo.12021086
# scp from local machine

# Renaming gene identifiers for ease
sed -i 's/FUN_/Acervicornis/g' Acervicornis.fasta

# transcriptome statistics
conda activate bioperl
echo "seq_stats.pl Acervicornis.fasta > seqstats_Acervicornis.txt" > seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Acervicornis.fasta
-------------------------
36455 sequences.
1678 average length.
65308 maximum length.
61 minimum length.
N50 = 2747
61.2 Mb altogether (61163689 bp).
0 ambiguous Mb. (100 bp, 0%)
0 Mb of Ns. (100 bp, 0%)
-------------------------


#------------------------------
# uniprot annotations with blast

# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# unzipping
gunzip uniprot_sprot.fasta.gz &

# indexing the fasta database
module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch mdb.slurm

# splitting the transcriptome into 50 chunks, or however many is needed to keep the number of seqs per chunk under 1000
splitFasta.pl Acervicornis.fasta 50

# blasting all 50 chunks to uniprot in parallel, 4 cores per chunk
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
launcher_creator.py -j bl -n blast -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch blast.slurm

# watching progress:
grep "Query= " subset*.br | wc -l
# you should end up with the same number of queries as sequences from the seq_stats script

# combining all blast results
cat subset*br > myblast.br
rm subset*

# for trinity-assembled transcriptomes: annotating with isogroups
grep ">" Acervicornis.fasta | perl -pe 's/>Acervicornis(\d+)(\S+)\s.+/Acervicornis$1$2\tAcervicornis$1/'>Acervicornis_seq2iso.tab
cat Acervicornis.fasta | perl -pe 's/>Acervicornis(\d+)(\S+).+/>Acervicornis$1$2 gene=Acervicornis$1/'>Acervicornis_iso.fasta


#-------------------------
# extracting coding sequences and corresponding protein translations:
conda activate bioperl # if not active already
echo "perl ~/bin/CDS_extractor_v2.pl Acervicornis_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch cddd


#------------------------------
# GO annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

# selecting the longest contig per isogroup (also renames using isogroups based on Acervicornis and Cladocopium annotations):
fasta2SBH.pl Acervicornis_iso_PRO.fas >Acervicornis_out_PRO.fas

# scp your *_out_PRO.fas file to laptop, submit it to
http://eggnog-mapper.embl.de
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*_out_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# http://eggnog-mapper.embl.de/job_status?jobname=MM_f9vtfjaq

# once it is done, download results to HPC:
wget http://eggnog-mapper.embl.de/MM_f9vtfjaq/out.emapper.annotations

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Acervicornis_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Acervicornis_iso2geneName.tab


#------------------------------
# KOG annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
# this doesn't appear to be working, but the below line works awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "[,#S]" >Acervicornis_iso2kogClass1.tab
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "\tNA" >Acervicornis_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Acervicornis_iso2kogClass1.tab > Acervicornis_iso2kogClass.tab


#------------------------------
# KEGG annotations:

# selecting the longest contig per isogroup:
srun fasta2SBH.pl Acervicornis_iso.fasta >Acervicornis_4kegg.fasta

# scp *4kegg.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*4kegg.fasta .
# use web browser to submit 4kegg.fasta file to KEGG's KAAS server (http://www.genome.jp/kegg/kaas/)
# select SBH method, upload nucleotide query
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1703099308&key=XH4o3NDM

# Once it is done, download to HPC - it is named query.ko by default
wget https://www.genome.jp/tools/kaas/files/dl/1703099308/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Acervicornis_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc


#------------------------------
# file transfer

# copy all files to local machine
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/\* .
