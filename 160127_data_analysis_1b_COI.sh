## 16.01.27 - /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_1b_COI.sh
#  =================================================================================================================================
# shell analysis script COI data preparation from raw data:
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_1a_18S.sh
# and: 
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_1b_COI.sh
# (this file) lead to files processed in:
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_2.sh
# subsequently used in R analysis:
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/S9_Data_analysis.R

## performed steps:
#  ------------------
# Step 1: 454 read deconvolution, no denoising, RC of flipped reads
# Step 2: chimera detection
# Step 3: OTU picking
# Step 4: picking of a representative sequence set (rep. set)
# Step 5: taxonomy assignment of the various rep. sets using multiple tresholds
# Step 6: OTU table generation for all parameter combinations
# Step 7: Addition of metadata to biom files
# Step 8: Removal of unassigned reads
# Step 9: Subtraction of data contained in blank reactions
# Step 10: Subtraction of data contained in extraction controls
# Step 11: Retetion of invertebrate phylotypes
# Step 12: Retention of Austrlian and Antarctic samples specific to this sub-project

## file history 
# ---------------
# copy of 140811_COI_analysis.sh for analysis repetition after long break
# copied and modified from /Users/paul/Documents/140911_c3_analysis/1_analysis_shell_COI_seqdata/140916_C3_COI_analysis.sh
# latest version /Users/paul/Documents/140911_c3_analysis/1_analysis_shell_COI_seqdata/150121_C3_COI_analysis.sh


# ~~~~~~~~~~~~[ started including into methods section - below]~~~~~~~~~~
# re-read 7.1.2015

	# [COI workflow: use denoising approach as 17.7.14, but trim primers off before denoising]
	# ========================================================================================

	# split mapping files based on primer orientation "fwd" / "rev"
	# done and checked (see log file)

	# =======[demultiplexing]=======

	# demultiplexing (1/4)
	# 1st plate COI fwd, n= 1000000 set n higher for later plates
	split_libraries.py -m /Users/paul/Documents/140524_c3_analysis/140716_mapping_file/140714_COI_p1_mf_fwd.txt -f /Users/paul/Documents/140407_454_runs/PC_COI_1/1.TCA.454Reads.fna -q /Users/paul/Documents/140407_454_runs/PC_COI_1/1.TCA.454Reads.qual -l 350 -L 400 -b 12 -M 3 -o /Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p1_fwd/ -w 50 -z truncate_only -d  -n 1000000

	# demultiplexing (2/4)
	# 1st plate COI rev, n= 2000000 set n higher for later plates
	split_libraries.py -m /Users/paul/Documents/140524_c3_analysis/140716_mapping_file/140714_COI_p1_mf_rev.txt -f /Users/paul/Documents/140407_454_runs/PC_COI_1/1.TCA.454Reads.fna -q /Users/paul/Documents/140407_454_runs/PC_COI_1/1.TCA.454Reads.qual -l 350 -L 400 -b 12 -M 3 -o /Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p1_rev/ -w 50 -z truncate_only -d  -n 2000000

	# demultiplexing (3/4)
	# 2nd plate COI fwd, n= 3000000 set n higher for later plates
	split_libraries.py -m /Users/paul/Documents/140524_c3_analysis/140716_mapping_file/140714_COI_p2_mf_fwd.txt -f /Users/paul/Documents/140407_454_runs/PC_COI_2/2.TCA.454Reads.fna -q /Users/paul/Documents/140407_454_runs/PC_COI_2/2.TCA.454Reads.qual -l 350 -L 400 -b 12 -M 3 -o /Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p2_fwd/ -w 50 -z truncate_only -d  -n 3000000

	# demultiplexing (4/4)
	# 2nd plate COI rev, n= 4000000 set n higher for later plates
	split_libraries.py -m /Users/paul/Documents/140524_c3_analysis/140716_mapping_file/140714_COI_p2_mf_rev.txt -f /Users/paul/Documents/140407_454_runs/PC_COI_2/2.TCA.454Reads.fna -q /Users/paul/Documents/140407_454_runs/PC_COI_2/2.TCA.454Reads.qual -l 350 -L 400 -b 12 -M 3 -o /Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p2_rev/ -w 50 -z truncate_only -d  -n 4000000

	# =====[preparation of additional primer trimming]=====

	# copy files for additional primer trimming
	cp /Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p1_fwd/seqs.fna.gz \
	/Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p1_fwd/seqs_trim.fna.gz && \
	cp /Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p1_rev/seqs.fna.gz \
	/Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p1_rev/seqs_trim.fna.gz && \
	cp /Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p2_fwd/seqs.fna.gz \
	/Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p2_fwd/seqs_trim.fna.gz && \
	cp /Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p2_rev/seqs.fna.gz \
	/Users/paul/Documents/140524_c3_analysis/140717_splitlibs_p2_rev/seqs_trim.fna.gz

	## ====[additional primer trimming with Geneious]===
	# since counting didnt work this is done in Genious again - without sorting the sequences 
	# folder 140806_COI
	# Leray primers (fwd / rev) also looked for RC pair
	# 5 Mismatches
	# 3 Minimum match length
	# export files to wd, adding "trm" - indices checked and ok (head)
	/Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p1_fwd_trm.fna
	/Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p1_rev_trm.fna
	/Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p2_fwd_trm.fna
	/Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p2_rev_trm.fna

	# ===== denoising [doesn't work when .fnas are touched after splitting, go without denoising for now ] =====
	# notes of this step erased to remove clutter check old file is details are needed

	# ===== inflating [ not necessary without denoising, -i and -o removed for safety reasons ] ===== 
	# notes of this step erased to remove clutter check old file is details are needed

	# =======[flipping of reversed reads]=======

	# flipping (1/2) COI plate 1 rev to RC
	fastx_reverse_complement \
	-i /Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p1_rev_trm.fna \
	-o /Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p1_rev_trm_RC.fna

	# flipping (2/2) COI plate 2 rev to RC
	fastx_reverse_complement \
	-i /Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p2_rev_trm.fna \
	-o /Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p2_rev_trm_RC.fna


	# ========[concatenate fasta files (p1 fwd rev_RC, p2 fwd rev_RC)]========

	cat \
	/Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p1_fwd_trm.fna \
	/Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p1_rev_trm_RC.fna \
	/Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p2_fwd_trm.fna \
	/Users/paul/Documents/140524_c3_analysis/140811_cutadapt/140811_COI_p2_rev_trm_RC.fna \
	> \
	/Users/paul/Documents/140524_c3_analysis/140811_cated_fnas/140811_COI_reads_undenoised.fna


	# ========[cleaning up and getting file to cluster]===================

	# packing up to to copy to cluster
	pigz /Users/paul/Documents/140524_c3_analysis/140811_cated_fnas/140811_COI_reads_undenoised.fna

	# copy file to Tizard
	scp \
	/Users/paul/Documents/140524_c3_analysis/140811_cated_fnas/140811_COI_reads_undenoised.fna.gz \
	pczechowski@bigmem1024-1.tizard.ersa.edu.au:/home/users/pczechowski/140728_starcluster_data/

	# copy unckustered undenoised sequence file (and a ref db) to Stacluster (see /Users/paul/Documents/140524_c3_analysis/140715_analysis_notes/140808_COI_ref_db.sh)
	scp \
	-i /opt/local/shared_nectarkeys/ACAD/acad_keypair.pem pczechowski@130.220.209.109 \
	140811_COI_reads_undenoised.fna.gz \
	140811_COI_ref_db.tgz \

	# ==============[on starcluster - reorganize folder in COI data]========

	# reference dbs stored in - /mnt/paul_folder/140719_COI_data/140812_COI_ref_DBs
	mkdir 140812_COI_ref_DBs
	mv Q_COI_rdb_fasta.fna  ./140812_COI_ref_DBs/
	mv Q_COI_rdb_taxonomy.txt ./140812_COI_ref_DBs/

	mkdir 140719_raw_sequences

	# raw sequences to here
	scp \
	-i /opt/local/shared_nectarkeys/ACAD/acad_keypair.pem \
	140811_COI_reads_undenoised.fna.gz \
	pczechowski@130.220.209.109:/mnt/paul_folder/140719_COI_data/140719_raw_sequences/

	#db goes here 
	scp \
	-i /opt/local/shared_nectarkeys/ACAD/acad_keypair.pem \
	140811_COI_ref_db.tgz \
	pczechowski@130.220.209.109:/mnt/paul_folder/140719_COI_data/140812_COI_ref_DBs/

	# more moving around
	mv -t ./140812_COI_ref_DBs/ 140725_COI_NCBI_invert.zip 140723_COI_refseqs.zip 140723_blast.zip

	# unpacking files
	tar -zxvf /mnt/paul_folder/140719_COI_data/140812_COI_ref_DBs/140811_COI_ref_db.tgz
	gzip -d /mnt/paul_folder/140719_COI_data/140719_raw_sequences/140811_COI_reads_undenoised.fna.gz


	# ===== [Chimera detection in COI data] =====

	# abundance based, without reference, with usearch61
	identify_chimeric_seqs.py \
	-i /mnt/paul_folder/140719_COI_data/140719_raw_sequences/140811_COI_reads_undenoised.fna \
	-o /mnt/paul_folder/140719_COI_data/140812_chim_search \
	-m usearch61 \
	--suppress_usearch61_ref

	# checking log in created dir
	cat /mnt/paul_folder/140719_COI_data/140812_chim_search/identify_chimeric_seqs.log
	# denovo_chimeras		37
	# denovo_non_chimeras	90202

	# =====[ Chimera removal on COI data (on Starcluster) ]==== 
	#
	# line count of fasta file prior to removal of chimeras
	grep ">" /mnt/paul_folder/140719_COI_data/140719_raw_sequences/140811_COI_reads_undenoised.fna | wc -l 		# 90239 - ok

	#removal of chimeras
	filter_fasta.py \
	-f /mnt/paul_folder/140719_COI_data/140719_raw_sequences/140811_COI_reads_undenoised.fna \
	-o /mnt/paul_folder/140719_COI_data/140719_raw_sequences/140811_COI_reads_undenoised_nChm.fna \
	-s /mnt/paul_folder/140719_COI_data/140812_chim_search/chimeras.txt -n

	# line count of fasta file prior to removal of chimeras
	grep ">" /mnt/paul_folder/140719_COI_data/140719_raw_sequences/140811_COI_reads_undenoised_nChm.fna | wc -l # 90202 -ok 37 removed


	# =====[ 1st OTU picking: OTU picking on cluster - at 97% ]===== 
	# with reverse strand matching this time, optimal and exact
	# tested w/o nohup then started Tue Aug 12 04:21:24 EDT 2014
	# results verified 7.1.15
	
	nohup \
	pick_otus.py \
	-i /mnt/paul_folder/140719_COI_data/140719_raw_sequences/140811_COI_reads_undenoised_nChm.fna \
	-o /mnt/paul_folder/140719_COI_data/140812_uclust_otus/ \
	-s 0.97 \
	-m uclust \
	--optimal \
	--exact_uclust \
	-z \
	&> 140812_uclust_otus_nohup.out

	# log file: /mnt/paul_folder/140719_COI_data/140812_uclust_otus/140811_COI_reads_undenoised_nChm_otus.log
	# Num OTUs:1641
	# Result path: /mnt/paul_folder/140719_COI_data/140812_uclust_otus/140811_COI_reads_undenoised_nChm_otus.txt

# re-read 7.1.2015 - end
# ~~~~~[not included into methods section - below] ~~~~~~

# =====[ 2nd OTU picking: OTU picking on cluster - at 100 %]===== 
# with reverse strand matching this time, optimal and exact
# tested w/o nohup then started Tue Aug 12 04:21:24 EDT 2014
# 7.1.15 - this was not really necessary, it results in a lot of OTUs, these results will not be used probably
nohup \
pick_otus.py \
-i /mnt/paul_folder/140719_COI_data/140719_raw_sequences/140811_COI_reads_undenoised_nChm.fna \
-o /mnt/paul_folder/140719_COI_data/140824_uclust_otus/ \
-s 1.00 \
-m uclust \
--optimal \
--exact_uclust \
-z \
&> 140824_uclust_otus_nohup.out

# Num OTUs:15105

# =====[ 3rd OTU picking: OTU picking on cluster - at 99 %]===== 
# with reverse strand matching this time, optimal and exact
# tested w/o nohup then started Sat Sep 27 03:34:04 EDT 2014
# verified result 7.1.2015
nohup \
pick_otus.py \
-i /mnt/paul_folder/140719_COI_data/140719_raw_sequences/140811_COI_reads_undenoised_nChm.fna \
-o /mnt/paul_folder/140719_COI_data/140825_uclust_otus/ \
-s 0.99 \
-m uclust \
--optimal \
--exact_uclust \
-z \
&> 140825_uclust_otus_nohup.out


# log file: /mnt/paul_folder/140719_COI_data/140825_uclust_otus/140811_COI_reads_undenoised_nChm_otus.log
# Num OTUs:2957 
# Result path: /mnt/paul_folder/140719_COI_data/140825_uclust_otus/140811_COI_reads_undenoised_nChm_otus.txt

# ~~~~~[not included into methods section - above] ~~~~~~



# ======[pick rep sets COI based on 97% and 99% and 100% OTU clustering - start]====
# verified 7.1.2015

# ~~~~~~ [needs to be adopted in manuscript - below]~~~~~~~~~~~~~~

# = run the script 
# = rep set picking of COI data C3 30.9.14 

# define fixed paths
repsetdir="/mnt/paul_folder/140719_COI_data/140925_rep_set"
unclustfna="/mnt/paul_folder/140719_COI_data/140719_raw_sequences/140811_COI_reads_undenoised_nChm.fna"

#define OTU maps
otu_map[1]="/mnt/paul_folder/140719_COI_data/140824_uclust_otus/140811_COI_reads_undenoised_nChm_otus.txt"
otu_map[2]="/mnt/paul_folder/140719_COI_data/140812_uclust_otus/140811_COI_reads_undenoised_nChm_otus.txt"
otu_map[3]="/mnt/paul_folder/140719_COI_data/140825_uclust_otus/140811_COI_reads_undenoised_nChm_otus.txt"

#define output file prefix names
prefix[1]="140925_rep_set_100"
prefix[2]="140825_rep_set_097"
prefix[3]="140825_rep_set_099"

# create output directory"
mkdir -p "$repsetdir"

# run loop over files 
for ((i = 1; i <= "${#otu_map[@]}"; i++)) ; do

# the real command, remove echo if script runs fine
pick_rep_set.py -i "${otu_map[$i]}" -f "$unclustfna" -o "$repsetdir"/"${prefix[$i]}".fasta

done 

# = count sequences in files
# verified 7.1.2015
grep ">" /mnt/paul_folder/140719_COI_data/140925_rep_set/140825_rep_set_097.fasta | wc -l # 1641
grep ">" /mnt/paul_folder/140719_COI_data/140925_rep_set/140825_rep_set_099.fasta | wc -l # 2957
grep ">" /mnt/paul_folder/140719_COI_data/140925_rep_set/140925_rep_set_100.fasta | wc -l # 15105

 
# ~~~~~~~~~~~~~[ needs to be adopted in manuscript  - above ]~~~~~~~~~~~~~~

# ======[taxonomy assignment of COI rep. sets ]====
# re-edited and verified 21.1.2015 - start

# use Uclust because it is the fastest method (suitable for large data sets)
# http://qiime.org/scripts/parallel_assign_taxonomy_uclust.html

# use 3 different assignment treshholds to compare results:
# ---------------------------------------------------------
# --uclust_similarity 0.90
# --uclust_similarity 0.95
# --uclust_similarity 0.99


# load Qiime module on Starcluster
# ---------------------------------
module load Qiime/1.8.0

# declare variables
# -----------------

# 21.1.2015 new ref db (unclustered, dereplicated by clustering at 100%, inkl AVC data - 100 sequences)
refdb100tax="/mnt/paul_folder/140719_COI_data/150120_COI_refdb/150120_COI_ref_TAX.txt"
refdb100seq="/mnt/paul_folder/140719_COI_data/150120_COI_refdb/150120_COI_ref_SEQ.fna"

#rep sets (omitting 100% clustered rep set, using only 97%, 99%):
repset[1]="/mnt/paul_folder/140719_COI_data/140925_rep_set/140825_rep_set_097.fasta"
repset[2]="/mnt/paul_folder/140719_COI_data/140925_rep_set/140825_rep_set_099.fasta"
# vars checked

# define taxonomy assignment thresholds in array (testing differnt assignment treshholds)
threshold[1]="0.90"
threshold[2]="0.95"
threshold[3]="0.99"
threshold[4]="0.70"
threshold[5]="0.75"
threshold[6]="0.80"
threshold[7]="0.85"
# should be ok

# define output destinations matching array position in "threshhold"
destination[1]="/mnt/paul_folder/140719_COI_data/150121_tax_assignment_90"
destination[2]="/mnt/paul_folder/140719_COI_data/150121_tax_assignment_95"
destination[3]="/mnt/paul_folder/140719_COI_data/150121_tax_assignment_99"
destination[4]="/mnt/paul_folder/140719_COI_data/150121_tax_assignment_70"
destination[5]="/mnt/paul_folder/140719_COI_data/150121_tax_assignment_75"
destination[6]="/mnt/paul_folder/140719_COI_data/150121_tax_assignment_80"
destination[7]="/mnt/paul_folder/140719_COI_data/150121_tax_assignment_85"
# new output destinations created

# define output destinations 
# new output destinations created
dir[1]="150121_rep_set_097_vs_db_100"
dir[2]="150121_rep_set_099_vs_db_100"

# run nested loops 
# ---------------------------------
# restarted Fri Jan  9 18:50:29 CST 2015
# with normal and flipped reference database
# with all seven treshholds

# loop through rep sets (2)
for ((i = 1; i <= "${#repset[@]}"; i++))
do
# loop through threshholds (3)
for ((k = 1; k <= "${#threshold[@]}"; k++))
do
# for debugging only 
date
echo --+ repset: "${repset[$i]}" 
echo --+ threshhold : "${threshold[$k]}"

# the real command, remove echo once ready to run 
parallel_assign_taxonomy_uclust.py -i "${repset[$i]}" -o "${destination[$k]}"/"${dir[$i]}" -r "$refdb100seq" -t "$refdb100tax" --uclust_similarity "${threshold[$k]}" -O 8

# closing inner loop
done 
# closing outer loop 
done

# taxonomy assignment script for Starcluster - end
# ------------------------------------------------

# 21.1.14: verify tax assignment results:
# --------------------------------------

# check all log files
less `ls 150121_tax_assignment_*/150121_rep_set_*/*.log`


# ======[taxonomy assignment of COI rep. sets - end]====


# 7.1.14 OTU table generation - start
# ==================================
## re-done 10.01.2015
# make output dir
# parse input filenames and generate output filenames
# run OTU table command
# should generate 28 biom tables, of which "FLP" are stemming from the RCd reference database and can be aerased later

mkdir -p /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/

for file in /mnt/paul_folder/140719_COI_data/150121_tax_assignment_*/150121_rep_set_*/*.txt; do	
case "$file" in
*_rep_set_097_vs_db_100* ) otumap="/mnt/paul_folder/140719_COI_data/140812_uclust_otus/140811_COI_reads_undenoised_nChm_otus.txt" && clustering="_clust97" ;;&
*_rep_set_099_vs_db_100* ) otumap="/mnt/paul_folder/140719_COI_data/140825_uclust_otus/140811_COI_reads_undenoised_nChm_otus.txt" && clustering="_clust99" ;;&
*150121_tax_assignment_70* ) taxassgn="_tassgn70" ;;&
*150121_tax_assignment_75* ) taxassgn="_tassgn75" ;;&
*150121_tax_assignment_80* ) taxassgn="_tassgn80" ;;&
*150121_tax_assignment_85* ) taxassgn="_tassgn85" ;;&
*150121_tax_assignment_90* ) taxassgn="_tassgn90" ;;&
*150121_tax_assignment_95* ) taxassgn="_tassgn95" ;;&
*150121_tax_assignment_99* ) taxassgn="_tassgn99" ;;&
esac
echo -t "$file" # print input taxonomy assignment file
echo -i "$otumap" # print input OTU table
biom_path="/mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/150121_COI_OTUs$clustering$taxassgn.biom"
echo -o "$biom_path"
make_otu_table.py -v -i "$otumap" -t "$file" -o "$biom_path" 
done

# 14 biom tables present? -ok 
ls /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/*.biom | wc -l


# 21.1.14 Metadata addition - start
# ==================================
# define input path names
mf="/mnt/paul_folder/140719_COI_data/150108_COI_mf_w_md/150108_COI_MF.txt"

# looping over *.biom files in a good way
# ading metadata from mapping file 
for file in /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/*.biom; do
[[ -e $file ]] || continue

# filename manipulations
targetdir=$(dirname $file) 
filename=$(basename "$file")
suffix="${filename##*.}"
prefix="${filename%.*}"

# command
biom add-metadata \
-i "$file" \
-o "$targetdir"/"$prefix"_md."$suffix" \
-m "$mf"

# generate summary files
biom summarize-table -i  "$targetdir"/"$prefix"_md."$suffix" -o "$targetdir"/"$prefix"_md.sum_qual.txt --qualitative
biom summarize-table -i  "$targetdir"/"$prefix"_md."$suffix" -o "$targetdir"/"$prefix"_md.sum_quan.txt
done

# read summaries on Starcluster
# - usable tables no yet created -
# for now all summaries stored in /mnt/paul_folder/140719_COI_data/151001_COI_OTU_tables/150110_all_summaries.txt
for file in /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/*.txt; do
echo "$file" >> /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/150121_all_summaries.txt
cat "$file" >> /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/150121_all_summaries.txt
done

# 21.1.14 Metadata addition - end
# ==================================


# 21.1.14 COI table filetering - remove unassigned - start
# =======================================================

## 21.1.14 remove unassigned
# --------------------------
for file in /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/*_md.biom; do
[[ -e $file ]] || continue

# filename manipulations
targetdir=$(dirname $file) 
filename=$(basename "$file")
suffix="${filename##*.}"
prefix="${filename%.*}"

# command
filter_taxa_from_otu_table.py \
-i "$file" \
-o "$targetdir"/"$prefix"_assigned_only."$suffix" \
-n Unassigned

# generate summary files
biom summarize-table -i  "$targetdir"/"$prefix"_assigned_only."$suffix" -o "$targetdir"/"$prefix"_assigned_only.sum_qual.txt --qualitative
biom summarize-table -i  "$targetdir"/"$prefix"_assigned_only."$suffix" -o "$targetdir"/"$prefix"_assigned_only.sum_quan.txt

done

# read summaries on Starcluster
# - usable tables no yet created -
# for now all summaries stored in /mnt/paul_folder/140719_COI_data/151001_COI_OTU_tables/150110_all_summaries_md_assigned_only.txt
for file in /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/*md_assigned_only.sum_*.txt; do
echo "$file" >> /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/150121_all_summaries_md_assigned_only.txt
cat  "$file" >> /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/150121_all_summaries_md_assigned_only.txt
done

# 8.1.15 COI table filetering - remove unassigned - end
# =======================================================




# 10.1.15 COI table filetering - substract blanks - start
# =======================================================
## this will not work as there are no blank controls lefet to substarct in the COI data:
# this is can be seen as the following command does not return anything:
grep "PCR" /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/150121_all_summaries_md_assigned_only.txt
echo "This directory can't be filled as there are no blanks left to process in COI data" >> /mnt/paul_folder/140719_COI_data/150110_COI_OTU_tables_substracted_blanks/150110_README.txt

## the following code is untested but should be more or less functional 

# substarct blanks in biom tables that have metadata associated with them i.e. these 
ls /mnt/paul_folder/140719_COI_data/151001_COI_OTU_tables/*md_assigned_only.biom

# define output directory and create it 
targetdir="/mnt/paul_folder/140719_COI_data/150110_COI_OTU_tables_substracted_blanks/"
mkdir -p  $targetdir

# define mapping file path
mf="/mnt/paul_folder/140719_COI_data/150108_COI_mf_w_md/150108_COI_MF.txt"

# Print 9th field
awk -F '\t' '{print $9}' /mnt/paul_folder/140719_COI_data/150108_COI_mf_w_md/150108_COI_MF.txt

# the filetring commands in the loop 
for file in /mnt/paul_folder/140719_COI_data/151001_COI_OTU_tables/*md_assigned_only.biom; do
[[ -e $file ]] || continue

# debugging
echo "--> $file"
echo ""

# Create an OTU table with the just the blank control samples:
filter_samples_from_otu_table.py -i "$file" -o "$targetdir"$(basename "$file" .biom).biom_tmp -m "$mf" -s "XtrOri:PCR_blk_h20"
echo "1/5"

# Filter out OTU ids that have zero counts, as we only want the OTUs with positive counts from the Control_Blank samples:
filter_otus_from_otu_table.py \
-i "$targetdir"$(basename "$file" .biom).biom_tmp \
-o "$targetdir"$(basename "$file" .biom).biom_tmp2 \
-n 1
echo "2/5"

# Then create a tab separated version of this OTU table:
biom convert -b \
-i "$targetdir"$(basename "$file" .biom).biom_tmp2 \
-o "$targetdir"$(basename "$file" .biom).txt_tmp3
echo "3/5"

# Filter out OTU ids from the run 1 OTU table that were determined to be present in the Control_Blank samples:
filter_otus_from_otu_table.py \
-i "$file" -o "$targetdir"$(basename "$file" .biom).biom_tmp4 \
-e "$targetdir"$(basename file .biom).txt_tmp3
echo "4/5"

# The otu_table_1_minus_contaminants.biom file now has two samples with zero sequences associated with it. These can be removed to get a final OTU table:
filter_samples_from_otu_table.py \
-i "$targetdir"$(basename "$file" .biom).biom_tmp4 \
-o "$targetdir"$(basename "$file" .biom)_no_contaminants.biom -n 1
echo "5/5"

done

# 10.1.15 COI table filetering - substract blanks - end
# ======================================================



# 10.1.15 COI table filetering - substract extraction controls - start
# ====================================================================
# checking if samples for filtering are  contained in the COI data - nope - already filtered out (or perhaps not sequenced due to low coverage 
grep "XTRneg" /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/150121_all_summaries_md_assigned_only.txt

# 10.1.15 COI table filetering - substract extraction controls - end
# ====================================================================



# 21.1.15 COI table filetering - retain invertebrates - start 
# ===========================================================


# print all possible terms of the first taxonomy field
infile="/mnt/paul_folder/140719_COI_data/150120_COI_refdb/150120_COI_ref_TAX.txt"
awk -F '[\t;]' '{print $5}' "$infile" | sort -d | uniq -c
# 74603  Arthropoda
#  2794  Nematoda
#  1990  Rotifera
#   519  Tardigrada

# create output directory
targetdir="/mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables_inverts_only"
mkdir -p $targetdir

# check input file paths - choose only files with metadata and taxobnomy assignments
ls /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/*md_assigned_only.biom

# load qiime
module load Qiime/1.8.0

# loop
for file in /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/*md_assigned_only.biom; do
[[ -e $file ]] || continue

# pathname manipulation and verification - input file 
echo @ processed file: "$file"

# pathname manipulation and verification - tmp file 1
out_tmp1=$targetdir/$(basename $file .biom).tmp1
# echo @ temp file 1 "$out_tmp1"

# pathname manipulation and verification - tmp file 2
out_tmp2=$targetdir/$(basename $file .biom).tmp2
# echo @ temp file 2 "$out_tmp2"

# pathname manipulation and verification - biom file
out_biom1=$targetdir/$(basename $file .biom)_above5.biom
# echo @ biom file 1 "$out_biom1"

# pathname manipulation and verification - biom file
out_biom2=$targetdir/$(basename $file .biom)_above5_inverts.biom
# echo @ biom file 2 "$out_biom2"

# positive filtering of OTUs - discard anything less then 5
filter_otus_from_otu_table.py -i  "$file" -o  "$out_tmp1" -n 5
echo 1/4 - filtering low count OTUs

# positive filtering of taxa
filter_taxa_from_otu_table.py -i  "$out_tmp1" -o  "$out_tmp2" -p Arthropoda,Nematoda,Rotifera,Tardigrada

echo 2/4 - retaining invertebrates

# filter out zero count samples
filter_samples_from_otu_table.py -i "$out_tmp1" -o "$out_biom1" -n 1
echo 3/4 - filter 0s from all taxa

# filter out zero count sumples
filter_samples_from_otu_table.py -i "$out_tmp2" -o "$out_biom2" -n 1
echo 4/4 - filter 0s from invertebrates

done

# (erased tmp file in dir)

cd /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables_inverts_only

for file in /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables_inverts_only/*.biom; do
[[ -e $file ]] || continue

# generate summary files
echo "$file"
biom summarize-table -i  "$file" -o $(basename "$file" .biom).sum_qual.txt --qualitative
biom summarize-table -i  "$file" -o $(basename "$file" .biom).sum_quan.txt

done

# 21.1.15 COI table filtering - retain invertebrates - end 
# ==========================================================



# 21.1.15 COI table filtering - retain controls - start 
# ==========================================================

## check mapping file - which categories to filter?
Location:E_Ant_coast # AVC controls
Location:Australia # soil controls and mock mixes

# load qiime
module load Qiime/1.8.0

# file list of relavant biom tables
ls /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables_inverts_only/*_above5_inverts.biom

# define output directory and create it 
targetdir="/mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables_inverts_only_cntrls"
mkdir -p  $targetdir

# define mapping file path
mf="/mnt/paul_folder/140719_COI_data/150108_COI_mf_w_md/150108_COI_MF.txt"

# the filetring commands in the loop  - did not run weel 1st time hence commented sections
for file in for file in /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables_inverts_only/*_above5_inverts.biom; do
[[ -e $file ]] || continue

# pathname manipulation and verification - input file 
echo - 0/6 - processed file: "$file"

# pathname manipulation and verification - biom file
out_biom=$targetdir/$(basename $file .biom)_cntrl_only.biom
# echo - 1/6 - biom file "$out_biom"

# pathname manipulation and verification - tmp file 1
out_tmp1=$targetdir/$(basename $file .biom)_cntrl_only.tmp1
# echo - 2/6 - temp file "$out_tmp1"

# pathname manipulation and verification - tmp file 2
out_tmp2=$targetdir/$(basename $file .biom)_cntrl_only.tmp2
# echo - 3/6 - temp file "$out_tmp2"

# pathname manipulation and verification - tmp file 3
out_tmp3=$targetdir/$(basename $file .biom)_cntrl_only.tmp3
# echo - 4/6 - temp file "$out_tmp3"

# pathname manipulation and verification - tmp file 4
out_tmp4=$targetdir/$(basename $file .biom)_cntrl_only.tmp4
# echo - 5/6 - temp file "$out_tmp4"

# pathname manipulation and verification - tmp file 5
out_tmp5=$targetdir/$(basename $file .biom)_cntrl_only.tmp5
# echo - 6/6 - temp file "$out_tmp5"

# Remove Antarctic control samples from input table:
echo "1/3"
filter_samples_from_otu_table.py -i "$file" -o "$out_tmp1" -m "$mf" -s "Location:E_Ant_coast"

# Remove Australian control samples from input table:
echo "2/3"
filter_samples_from_otu_table.py -i "$file" -o "$out_tmp2" -m "$mf" -s "Location:Australia"

# merge both tables
echo "3/3"
merge_otu_tables.py -i "$out_tmp1","$out_tmp2" -o "$out_biom"

# generate summary files
echo "writing summaries"
biom summarize-table -i  "$out_biom" -o $targetdir/$(basename $"$out_biom" .biom).sum_qual.txt --qualitative
biom summarize-table -i  "$out_biom" -o $targetdir/$(basename $"$out_biom" .biom).sum_quan.txt

# erase temp files
rm -v /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables_inverts_only_cntrls/*.tmp?

done

# concatenate summary files
# - usable tables no yet created -
# for now all summaries stored in /mnt/paul_folder/140719_COI_data/150112_COI_OTU_tables_inv_only_controls/150112_all_summaries_md_assigned_only_inverts_controls.txt
for file in /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables_inverts_only_cntrls/*.sum_*.txt; do
echo "$file" >> /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables_inverts_only_cntrls/150121_all_summaries_md_assigned_only_inverts_controls.txt
cat  "$file" >> /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables_inverts_only_cntrls/150121_all_summaries_md_assigned_only_inverts_controls.txt
done

# 12.1.15 COI table filtering - retain controls - end 
# ==========================================================


# backup: Wed Jan 21 16:49:19 CST 2015; continued in 
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_2.sh
