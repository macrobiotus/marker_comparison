## 16.01.27 - /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_1a_18S.sh
#  =================================================================================================================================
# shell analysis script 18S data preparation form raw data:
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_1a_18S.sh
# (this file) and: 
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_1b_COI.sh
# lead to files processed in:
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_2.sh
# subsequently used in R analysis:
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/S9_Data_analysis.R


## performed steps:
#  ------------------
# Step  1: Splitting of 18S lanes form Illumina data
# Step  2: Merging of fasta files to create complete fna file with all 18S data
# Step  3: Phylotype clustering at 97% and 99% 
# Step  4: Pick rep. set for 97% and 99% clustered data
# Step  5: Taxonopmy assignments for 97% and 99% data using multiple treshholds
# Step  6: Verify rep set picking and assignment results
# Step  7: OTU table generation
# Step  8: Addition of metadata to OTU tables (.biom files)
# Step  9: Removal of unassigned phylotypes
# Step 10: Subtraction of data contained in PCR / extraction blanks blanks
# Step 11: Subtraction of phylotypes fewer then 5, retention of Invertebrates
# Step 12: Retention of Austrlian and Antarctic samples specific to this sub-project


# ===================[local mapping file versions]=====
# 12.11.14: 18S mapping file copied to current analysis folder
cp /Users/paul/Documents/140524_c3_analysis/140716_mapping_file/140714_18S_mf.txt /Users/paul/Documents/140911_c3_analysis/mapping_files
# validation of 18S mapping file in Qiime
validate_mapping_file.py -m 140714_18S_mf.txt -o 140714_18S_mf_check_141112
# for analysis after taxonomy assignment additional metadata might be added to further version of this file
 
# ===================[Mapping file and 18S data to Starcluster in Tizard]=====
 
# [pczechowski@tizard2 ~/140730_18S_data]$ pwd /home/users/pczechowski/140730_18S_data
ls -lh
# -rw-rw-r-- 1 pczechowski pczechowski 129M Jul 30 16:45 18S3_PC_merged.Index.fastq.gz
# -rw-rw-r-- 1 pczechowski pczechowski 328M Jul 30 16:45 18S3_PC_merged.Read.fastq.gz
# -rw-rw-r-- 1 pczechowski pczechowski 129M Jul 30 16:45 18S4_PC_merged.Index.fastq.gz
# -rw-rw-r-- 1 pczechowski pczechowski 324M Jul 30 16:45 18S4_PC_merged.Read.fastq.gz
# zipping up 
pwd #/home/users/pczechowski
zip -r0 140730_18S_data.zip 140730_18S_data
# zip file:
/home/users/pczechowski140730_18S_data.zip
# Starcluster
pwd # /mnt/paul_folder
mkdir 140730_18S_data
#getting data to starcluster - on tizard:
scp -i /opt/local/shared_nectarkeys/ACAD/acad_keypair.pem \
/home/users/pczechowski/140730_18S_data.zip \
pczechowski@130.220.209.109:/mnt/paul_folder/140730_18S_data/
# Starcluster (unpacking steps not verbosed)
pwd #/mnt/paul_folder/140730_18S_data
gunzip *.gz
# copying mapping file to tizard
# [17:50:14] paul@dsl3a5885:~ $ scp /Users/paul/Documents/140524_c3_analysis/140716_mapping_file/140714_18S_mf.txt   pczechowski@tizard2.ersa.edu.au:/home/users/pczechowski/140728_starcluster_data/
scp -i /opt/local/shared_nectarkeys/ACAD/acad_keypair.pem \
/home/users/pczechowski/140728_starcluster_data/140714_18S_mf.txt \
pczechowski@130.220.209.109:/mnt/paul_folder/140730_18S_data/

# ===================[merging of  18S reads]=====

# data is obviously already merged by Jimmy
join_paired_ends.py # request parameters from Jimmy 

# ==== split mapping file and check =====

# don't know if this is necacssary but might help with barcode recognition issues 
sed  '/18S.1/d' 140714_18S_mf.txt > 140714_18S_mf_p2.txt
sed  '/18S.2/d' 140714_18S_mf.txt > 140714_18S_mf_p1.txt
wc -l 140714_18S_mf_p2.txt 98
wc -l 140714_18S_mf_p1.txt 98
validate_mapping_file.py \
-m 140714_18S_mf_p1.txt \
-o 140714_18S_mf_p1_CHECK # ok
validate_mapping_file.py \
-m 140714_18S_mf_p2.txt \
-o 140714_18S_mf_p2_CHECK # ok

# ===================[Splitting 18S lanes]=====

# 1.8.2014: several things tried beforehand, plates are switched in Jimmys data and usage of the complete mapping file might
# give trouble because recognition isn't perfect and samples get recognized on a run were they
# weren't - so the following options are the best.

nohup \
split_libraries_fastq.py \
-i 18S4_PC_merged.Read.fastq \
-b 18S4_PC_merged.Index.fastq \
-o /mnt/paul_folder/140730_18S_data/140801_1st_plate_splitted \
-m 140714_18S_mf_p1.txt \
-q 10 \
--store_qual_scores \
--rev_comp_mapping_barcodes \
-s 100000000 \
&> 140731_split_libraries_fastq_nohup_out_plate1.txt

nohup \
split_libraries_fastq.py \
-i 18S3_PC_merged.Read.fastq \
-b 18S3_PC_merged.Index.fastq \
-o /mnt/paul_folder/140730_18S_data/140801_2nd_plate_splitted \
-m 140714_18S_mf_p2.txt \
-q 10 \
--store_qual_scores \
--rev_comp_mapping_barcodes \
-s 200000000 \
&> 140731_split_libraries_fastq_nohup_out_plate2.txt
# checked log files now - all good
# 12.11.14: check library splitting success
# see file /Users/paul/Documents/140911_c3_analysis/3_analysis_troubleshooting/141112_validate_demultiplexing.sh/141112_validate_demultiplexing.sh
# for log - not done checked with Andrew seems to be ok

# ===================[merging fasta files]=====

cat \
/mnt/paul_folder/140730_18S_data/140801_1st_plate_splitted/seqs.fna \
/mnt/paul_folder/140730_18S_data/140801_2nd_plate_splitted/seqs.fna \
> /mnt/paul_folder/140730_18S_data/140801_unclustered/140801_18S_unclustered.fna

# ===================[preparing for ACAD cluster ] ===== 
# section deleted to minimize clutter - refer to /Users/paul/Documents/140911_c3_analysis/1_analysis_shell_18S_seqdata/140730_18S_analysis.sh for more details


## repetition necessary on 15.01.12 - move old stuff elsewhere - start
# ====================================================================


# ~~~~~~[ OTU picking, rep set picking, taxonomy assignment - below ]~~~~~~

# ~~~~~~[ Needs to be included / updated  in methods section - below - check section "9.12.14-2" (de-novo approach) and "9.12.14-3" (open-ref approach)]~~~~~~~~~~~~~


## 18S OTU picking on Starcluster 15.01.12 - move old stuff elsewhere - start
# ====================================================================

# generate script file for 97% clustering
# ---------------------------------------
echo '
#!/bin/bash -x

# define treshholds
# -----------------
clustering[1]="0.97"
clustering[2]="0.99"

# define paths:
# --------------------

# path to unclustered 18S data
otupickin="/mnt/paul_folder/140730_18S_data/140801_unclustered/140801_18S_unclustered.fna"

# output paths
otupickout[1]="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_denovo097"
otupickout[2]="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_denovo099" 


# load module
# ------------
module load Qiime/1.8.0

echo pick_otus.py -i "$otupickin" -o "${otupickout[1]}" -s "${clustering[1]}" -m uclust -z
' > /mnt/paul_folder/140730_18S_data/150112_18S_denovo_OTU_picking_97.sh

# generate script file for 99% clustering
# ---------------------------------------
echo '
#!/bin/bash -x

# define treshholds
# -----------------
clustering[1]="0.97"
clustering[2]="0.99"

# define paths:
# --------------------

# path to unclustered 18S data
otupickin="/mnt/paul_folder/140730_18S_data/140801_unclustered/140801_18S_unclustered.fna"

# output paths
otupickout[1]="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_denovo097"
otupickout[2]="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_denovo099" 


# load module
# ------------
module load Qiime/1.8.0

echo pick_otus.py -i "$otupickin" -o "${otupickout[2]}" -s "${clustering[2]}" -m uclust -z
' > /mnt/paul_folder/140730_18S_data/150112_18S_denovo_OTU_picking_99.sh

# make executable
chmod +x 150112_18S_denovo_OTU_picking_97.sh && \
chmod +x 150112_18S_denovo_OTU_picking_99.sh

# test scripts and removing echos
nohup ./150112_18S_denovo_OTU_picking_97.sh &> 150112_18S_denovo_OTU_picking_97.nohup.txt
nohup ./150112_18S_denovo_OTU_picking_99.sh &> 150112_18S_denovo_OTU_picking_99.nohup.txt
# started in 2 subshells Mon Jan 12 17:15:21 CST 2015

# = check after completion: 
cat `ls 150112_18S_OTU_denovo*/*.log` | grep "Num OTUs"
# Num OTUs:246075
# Num OTUs:821838

## 18S rep set picking on Starcluster 15.01.13 - at 97% and 99%
# ==============================================================

module load Qiime/1.8.0

# define fixed paths
repsetdir="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_rep_sets"
unclustfna="/mnt/paul_folder/140730_18S_data/140801_unclustered/140801_18S_unclustered.fna"

# define OTU maps
otu_map[1]="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_denovo097/140801_18S_unclustered_otus.txt"
otu_map[2]="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_denovo099/140801_18S_unclustered_otus.txt"

# define output file prefix names
prefix[1]="150113_rep_set_097"
prefix[2]="150113_rep_set_099"

# create output directory"
mkdir -p "$repsetdir"

# run loop over files 
for ((i = 1; i <= "${#otu_map[@]}"; i++)) ; do

# the real command, remove echo if script runs fine
pick_rep_set.py -i "${otu_map[$i]}" -f "$unclustfna" -o "$repsetdir"/"${prefix[$i]}".fasta

done

# = count sequences in files
for i in /mnt/paul_folder/140730_18S_data/150112_18S_OTU_rep_sets/*fasta
do
grep -c ">" "$i"
done

# 246075
# 821838

## 18S rep set alignment for ChimeraSlayer
# ========================================
# if needed, for later, aborted (and crashed)
echo '
#!/bin/bash -x
# 97
module load Qiime/1.8.0
# define file paths
template="/mnt/paul_folder/140730_18S_data/Silva_111_post/rep_set_aligned/97_Silva_111_rep_set.fasta"
infna="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_rep_sets/150113_rep_set_097.fasta"
odir="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_rep_sets_aligned_97"

echo "for debugging only"
echo "$template" 
echo "$infna"  
echo "$odir" 

echo parallel_align_seqs_pynast.py -t "$template" -i "$infna"  -o "$odir" -O 4
' > /mnt/paul_folder/140730_18S_data/150113_repset_align_97.sh

echo '
#!/bin/bash -x
# 99
module load Qiime/1.8.0
# define file paths
template="/mnt/paul_folder/140730_18S_data/Silva_111_post/rep_set_aligned/99_Silva_111_rep_set.fasta"
infna="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_rep_sets/150113_rep_set_099.fasta"
odir="/mnt/paul_folder/140730_18S_data/150112_18S_OTU_rep_sets_aligned_99"

echo "for debugging only"
echo "$template" 
echo "$infna"  
echo "$odir" 

echo parallel_align_seqs_pynast.py -t "$template" -i "$infna"  -o "$odir" -O 4
' > /mnt/paul_folder/140730_18S_data/150113_repset_align_99.sh

# execute Tue Jan 13 11:37:38 CST 2015
nohup ./150113_repset_align_97.sh &> 150113_repset_align_97.nohup.txt
nohup ./150113_repset_align_99.sh $> 150113_repset_align_99.nohup.txt

## 13.01.2015
# Jimmy already conducted rep set picking and taxonomy assignment and now the files are in and will be used instead 
# of clustering results above 

## repetition necessary on 15.01.12 - move old stuff elsewhere - end
# ====================================================================
# older contens of this file have been stored in /Users/paul/Documents/140911_c3_analysis/1_analysis_shell_18S_seqdata/140916_C3_18S_analysis_contents_150112.sh

## 13.01.15 copying older result files from Jimmy to Starcluster
# ==============================================================
local_file="/Users/paul/Desktop/PaulNew.tar.gz" # use full paths with leading slash
remote_file="/mnt/paul_folder/140730_18S_data" # use full paths with leading slash
scp -i /Users/paul/Documents/06_shell/acad_keypair.pem  "$local_file" pczechowski@130.220.209.109:"$remote_file"
tar zxvf PaulNew.tar.gz -C 150113_18S_results_ACAD_cluster_JB
# verification of OTU counts
# = count sequences in files
for i in /mnt/paul_folder/140730_18S_data/150113_18S_results_ACAD_cluster_JB/*/*.log
do
grep  "Num OTUs" "$i"
done
# Num OTUs:246 075
# Num OTUs:821 838
# numbers match with OTU picking conducted on cluster

## 13.01.15 taxonomy assignment of 18S rep. set
# =============================================
# done by JB on ACAD cluster not conducted here, older results are used

echo '
#!/bin/bash -x

# =====[ SCRIPT start ]======

# taxonomy assignemnt of Illumina 18S sequence data (C3) to 100% SILVA refernce DB

# modified from source file of script:
# /Users/paul/Documents/140911_c3_analysis/1_analysis_shell_18S_seqdata/140916_18S_analysis.sh

# use Uclust because it is the fastest method (suitable for large data sets)
# http://qiime.org/scripts/parallel_assign_taxonomy_uclust.html

# use 3 different assignment treshholds to compare results:
# --uclust_similarity 0.90
# --uclust_similarity 0.95
# --uclust_similarity 0.99

# use two different rep sets:
# ----------------------------
# clustered at 97%: /home/users/jbreen/2014PC/140930_18S_rep_set/140930_rep_set_097.fasta
# clustered at 99%: /home/users/jbreen/2014PC/140930_18S_rep_set/140930_rep_set_097.fasta
# (100%: too many OTUs - too much noise in data)

# use unclusted reference database
# --------------------------------
# database creation described in:
# /mnt/paul_folder/140730_18S_data/Silva_111_post/SILVA_111_QIIME_format_notes.txt


# 18S ref db files (mt, animal, unclustered):
# -------------------------------------
# /home/users/jbreen/2014PC/140930_18S_reference_data/Silva_111_taxa_map_no_ambiguous.txt
# /home/users/jbreen/2014PC/140930_18S_reference_data/Silva_111_full_unaligned_no_ambig.fasta 
# 13.01.2015: or equivalent paths 

# declare variables
# -----------------

# ref db (unclustered):
# 13.01.2015: or equivalent paths
refdb100tax="/mnt/paul_folder/140730_18S_data/141009_ACAD_cluster_test/141009_18S_reference_data/Silva_111_taxa_map_no_ambiguous.txt"
refdb100seq="/mnt/paul_folder/140730_18S_data/141009_ACAD_cluster_test/141009_18S_reference_data/Silva_111_full_unaligned_no_ambig.fasta"

# rep sets (omitting 100% clustered rep set, using only 97%, 99%):
# 13.01.2015: or equivalent paths
repset[1]="/mnt/paul_folder/140730_18S_data/141009_ACAD_cluster_test/141009_18S_rep_set/141009_rep_set_097.fasta"
repset[2]="/mnt/paul_folder/140730_18S_data/141009_ACAD_cluster_test/141009_18S_rep_set/141009_rep_set_099.fasta"

# define taxonomy assignment thresholds in array (testing differnt assignment treshholds)
threshold[1]="0.90"
threshold[2]="0.95"
threshold[3]="0.99"

# define output destinations matching array position in "threshhold"
# 13.01.2015: or equivalent paths
destination[1]="/mnt/paul_folder/140730_18S_data/141009_ACAD_cluster_test/141009_18S_tax_assignment_90"
destination[2]="/mnt/paul_folder/140730_18S_data/141009_ACAD_cluster_test/141009_18S_tax_assignment_95"
destination[3]="/mnt/paul_folder/140730_18S_data/141009_ACAD_cluster_test/141009_18S_tax_assignment_99"

# define output destinations 
# 13.01.2015: or equivalent paths
dir[1]="141009_18S_rep_set_097_vs_db_100"
dir[2]="141009_18S_rep_set_099_vs_db_100"

# run nested loops 

# loop through rep sets (2)

for ((i = 1; i <= "${#repset[@]}"; i++))
	do
		# loop through threshholds (3)
		for ((k = 1; k <= "${#threshold[@]}"; k++))
			do
			# for debugging only 
			date
			echo  ~~~~~~~~~~~~~
			echo  - outer loop: "${repset[$i]}"
			echo  - inner loop: "${threshold[$k]}"
			echo  ~~~~~~~~~~~~~
			# ADJUST threads (-o integer) - each thread used about 12 GB Ram for this step 
			parallel_assign_taxonomy_uclust.py -i "${repset[$i]}" -o "${destination[$k]}"/"${dir[$i]}" -r "$refdb100seq" -t "$refdb100tax" --uclust_similarity "${threshold[$k]}" -O 1
			# closing inner loop
		done 
			# closing outer loop 
done


# =====[ SCRIPT end ]======
' > /Users/paul/Documents/140523_c2_analysis_re/141009_ACAD_cluster_scripts/141009_18S_tax_assignment_checked.sh

## 13.1.15: verify rep set picking and tax assignment results
# ===========================================================

# rep set picking of de-novo workflow
# -----------------------------------
# count rep-set sequences of JBs workflow
grep ">" -c /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_rep_set/141211_rep_set_097.fasta
# 246 075
grep ">" -c /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_rep_set/141211_rep_set_099.fasta
# 821 838

# count rep-set sequences of workflow I created independently 

grep ">" -c  /mnt/paul_folder/140730_18S_data/150112_18S_OTU_rep_sets/150113_rep_set_097.fasta
# 246 075
grep ">" -c  /mnt/paul_folder/140730_18S_data/150112_18S_OTU_rep_sets/150113_rep_set_099.fasta
# 821 838

# tax assignemnt de-novo workflow
# -------------------------------
# count taxonomy assignments:

# http://stackoverflow.com/questions/9612090/how-to-loop-list-of-file-names-returned-by-find
# only without whitespaces - careful!
for i in $(find /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo -type f -name "*.txt"); do
    echo -e  '\n-Line count:'
    wc -l   "$i"
    echo -e  '\n-Unassigned:'
    cut -f2 "$i" | sort -d | uniq -c | tail -n1
done

# -Line count:
# 246075 /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_tax_assignment_90/141211_18S_rep_set_097_vs_db_100/141211_rep_set_097_tax_assignments.txt
# -Unassigned:
#  42 098 Unassigned,
# 203 977 Assigned
# 
# -Line count:
# 821838 /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_tax_assignment_90/141211_18S_rep_set_099_vs_db_100/141211_rep_set_099_tax_assignments.txt
# -Unassigned:
# 129 654 Unassigned, 
# 692 184 Assigned
#
# -Line count:
# 246075 /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_tax_assignment_95/141211_18S_rep_set_097_vs_db_100/141211_rep_set_097_tax_assignments.txt
# -Unassigned:
# 121 450 Unassigned, 
# 124 625 Assigned
# 
# -Line count:
# 821838 /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_tax_assignment_95/141211_18S_rep_set_099_vs_db_100/141211_rep_set_099_tax_assignments.txt
# 294 317 Unassigned, 
# 526 824 Assigned
#
# -Line count:
# 246075 /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_tax_assignment_99/141211_18S_rep_set_097_vs_db_100/141211_rep_set_097_tax_assignments.txt
# 245 014 Unassigned
#   1 061 Assigned
# 
# -Line count:
#  /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_tax_assignment_99/141211_18S_rep_set_099_vs_db_100/141211_rep_set_099_tax_assignments.txt
#  811 007 Unassigned
#   10 831 Assigned

# 7.1.14: verify rep set picking and  tax assignment results - end
# ===================================================================


# 13.1.14 OTU table generation - start
# ==================================

# get list with taxonomy assignments ( resulting txt file from assign_taxonomy.py for -t TAXONOMY_FNAME)
find /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo -type f -name "*.txt"

module load Qiime/1.8.0

# loop
for file in /mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_tax_assignment_*/141211_18S_rep_set_*_vs_db_100/141211_rep_set_*_tax_assignments.txt; do	
	case "$file" in
		*_rep_set_097_vs_db_100* ) otumap="/mnt/paul_folder/140730_18S_data/150113_18S_results_ACAD_cluster_JB/140930_18S_OTU_097/141009_18S_unclustered_otus.txt" && clustering="_clust97" ;;&
		*_rep_set_099_vs_db_100* ) otumap="/mnt/paul_folder/140730_18S_data/150113_18S_results_ACAD_cluster_JB/140930_18S_OTU_099/141009_18S_unclustered_otus.txt" && clustering="_clust99" ;;&
		*141211_18S_tax_assignment_90* ) taxassgn="_tassgn90" ;;&
		*141211_18S_tax_assignment_95* ) taxassgn="_tassgn95" ;;&
		*141211_18S_tax_assignment_99* ) taxassgn="_tassgn99" ;;&
	esac
	# echo -t "$file" # print input taxonomy assignment file
	# echo -i "$otumap" # print input OTU table
	biom_path="/mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/150113_18_OTUs$clustering$taxassgn.biom"
	# echo -o "$biom_path"
	make_otu_table.py -v -i "$otumap" -t "$file" -o "$biom_path"
	echo 
done

# 13.1.14 OTU table generation - end
# ==================================


# 13.1.14 Metadata addition - start
# ==================================
## edited 15.01.10

# define input path names
mf="/mnt/paul_folder/140730_18S_data/150108_18S_mf_w_md/150108_18S_MF.txt"

# ading metadata from mapping file 
for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/*.biom; do
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
for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/*.txt; do
echo "$file" >> /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/150113_18S_all_summaries.txt
cat  "$file" >> /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/150113_18S_all_summaries.txt
done

# 13.1.14 Metadata addition - end
# ==================================


# 13.1.14 18S table filetering - remove unassigned - start
# =======================================================

for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/*_md.biom; do
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
# for now all summaries stored in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/150113_18S_all_summaries_assigned.txt
for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/*md.sum_*.txt; do
echo "$file" >> /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/150113_18S_all_summaries_assigned.txt
cat  "$file" >> /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/150113_18S_all_summaries_assigned.txt
done

# verifification of decreasing file sizes - ok:
cd /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables
ls *tassgn??_md.biom -lh
ls *tassgn??_md_*.biom -lh

# 13.1.15 18S table filetering - remove unassigned - end
# =======================================================


# 13.1.15 18S table filetering - substract blanks - start
# =======================================================

# substract blanks in biom tables that have metadata associated with them i.e. these 
ls /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/*md_assigned_only.biom

# define output directory and create it 
targetdir="/mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks/"
mkdir -p  $targetdir

# define mapping file path
mf="/mnt/paul_folder/140730_18S_data/150108_18S_mf_w_md/150108_18S_MF.txt"

module load Qiime/1.8.0

# the filetring commands in the loop 

# testing:  for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/150113_18_OTUs_clust97_tassgn99_md_assigned_only.biom; do

for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables/*md_assigned_only.biom; do
[[ -e $file ]] || continue

echo "$file"

# Create an OTU table with the just the blank control samples:

filter_samples_from_otu_table.py -i "$file" -o "$targetdir"$(basename "$file" .biom).biom_tmp_a1 -m "$mf" -s "XtrOri:PCR_blk_h20"
echo "1/9"

filter_samples_from_otu_table.py -i "$file" -o "$targetdir"$(basename "$file" .biom).biom_tmp_b1 -m "$mf" -s "XtrOri:XTR_blk_multiple"
echo "2/9"

# Filter out OTU ids that have zero counts, as we only want the OTUs with positive counts from the Control_Blank samples:

filter_otus_from_otu_table.py -i "$targetdir"$(basename "$file" .biom).biom_tmp_a1 -o "$targetdir"$(basename "$file" .biom).biom_tmp_a2 -n 1
echo "3/9"

filter_otus_from_otu_table.py -i "$targetdir"$(basename "$file" .biom).biom_tmp_b1 -o "$targetdir"$(basename "$file" .biom).biom_tmp_b2 -n 1
echo "4/9"


# Then create a tab separated version of this OTU table:

biom convert -b -i "$targetdir"$(basename "$file" .biom).biom_tmp_a2 -o "$targetdir"$(basename "$file" .biom).txt_tmp_a3
echo "5/9"

biom convert -b -i "$targetdir"$(basename "$file" .biom).biom_tmp_b2 -o "$targetdir"$(basename "$file" .biom).txt_tmp_b3
echo "6/9"

# Filter out OTU ids from the run 1 OTU table that were determined to be present in the Control_Blank samples:

filter_otus_from_otu_table.py -i "$file" -o "$targetdir"$(basename "$file" .biom).biom_tmp_4 -e "$targetdir"$(basename "$file" .biom).txt_tmp_a3
echo "7/9"

filter_otus_from_otu_table.py -i "$targetdir"$(basename "$file" .biom).biom_tmp_4 -o "$targetdir"$(basename "$file" .biom).biom_tmp_5 -e "$targetdir"$(basename "$file" .biom).txt_tmp_b3
echo "8/9"

# The otu_table_1_minus_contaminants.biom file now has two samples with zero sequences associated with it. These can be removed to get a final OTU table:
filter_samples_from_otu_table.py \
-i "$targetdir"$(basename "$file" .biom).biom_tmp_5 \
-o "$targetdir"$(basename "$file" .biom)_no_contaminants.biom -n 1
echo "9/9"

done

# 13.01.15: generate and read summaries on Starcluster
#
cd /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks
for file in ./*.biom; do
[[ -e $file ]] || continue
biom summarize-table -i  "$file" -o $(basename "$file" .biom).sum_qual.txt --qualitative
biom summarize-table -i  "$file" -o $(basename "$file" .biom).sum_quan.txt
done


# for now all summaries stored in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks/150113_18S_all_summaries_assigned_no_blanks.txt
for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks/*.txt; do
echo "$file" >> /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks/150113_18S_all_summaries_assigned_no_blanks.txt
cat  "$file" >> /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks/150113_18S_all_summaries_assigned_no_blanks.txt
done



# validating size decrease
# 13.01.15: ok
cd /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables && ls *tassgn??_md_*.biom -lh
cd /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks && ls *tassgn??_md_*.biom -lh

# removing temp files
# 13.01.15: ok

rm -v /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks/*_tmp_*


# 13.1.15 18S table filetering - substract blanks - end
# =======================================================


# 14.1.15 18S table filetering - rmeove OTUs less then 5 clusters, retain invertebrates - start 
# =========================================================================================


## find all  unique taxonomy strings in one assignemnet file - most OTUs, lowest treshholds, should be best IDs in here
target="/mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_tax_assignment_90/141211_18S_rep_set_099_vs_db_100/141211_rep_set_099_tax_assignments.txt"
grep "denovo"  "$target" | awk -F '[\t]' '{print $2 }' | sort -d | uniq -c

# Eukaryota; __Opisthokonta; __Metazoa; __Nematoda; __Chromadorea; __Monhysteridae; __Halomonhystera; __Halomonhystera_disjunct
# Eukaryota; __Opisthokonta; __Metazoa; __Rotifera; __Philodinidae; __Adineta; __Adineta_vaga
# Eukaryota; __Opisthokonta; __Metazoa; __Tardigrada; __Echiniscidae; __Echiniscus; __Echiniscus_canadensis
# Eukaryota  __Opisthokonta  __Metazoa  __Arthropoda  __Hexapoda

# find all Metazoans
awk -F '[\t;]' '$5 == " __Metazoa" {print $2,$3,$4,$5,$6 }' "$target" | sort -d | uniq -c

# find all relavant invertebrates
awk -F '[\t;]' '$5 == " __Tardigrada" || $5 == " __Rotifera" || $5 == " __Arthropoda" || $5 == " __Nematoda" || $5 == " __Rotifera"  {print $2,$3,$4,$5,$6 }' "$target" | sort -d | uniq -c

# create output directory
targetdir="/mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only"
mkdir -p $targetdir

#check input file paths - choose only files with metadata and taxobnomy assignments
ls /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks/*.biom

# load qiime
module load Qiime/1.8.0

# loop
for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks/*.biom; do
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
filter_taxa_from_otu_table.py -i  "$out_tmp1" -o  "$out_tmp2" -p __Nematoda,__Rotifera,__Tardigrada,__Arthropoda
echo 2/4 - retaining invertebrates

# filter out zero count samples
filter_samples_from_otu_table.py -i "$out_tmp1" -o "$out_biom1" -n 1
echo 3/4 - filter 0s from all taxa

# filter out zero count sumples
filter_samples_from_otu_table.py -i "$out_tmp2" -o "$out_biom2" -n 1
echo 4/4 - filter 0s from invertebrates

done

for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only/*.biom; do
[[ -e $file ]] || continue

# generate summary files
echo "$file"
biom summarize-table -i  "$file" -o $(basename "$file" .biom).sum_qual.txt --qualitative
biom summarize-table -i  "$file" -o $(basename "$file" .biom).sum_quan.txt

done

# concatenate invert summary files
out="/mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only/150114_inverts_sum_quan_all.txt"
for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only/*inverts.sum_quan.txt; do
[[ -e $file ]] || continue
echo "$file" >> "$out"
cat "$file" >> "$out"
done

out="/mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only/150114_inverts_sum_qual_all.txt"
for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only/*inverts.sum_qual.txt; do
[[ -e $file ]] || continue
echo "$file" >> "$out"
cat "$file" >> "$out"
done

# erase temp files
rm -v /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only/*.tmp?

# check size decrease
# 13.01.15: ok
cd /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks && ls *tassgn??_md_*.biom -lh
cd /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only && ls *_above5.biom -lh
cd /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only && ls *_above5_inverts.biom -lh


# 14.1.15 18S table filetering - rmeove OTUs less then 5 clusters, retain invertebrates - end 
# =========================================================================================


## 14.1.15 18S table filtering - retain Antarctic retain controls - start 
# =======================================================================

## check mapping file - which categories to filter?

# convert line endings
mfmac="/mnt/paul_folder/140730_18S_data/150108_18S_mf_w_md/150108_18S_MF.txt"
mfunix="/mnt/paul_folder/140730_18S_data/150108_18S_mf_w_md/150108_18S_MF_unix.txt"
awk '{ gsub("\r", "\n"); print $0;}' "$mfmac" > "$mfunix"

# Print unique relevant fields
awk -F '\t' '{print $9}' "$mfunix" | sort -d | uniq -d

#retain 
Location:E_Ant_coast # AVC controls
Location:Australia # soil controls and mock mixes

# load qiime
module load Qiime/1.8.0

# file list of relavant biom tables
ls /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only/*_above5_inverts.biom

# define output directory and create it 
targetdir="/mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only_controls"
mkdir -p  $targetdir

# define mapping file path
mf="/mnt/paul_folder/140730_18S_data/150108_18S_mf_w_md/150108_18S_MF_unix.txt"

# the filetring commands in the loop  - did not run weel 1st time hence commented sections
for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only/*_above5_inverts.biom; do
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
rm -vi /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only_controls/*.tmp?

done

# concatenate summarie files

# for now all summaries stored in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only_controls/150114_all_summaries_md_assigned_only_inverts_controls.txt

cd /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_inverts_only_controls
for file in ./*.sum_*.txt; do
echo "$file" >> 150114_all_summaries_md_assigned_only_inverts_controls.txt
cat  "$file" >> 150114_all_summaries_md_assigned_only_inverts_controls.txt
done

# 12.1.15 18S table filtering - retain controls - end 
# ==========================================================

# backup; continued in /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_2.sh