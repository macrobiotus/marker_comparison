## 16.01.27 - /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_2.sh
#  ===========================================================================================================================
# shell analysis script merging results from:
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_1a_18S.sh
# and: 
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/160127_data_analysis_1b_COI.sh
# for R analysis:
# /Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/S9_Data_analysis.R
# to be located in upload folder 

## performed steps:
#  ------------------
#  Step 1: Get Invertebrate OTU tables that are not filtered on abundance
#  Step 2: filter by different abundances and sort data (these creates the biom tables provided in the upload folder)
#  Step 3: Do all other analyses in R script


## 15.04.29 - Re-analysis for Chapter 3 - Shell part
## =================================================
# modified from /Users/paul/Documents/140911_c3_analysis/1_analysis_shell_COI-18S-AVC/150504_C3_old_scripts/150302_C3_metagenetic_approach.txt
# see R part: /Users/paul/Documents/140911_c3_analysis/2_analysis_R_scripts/150504_C3_revisions.r

## Step 1: Get Invertebrate OTU tables that are not filtered on abundance
## ======================================================================
## This step is okay 
# filter for inverts
# copy to local

## 18S:
## ----
# find all Metazoans
target="/mnt/paul_folder/140730_18S_data/141211_18S_DeNovo/141211_18S_tax_assignment_90/141211_18S_rep_set_099_vs_db_100/141211_rep_set_099_tax_assignments.txt"
awk -F '[\t;]' '$4 == " __Metazoa" {print $5,$6 }' "$target" | sort -d | uniq -c
# create output directory
targetdir="/mnt/paul_folder/140730_18S_data/150302_18S_phylotypes_invertebrates"
mkdir -p $targetdir
# check input file paths - choose only files with metadata and taxobnomy assignments
ls /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks/*.biom
# load qiime
module load Qiime/1.8.0
# loop
for file in /mnt/paul_folder/140730_18S_data/150113_18S_OTU_tables_substracted_blanks/*.biom; do
	[[ -e "$file" ]] || continue
	echo @ processed file: "$file"
	out_tmp1=$targetdir/$(basename $file .biom).tmp1
	echo @ temp file 1 "$out_tmp1"
	out_tmp2=$targetdir/$(basename $file .biom).tmp2
	echo @ temp file 2 "$out_tmp2"
	out_biom="$targetdir"/$(basename "$file" .biom)_invertebrates.biom
	echo @ biom file 2 "$out_biom"
	echo @ retaining invertebrates
	# filter_taxa_from_otu_table.py -i  "$file" -o  "$out_tmp1" -p __Nematoda,__Rotifera,__Tardigrada,__Chelicerata
	filter_taxa_from_otu_table.py -i  "$file" -o  "$out_tmp1" -p __Nematoda,__Rotifera,__Tardigrada,__Arthropoda
	echo @ filter empty samples
	filter_samples_from_otu_table.py -i "$out_tmp1" -o "$out_tmp2" -n 1
	echo @ filter 5 count OTUs
	filter_otus_from_otu_table.py -i "$out_tmp2" -o "$out_biom" -n 5
	echo @ generating summary files 
	biom summarize-table -i  "$out_biom" -o "$targetdir"/$(basename "$out_biom" .biom).sum_qual.txt --qualitative
	biom summarize-table -i  "$out_biom" -o "$targetdir"/$(basename "$out_biom" .biom).sum_quan.txt
	echo @ erasing temp files
	rm -v /mnt/paul_folder/140730_18S_data/150302_18S_phylotypes_invertebrates/*.tmp?
done

## COI:
## ----
# create output directory
targetdir="/mnt/paul_folder/140719_COI_data/150302_COI_phylotypes_invertebrates"
mkdir -p $targetdir
# check input file paths - choose only files with metadata and taxobnomy assignments
ls /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/*md_assigned_only.biom
# load qiime
module load Qiime/1.8.0
# loop
for file in /mnt/paul_folder/140719_COI_data/151021_COI_OTU_tables/*md_assigned_only.biom; do
	[[ -e "$file" ]] || continue
	echo @ processed file: "$file"
	out_tmp1=$targetdir/$(basename $file .biom).tmp1
	echo @ temp file 1 "$out_tmp1"
	out_tmp2=$targetdir/$(basename $file .biom).tmp2
	echo @ temp file 2 "$out_tmp2"
	out_biom="$targetdir"/$(basename "$file" .biom)_invertebrates.biom
	echo @ biom file 2 "$out_biom"
	echo @ retaining invertebrates
	# filter_taxa_from_otu_table.py -i  "$file" -o  "$out_tmp1" -p Chelicerata,Nematoda,Rotifera,Tardigrada
	filter_taxa_from_otu_table.py -i  "$file" -o  "$out_tmp1" -p Arthropoda,Nematoda,Rotifera,Tardigrada
	echo @ filter empty samples
	filter_samples_from_otu_table.py -i "$out_tmp1" -o "$out_tmp2" -n 1
	echo @ filter 5 count OTUs
	filter_otus_from_otu_table.py -i "$out_tmp2" -o "$out_biom" -n 5
	echo @ generating summary files # update pathnames
	biom summarize-table -i  "$out_biom" -o "$targetdir"/$(basename "$out_biom" .biom).sum_qual.txt --qualitative
	biom summarize-table -i  "$out_biom" -o "$targetdir"/$(basename "$out_biom" .biom).sum_quan.txt
	echo @ erasing temp files
	rm -v /mnt/paul_folder/140719_COI_data/150302_COI_phylotypes_invertebrates/*.tmp?
done

## copy to local
# 18S:
local_file="/Users/paul/Documents/140911_c3_analysis/150302_Invertebrates/150302_18S" 
remote_file="/mnt/paul_folder/140730_18S_data/150302_18S_phylotypes_invertebrates/*" 
scp -r -C -i /Users/paul/Documents/06_shell/acad_keypair.pem  pczechowski@130.220.209.109:"$remote_file" "$local_file"
# COI:
local_file="/Users/paul/Documents/140911_c3_analysis/150302_Invertebrates/150302_COI"
remote_file="/mnt/paul_folder/140719_COI_data/150302_COI_phylotypes_invertebrates/*"
scp -r -C -i /Users/paul/Documents/06_shell/acad_keypair.pem  pczechowski@130.220.209.109:"$remote_file" "$local_file"


## Step 2: filter by different abundances and sort data
## ======================================================================
## This step is okay

# filtering percentages
percentage[1]="0.001"
percentage[2]="0.002"
percentage[3]="0.003"
percentage[4]="0.005"
# input file lists
list[1]="/Users/paul/Documents/140911_c3_analysis/150302_Invertebrates/150302_18S/*.biom"
list[2]="/Users/paul/Documents/140911_c3_analysis/150302_Invertebrates/150302_COI/*.biom"
# mapping files
mf[1]="/Users/paul/Documents/140911_c3_analysis/mapping_files/150107_mf_metadata/150108_18S_MF.txt"
mf[2]="/Users/paul/Documents/140911_c3_analysis/mapping_files/150107_mf_metadata/150108_COI_MF.txt"
# define path and filenames
outdir[1]="/Users/paul/Documents/140911_c3_analysis/150302_Invertebrates/150304_18S_abundance"
outdir[2]="/Users/paul/Documents/140911_c3_analysis/150302_Invertebrates/150304_COI_abundance"
# parameter file
pf="/Users/paul/Documents/140911_c3_analysis/2_analysis_Qiime_PFs/150215_plot_pf.txt"
# loop over percentages
for ((k = 1; k <= "${#percentage[@]}"; k++)); do
	echo "@ enter loop filtering out fraction: ${percentage[k]}"
	# loop over lists
	for ((j = 1; j <= "${#list[@]}"; j++)); do
		echo "@ enter loop list: ${list[j]}"
		# loop over files
		for file in ${list[j]}; do
			[[ -e $file ]] || continue
			echo "@ enter loop table: $file"
			# filter by abundace
			abdir="${outdir[j]}"/150304_"${percentage[k]}"_discarded && mkdir -p "$abdir"
			echo "@ defined and created abundance directory: $abdir"
			abundance="$abdir"/$(basename "$file" .biom).tmp
			echo "@ defined filename: $abundance"
			echo "@ retaining percentage"
			filter_otus_from_otu_table.py --min_count_fraction "${percentage[k]}" -i "$file" -o "$abundance"
			# retain controls in single file
			antarctic="$abdir"/$(basename "$abundance" .tmp)_ANT_"${percentage[k]}".tmp
			echo "@ defined filename: $antarctic"
			australia="$abdir"/$(basename "$abundance" .tmp)_AUST_"${percentage[k]}".tmp
			echo "@ defined filename: $australia"
			all_data="$abdir"/$(basename "$abundance" .tmp)_ALL_"${percentage[k]}".tmp
			echo "@ defined filename: $all_data"
			echo "@ retaining Antarctic and Austrlian controls"
			filter_samples_from_otu_table.py -i "$abundance" -o "$antarctic" -m "${mf[j]}" -s "Location:E_Ant_coast"
			filter_samples_from_otu_table.py -i "$abundance" -o "$australia" -m "${mf[j]}" -s "Location:Australia"
			merge_otu_tables.py -i "$australia","$antarctic" -o "$all_data"
			# isolate soil controls from Australian samples
			aust_soil="$abdir"/$(basename "$australia" .tmp)_SOIL.tmp
			echo "@ defined filename: $aust_soil"
			echo "@ isolating Austrlian soil controls"
			filter_samples_from_otu_table.py -i "$australia" -o "$aust_soil" -m "${mf[j]}" -s "XtrOri:XTR_pos"
			# isolate mock controls from Australian samples
			aust_mock="$abdir"/$(basename "$australia" .tmp)_MOCK.tmp
			echo "@ defined filename: $aust_mock"
			echo "@ isolating Austrlian mock controls"
			filter_samples_from_otu_table.py -i "$australia" -o "$aust_mock" -m "${mf[j]}" -s "XtrOri:PCR_pos"
			# defining output filenames
			all_data_b="$abdir"/$(basename "$all_data" .tmp).biom
			echo "@ defined output filename: $all_data_b"
			antarctic_b="$abdir"/$(basename "$antarctic" .tmp).biom
			echo "@ defined output filename: $antarctic_b"
			aust_soil_b="$abdir"/$(basename "$aust_soil" .tmp).biom
			echo "@ defined output filename: $aust_soil_b"
			aust_mock_b="$abdir"/$(basename "$aust_mock" .tmp).biom
			echo "@ defined output filename: $aust_mock_b"
			# final filtering
			echo "@ removing 0 count OTUs"
			filter_otus_from_otu_table.py -i "$all_data" -o "$all_data_b" -n 1
			filter_otus_from_otu_table.py -i "$antarctic" -o "$antarctic_b" -n 1
			filter_otus_from_otu_table.py -i "$aust_soil" -o "$aust_soil_b" -n 1
			filter_otus_from_otu_table.py -i "$aust_mock" -o "$aust_mock_b" -n 1
			# convert to text
			echo "@ converting to text"
			all_data_txt="$abdir"/$(basename "$all_data_b" .biom).txt
			antarctic_txt="$abdir"/$(basename "$antarctic_b" .biom).txt
			aust_soil_txt="$abdir"/$(basename "$aust_soil_b" .biom).txt
			aust_mock_txt="$abdir"/$(basename "$aust_mock_b" .biom).txt
			biom convert -i "$all_data_b" -o "$all_data_txt" -b --header-key taxonomy --output-metadata-id "Taxonomy"
			biom convert -i "$antarctic_b" -o "$antarctic_txt" -b --header-key taxonomy --output-metadata-id "Taxonomy"
			biom convert -i "$aust_soil_b" -o "$aust_soil_txt" -b --header-key taxonomy --output-metadata-id "Taxonomy"
			biom convert -i "$aust_mock"  -o "$aust_mock_txt" -b --header-key taxonomy --output-metadata-id "Taxonomy"
			# summaries
			echo "@ getting summaries"
			biom summarize-table -i "$all_data_b" -o "$abdir"/$(basename "$all_data_b" .biom).sum_qual.txt --qualitative
			biom summarize-table -i "$all_data_b" -o "$abdir"/$(basename "$all_data_b" .biom).sum_quan.txt
			biom summarize-table -i "$antarctic_b" -o "$abdir"/$(basename "$antarctic_b" .biom).sum_qual.txt --qualitative
			biom summarize-table -i "$antarctic_b" -o "$abdir"/$(basename "$antarctic_b" .biom).sum_quan.txt
			biom summarize-table -i "$aust_soil_b" -o "$abdir"/$(basename "$aust_soil_b" .biom).sum_qual.txt --qualitative
			biom summarize-table -i "$aust_soil_b" -o "$abdir"/$(basename "$aust_soil_b" .biom).sum_quan.txt
			biom summarize-table -i "$aust_mock_b" -o "$abdir"/$(basename "$aust_mock_b" .biom).sum_qual.txt --qualitative
			biom summarize-table -i "$aust_mock_b" -o "$abdir"/$(basename "$aust_mock_b" .biom).sum_quan.txt
			# erasing temp files
			echo "@ erasing temp files"
			if [ -d "$abdir" ]; then
			rm -v "$abdir"/*.tmp
			fi
		done
	done
done

## Step 3: Do all other analyses
## ======================================================================================

## R script at: 
## /Users/paul/Documents/140911_c3_analysis/2_analysis_R_scripts/150504_C3_revisions.r
