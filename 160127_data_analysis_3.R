## 16.02.09 - Data analysis using R
## ================================

# Paul Czechowski, paul.czechowski(at)adelaide(dot)edu(dot)au

# 1: Import all biom files derived from soil data

# - Evaluation of Parameters:
# 2: Plot effect of filtering on 18S and COI soil data
# 3: Get tables with highest numbers of Phylotype counts into new list
# 4: Get Similarity values between replicates of two soil controls 
# 5: Combine count and similarity values
# 6: * Create display Items: Plot Phylotype counts and Similarity between replicates

# - Evaluation of Insect controls:
# 7: Generate comparale objects of insect controls for plotting
# 7a: Plot composition of "Artificial blends" reference data
# 8: Compare similarities between insect reference and insect controls
# 9:  * Create display Items: Heatmap as comparsion between samples
# 10: * Create display Items: Barplot  as comparsion between samples

# - Evaluation of Antarctic samples:
# 11: Import Antarctic morphology data and metagenetic samples
# 11a: Plot composition of "Antarctic soil" reference data
# 12: compare assignemnts
# 13: plot taxonomic composition
# 14: Format Antarctic data into a more useful data frame

# - Stats

# 15: Calculate ICC on Control and Antarctic data

### clear environment, set working directory
rm(list=ls()) # clear R environment
setwd("/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo") # working directory is here 

### load packages
require("phyloseq")
require("biom")
require("plyr")
require("dplyr")
require("vegan")
require("gplots")
require("ade4")
require("ggplot2")
require("reshape2")
require("gridExtra")
require("foreach")
require("irr")

### get packages citations, were available
capture.output(utils:::print.bibentry(citation("phyloseq"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("biom"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("plyr"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("dplyr"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("vegan"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("gplots"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("ade4"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("ggplot2"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("reshape2"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("gridExtra"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("foreach"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)
capture.output(utils:::print.bibentry(citation("irr"), style = "Bibtex"), file = "160621_package_citations.bib", append = TRUE)

########### 1: data import
# Find all .biom files and store in lists, combine list, name list itemes
# data sources updated to upload folder, compare file 160121_file_collation_for_upload.txt
files_18S <- list.files("/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_18S/", full.names=TRUE, pattern="*.biom")
files_COI <- list.files("/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_COI/", full.names=TRUE, pattern="*.biom")
# filenames become names of list items
names(files_18S) <- files_18S[] 
names(files_COI) <- files_COI[] 

# Make a function to process each file
# needs file list, returns biom objects in list
createPhyloseq <- function(list_item) {
	# import biom file with pathname in list 
	physeq_ob <- try(import_biom(list_item))
	# return biom object
	return(physeq_ob)
}

# read in biom files and generate a list of Phyloseq objects
psob_list_18S <- lapply(files_18S, createPhyloseq)
psob_list_COI <- lapply(files_COI, createPhyloseq)

# remove failed imports form list
psob_list_COI[which(names(psob_list_COI) %in% c("/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_COI//150121_COI_OTUs_clust97_tassgn75_md_assigned_only_invertebrates_AUST_0.005_SOIL.biom",
	"/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_COI//150121_COI_OTUs_clust97_tassgn80_md_assigned_only_invertebrates_AUST_0.002_SOIL.biom",
	"/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_COI//150121_COI_OTUs_clust97_tassgn80_md_assigned_only_invertebrates_AUST_0.003_SOIL.biom",
	"/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_COI//150121_COI_OTUs_clust97_tassgn80_md_assigned_only_invertebrates_AUST_0.005_SOIL.biom",
	"/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_COI//150121_COI_OTUs_clust97_tassgn90_md_assigned_only_invertebrates_AUST_0.001_SOIL.biom",
	"/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_COI//150121_COI_OTUs_clust99_tassgn75_md_assigned_only_invertebrates_AUST_0.003_SOIL.biom",
	"/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_COI//150121_COI_OTUs_clust99_tassgn75_md_assigned_only_invertebrates_AUST_0.005_SOIL.biom",
	"/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_COI//150121_COI_OTUs_clust99_tassgn80_md_assigned_only_invertebrates_AUST_0.002_SOIL.biom",
	"/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/0_AU_SOIL_biom_files_COI//150121_COI_OTUs_clust99_tassgn90_md_assigned_only_invertebrates_AUST_0.001_SOIL.biom"))] <- NULL

# finished data - lists of PhyloSeq objects:
psob_list_18S
psob_list_COI

########### 2: Plot effect of filtering
# get phylotype count in each table and store in integer (vector?)
ets_counts <- sapply(psob_list_18S, function(x){nrow(otu_table(x))})
coi_counts <- sapply(psob_list_COI, function(x){nrow(otu_table(x))})

# function: simplify the names of each element - 18S 
ShortenNames18S <- function(ct) { 
	# substring command needs to be adopted for path length
	# everything before "18S_OTUs" or "COI_OTUs" needs to cut out
	# otherwise the graph labeling will be stuffed 
	names(ct) <- gsub("_md_assigned_only_no_contaminants_invertebrates_", " ", substr(names(ct), 121, 202))
	names(ct) <- gsub("_OTUs_clust", "S c", names(ct))
	names(ct) <- gsub("_tassgn", " t", names(ct))
	names(ct) <- gsub("_tassgn", " t", names(ct))
	names(ct) <- gsub("AUST_", " f", names(ct))
	return(ct)
	}

# function: simplify the names of each element - COI - needs fixing for new data locations
ShortenNamesCOI <- function(ct) { 
	# substring command needs to be adopted for path length
	# everything before "18S_OTUs" or "COI_OTUs" needs to cut out
	# otherwise the graph labeling will be stuffed
	names(ct) <- gsub("__md_assigned_only_invertebrates_", " ", substr(names(ct), 121, 187))
	names(ct) <- gsub("_OTUs_clust", " c", names(ct))
	names(ct) <- gsub("_tassgn", " t", names(ct))
	names(ct) <- gsub("_tassgn", " t", names(ct))
	names(ct) <- gsub("AUST_", " f", names(ct))
	names(ct) <- gsub("_md_assigned_only_invertebrates_", " ", names(ct))
	return(ct)
	}

# apply name shortening commands
ets_counts <- ShortenNames18S(ets_counts)
coi_counts <- ShortenNamesCOI(coi_counts)

# function: convert count per setting to data frame suitable for plottig
MakeSortedDF <- function(ct){
	ct <- data.frame(ct)
	colnames(ct) <- "count"
	ct[2] <- as.character(rownames(ct)) 	
	colnames(ct) <- c("count", "parameters")
	# sort OTU counts per table descending 
	ct <-  dplyr::arrange(ct, desc(count))
	return(ct)
}

# create sorted data frames
ets_counts_df <- MakeSortedDF(ets_counts)
coi_counts_df <- MakeSortedDF(coi_counts)

# function: create plot counts versus parameters
PlotCountsParams <- function(cdf, data_type = c("","18S", "COI")){
	# melt data frame for plotting 
	melted <- melt(cdf, id.vars = "parameters") 
	# do basic barplot
	p <- ggplot(melted, aes(x = reorder(parameters, -value), y = value)) + geom_bar(stat = "identity", fill = "grey")
	last_plot() + theme_minimal()
	last_plot() + theme(axis.text.y = element_text(hjust=0, angle=0), axis.text.x = element_text(hjust=0, angle=90))
	last_plot() + ggtitle(paste("Cumulative effect of filtering parameters on", data_type, "data"))
	last_plot() + xlab("Data sets with variable parameters")
	last_plot() + ylab("Total phylotype count")
	}

# plot both counts
# this generated "160121_Rplot_1.pdf" in upload folder
# (Abundance sorted effect of processing parameters)
plot1 <- PlotCountsParams(ets_counts_df, "18S")
plot2 <- PlotCountsParams(coi_counts_df, "COI")
grid.arrange(plot1, plot2, ncol=2)


########### 3: get tables with highest numbers of Phylotype counts into new list 

# finished data - lists of PhyloSeq objects:
psob_list_18S
psob_list_COI


########### 4: Get replicate similarities during filtering

GetReplicateSimilarities <- function(list_item){
	
	# for testing
	# cntrls <- psob_list_18S_sel[[1]]
	
	# get item from list
	cntrls <- list_item
	
	# store OTU table in matrix
	cntrls.t <- otu_table(cntrls)
 
	# isolate replicates via dplyr, via temporary data frame
	soil_a <- as.matrix(dplyr::select(data.frame(cntrls.t),contains("10.H")))
	soil_b <- as.matrix(dplyr::select(data.frame(cntrls.t),contains("11.H")))

	# drop empty observations and transpose
	soil_a <- t(soil_a[rowSums(soil_a) > 0,,  drop=FALSE])
	soil_b <- t(soil_b[rowSums(soil_b) > 0,,  drop=FALSE])

	# calculate distance between replicates
	soil_a.sim <- 1-(vegdist(soil_a, method = "jaccard", na.rm = TRUE))
	soil_b.sim <- 1-(vegdist(soil_b, method = "jaccard", na.rm = TRUE))
	mean.sim <- mean(c(soil_a.sim, soil_b.sim))
	sd.sim <- sd(c(soil_a.sim, soil_b.sim))

	
	soil_ab.mean.sim <- c(soil_a.sim, soil_b.sim, mean.sim, sd.sim)
	
	# return stress items for list items
	return(soil_ab.mean.sim)
	}

# apply function
results_a <- lapply(psob_list_18S, GetReplicateSimilarities)
results_b <- lapply(psob_list_COI, GetReplicateSimilarities)

results_a
results_b

# apply name shortening commands
simis_ets <- ShortenNames18S(results_a)
simis_coi <- ShortenNamesCOI(results_b)

########### 5: Combine count and similarity values

# function: combine results
# from counts and similarity comparsion for plotting 
# uses count data and results from similarity comparsion
# generates data frame
CombineResults <- function(res_count,res_sim){
	# copy input objects, in case testing is needed
	counts.df <- res_count
	simis.list <- res_sim
	# covert similarity results to dataframe 
	simis.df <- t(data.frame(simis.list))
	# add column names to dataframe
	colnames(simis.df) <- c("soil_a", "soil_b", "mean", "sd")
	# create combined dataframe
	comb.df <- cbind(counts.df, simis.df)
	# return results
	return(comb.df)
	}

# apply function to both data sets
df.csp.ets <- CombineResults(ets_counts_df,simis_ets)
df.csp.coi <- CombineResults(coi_counts_df,simis_coi)

df.csp.ets
df.csp.coi


########### 6: Plot Phylotype counts and Similarity between replicates

# function: create plot counts versus parameters
PlotCPS <- function(dfcps, string = c("","18S", "COI")){
	#store paramters(for testing)
	cdf <- dfcps
	data_type <- string
	
	# rename one dataframe column (for plotting of standard deviation) - for now removed
	names(cdf)[names(cdf)=="sd"] <- "stdv"
	cdf$stdv <- NULL
	
	# melt data frame for plotting 
	melted <- melt(cdf, id.vars = "parameters") 
	
	# add variable to highlight maxima
	melted <- mutate(melted, max = "notMax")
  	maxs <- melted %>% group_by(variable) %>% summarise(max = max(value))
 	melted$max[which(melted$value %in% maxs$max)] <- "max"
 	
 	# print
 	print(melted)
	
	# store facet labelling alterations in lis 
	f_labs <- list(
	'count'="Count",
	'soil_a'="Soil 1",
	'soil_b'="Soil 2",
	'mean'="Soils 1 & 2")
	
	# Function to rename facet lables
	f_labeller <- function(variable,value){return(f_labs[value])}
	
	# barplot - count data
	p <- ggplot(melted, aes(x = reorder(parameters, -value), y = value, fill = max, group = variable))#, colour=c("red","grey")))
	last_plot() + facet_grid(variable~., scales = "free", labeller=f_labeller)
	last_plot() + geom_bar(data = subset(melted,variable=="count"),stat = "identity")
	last_plot() + geom_bar(data = subset(melted,variable=="soil_a"),stat = "identity")
	last_plot() + geom_bar(data = subset(melted,variable=="soil_b"),stat = "identity")
	last_plot() + geom_bar(data = subset(melted,variable=="mean"),stat = "identity") 
	# + geom_errorbar(aes(ymax = "stdv", ymin= "stdv"))
	last_plot() + theme_bw()
	last_plot() + theme(axis.text.y = element_text(hjust=0, angle=0), axis.text.x = element_text(hjust=0, angle=90))
	last_plot() + ggtitle(paste("Cumulative effect of processing parameters"))
	last_plot() + xlab(paste("Parameter combinations applied to", data_type, "data"))
	last_plot() + ylab("")
	last_plot() + guides(fill=FALSE) + scale_fill_manual(values=c("dimgray","darkgrey"))
	}
	
# this plotting command generated "160121_Rplot_2.pdf" in upload folder
# ("Cumulative effect of processing parameters")
plot1 <- PlotCPS(df.csp.ets, "18S")
plot2 <- PlotCPS(df.csp.coi, "COI")
grid.arrange(plot1, plot2, ncol=2)

# settings with highest mean similarity between the two soils:
# 18S: c99 t95 f0.001
# COI: c97 t80 f0.001

########### 7: Generate comparable objects of insect controls for plotting

# define import pathnames - updated for data upload - see 160121_file_collation_for_upload.txt
ets_path = "/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/1_AU_MOCK_biom_files_18S/150113_18_OTUs_clust97_tassgn99_md_assigned_only_no_contaminants_invertebrates_AUST_0.001_MOCK.biom"
coi_path = "/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/1_AU_MOCK_biom_files_COI/150121_COI_OTUs_clust97_tassgn80_md_assigned_only_invertebrates_AUST_0.001_MOCK.biom"

# Import function to get Phyloseq objects
# needs file list, returns biom objects in list
createPhyloseq <- function(path) {
	
	# import biom file with pathname in list 
	ps_object <- try(import_biom(path))
	# clean out strings in tax table for further handling
	
	# remove OTUs that do not occur in any samples
	ps_object <- prune_taxa(taxa_sums(ps_object) > 0, ps_object)
	
	# remove crap from taxonomy strings
	tax_table(ps_object) <- gsub("__", "", tax_table(ps_object))
	tax_table(ps_object) <- gsub("_", " ", tax_table(ps_object))
	tax_table(ps_object) <- gsub("-", " ", tax_table(ps_object))
	tax_table(ps_object) <- gsub("  ", " ", tax_table(ps_object))
	
	# rename taxonomic ranks according to design strcuture
	colnames(tax_table(ps_object))[1] = "Superphylum"
	colnames(tax_table(ps_object))[2] = "Phylum"
	colnames(tax_table(ps_object))[3] = "Class" 
	colnames(tax_table(ps_object))[4] = "Order"
	colnames(tax_table(ps_object))[5] = "Family"
	colnames(tax_table(ps_object))[6] = "Genus"
	colnames(tax_table(ps_object))[7] = "Species"
	
	# print a txonomy table to see what needs correction
	print(tax_table(ps_object))
	# return biom object
	return(ps_object)
}

# read in biom files and generate a list of Phyloseq objects
ets <- createPhyloseq(ets_path)
coi <- createPhyloseq(coi_path)

## function correct taxonomy strings in 18S and COI
RenameTaxa18S <- function(ets){
	tax_table(ets)
	tax_table(ets)[,1] <- c("Ecdysozoa")
	tax_table(ets)[,2] <- tax_table(ets)[,4] 
	tax_table(ets)[,3] <- tax_table(ets)[,6] 
	tax_table(ets)[2,3] <- c("Eutardigrada")
	tax_table(ets)[,4] <- c("Blattodea", NA, NA, "Lepidoptera", NA, "Hymenoptera" , NA, "Odonata")
	tax_table(ets)[,5] <- c("Blattidae", NA, NA, "Zygaenidae", NA, "Ichneumonidae", NA, "Coenagrionidae")
	tax_table(ets)[,6] <- tax_table(ets)[,7]
	tax_table(ets)[,7] <- tax_table(ets)[,8]
	tax_table(ets)[,7] <- gsub(" red eyed damselfly", "", tax_table(ets)[,7])
	tax_table(ets)[,7] <- gsub(" DJGI 2006", "", tax_table(ets)[,7])
	tax_table(ets) <- tax_table(ets)[,-8]
	print(tax_table(ets))
	return(ets)
	}
RenameTaxaCOI <- function(coi){

	tax_table(coi)	
	tax_table(coi)[,1] <- tax_table(coi)[,3] 	
	tax_table(coi)[,2] <- tax_table(coi)[,4] 	
	tax_table(coi)[,3] <- tax_table(coi)[,6] 	
	tax_table(coi)[,4] <- tax_table(coi)[,10] 	
	tax_table(coi)[,5] <- tax_table(coi)[,11] 	
	tax_table(coi)[,6] <- tax_table(coi)[,16] 	
	tax_table(coi)[,7] <- tax_table(coi)[,18] 	
	tax_table(coi)[,-8] -> tax_table(coi)  
	tax_table(coi)[,-8] -> tax_table(coi)  
	tax_table(coi)[,-8] -> tax_table(coi)  
	tax_table(coi)[,-8] -> tax_table(coi)  
	tax_table(coi)[,-8] -> tax_table(coi)  
	tax_table(coi)[,-8] -> tax_table(coi)  
	tax_table(coi)[,-8] -> tax_table(coi)  
	tax_table(coi)[,7] <- c(NA)
	tax_table(coi)[10,7] <- tax_table(coi)[10,9]
	tax_table(coi)[10,6] <- c("Acyphas")
	tax_table(coi)[11,6] <- c("Cryptinae")
	tax_table(coi)[6,6] <- c("Houghia")
	tax_table(coi)[,-8] -> tax_table(coi)
	tax_table(coi)[,-8] -> tax_table(coi)
	tax_table(coi)[,-8] -> tax_table(coi)
	tax_table(coi)[,-8] -> tax_table(coi)
	tax_table(coi)[,4] <- gsub("Zygoptera", "Odonata", tax_table(coi)[,4])
	tax_table(coi)[,4] <- gsub("Dionycha", "Araneae", tax_table(coi)[,4])
	tax_table(coi)[,5] <- gsub("Apocrita", "Ichneumonidae", tax_table(coi)[,5])
	tax_table(coi)[3,5] <- "Formicidae"
	tax_table(coi)[,5] <- gsub("Glossata", "Erebidae", tax_table(coi)[,5])
	tax_table(coi) <- gsub("NA", NA, tax_table(coi))

	print(tax_table(coi))
	return(coi)
	}

# get clean taxonomy tables
ets <- RenameTaxa18S(ets)
coi <- RenameTaxaCOI(coi)

# 22.01.2016 checking tax tables for errors
tax_table(coi)
# "Ecdysozoa" "Arthropoda" "Insecta"   "Diptera"     "Brachycera"     "Houghia"     NA
# needs to be 
# "Ecdysozoa" "Arthropoda" "Insecta"   "Diptera"     "Tachinidae"     "Houghia"     NA
tax_table(coi)[,5] <- gsub("Brachycera", "Tachinidae", tax_table(coi)[,5])
# "Ecdysozoa" "Arthropoda" "Insecta"   "Hemiptera"   "Euhemiptera"    NA            NA
# needs to be  
# "Ecdysozoa" "Arthropoda" "Insecta"   "Hemiptera"   NA    NA            NA
tax_table(coi)[,5] <- gsub("Euhemiptera", NA, tax_table(coi)[,5])


# function: Create object from Laurences controls:
newPsPb = function(physeq){
	# Create an OTU table with one obsrevation per species, as there are no abundance reported
	otumat = matrix(1,nrow = 14, ncol = 1)
	otumat

	# name matrix elements
	rownames(otumat) <- paste0("morpho", 1:nrow(otumat))
	colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
	otumat

	# generate a taxon matrix - translating only to available genetic da 
	taxmat = matrix(nrow = nrow(otumat), ncol = 7)

	# name taxa according to LCs manuscript
	taxmat[1,] <- c("Ecdysozoa","Arthropoda", "Arachnida", "Araneae", "Sparassidae", NA, NA)
	taxmat[2,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Blattodea", "Blattidae", "Drymaplaneta", "Drymaplaneta communis")
	taxmat[3,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Coleoptera", "Scarabaeidae", NA, NA)
	taxmat[4,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Dermaptera", "Forficulidae", "Forficula", "Forficula auricularia")
	taxmat[5,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Diptera", "Lauxaniidae", NA, NA)
	taxmat[6,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Hemiptera", "Eurybrachyidae", NA, NA)
	taxmat[7,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Hymenoptera", "Formicidae", "Myrmecia", NA)
	taxmat[8,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Hymenoptera", "Ichneumonidae", NA, NA)
	taxmat[9,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Isoptera", "Rhinotermitidae", "Coptotermes", NA)
	taxmat[10,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Neuroptera", "Chrysopidae", NA, NA)
	taxmat[11,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Lepidoptera", "Lymantriidae", NA, NA)
	taxmat[12,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Odonata", "Coenagrionidae", "Ischnura", "Ischnura heterosticta")
	taxmat[13,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Orthoptera", "Acrididae", NA, NA)
	taxmat[14,] <- c("Ecdysozoa","Arthropoda", "Insecta", "Orthoptera", "Tettigoniidae", NA, NA)

	# name matrix elements
	rownames(taxmat) <- rownames(otumat)
	colnames(taxmat) <- c("Superphylum", "Phylum", "Class", "Order", "Family", "Genus", "Species")
	taxmat

	# convert to Phyloseq
	OTU = otu_table(otumat, taxa_are_rows = TRUE)
	TAX = tax_table(taxmat)
	OTU

	physeq = phyloseq(OTU, TAX)
	physeq
	
	#return object
	return(physeq)
}

# create Phyloseq object from Laurences Paper
ins <- newPsPb(ins) 

tax_table(ets)
tax_table(coi)

otu_table(ins)
tax_table(ins)


########### 7a: Plot composition of "Artificial blends"

# copy object for plotting
ins.p <-ins
otu_table(ins.p)
tax_table(ins.p)

# plot out different levels (Family, Genus, Species)
# this generated "160121_Rplot_3.pdf" in upload folder
# ("Reference data composition of Australian blend")
plot_bar(ins.p, "Species", fill = "Genus", facet_grid="Family~.") + 
theme_bw()  + 
theme(axis.title.y = element_blank(), axis.text.y  = element_blank(), axis.ticks.y  = element_blank()) + 
theme(strip.text.y = element_text(size=8, angle=0))


########### 8: Compare similarities between insect reference and insect controls

# get info from objects for function writing
ets
coi
ins
sample_names(ets)
sample_names(coi)
sample_names(ins)


# isolate mixed sample from both samples
ets.mix <- prune_samples("18S.2.8.H.Inscntrl",ets)
coi.mix <- prune_samples("COI.2.8.H.Inscntrl",coi)

# remove 0 count OTUs
ets.mix <- prune_taxa(taxa_sums(ets.mix) > 0, ets.mix)
coi.mix <- prune_taxa(taxa_sums(coi.mix) > 0, coi.mix)

# rename mixed samples and reference data
sample_names(ets.mix) <- "18S Cntrl"
sample_names(coi.mix) <- "COI Cntrl"
sample_names(ins) <- "Ref"

# finished Phyloseq objects for testing
ets.mix
coi.mix
ins

sample_names(ets.mix)
sample_names(coi.mix)
sample_names(ins)

tax_table(ins)
tax_table(coi.mix)
tax_table(ets.mix)

# function: Compare reference data with otu data
GetMatchPerLevel <- function(refr, quer, vec_ind ){

	# requires:
	# - reference data as phyloseq object - subset to correct sample
	# - query data as phyloseq object - subset to correct sample
	# - vector index for taxonomic level
	# returns:
	# - fraction of TRUE strings contained in query when compared to reference at specidic taxonomic level


	# get input phyloseq objects for function, reference and query data
	# - these will be function arguments
	refr <- refr
	quer <- quer
	r <- vec_ind

	# get names from refernce OTU data, for level x
	#r.str.lev <- na.omit(get_taxa_unique(refr, rank_names(refr)[r]))
	#q.str.lev <- na.omit(get_taxa_unique(quer, rank_names(quer)[r]))
	
	# get names from refernce OTU data, for level x
	r.str.lev <- get_taxa_unique(refr, rank_names(refr)[r])
	q.str.lev <- get_taxa_unique(quer, rank_names(quer)[r])
	
	
	print(paste("Evaluating taxonomic level: ", rank_names(refr)[r]))
	print(paste("Referernce:", sample_names(refr)))
	print(paste("Unique taxa in reference: ", r.str.lev))
	print(paste("Query:", sample_names(quer)))
	print(paste("Unique taxa in query: ", q.str.lev))

	# refernce detected in query? - store value
	ref.was.detected.vec <- r.str.lev %in% q.str.lev
	
	print(paste("Query string matches to reference string: ", ref.was.detected.vec))

	# calculate fraction of detected entities
	ref.was.detected.frc <- sum(ref.was.detected.vec)/ length(ref.was.detected.vec) # length(r.str.lev)
	print(paste("Fraction of detected entities: ", ref.was.detected.frc))

	# return similarity values
	return(ref.was.detected.frc)

}

# GetMatchPerLevel: function testing
GetMatchPerLevel(ins,ets.mix,1)
GetMatchPerLevel(ins,ets.mix,2)
GetMatchPerLevel(ins,ets.mix,3)
GetMatchPerLevel(ins,ets.mix,4)
GetMatchPerLevel(ins,ets.mix,5)
GetMatchPerLevel(ins,ets.mix,6)
GetMatchPerLevel(ins,ets.mix,7)

# function: GetMatchMultipleLevels
GetMatchMultipleLevels <- function(refr, quer){

	# requires:
	# - reference data as phyloseq object - subset to correct sample
	# - query data as phyloseq object - subset to correct sample
	# - function GetMatchPerLevel()
	# returns
	# - vector with concordance as calculated by GetMatchPerLevel()

	# get input phyloseq objects for function, reference and query data
	# - these will be function arguments
	refr <- refr
	quer <- quer

	# loop through rank names via vector index
	foreach(i=1:(length(rank_names(refr))), .combine= 'c') %do% {
		# return vector index
		GetMatchPerLevel(refr,quer,(i))
		}	
	}

# Compare and store in Matrix:
cntrl.vs.mg <- rbind(
	GetMatchMultipleLevels(ins, ets.mix),
	# GetMatchMultipleLevels(ets.mix, ins),
	GetMatchMultipleLevels(ins, coi.mix)) # ,
	# GetMatchMultipleLevels(coi.mix, ins),
	# GetMatchMultipleLevels(coi.mix, ets.mix),
	# GetMatchMultipleLevels(ets.mix, coi.mix))

# Name columns and rows
colnames(cntrl.vs.mg) <- rank_names(ins)
# rownames(cntrl.vs.mg) <-c("18S to Cntrl","Cntrl to 18S","COI to Cntrl","Cntrl to COI", "18S to COI", "COI to 18S")
rownames(cntrl.vs.mg) <-c("18S to Cntrl", "COI to Cntrl") # , "18S to COI", "COI to 18S")

cntrl.vs.mg


########### 9: Create display Items: Insect vs Controls

## heatmap - unused in mansucript
plot.new()
heatmap.2(cntrl.vs.mg,
	Rowv = FALSE,
	Colv = FALSE, 
	dendrogram = "none",
	trace = "row",
	tracecol = "blue",
	# scale = "row",
	na.rm=TRUE,
	# sepcolor="black",
	sepwidth=c(0.0001,0.0001),
	colsep=1:ncol(cntrl.vs.mg),
        rowsep=1:nrow(cntrl.vs.mg),
        cexRow = 1.2,
	cexCol = 1.4,
	density.info=c("none"),
	lmat = rbind(c(1,2),c(4,3)),
	lwid = c(4,0.5),
	lhei = c(2,.5),
	keysize = 1.6,
	margins = c(5,1))

## barplot

# create data frame
cntrl.vs.mg.df <- data.frame(cntrl.vs.mg)
cntrl.vs.mg.df$Comparison <- factor(c("18S to Cntrl.", "COI to Cntrl."))
cntrl.vs.mg.df

# melt data frame
melted <- melt(cntrl.vs.mg.df, id.vars=c("Comparison"))
str(melted)

# Define palette
cbPalette <- c( "#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC", "#E69F00")

# create plot - this generated - "160121_Rplot_4.pdf" in upload folder
# (Australien blend - reference data generated from paper vs sequenced extract)
plot3 <-  ggplot(melted, aes(variable, value)) + geom_bar(aes(fill = Comparison), position = "dodge", stat="identity")
plot3 <- last_plot() + theme_bw()
plot3 <- last_plot() + xlab("") + ylab("Detected fraction of taxa contained in reference")
plot3 <- last_plot() + theme(axis.text=element_text(size=13))
plot3 <- last_plot() + theme(axis.title=element_text(size=12))
plot3 <- last_plot() + scale_fill_manual(values=cbPalette) +   scale_colour_manual(values=cbPalette)
plot3


########### 10: Create display Items: Barplot as comparsion between samples

# function: do fecetted barplot
# works only with formated input Phyloseq objects
GetBarplotsInsCntrl <- function(ins,ets.mix,coi.mix){

	# Define palette
	cbPalette <- c( "#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC", "#E69F00")

	plot1 <- plot_bar(ins, fill="Family")
	plot1 <- last_plot() + theme_bw()
	plot1 <- last_plot() + facet_grid(Order~., scales = "free")
	plot1 <- last_plot() + xlab("12 Orders, 14 Families")
	plot1 <- last_plot() + ylab("Count")
	plot1 <- last_plot() + theme(axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank())
	plot1

	plot2 <- plot_bar(ets.mix, fill="Family")
	plot2 <- last_plot() + theme_bw()
	plot2 <- last_plot() + facet_grid(Order~., scales = "free")
	plot2 <- last_plot() + xlab("18S: 4 Orders, 4 Families")
	plot2 <- last_plot() + ylab("")
	plot2 <- last_plot() + theme(axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank())
	plot2

	plot3 <- plot_bar(coi.mix, fill="Family")
	plot3 <- last_plot() + theme_bw()
	plot3 <- last_plot() + facet_grid(Order~., scales = "free")
	plot3 <- last_plot() + xlab("COI: 6 Orders, 8 Families")
	plot3 <- last_plot() + ylab("")
	plot3 <- last_plot() + theme(axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank())
	plot3
	
	grid.arrange(plot1, plot2, plot3, ncol=3)
	}

# create display Items: Composition of 18S and COI Insect controls vs Reference data
# this generated "160121_Rplot_5.pdf"

# (comparsion of reference data to sequenced data in Austrlian soils)
GetBarplotsInsCntrl(ins,ets.mix,coi.mix)

############ 11: Import Antarctic morphology data and metagenetic samples

# define import pathnames - updated for upload folder
ets_path = "/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/2_ANT_SOIL_biom_18S/150113_18_OTUs_clust97_tassgn99_md_assigned_only_no_contaminants_invertebrates_ANT_0.001.biom"
coi_path = "/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/2_ANT_SOIL_biom_COI/150121_COI_OTUs_clust97_tassgn80_md_assigned_only_invertebrates_ANT_0.001.biom"

# read in biom files and generate a list of Phyloseq objects
ets.ant <- createPhyloseq(ets_path)
coi.ant <- createPhyloseq(coi_path)

# insect mix
ets.mix
coi.mix

# soil samples
ets.soil <- psob_list_18S[[17]]
coi.soil <- psob_list_COI[[8]]

# function: define new operator "not in"
'%!in%' <- function(x,y)!('%in%'(x,y))

# function: filter OTUs from controls from Antarctic samples
# input: Phyloseq objcet to be filtered followed by the 2 ps object that need to be exluded
# returns: Filtered Phyloseq object:
KeepAntPS <- function(ant_ps, mix_ps, soil_ps){

	# copy object 
	filtered_ps <- ant_ps
	
	# filter samples contained in insect mix
	otu_table(filtered_ps) <- subset(otu_table(ant_ps), subset = rownames(otu_table(ant_ps)) %!in% rownames(otu_table(mix_ps)), drop = FALSE)

	# filter samples contained in soil amplicons
	otu_table(filtered_ps) <- subset(otu_table(filtered_ps), subset = rownames(otu_table(filtered_ps)) %!in% rownames(otu_table(soil_ps)), drop = FALSE)

	# return filtered tables
	return(filtered_ps)

	}

# Filter Antarctic metagenetic samples to contain only Phylotypes not contained in any of the controls
ets.ant.fil <- KeepAntPS(ets.ant, ets.mix, ets.soil)
coi.ant.fil <- KeepAntPS(coi.ant, coi.mix, coi.soil)

ets.ant.fil
coi.ant.fil

# functions: correct taxonomy strings in 18S and COI
RenameAntTaxa18S <- function(ets.f){

	# copy for function testing 
	# ets.f <- ets.ant.fil 
	
	# rename items
	tax_table(ets.f)
	tax_table(ets.f)[1,] <- c("Ecdysozoa", "Arthropoda", "Arachnida", "Oribatida", "Phenopelopidae", "Eupelops", "Eupelops plicatus", NA)
	tax_table(ets.f)[2,] <- c("Lophotrochozoa", "Rotifera", NA, NA, NA, NA, NA, NA)
	tax_table(ets.f)[3,] <- c("Ecdysozoa", "Nematoda", "Chromadorea", "Araeolaimida", "Plectidae", NA, NA, NA)
	tax_table(ets.f)[4,] <- c("Ecdysozoa", "Nematoda", "Chromadorea", "Monhysterida", "Monhysteridae", "Halomonhystera", "Halomonhystera disjuncta", NA)
	tax_table(ets.f)[5,] <- c("Lophotrochozoa", "Rotifera", NA, NA, NA, NA, NA, NA)
	tax_table(ets.f) <- tax_table(ets.f)[,-8]

	# print and return items
	print(tax_table(ets.f))
	return(ets.f)
	
	}
RenameAntTaxaCOI <- function(coi.f){

	# rename items
	tax_table(coi.f)
	
	tax_table(coi.f)[,1] <- tax_table(coi.f)[,3]
	tax_table(coi.f)[,2] <- tax_table(coi.f)[,4]
	tax_table(coi.f)[,3] <- tax_table(coi.f)[,5]
	tax_table(coi.f)[,3] <- gsub("Hexapoda", "Insecta", tax_table(coi.f)[,3])
	tax_table(coi.f)[,4] <- tax_table(coi.f)[,10]
	tax_table(coi.f)[,5] <- tax_table(coi.f)[,13]
	tax_table(coi.f)[,6] <- tax_table(coi.f)[,16]
	tax_table(coi.f)[,7] <- tax_table(coi.f)[,17] 
	
	tax_table(coi.f)[1,4:7] <- c("Adinetida",  "Adinetidae", "Adinata", NA)
	tax_table(coi.f)[2,4:7] <- c("Araeolaimida", "Plectidae", "Plectus", "Plectus murrayi")
	tax_table(coi.f)[4,5:7] <- c("Geometridae","Pleurolopha", NA)
	tax_table(coi.f)[6,7] <- c("Pterostichus cristatus")
	tax_table(coi.f)[9,5:7] <- c("Muscidae", "Hydrotaea", "Hydrotaea aenescens")
	tax_table(coi.f)[10,5:7] <- c("Carabidae", "Pterostichus", "Pterostichus cristatus")
	tax_table(coi.f)[11,5:7] <- c("Ulidiidae", "Timia", "Timia nigripes")
	tax_table(coi.f)[12, 5:7] <- c("Sphingidae", "Cerberonoton", "Cerberonoton rubescens")
	tax_table(coi.f)[13, 5] <- c("Geometridae")

	tax_table(coi.f) <- tax_table(coi.f)[,-c(8:18)]

	# print and return items
	print(tax_table(coi.f))
	return(coi.f)
	
	}

# correct 18S and COI tax strings
ets.ant.fil.rn <- RenameAntTaxa18S(ets.ant.fil)
coi.ant.fil.rn <- RenameAntTaxaCOI(coi.ant.fil)

# function: read in morphotype data and correct taxon and sample names 
# works only for spcied object!
# return cleaned phyloseq object

GetAntMorph <- function(morph_path){

	# import phyloseq object 
	path <- morph_path
	morph.ant <- import_biom(path)

	# create 7th column
	tax_table(morph.ant) <- cbind(tax_table(morph.ant), NA) 

	# rename ranks
	colnames(tax_table(morph.ant)) <- c("Superphylum", "Phylum", "Class", "Order", "Family", "Genus", "Species" )
	rank_names(morph.ant)

	#correct OTU table
	tax_table(morph.ant)
	tax_table(morph.ant)[4,4:7] <- c("Parachela","Hypsibiidae","Acutuncus","Acutuncus antarcticus")
	tax_table(morph.ant)[5,4:7] <- c("Parachela","Macrobiotidae","Macrobiotus",NA)
	
	# 22.01.2056 correction of a mistake in the morphology table
	#   morph2  "Ecdysozoa"      "Nematoda"   "Nematoda"     "Enoplea"      "Qudsianematidae" "Eudorylaimus" NA                     
	#   morph3  "Ecdysozoa"      "Nematoda"   "Nematoda"     "Chromadorea"  "Rhabditidae"     "Scottnema"    NA  
	#                                 Class        Order          Family             Genus           Species
	tax_table(morph.ant)[7,3:7] <- c("Enoplea",    "Dorylaimida", "Qudsianematidae", "Eudorylaimus", NA)
	tax_table(morph.ant)[8,3:7] <- c("Chromadorea","Rhabditidae", "Cephalobidae",     "Scottnema",   NA)
	
	# remove OTUs that do not occur in any samples
	morph.ant <- prune_taxa(taxa_sums(morph.ant) > 0, morph.ant)
	
	# print object info, corrected
	print(otu_table(morph.ant))
	print(tax_table(morph.ant))

	# return corrected phyloseq object
	return(morph.ant)

	}

# import AVCs morpho dat and clean up  
morph_path <- "/Users/paul/Documents/14_Ph_D_thesis_C3_Invertebrates_Markers/160120_submission/zenodo/3_ANT_MORPH/150307_AVC_OTU_table.biom"
morph.ant.rn <- GetAntMorph(morph_path) 

# rename sample in all sets
sample_names(ets.ant.fil.rn) <- c("LH-1", "MS-1", "VH-2", "mix", "HI-1", "CS-2", "VH-1", "LH-2", "CS-1")
sample_names(coi.ant.fil.rn) <- c("HI-1", "CS-1", "LH-2", "mix", "CS-2")
sample_names(morph.ant.rn) <- c("CS-1", "CS-2", "VH-1", "LH-1", "LH-2", "HI-1", "VH-2")

sample_names(ets.ant.fil.rn)
sample_names(coi.ant.fil.rn)
sample_names(morph.ant.rn) 

tax_table(morph.ant.rn)

############ 11a: Plot composition of "Antarctic soil" reference data

# copy object for plotting
morph.ant.rn -> p.ant

otu_table(p.ant)
tax_table(p.ant)
sample_names(p.ant)

# plot out different levels (Family, Genus, Species)
# this generated "160122_Rplot_6.pdf" in upload folder - updated after bug chasing 22.01.16
plot_bar(p.ant, "Order", fill = "Class", facet_grid="Phylum~Sample") + 
theme_bw()  + 
theme(axis.title.y = element_blank(), axis.text.y  = element_blank(), axis.ticks.y  = element_blank()) + 
theme(strip.text.y = element_text(size=8, angle=0)) +
theme(axis.text.x  = element_text(angle = 90, hjust = 1))


############ 12: compare assignemnts

# function: subset matching samples 18S or COI vs morphtype data
# do pairwise comparsion
# requires
# - three phyloseq objects with matching sample names:
# -- morphologic reference
# -- metagenetic data 1 (18S)
# -- metagenetic data 2 (COI) 
# returns
# - similarity matrix
CompareAntSamples <- function(ref,quer_a, quer_b){ 

	# copy variables so that variable in function do not need renaming after design:
	morph.ant.rn <- ref
	ets.ant.fil.rn <- quer_a
	coi.ant.fil.rn <- quer_b

	# get samples that are available in both sets
	ets.morph.sample.names <- sort(sample_names(morph.ant.rn)[which(sample_names(morph.ant.rn) %in% sample_names(ets.ant.fil.rn))]) 
	coi.morph.sample.names <- sort(sample_names(morph.ant.rn)[which(sample_names(morph.ant.rn) %in% sample_names(coi.ant.fil.rn))])

	# Compare 18S samples
	# create results list
	ResultListEts <- list()
	# loop over vector with sample names and do comparison
	for(i in ets.morph.sample.names){

		# print(paste("Using sample: ",i))
	
		# create new phyloseq object 1  - pruned to sample specified at vector position
		ets.ant.fil.rn.iso <- prune_samples(i,ets.ant.fil.rn)
		# create new phyloseq object 2  - pruned to sample specified at vector position
		morph.ant.rn.iso <-  prune_samples(i,morph.ant.rn)
	
		# print(paste("Subset to sample: ", sample_names(ets.ant.fil.rn.iso)))
		# print(paste("Subset to sample: ", sample_names(morph.ant.rn.iso)))
	
		# remove OTUs that do not occur in any samples
		ets.ant.fil.rn.iso <- prune_taxa(taxa_sums(ets.ant.fil.rn.iso) > 0, ets.ant.fil.rn.iso)
		morph.ant.rn.iso <- prune_taxa(taxa_sums(morph.ant.rn.iso) > 0, morph.ant.rn.iso)
	
		# print(otu_table(ets.ant.fil.rn.iso))
		# print(otu_table(morph.ant.rn.iso))
		
		# Compare composition
		comp_f <- GetMatchMultipleLevels(ets.ant.fil.rn.iso, morph.ant.rn.iso)
		# comp_r <- GetMatchMultipleLevels(morph.ant.rn.iso, ets.ant.fil.rn.iso)
	
		# create sub-result matrix
		# s_res <- rbind(comp_f, comp_r)
		# rownames(s_res) <- c(paste("18S, MvsG,", i),paste("18S, GvsM,",i)) 
		# colnames(s_res) <- c("Superphylum","Phylum", "Class", "Order", "Family", "Genus", "Species")
		
		s_res <- rbind(comp_f)
		rownames(s_res) <- c(paste(i, "18S")) 
		colnames(s_res) <- c("Superphylum","Phylum", "Class", "Order", "Family", "Genus", "Species")
		
	
		# store results in list 
		ResultListEts[[i]] <- s_res
  
		}
	# combine 18S results
	ResultMatrixEts <-  do.call(rbind, ResultListEts)
	ResultMatrixEts

	# Compare COI samples
	# create results list
	ResultListCoi <- list()
	# loop over vector with sample names
	for(i in coi.morph.sample.names){

		# print(paste("Using sample: ",i))
	
		# create new phyloseq object 1  - pruned to sample specified at vector position
		coi.ant.fil.rn.iso <- prune_samples(i,coi.ant.fil.rn)
		# create new phyloseq object 2  - pruned to sample specified at vector position
		morph.ant.rn.iso <-  prune_samples(i,morph.ant.rn)
	
		# print(paste("Subset to sample: ", sample_names(coi.ant.fil.rn.iso)))
		# print(paste("Subset to sample: ", sample_names(morph.ant.rn.iso)))
	
		# remove OTUs that do not occur in any samples
		coi.ant.fil.rn.iso <- prune_taxa(taxa_sums(coi.ant.fil.rn.iso) > 0, coi.ant.fil.rn.iso)
		morph.ant.rn.iso <- prune_taxa(taxa_sums(morph.ant.rn.iso) > 0, morph.ant.rn.iso)
	
		# print(otu_table(coi.ant.fil.rn.iso))
		# print(otu_table(morph.ant.rn.iso))
		
		# Compare composition
		comp_f <- GetMatchMultipleLevels(coi.ant.fil.rn.iso, morph.ant.rn.iso)
		# comp_r <- GetMatchMultipleLevels(morph.ant.rn.iso, coi.ant.fil.rn.iso)
	
		# create sub-result matrix
		# s_res <- rbind(comp_f, comp_r)
		# rownames(s_res) <- c(paste("COI, MvsG,", i),paste("COI, GvsM,",i)) 
		# colnames(s_res) <- c("Superphylum","Phylum", "Class", "Order", "Family", "Genus", "Species")
		
		# create sub-result matrix
		s_res <- rbind(comp_f)
		rownames(s_res) <- c(paste(i, "COI")) 
		colnames(s_res) <- c("Superphylum","Phylum", "Class", "Order", "Family", "Genus", "Species")
	
	
		# store results in list 
		ResultListCoi[[i]] <- s_res
  
		}

	# combine Coi results
	ResultMatrixCoi <-  do.call(rbind, ResultListCoi)
	ResultMatrixCoi

	# combine both results
	ResultsComp <- rbind(ResultMatrixEts, ResultMatrixCoi)
	ResultsComp

	return(ResultsComp)
	}

# calculate similarities at taxonomic level
res.mat <- CompareAntSamples(morph.ant.rn,ets.ant.fil.rn, coi.ant.fil.rn)
res.mat

# sort samples based on location
res.mat <- res.mat[order(substring( rownames(res.mat) ,1,4)),]

##  display items: generate heatmap
# unused in mansucript
plot.new()
heatmap.2(res.mat,
	Rowv = FALSE,
	Colv = FALSE, 
	dendrogram = "none",
	trace = "row",
	tracecol = "blue", 
	#sepcolor="black",
	sepwidth=c(0.0001,0.0001),
	colsep=1:ncol(res.mat),
        rowsep=1:nrow(res.mat),
        cexRow = 1.4,
	cexCol = 1.4,
	density.info=c("none"),
	lmat = rbind(c(1,2),c(4,3)),
	lwid = c(4,0.5),
	lhei = c(2,.5),
	keysize = 1.6,
	margins = c(5,1),
	na.color = "grey")

##  display items: barplot

# store matrix in data frame
res.mat.df <- data.frame(res.mat)

# name columns
res.mat.df$Marker <- as.factor(substring(rownames(res.mat),6,8)) 
res.mat.df$Location <- as.factor(substring(rownames(res.mat),1,4)) 

# filter out locations for which data is not available for both locations
res.mat.df <-  dplyr::filter(res.mat.df, Location != "LH-1" & Location != "VH-1" & Location != "VH-2" )
str(res.mat.df)

# redefine factor levels
res.mat.df$Location <- factor(res.mat.df$Location)
res.mat.df$Marker <- factor(res.mat.df$Marker)
res.mat.df
str(res.mat.df)

# melt data frame
melted <- melt(res.mat.df, id.vars=c("Location", "Marker"))
melted
str(melted)

# Define palette
cbPalette <- c( "#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC", "#E69F00")

# barplot - count data
plot4 <- ggplot(melted, aes(x = variable, y = value, fill = Marker))

# base plot
plot4 <- last_plot() + facet_grid(Location~., scales = "free")
plot4 <- last_plot() + geom_bar(data = melted ,stat = "identity", position = "dodge")

# make it pretty
# this generated "160122_Rplot_7.pdf"
plot4 <- last_plot() + theme_bw()
plot4 <- last_plot() + theme(axis.text=element_text(size=12))
plot4 <- last_plot() + theme(axis.title=element_text(size=11))
plot4 <- last_plot() + scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette)
plot4 <- last_plot() + xlab("") + ylab("Detected fraction of taxa contained in reference")
plot4

## do comosite figure

# define plots
upper <- plot3
lower <- plot4

# reposition legend
upper <- upper + theme(legend.justification=c(1,1), legend.position=c(1,1))
lower <- lower + theme(legend.position="none")

# do axis titles
upper <- upper + xlab("") + ylab("Matching taxonomic assignments")
upper <- upper + xlab("") + ylab("Matching taxonomic assignments")
theme(axis.title.x=element_text(size=12))
lower <- lower + xlab("") + ylab("")

# re-scale x asis lables
upper <- upper + theme(axis.text.x=element_text(size=13))
lower <- lower + theme(axis.text.x=element_text(size=13))

# re-scale y asis lables
upper <- upper + theme(axis.text.y=element_text(size=13))
lower <- lower + theme(axis.text.y=element_text(size=9))

# finished plot - "160122_Rplot_8.pdf" in upload folder
grid.arrange(upper, lower, nrow=2)

############ 13: plot taxonomic composition

# Barplot function
GetBarplotsAnt <- function(ant,ets,coi){

	print(paste("Omitting: ", sample_names(ets)[4]))
	print(paste("Omitting: ", sample_names(ets)[2]))
	print(paste("Omitting: ", sample_names(coi)[4]))
	
	ets <- prune_samples(sample_names(ets)[-4], ets)
	ets <- prune_samples(sample_names(ets)[-2], ets)
	coi <- prune_samples(sample_names(coi)[-4], coi)

	plot1 <- plot_bar(ant, fill="Family")
	plot1 <- last_plot() + theme_bw()
	plot1 <- last_plot() + facet_grid(Order~., scales = "free")
	plot1 <- last_plot() + xlab("6 Orders, 7 Families")
	plot1 <- last_plot() + ylab("Count")
	plot1 <- last_plot() + theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 6))
	plot1

	plot2 <- plot_bar(ets, fill="Family")
	plot2 <- last_plot() + theme_bw()
	plot2 <- last_plot() + facet_grid(Order~., scales = "free")
	plot2 <- last_plot() + xlab("18S: 3 Orders, 3 Families")
	plot2 <- last_plot() + ylab("")
	plot2 <- last_plot() + theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 6))
	plot2

	plot3 <- plot_bar(coi, fill="Family", facet_grid=~Order)
	plot3 <- last_plot() + theme_bw()
	plot3 <- last_plot() + facet_grid(Order~., scales = "free")
	plot3 <- last_plot() + xlab("COI: 5 Orders, 7 Families")
	plot3 <- last_plot() + ylab("")
	plot3 <- last_plot() + theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 8))
	plot3
	
	grid.arrange(plot1, plot2, plot3, ncol=3)
	}

# do Barplot - this generates "160122_Rplot_9.pdf"
GetBarplotsAnt(morph.ant.rn, ets.ant.fil.rn, coi.ant.fil.rn)


############ 14: Format Antarctic data into a more useful data frame

# store matrix in data frame
res.mat.df <- data.frame(res.mat)

# name columns
res.mat.df$Gene <- as.factor(substring(rownames(res.mat),6,8)) 
res.mat.df$Location <- as.factor(substring(rownames(res.mat),1,4)) 

# filter out locations for which data is not available for both locations
res.mat.df <-  dplyr::filter(res.mat.df, Location != "LH-1" & Location != "VH-1" & Location != "VH-2" )
# redefine factor levels
res.mat.df$Location <- factor(res.mat.df$Location)
res.mat.df

############ 15: Calculate Intraclass correlation coefficient on Control and Antarctic data

# function: Calculate Intraclass correlation coefficient (oneway, consistency)
GetICC <- function(two_df_rows){

	test.data <- t(as.matrix(two_df_rows))
	print(test.data[1,1:2])
	print(test.data[2,1:2])
	print(test.data[3,1:2])
	print(test.data[4,1:2])
	print(test.data[5,1:2])
	print(test.data[6,1:2])
	print(test.data[7,1:2])

	icc.result <- icc(test.data, model="oneway", type="consistency")
	print(icc.result)

	# print(icc.result)
	# print(paste("value: ", icc.a$value))
	# print(paste("p: ", icc.a$p.value))

	# print(paste("Genes: ",icc.a$raters))
	# print(paste("Ranks: ",icc.a$subjects))
	
	}


## Insect controls
GetICC(cntrl.vs.mg)

## Antarctic CS - 1
GetICC(res.mat.df[1:2,1:7])

## Antarctic CS - 2
GetICC(res.mat.df[3:4,1:7])

## Antarctic HI - 1
GetICC(res.mat.df[5:6,1:7])

## Antarctic LH - 2
GetICC(res.mat.df[7:8,1:7])
