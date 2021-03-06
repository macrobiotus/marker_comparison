# Ground-truthing Phylotype Assignments for Antarctic Invertebrates

## Introduction
The following code represents analyses conducted for the manuscript titled
**Ground-truthing Phylotype Assignments for Antarctic Invertebrates, DNA Barcodes 2017, 4: 26–38, DOI 10.1515/dna-2017-0001**.
Use `git clone` or your favourite Git manager to download these files to your system. **Please do not
execute the scripts on your system without looking at them first**. Please contact me should you have questions, I will attempt to help you if time permits.

## Disclaimer
**THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.**

## System requirements
1. R environment (currently using R 3.2.0 - R packages employed for analysis,
are listed in file `./160621_package_citations.bib`. This file can be opened
with most common reference managers)
2. Bash shell (OSX 10.10.5) - Should you work in a different environment consider
adopting the line endings of files as well as syntax differences in tools such as
`sed` or `grep`.
4. MacQiime (Qiime 1.8., later versions of Qiime may have slight syntax differences).


## Contents

### Directories
- `0_AU_SOIL_biom_files_18S/`: `.biom` files for 18S data from Australian soils
generated with variable thresholds as described in the manuscript.
- `0_AU_SOIL_biom_files_COI/`: `.biom` files for COI data from Australian soils
generated with variable thresholds as described in the manuscript.
- `1_AU_MOCK_biom_files_18S/`: `.biom` files for 18S data from Insect blends
generated with variable thresholds as described in the manuscript.
- `1_AU_MOCK_biom_files_COI/`: `.biom` files for COI data from Insect blends
generated with variable thresholds as described in the manuscript.
- `2_ANT_SOIL_biom_18S/`: `.biom` files for 18S data from Antarctic soils
generated with variable thresholds as described in the manuscript.
- `2_ANT_SOIL_biom_COI/`: `.biom` files for COI data from Antarctic soils
generated with variable thresholds as described in the manuscript.
- `3_ANT_MORPH/`: `.biom` files containing information from morphological assessments of
Antarctic soil samples.
- `4_rep_sets/`: compressed `.fasta` files containing relevant sequences for 18S and COI data sets
(as indicated in file `./160121_file_collation_for_upload.txt`).
- `5_mapping_files/`: mapping files for raw data generation using Qiime. **Note:** Samples
are filtered in shell scripts (`.sh`) and the main analysis file (`160127_data_analysis_3.R`).

### Files
- `./160121_comparison_values.txt`:
- `./160121_file_collation_for_upload.txt`: extract of output from `./160127_data_analysis_3.Rdata`; for writing.
- `./160121_Rplot_1.pdf`: Raw R plot, used in manuscript or supplemental file, generated by`./160127_data_analysis_3.R`.
- `./160121_Rplot_2.pdf`: Raw R plot, used in manuscript or supplemental file, generated by`./160127_data_analysis_3.R`.
- `./160121_Rplot_3.pdf`: Raw R plot, used in manuscript or supplemental file, generated by`./160127_data_analysis_3.R`.
- `./160121_Rplot_4.pdf`: Raw R plot, used in manuscript or supplemental file, generated by`./160127_data_analysis_3.R`.
- `./160122_Rplot_5.pdf`: Raw R plot, used in manuscript or supplemental file, generated by`./160127_data_analysis_3.R`.
- `./160122_Rplot_6.pdf`: Raw R plot, used in manuscript or supplemental file, generated by`./160127_data_analysis_3.R`.
- `./160122_Rplot_7.pdf`: Raw R plot, used in manuscript or supplemental file, generated by`./160127_data_analysis_3.R`.
- `./160122_Rplot_8.pdf`: Raw R plot, used in manuscript or supplemental file, generated by`./160127_data_analysis_3.R`.
- `./160122_Rplot_9.pdf`: Raw R plot, used in manuscript or supplemental file, generated by`./160127_data_analysis_3.R`.
- `./160127_data_analysis_1a_18S.sh`: Data preparation using `bash` and `macqiime` (18S). **Do not execute this script without looking at it.**
- `./160127_data_analysis_1b_COI.sh`: Data preparation using `bash` and `macqiime` (COI). **Do not execute this script without looking at it.**
- `./160127_data_analysis_2.sh`: Data preparation using `bash` and `macqiime` (18S and COI). **Do not execute this script without looking at it.**
- `./160127_data_analysis_3.R`: Main analysis. Uses `.biom` files in this directory structure. **Do not execute this script without looking at it.**
- `./160127_data_analysis_3.Rdata`: (May not be present). Copy of R environment after successful run of R analysis.
- `./160621_package_citations.bib`: Citations for R packages, generated by`./160127_data_analysis_3.R`.
- `./README.md`: This file.
