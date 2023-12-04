library(dplyr)
library(tibble)
library(tidyverse)
library(ggrepel) #for genes points on plot
library(gridExtra)


# script to produce plot of (zonal) positions of teto-Olfr17 expressing cells along a curve of the MOE (Figure 7E) 
# the curve is drawn to start medially along the septum (in zone1) and move distally (towards zone 5)

# input is TSV file generated with get_zonal_position_along_curve.py script with normalized distance along the curve for each cell


## read in tetoP2*/+ (rep1) and tetoP2*/*(rep2) imaging data. 15um sections. ~8 week old mouse.
# total of 29 sections, 11 from replicate1, and 18 from replicate2
tetoP2_1_1a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto006.tsv") %>%
  mutate(replicate=1, section="1a", condition="no_tTA")
tetoP2_1_1b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto006_2.tsv") %>%
  mutate(replicate=1, section="1b", condition="no_tTA")
tetoP2_1_2a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto010.tsv") %>%
  mutate(replicate=1, section="2a", condition="no_tTA")
tetoP2_1_2b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto010_2.tsv") %>%
  mutate(replicate=1, section="2b", condition="no_tTA")
tetoP2_1_3_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto014.tsv") %>%
  mutate(replicate=1, section="3", condition="no_tTA")
tetoP2_1_4_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto015.tsv") %>%
  mutate(replicate=1, section="4", condition="no_tTA")
tetoP2_1_5_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto019.tsv") %>%
  mutate(replicate=1, section="5", condition="no_tTA")
tetoP2_1_6_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto021.tsv") %>%
  mutate(replicate=1, section="6", condition="no_tTA")
tetoP2_1_7_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto022.tsv") %>%
  mutate(replicate=1, section="7", condition="no_tTA")
tetoP2_1_8_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto023.tsv") %>%
  mutate(replicate=1, section="8", condition="no_tTA")
tetoP2_1_9_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_10x_notAuto024.tsv") %>%
  mutate(replicate=1, section="9", condition="no_tTA")
tetoP2_2_1a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus003.tsv") %>%
  mutate(replicate=2, section="1a", condition="no_tTA")
tetoP2_2_1b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus003_2.tsv") %>%
  mutate(replicate=2, section="1b", condition="no_tTA")
tetoP2_2_2_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus004.tsv") %>%
  mutate(replicate=2, section="2", condition="no_tTA")
tetoP2_2_3_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus006.tsv") %>%
  mutate(replicate=2, section="3", condition="no_tTA")
tetoP2_2_4a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus014.tsv") %>%
  mutate(replicate=2, section="4a", condition="no_tTA")
tetoP2_2_4b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus014_2.tsv") %>%
  mutate(replicate=2, section="4b", condition="no_tTA")
tetoP2_2_5a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus015.tsv") %>%
  mutate(replicate=2, section="5a", condition="no_tTA")
tetoP2_2_5b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus015_2.tsv") %>%
  mutate(replicate=2, section="5b", condition="no_tTA")
tetoP2_2_6a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus017.tsv") %>%
  mutate(replicate=2, section="6a", condition="no_tTA")
tetoP2_2_6b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus017_2.tsv") %>%
  mutate(replicate=2, section="6b", condition="no_tTA")
tetoP2_2_7a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus020.tsv") %>%
  mutate(replicate=2, section="7a", condition="no_tTA")
tetoP2_2_7b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus020_2.tsv") %>%
  mutate(replicate=2, section="7b", condition="no_tTA")
tetoP2_2_8a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus022.tsv") %>%
  mutate(replicate=2, section="8a", condition="no_tTA")
tetoP2_2_8b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus022_2.tsv") %>%
  mutate(replicate=2, section="8b", condition="no_tTA")
tetoP2_2_9a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus023.tsv") %>%
  mutate(replicate=2, section="9a", condition="no_tTA")
tetoP2_2_9b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus023_2.tsv") %>%
  mutate(replicate=2, section="9b", condition="no_tTA")
tetoP2_2_10a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus025.tsv") %>%
  mutate(replicate=2, section="10a", condition="no_tTA")
tetoP2_2_10b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2/P2_homo_r4_10x_notAutofocus025_2.tsv") %>%
  mutate(replicate=2, section="10b", condition="no_tTA")

## read in tetoP2*/+; gg8tTA*/+ 5ug Dox imaging data (new OE embeded at ~8 weeks old). 14um sections
# 16 sections replicate1, 19 sections replicate2
lowDOX_5ug_1_1a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto005.tsv") %>%
  mutate(replicate=1, section="1a", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_1b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto005_2.tsv") %>%
  mutate(replicate=1, section="1b", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_2a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto006.tsv") %>%
  mutate(replicate=1, section="2a", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_2b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto006_2.tsv") %>%
  mutate(replicate=1, section="2b", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_3a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto009.tsv") %>%
  mutate(replicate=1, section="3a", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_3b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto009_2.tsv") %>%
  mutate(replicate=1, section="3b", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_4a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto010.tsv") %>%
  mutate(replicate=1, section="4a", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_4b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto010_2.tsv") %>%
  mutate(replicate=1, section="4b", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_5a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto011.tsv") %>%
  mutate(replicate=1, section="5a", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_5b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto011_2.tsv") %>%
  mutate(replicate=1, section="5b", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_6a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto012.tsv") %>%
  mutate(replicate=1, section="6a", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_6b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto012_2.tsv") %>%
  mutate(replicate=1, section="6b", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_7a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto013.tsv") %>%
  mutate(replicate=1, section="7a", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_7b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto013_2.tsv") %>%
  mutate(replicate=1, section="7b", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_8a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto015.tsv") %>%
  mutate(replicate=1, section="8a", condition="gg8tTA_dox5ug")
lowDOX_5ug_1_8b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_10x_notAuto015_2.tsv") %>%
  mutate(replicate=1, section="8b", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_1a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto010.tsv") %>%
  mutate(replicate=2, section="1a", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_1b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto010_2.tsv") %>%
  mutate(replicate=2, section="1b", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_2a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto011.tsv") %>%
  mutate(replicate=2, section="2a", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_2b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto011_2.tsv") %>%
  mutate(replicate=2, section="2b", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_3a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto012.tsv") %>%
  mutate(replicate=2, section="3a", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_3b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto012_2.tsv") %>%
  mutate(replicate=2, section="3b", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_4a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto013.tsv") %>%
  mutate(replicate=2, section="4a", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_4b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto013_2.tsv") %>%
  mutate(replicate=2, section="4b", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_5_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto015.tsv") %>%
  mutate(replicate=2, section="5", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_6a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto016.tsv") %>%
  mutate(replicate=2, section="6a", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_6b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto016_2.tsv") %>%
  mutate(replicate=2, section="6b", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_7a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto017.tsv") %>%
  mutate(replicate=2, section="7a", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_7b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto017_2.tsv") %>%
  mutate(replicate=2, section="7b", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_8a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto020.tsv") %>%
  mutate(replicate=2, section="8a", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_8b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto020_2.tsv") %>%
  mutate(replicate=2, section="8b", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_9a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto021.tsv") %>%
  mutate(replicate=2, section="9a", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_9b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto021_2.tsv") %>%
  mutate(replicate=2, section="9b", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_10a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto023.tsv") %>%
  mutate(replicate=2, section="10a", condition="gg8tTA_dox5ug")
lowDOX_5ug_2_10b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_5ug/P2_gg8_5ug_rep2_10x_notAuto023_2.tsv") %>%
  mutate(replicate=2, section="10b", condition="gg8tTA_dox5ug")

## read in tetoP2*/+; gg8tTA*/+ 1ug Dox imaging data (new OE embeded at ~8 weeks old). 14um sections
# 6 sections replicate1, 3 sections replicate2
lowDOX_1ug_1_1_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_1ug/P2_gg8_1ug_10x_notAuto021.tsv") %>%
  mutate(replicate=1, section="1", condition="gg8tTA_dox1ug")
lowDOX_1ug_1_2a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_1ug/P2_gg8_1ug_10x_notAuto024.tsv") %>%
  mutate(replicate=1, section="2a", condition="gg8tTA_dox1ug")
lowDOX_1ug_1_2b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_1ug/P2_gg8_1ug_10x_notAuto024.tsv") %>%
  mutate(replicate=1, section="2b", condition="gg8tTA_dox1ug")
lowDOX_1ug_1_3a_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_1ug/P2_gg8_1ug_10x_notAuto025.tsv") %>%
  mutate(replicate=1, section="3a", condition="gg8tTA_dox1ug")
lowDOX_1ug_1_3b_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_1ug/P2_gg8_1ug_10x_notAuto025.tsv") %>%
  mutate(replicate=1, section="3b", condition="gg8tTA_dox1ug")
lowDOX_1ug_1_4_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_1ug/P2_gg8_1ug_10x_notAuto028.tsv") %>%
  mutate(replicate=1, section="4", condition="gg8tTA_dox1ug")
lowDOX_1ug_2_1_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_1ug/P2_gg8_1ug_rep2_10x_notAuto006.tsv") %>%
  mutate(replicate=2, section="1", condition="gg8tTA_dox1ug")
lowDOX_1ug_2_2_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_1ug/P2_gg8_1ug_rep2_10x_notAuto010.tsv") %>%
  mutate(replicate=2, section="2", condition="gg8tTA_dox1ug")
lowDOX_1ug_2_3_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8_1ug/P2_gg8_1ug_rep2_10x_notAuto013.tsv") %>%
  mutate(replicate=2, section="3", condition="gg8tTA_dox1ug")

## read in tetoP2*/+; gg8tTA*/+ imaging data (new OE embeded at ~8 weeks old). 14um sections
# 4 sections replicate1, 2 sections replicate2
tetoP2_gg8_1_1_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8/P2_gg8_r2_10x_notAuto003_2.tsv") %>%
  mutate(replicate=1, section="1", condition="gg8tTA")
tetoP2_gg8_1_2_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8/P2_gg8_r2_10x_notAuto006_2.tsv") %>%
  mutate(replicate=1, section="2", condition="gg8tTA")
tetoP2_gg8_1_3_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8/P2_gg8_r2_10x_notAuto009_2.tsv") %>%
  mutate(replicate=1, section="3", condition="gg8tTA")
tetoP2_gg8_1_4_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8/P2_gg8_r2_10x_notAuto011_2.tsv") %>%
  mutate(replicate=1, section="4", condition="gg8tTA")
tetoP2_gg8_2_1_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8/P2_gg8_r3_10x_notAutofocus010_2.tsv") %>%
  mutate(replicate=2, section="1", condition="gg8tTA")
tetoP2_gg8_2_2_table<-read_tsv("/media/storageA/lisa/2021_zonal_paper_resubmission/analysis/imaging_tetoP2_expression/221101_lowDOX_final_imaging/P2_gg8/P2_gg8_r3_10x_notAutofocus011_2.tsv") %>%
  mutate(replicate=2, section="2", condition="gg8tTA")

#combine data
tetoP2_combined_table<-rbind(tetoP2_1_1a_table,tetoP2_1_1b_table,tetoP2_1_2a_table,tetoP2_1_2b_table,tetoP2_1_3_table,tetoP2_1_4_table,tetoP2_1_5_table,tetoP2_1_6_table,tetoP2_1_7_table,tetoP2_1_8_table,tetoP2_1_9_table,
                             tetoP2_2_1a_table,tetoP2_2_1b_table,tetoP2_2_2_table,tetoP2_2_3_table,tetoP2_2_4a_table,tetoP2_2_4b_table,tetoP2_2_5a_table,tetoP2_2_5b_table,tetoP2_2_6a_table,tetoP2_2_6b_table,tetoP2_2_7a_table,tetoP2_2_7b_table,
                             tetoP2_2_6a_table,tetoP2_2_6b_table,tetoP2_2_7a_table,tetoP2_2_7b_table,tetoP2_2_8a_table,tetoP2_2_8b_table,tetoP2_2_9a_table,tetoP2_2_9b_table,tetoP2_2_10a_table,tetoP2_2_10b_table)
dim(tetoP2_combined_table) #58 rows(data points)

lowDOX_5ug_combined_table<-rbind(lowDOX_5ug_1_1a_table,lowDOX_5ug_1_1b_table,lowDOX_5ug_1_2a_table,lowDOX_5ug_1_2b_table,lowDOX_5ug_1_3a_table,lowDOX_5ug_1_3b_table,lowDOX_5ug_1_4a_table,lowDOX_5ug_1_4b_table,lowDOX_5ug_1_5a_table,lowDOX_5ug_1_5b_table,
                                 lowDOX_5ug_1_5a_table,lowDOX_5ug_1_5b_table,lowDOX_5ug_1_6a_table,lowDOX_5ug_1_6b_table,lowDOX_5ug_1_7a_table,lowDOX_5ug_1_7b_table,lowDOX_5ug_1_8a_table,lowDOX_5ug_1_8b_table,
                                 lowDOX_5ug_2_1a_table,lowDOX_5ug_2_1b_table,lowDOX_5ug_2_2a_table,lowDOX_5ug_2_2b_table,lowDOX_5ug_2_3a_table,lowDOX_5ug_2_3b_table,lowDOX_5ug_2_4a_table,lowDOX_5ug_2_4b_table,lowDOX_5ug_2_5_table,
                                 lowDOX_5ug_2_6a_table,lowDOX_5ug_2_6b_table,lowDOX_5ug_2_7a_table,lowDOX_5ug_2_7b_table,lowDOX_5ug_2_8a_table,lowDOX_5ug_2_8b_table,lowDOX_5ug_2_9a_table,lowDOX_5ug_2_9b_table,lowDOX_5ug_2_10a_table,lowDOX_5ug_2_10b_table)
dim(lowDOX_5ug_combined_table) #803 rows(data points)

lowDOX_1ug_combined_table<-rbind(lowDOX_1ug_1_1_table,lowDOX_1ug_1_2a_table,lowDOX_1ug_1_2b_table,lowDOX_1ug_1_3a_table,lowDOX_1ug_1_3b_table,lowDOX_1ug_1_4_table,
                                 lowDOX_1ug_2_1_table,lowDOX_1ug_2_2_table,lowDOX_1ug_2_3_table)
dim(lowDOX_1ug_combined_table) #5051 rows (data points)

tetoP2_gg8_combined_table<-rbind(tetoP2_gg8_1_1_table,tetoP2_gg8_1_2_table,tetoP2_gg8_1_3_table,tetoP2_gg8_1_4_table,tetoP2_gg8_2_1_table,tetoP2_gg8_2_2_table)
dim(tetoP2_gg8_combined_table) #10497 rows(data points)

# take a random subset of data points for P2;gg8 and P2;gg8 1ug DOX first before combining genotypes into data frame
# take 1000 data points 
lowDOX_1ug_combined_table_1000<-lowDOX_1ug_combined_table[sample(nrow(lowDOX_1ug_combined_table),1000),]
tetoP2_gg8_combined_table_1000<-tetoP2_gg8_combined_table[sample(nrow(tetoP2_gg8_combined_table),1000),]
experiments_combined_subsetted<-rbind(lowDOX_5ug_combined_table,tetoP2_combined_table,lowDOX_1ug_combined_table_1000,tetoP2_gg8_combined_table_1000)

p.tetoP2gg8<-ggplot(filter(experiments_combined_subsetted,condition == "gg8tTA"), aes(x=normalized_position,y=0.5,col = factor(replicate))) + 
  #geom_point(size=3, alpha=0.7) + 
  geom_jitter(width = 0, height = 0.5, alpha=0.5,size=2)+
  xlim(0,1)+ #ylim (0,30)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  ylab("random jitter") + 
  xlab("normalized zonal position")+
  scale_color_manual(values = c("#006600","#009900","#33CC33","#009966"))+ 
  ggtitle("tetoP2;gg8tTA")
p.tetoP2gg8

p.1ugdox<-ggplot(filter(experiments_combined_subsetted,condition == "gg8tTA_dox1ug"), aes(x=normalized_position,y=0.5,col = factor(replicate))) + 
  #geom_point(size=3, alpha=0.7) + 
  geom_jitter(width = 0, height = 0.5, alpha=0.5,size=2)+
  xlim(0,1)+ #ylim (0,30)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  ylab("random jitter") + 
  xlab("normalized zonal position")+
  scale_color_manual(values = c("#006600","#009900","#33CC33","#009966"))+ 
  ggtitle("tetoP2;gg8tTA 1ug DOX")
p.1ugdox

p.5ugdox<-ggplot(filter(experiments_combined_subsetted,condition == "gg8tTA_dox5ug"), aes(x=normalized_position,y=0.5,col = factor(replicate))) + 
  #geom_point(size=3, alpha=0.7) + 
  geom_jitter(width = 0, height = 0.5, alpha=0.5,size=2)+
  xlim(0,1)+ #ylim (0,30)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  ylab("random jitter") + 
  xlab("normalized zonal position")+
  scale_color_manual(values = c("#006600","#009900","#33CC33","#009966"))+ 
  ggtitle("tetoP2;gg8tTA 5ug DOX")
p.5ugdox

p.tetoP2<-ggplot(filter(experiments_combined_subsetted,condition == "no_tTA"), aes(x=normalized_position,y=0.5,col = factor(replicate))) + 
  #geom_point(size=3, alpha=0.7) + 
  geom_jitter(width = 0, height = 0.5, alpha=0.5,size=2)+
  xlim(0,1)+ #ylim (0,30)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  ylab("random jitter") + 
  xlab("normalized zonal position")+
  scale_color_manual(values = c("#006600","#009900","#33CC33","#009966"))+ 
  ggtitle("tetoP2")
p.tetoP2

grid.arrange(p.tetoP2gg8,p.1ugdox,p.5ugdox,p.tetoP2, ncol = 1,nrow=4)
g<-arrangeGrob(p.tetoP2gg8,p.1ugdox,p.5ugdox,p.tetoP2, ncol = 1,nrow=4)
ggsave(file="scatterplot.tetoP2_lowDOX_221101_replicates_shown.pdf",g, width = 12, height = 8, units = "in", useDingbats = FALSE )

# do facet grid
p.combined<-ggplot(experiments_combined_subsetted, aes(x=normalized_position,y=0.5,col = factor(replicate))) + 
  geom_jitter(width = 0, height = 0.5, alpha=0.5,size=2)+
  xlim(0,1)+ #ylim (0,30)+
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10),strip.text.y = element_text(color="black",size = 12))+ 
  ylab("random scatter") + xlab("normalized zonal position")+
  scale_color_manual(values = c("#006600","#009900","#33CC33","#009966", "#006633","#66CC00"))+ 
  ggtitle("tetoP2 expression")+
  facet_grid(rows = "condition")+
  #scale_x_continuous(breaks = round(seq(min(dat$x), max(dat$x), by = 0.5),1)) +
  scale_y_continuous(breaks = round(seq(0, 1, by = 0.5),1),minor_breaks = 0)
p.combined
ggsave(file="scatterplot.tetoP2_lowDOX_221101_facet_grid_show_replicates.pdf",p.combined, width = 12, height = 8, units = "in" , useDingbats = FALSE )


## add density plots
p.P2_gg8_density<-ggplot(filter(experiments_combined_subsetted,condition == "gg8tTA")) + 
  #geom_point(size=3, alpha=0.7) + 
  geom_density(aes(x=normalized_position,..scaled..),alpha=.2, fill="#33CC33", color="#33CC33")+
  #geom_jitter(aes(x=normalized_position,y=0.5,col = factor(replicate)),width = 0, height = 0.5, alpha=0.5,size=2)+
  geom_jitter(aes(x=normalized_position,y=0.5),col = "#006600",width = 0, height = 0.5, alpha=0.5,size=2)+
  xlim(0,1)+ #ylim (0,30)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  ylab("random jitter") + 
  xlab("normalized zonal position")+
  scale_color_manual(values = c("#006600","#009900","#33CC33","#009966"))+ 
  ggtitle("tetoP2;gg8tTA")
p.P2_gg8_density

p.1ugdox_density<-ggplot(filter(experiments_combined_subsetted,condition == "gg8tTA_dox1ug")) + 
  #geom_point(size=3, alpha=0.7) + 
  geom_density(aes(x=normalized_position,..scaled..),alpha=.2, fill="#33CC33", color="#33CC33")+
  #geom_jitter(aes(x=normalized_position,y=0.5,col = factor(replicate)),width = 0, height = 0.5, alpha=0.5,size=2)+
  geom_jitter(aes(x=normalized_position,y=0.5),col = "#006600",width = 0, height = 0.5, alpha=0.5,size=2)+
  xlim(0,1)+ #ylim (0,30)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  ylab("random jitter") + 
  xlab("normalized zonal position")+
  scale_color_manual(values = c("#006600","#009900","#33CC33","#009966"))+ 
  ggtitle("tetoP2;gg8tTA 1ug DOX")
p.1ugdox_density

p.5ugdox_density<-ggplot(filter(experiments_combined_subsetted,condition == "gg8tTA_dox5ug")) + 
  #geom_point(size=3, alpha=0.7) + 
  geom_density(aes(x=normalized_position,..scaled..),alpha=.2, fill="#33CC33", color="#33CC33")+
  #geom_jitter(aes(x=normalized_position,y=0.5,col = factor(replicate)),width = 0, height = 0.5, alpha=0.5,size=2)+
  geom_jitter(aes(x=normalized_position,y=0.5),col = "#006600",width = 0, height = 0.5, alpha=0.5,size=2)+
  xlim(0,1)+ #ylim (0,30)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  ylab("random jitter") + 
  xlab("normalized zonal position")+
  scale_color_manual(values = c("#006600","#009900","#33CC33","#009966"))+ 
  ggtitle("tetoP2;gg8tTA 5ug DOX")
p.5ugdox_density

p.P2_density<-ggplot(filter(experiments_combined_subsetted,condition == "no_tTA")) + 
  #geom_point(size=3, alpha=0.7) + 
  geom_density(aes(x=normalized_position,..scaled..),alpha=.2, fill="#33CC33", color="#33CC33")+
  #geom_jitter(aes(x=normalized_position,y=0.5,col = factor(replicate)),width = 0, height = 0.5, alpha=0.5,size=2)+
  geom_jitter(aes(x=normalized_position,y=0.5),col = "#006600",width = 0, height = 0.5, alpha=0.5,size=2)+
  xlim(0,1)+ #ylim (0,30)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  ylab("random jitter") + 
  xlab("normalized zonal position")+
  scale_color_manual(values = c("#006600","#009900","#33CC33","#009966"))+ 
  ggtitle("tetoP2")
p.P2_density

grid.arrange(p.P2_gg8_density,p.1ugdox_density,p.5ugdox_density,p.P2_density, ncol = 1,nrow=4)
g<-arrangeGrob(p.P2_gg8_density,p.1ugdox_density,p.5ugdox_density,p.P2_density, ncol = 1,nrow=4)
ggsave(file="scatterplot.tetoP2_lowDOX_221101_density_jitter_combined_reps.pdf",g, width = 12, height = 8, units = "in" , useDingbats = FALSE )

# do facet grid
p.combined<-ggplot(experiments_combined_subsetted, aes(x=normalized_position,y=0.5)) + 
  geom_density(aes(x=normalized_position,..scaled..),alpha=.2, fill="#33CC33", color="#33CC33")+
  geom_jitter(aes(x=normalized_position,y=0.5),col = "#006600",width = 0, height = 0.5, alpha=0.5,size=2)+
  xlim(0,1)+ #ylim (0,30)+
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10),strip.text.y = element_text(color="black",size = 12))+ 
  ylab("random scatter") + xlab("normalized zonal position")+
  scale_color_manual(values = c("#006600","#009900","#33CC33","#009966", "#006633","#66CC00"))+ 
  ggtitle("tetoP2 expression")+
  facet_grid(rows = "condition")+
  #scale_x_continuous(breaks = round(seq(min(dat$x), max(dat$x), by = 0.5),1)) +
  scale_y_continuous(breaks = round(seq(0, 1, by = 0.5),1),minor_breaks = 0)
p.combined
ggsave(file="scatterplot.tetoP2_lowDOX_221101_density_jitter_combined_reps_facet_grid2.pdf",p.combined, width = 10, height = 8, units = "in", useDingbats = FALSE  )


# combined density
p.combined_density<-ggplot(experiments_combined_subsetted) + 
  geom_density(aes(x=normalized_position,color=condition, fill= condition,..scaled..),alpha=.2)+
  xlim(0,1)+ #ylim (0,30)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =15), axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  ylab("random jitter") + 
  xlab("normalized zonal position")+
  #scale_color_manual(values = c("#006600","#009900","#33CC33","#009966"))+ 
  ggtitle("tetoP2 expression")
p.combined_density

ggsave(file="densityplot.tetoP2_lowDOX_221101.pdf",p.combined_density, width = 12, height = 3, units = "in" , useDingbats = FALSE )


