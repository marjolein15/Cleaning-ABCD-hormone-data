#Script to clean and categorize medication use data in ABCD baseline dataset. 
#Focusing on medication use in last 2 weeks (medsy01.txt)
#Created by MEAB on 2021-10-29
#last edited 2022-01-31

#Libraries and set folders
library(dplyr)
library(psych)
library(tidyr)
library(stringr)
origdatadir <- 'C:/path/to/ABCD/package/4.0 - Package_1194258/'
scriptsdir <- 'C:/path/to/folder/with/this/script/'

#Read in Medication codes #These contain RXCUIs and corresponding MESHPA codes matched using RxMix (https://mor.nlm.nih.gov/RxMix/#)
RXCUI_to_MESHPA <- read.csv(paste0(scriptsdir,"rxmix_output_meshpa_20211101.txt"),header=T,sep="\t")

#Read in the ABCD med use data
header_meduse_raw <- read.table(paste0(origdatadir,"medsy01.txt", sep=''), header = F, nrows = 1, as.is = T) 
meduse_raw <- read.csv(paste0(origdatadir,"medsy01.txt"), skip = 2, header = F, sep="\t",stringsAsFactors = F) 
colnames(meduse_raw)<- header_meduse_raw

#This finds the MESHPA codes (classId) that match the RXnorm code for each of medications the parent reported on 
#There can be multiple codes per RXnorm code
#MESHPA=Medical Subject Headings (MeSH) pharmacological actions

#RX = prescription drugs. Parents could report up to 15 rx meds
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx1_MESHPA_codes[i]=ifelse(meduse_raw$med1_rxnorm_p[i]=="", NA, 
                                      do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                        paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med1_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx2_MESHPA_codes[i]=ifelse(meduse_raw$med2_rxnorm_p[i]=="", NA, 
                                      do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                    paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med2_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx3_MESHPA_codes[i]=ifelse(meduse_raw$med3_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                         paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med3_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx4_MESHPA_codes[i]=ifelse(meduse_raw$med4_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                           paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med4_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx5_MESHPA_codes[i]=ifelse(meduse_raw$med5_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                           paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med5_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx6_MESHPA_codes[i]=ifelse(meduse_raw$med6_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                             paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med6_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx7_MESHPA_codes[i]=ifelse(meduse_raw$med7_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                             paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med7_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx8_MESHPA_codes[i]=ifelse(meduse_raw$med8_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                             paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med8_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx9_MESHPA_codes[i]=ifelse(meduse_raw$med9_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                            paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med9_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx10_MESHPA_codes[i]=ifelse(meduse_raw$med10_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                            paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med10_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx11_MESHPA_codes[i]=ifelse(meduse_raw$med11_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                           paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med11_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx12_MESHPA_codes[i]=ifelse(meduse_raw$med12_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                           paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med12_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx13_MESHPA_codes[i]=ifelse(meduse_raw$med13_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                            paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med13_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx14_MESHPA_codes[i]=ifelse(meduse_raw$med14_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                           paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med14_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$rx15_MESHPA_codes[i]=ifelse(meduse_raw$med15_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                          paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med15_rxnorm_p[i]))),"$")), "classId"])))) }


#OTC = Over the counter drugs. Parents could report up to 15 meds here
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc1_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_1_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                     paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_1_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc2_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_2_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                     paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_2_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc3_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_3_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                     paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_3_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc4_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_4_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                     paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_4_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc5_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_5_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                     paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_5_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc6_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_6_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                     paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_6_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc7_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_7_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                     paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_7_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc8_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_8_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                     paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_8_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc9_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_9_rxnorm_p[i]=="", NA, 
                                        do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                     paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_9_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc10_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_10_rxnorm_p[i]=="", NA, 
                                         do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                      paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_10_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc11_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_11_rxnorm_p[i]=="", NA, 
                                         do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                      paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_11_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc12_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_12_rxnorm_p[i]=="", NA, 
                                         do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                      paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_12_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc13_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_13_rxnorm_p[i]=="", NA, 
                                         do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                      paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_13_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc14_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_14_rxnorm_p[i]=="", NA, 
                                         do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                      paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_14_rxnorm_p[i]))),"$")), "classId"])))) }
for (i in 1:nrow(meduse_raw)) { 
  meduse_raw$otc15_MESHPA_codes[i]=ifelse(meduse_raw$med_otc_15_rxnorm_p[i]=="", NA, 
                                         do.call(paste,as.list(as.character(RXCUI_to_MESHPA[str_which(RXCUI_to_MESHPA$I0,
                                                                                                      paste0("^",as.character(as.numeric(gsub("([0-9]+).*$", "\\1", meduse_raw$med_otc_15_rxnorm_p[i]))),"$")), "classId"])))) }



#Classify into yes/no variables for drug classes
meduse_categorized <- meduse_raw 
# D003271,D003272,D000080066 = Contraceptive agents
meduse_categorized$contraceptive <- !!rowSums(sapply(meduse_categorized[grepl("_MESHPA_codes",names(meduse_categorized))], grepl, pattern = "D003271|D003272|D000080066"))
# D005938 = Glucocorticoids
meduse_categorized$glucocorticoid <- !!rowSums(sapply(meduse_categorized[grepl("_MESHPA_codes",names(meduse_categorized))], grepl, pattern = "D005938"))

#This filters to the necessary variables and timepoints
meduse_categorized_filtered <- meduse_categorized %>% 
  filter(eventname == "baseline_year_1_arm_1")  %>% 
  select(subjectkey,src_subject_id,eventname,interview_date,sex,contraceptive,glucocorticoid)
meduse_categorized_baseline_2year <- meduse_categorized %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1")  %>% 
  select(subjectkey,src_subject_id,eventname,interview_date,sex,contraceptive,glucocorticoid)

#Save
write.csv(meduse_categorized_filtered,row.names = F,paste0(scriptsdir,"meduse_covariates_forhormones.csv"))
write.csv(meduse_categorized_baseline_2year,row.names = F,paste0(scriptsdir,"meduse_covariates_forhormones_2tps.csv"))
