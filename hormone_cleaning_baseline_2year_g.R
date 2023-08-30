############################
#
# Script for cleaning ABCD baseline hormone data 
# Created by MEAB, last edited 2023-2-14
#
############################

#STEPS:
#Follow decision tree in Herting et al. 2021
#Clean confounds: Find nonsensical values in time confounds and set these to NA
#Add in medication use variables. Using glucocorticoid med use and contraceptive use 
#Calculate mentrual cycle phase. One variable indicating regularity of cycle, others indicating cycle phase 
#Log-transform
#Regress hormones on confounds and extract residuals (hours since wake, collection duration, recent caffeine use, recent physical activity,  med use)



###################
# Loading packages and data
###################

#loading packages
library(dplyr)
library(ggplot2)
library(chron)
library(mfp)
data_dir="C:/path/to/ABCD/package/4.0 - Package_1194258/"
meduse_dir="C:/path/to/folder/with/output/of/medication_use/script/"
hormones_dir="C:/path/to/output/folder/"


#loading data
header_abcd_hsss01 <- read.table(paste0(data_dir,"abcd_hsss01.txt", sep=''), header = F, nrows = 1, as.is = T) 
abcd_hsss01 <- read.csv(paste0(data_dir,"abcd_hsss01.txt"), skip = 2, header = F, sep="\t",stringsAsFactors=F,na.strings = c("NA",""),col.names=header_abcd_hsss01) 
abcd_hsss01 <- abcd_hsss01 %>% 
  filter(eventname=="baseline_year_1_arm_1"|eventname=="2_year_follow_up_y_arm_1") %>% 
  select(subjectkey,src_subject_id,eventname,hormone_scr_dhea_mean,hormone_scr_dhea_rep1,	hormone_scr_dhea_rep1_ll,hormone_scr_dhea_rep1_qns,	hormone_scr_dhea_rep1_nd,
         hormone_scr_dhea_rep2,	hormone_scr_dhea_rep2_ll,	hormone_scr_dhea_rep2_qns,	hormone_scr_dhea_rep2_nd,
         hormone_scr_hse_mean,	hormone_scr_hse_rep1,	hormone_scr_hse_rep1_ll,	hormone_scr_hse_rep1_qns,	hormone_scr_hse_rep1_nd,
         hormone_scr_hse_rep2,	hormone_scr_hse_rep2_ll,	hormone_scr_hse_rep2_qns,	hormone_scr_hse_rep2_nd,
         hormone_scr_ert_mean,	hormone_scr_ert_rep1, hormone_scr_ert_rep1_ll,	hormone_scr_ert_rep1_qns,	hormone_scr_ert_rep1_nd,
         hormone_scr_ert_rep2,	hormone_scr_ert_rep2_ll,	hormone_scr_ert_rep2_qns,	hormone_scr_ert_rep2_nd)

header_sph01 <- read.table(paste0(data_dir,"sph01.txt", sep=''), header = F, nrows = 1, as.is = T) 
sph01 <- read.csv(paste0(data_dir,"sph01.txt"), skip = 2, header = F, sep="\t",stringsAsFactors=F, na.strings = c("NA",""),col.names=header_sph01) 
sph01 <- sph01 %>% 
  filter(visit=="baseline_year_1_arm_1"|visit=="2_year_follow_up_y_arm_1") %>% 
  select(subjectkey,src_subject_id,visit, sex, hormone_sal_wake_y,hormone_sal_caff_y, hormone_sal_caff_mg_y, hormone_sal_active,
         hormone_sal_active_minutes_y, hormone_sal_start_y,hormone_sal_end_y, hormons_sal_freeztemp_y,hormone_sal_freezer_y,
         hormone_sal_sex,hormon_sal_notes_y___1,hormon_sal_notes_y___2,hormon_sal_notes_y___3,hormon_sal_notes_y___4,
         hormon_sal_notes_y___5,hormon_sal_notes_y___6) %>% 
  rename(eventname=visit)

header_abcd_lt01 <- read.table(paste0(data_dir,"abcd_lt01.txt", sep=''), header = F, nrows = 1, as.is = T) 
abcd_lt01 <- read.csv(paste0(data_dir,"abcd_lt01.txt"), skip = 2, header = F, sep="\t",stringsAsFactors=F, na.strings = c("NA",""),col.names=header_abcd_lt01) 
abcd_lt01 <- abcd_lt01 %>% filter(eventname=="baseline_year_1_arm_1"|eventname=="2_year_follow_up_y_arm_1")  %>% 
  select(subjectkey,src_subject_id,eventname, site_id_l)

#Create full dataset (N=11,878)
predictors_inclcrit_hormones <- merge(abcd_hsss01, sph01, by=c("subjectkey","src_subject_id","eventname")) 


###################
# filtering and cleaning following the decision tree in Herting et al 2021 (doi: 10.3389/fendo.2020.549928)
###################

#step1 and 2: saliva sex matches main sex var, and exclude if sample not completed
# saliva sex and main sex don't match n=40, no saliva sex n=2614, unable to complete sample n=76, refused sample n=43, other reason not completed n=287
# n=19,230 left. About 2500 missing hormone_sal_sex, no spit sample collected
predictors_inclcrit_hormones %>% group_by(hormone_sal_sex,sex) %>% summarize(ns=n())
filtereddata_hormones <- predictors_inclcrit_hormones %>% 
  mutate(sex_rec=ifelse(sex=="F",1,ifelse(sex=="M",2,NA)) ) %>%
           filter( hormone_sal_sex==sex_rec & (hormone_sal_sex==1 | hormone_sal_sex==2) )

#step3: NA (removed samples not yet assayed, only applicable if assays are incomplete)

#step4: exclude if above higher limits of assay AND 1 or more problems with sample quality
# 0 participants excluded
filtereddata_hormones <- filtereddata_hormones %>%
  mutate(dhea1=as.numeric(ifelse(as.numeric(hormone_scr_dhea_rep1)>1000 & (hormon_sal_notes_y___2==1 | hormon_sal_notes_y___3==1 | hormon_sal_notes_y___4==1 | hormon_sal_notes_y___5==1 | hormon_sal_notes_y___6==1),
                      NA,hormone_scr_dhea_rep1)),
         dhea2=as.numeric(ifelse(as.numeric(hormone_scr_dhea_rep2)>1000 & (hormon_sal_notes_y___2==1 | hormon_sal_notes_y___3==1 | hormon_sal_notes_y___4==1 | hormon_sal_notes_y___5==1 | hormon_sal_notes_y___6==1),
                      NA,hormone_scr_dhea_rep2)),
         estr1=as.numeric(ifelse(as.numeric(hormone_scr_hse_rep1)>32 & (hormon_sal_notes_y___2==1 | hormon_sal_notes_y___3==1 | hormon_sal_notes_y___4==1 | hormon_sal_notes_y___5==1 | hormon_sal_notes_y___6==1),
                      NA,hormone_scr_hse_rep1)),
         estr2=as.numeric(ifelse(as.numeric(hormone_scr_hse_rep2)>32 & (hormon_sal_notes_y___2==1 | hormon_sal_notes_y___3==1 | hormon_sal_notes_y___4==1 | hormon_sal_notes_y___5==1 | hormon_sal_notes_y___6==1),
                      NA,hormone_scr_hse_rep2)),
         test1=as.numeric(ifelse(as.numeric(hormone_scr_ert_rep1)>600 & (hormon_sal_notes_y___2==1 | hormon_sal_notes_y___3==1 | hormon_sal_notes_y___4==1 | hormon_sal_notes_y___5==1 | hormon_sal_notes_y___6==1),
                      NA,hormone_scr_ert_rep1)),
         test2=as.numeric(ifelse(as.numeric(hormone_scr_ert_rep2)>600 & (hormon_sal_notes_y___2==1 | hormon_sal_notes_y___3==1 | hormon_sal_notes_y___4==1 | hormon_sal_notes_y___5==1 | hormon_sal_notes_y___6==1),
                      NA,hormone_scr_ert_rep2)))

#step5: if too low for detection (nd), change replicate's hormone value to 0
# non-detectables baseline DHEA r1 n=85 , r2 n=111 ; T r1 n=1 , r2 n=7 ; E2 r1 n=120 , r2 n=127
# non-detectables 2y follow up DHEA r1 n=52 , r2 n=66 T r1 n=0, r2 n=0 ; E2 r1 n=2 , r2 n=2
# length(which(filtereddata_hormones$hormone_scr_hse_rep2_nd=="1"))
filtereddata_hormones <- filtereddata_hormones %>%
  mutate(dhea1=ifelse(as.numeric(hormone_scr_dhea_rep1_nd)==1,0,dhea1),
         dhea2=ifelse(as.numeric(hormone_scr_dhea_rep2_nd)==1,0,dhea2),
         estr1=ifelse(as.numeric(hormone_scr_hse_rep1_nd)==1,0,estr1),
         estr2=ifelse(as.numeric(hormone_scr_hse_rep2_nd)==1,0,estr2),
         test1=ifelse(as.numeric(hormone_scr_ert_rep1_nd)==1,0,test1),
         test2=ifelse(as.numeric(hormone_scr_ert_rep2_nd)==1,0,test2) )

#step6: obtain a single value per hormone by averaging the two replicates or using to only value available
filtereddata_hormones$dhea <- rowMeans(filtereddata_hormones[c("dhea1","dhea2")],na.rm = T)
filtereddata_hormones$estr <- rowMeans(filtereddata_hormones[c("estr1","estr2")],na.rm = T)
filtereddata_hormones$test <- rowMeans(filtereddata_hormones[c("test1","test2")],na.rm = T)


###############
# Add in medication use vars. Med use was categorized using the script categorize_medications.R
# Using glucocorticoid med use and contraceptives. 
###############

meduse <- read.csv(paste0(meduse_dir, "meduse_covariates_forhormones_2tps.csv"))
filtereddata_hormones <- merge(filtereddata_hormones, meduse, by=c("subjectkey","src_subject_id","sex","eventname")) 


###############
#Check missingness and remove those with no hormone values on any hormone
###############

nrow(filtereddata_hormones[which(is.na(filtereddata_hormones$dhea)&is.na(filtereddata_hormones$test)&is.na(filtereddata_hormones$estr)),]) #3724
nrow(filtereddata_hormones[which(is.na(filtereddata_hormones$dhea)),]) #3803
nrow(filtereddata_hormones[which(is.na(filtereddata_hormones$test)),]) #3821
nrow(filtereddata_hormones[which(filtereddata_hormones$hormone_sal_sex==1 & is.na(filtereddata_hormones$estr)),]) #1859
filtereddata_hormones <- filtereddata_hormones %>% 
  filter(!(is.na(filtereddata_hormones$dhea) & is.na(filtereddata_hormones$test) & is.na(filtereddata_hormones$estr)))
#remaining N=15,506



###############
#Cleaning the confound variables (e.g. time of day, start and end time)
###############

#Missings: start time n=10, end time n=23, wake time n=48, freeze time n=31, caffeine use n=28, physical act n=39
nrow(filtereddata_hormones[which(is.na(filtereddata_hormones$hormone_sal_start_y)),])
nrow(filtereddata_hormones[which(is.na(filtereddata_hormones$hormone_sal_end_y)),])
nrow(filtereddata_hormones[which(is.na(filtereddata_hormones$hormone_sal_wake_y)),])
nrow(filtereddata_hormones[which(is.na(filtereddata_hormones$hormone_sal_freezer_y)),])
nrow(filtereddata_hormones[which(is.na(filtereddata_hormones$hormone_sal_caff_y)),])
nrow(filtereddata_hormones[which(is.na(filtereddata_hormones$hormone_sal_active)),])

#make a time variable of the time confounds
hormones_confounds <- filtereddata_hormones %>% mutate(start_time=times(paste0(hormone_sal_start_y, ":00")),
                                                          end_time=times(paste0(hormone_sal_end_y, ":00")),
                                                          wake_time=times(paste0(hormone_sal_wake_y, ":00")),
                                                          freezer_time=times(paste0(hormone_sal_freezer_y, ":00")))
#find nonsensical values in waking time, start time
#waking before 5am or after 2pm
hormones_confounds[which(hormones_confounds$wake_time<"05:00:00" | hormones_confounds$wake_time>"13:59:00"),c("subjectkey","wake_time","start_time","end_time","freezer_time")]
#starting sample before 7am or after 9pm
hormones_confounds[which(hormones_confounds$start_time<"07:00:00" | hormones_confounds$start_time>"21:00:00"),c("subjectkey","wake_time","start_time","end_time","freezer_time")]

#set nonsensical values in waking time & start time to NA
hormones_confounds <- hormones_confounds %>% 
  mutate(wake_time=as.times(ifelse(wake_time>"04:59:00"&wake_time<"14:00:00",wake_time, NA)),
         start_time=as.times(ifelse(start_time>="07:00:00"&start_time<="21:00:00",start_time,NA)))

#create variables for time since waking, collection duration and time to freeze
hormones_confounds <- hormones_confounds %>% mutate(time_since_wake=start_time - wake_time,
                                                    collect_duration=end_time - start_time,
                                                    time_to_freeze=freezer_time - end_time)

#find nonsensical values in derived time variables
#negative collection duration or longer than 1hour 
hormones_confounds[which(hormones_confounds$collect_duration<0),c("subjectkey","wake_time","start_time","end_time","freezer_time")]
hormones_confounds[which(hormones_confounds$collect_duration>"01:00:00"),c("subjectkey","wake_time","start_time","end_time","freezer_time","collect_duration")]
#starting less than 30min after waking
hormones_confounds[which(hormones_confounds$time_since_wake<"00:30:00"),c("subjectkey","wake_time","start_time","end_time","freezer_time")]
#time to get sample to freezer more than 6 hours 
hormones_confounds[which(hormones_confounds$time_to_freeze>"06:00:00"),c("subjectkey","wake_time","start_time","end_time","freezer_time")]

#set nonsensical values in derived time variables to NA
filtered_hormones_confounds <- hormones_confounds %>% 
  mutate(time_since_wake=ifelse(time_since_wake>="00:30:00",time_since_wake, NA),
  collect_duration=ifelse(collect_duration>=0&collect_duration<="01:00:00",collect_duration, NA))


#create numerical variables of the time confounds
filtered_hormones_confounds$hours_since_wake <- as.numeric(filtered_hormones_confounds$time_since_wake)*24
filtered_hormones_confounds$minutes_to_collect <- as.numeric(filtered_hormones_confounds$collect_duration)*1440
#create factors of the caffeine and phys activity confounds
filtered_hormones_confounds$hormone_sal_caff_y <- as.factor(filtered_hormones_confounds$hormone_sal_caff_y)
filtered_hormones_confounds$hormone_sal_active <- as.factor(filtered_hormones_confounds$hormone_sal_active)


############
#Some checks and visualizations to see if data is behaving as expected
############

#plotting distribution of hormones
ggplot(data=filtered_hormones_confounds,aes(x=as.numeric(dhea))) + geom_histogram(bins=75)
ggplot(data=filtered_hormones_confounds,aes(x=as.numeric(test))) + geom_histogram(bins=75)
ggplot(data=filtered_hormones_confounds,aes(x=as.numeric(estr))) + geom_histogram(bins=75)

ggplot(data=filtered_hormones_confounds,aes(x=log(as.numeric(dhea)+1))) + geom_histogram(bins=75)
ggplot(data=filtered_hormones_confounds,aes(x=log(as.numeric(test)+1))) + geom_histogram(bins=75)
ggplot(data=filtered_hormones_confounds,aes(x=log(as.numeric(estr)+1))) + geom_histogram(bins=75)
#adding a constant because many values are <1, which will lead to negative log values

#plotting hormones by time since waking
ggplot(data=filtered_hormones_confounds,aes(x=hours_since_wake,y=log(as.numeric(dhea)+1))) + 
  geom_point() + geom_smooth() +
  labs(x="Hours since waking", y="log DHEA (pg/ml)")
ggplot(data=filtered_hormones_confounds,aes(x=hours_since_wake,y=log(as.numeric(test)+1))) + 
  geom_point() + geom_smooth()+
  labs(x="Hours since waking", y="log testosterone (pg/ml)")
ggplot(data=filtered_hormones_confounds,aes(x=hours_since_wake,y=log(as.numeric(estr)+1))) + 
  geom_point() + geom_smooth() +
  labs(x="Hours since waking", y="log estradiol (pg/ml)")

#plotting hormones by collection duration 
ggplot(data=filtered_hormones_confounds,aes(x=minutes_to_collect,y=log(as.numeric(dhea)+1))) + 
       geom_point() + geom_smooth(method="loess")+
       labs(x="Collection duration (min)", y="log DHEA (pg/ml)")
ggplot(data=filtered_hormones_confounds,aes(x=minutes_to_collect,y=log(as.numeric(test)+1))) + 
       geom_point() + geom_smooth(method="loess")+
       labs(x="Collection duration (min)", y="log testosterone (pg/ml)")
ggplot(data=filtered_hormones_confounds,aes(x=minutes_to_collect,y=log(as.numeric(estr)+1))) + 
       geom_point() + geom_smooth(method="loess")+
       labs(x="Collection duration (min)", y="log estradiol (pg/ml)")


#Check if confounds vary by site
site_check <- merge(filtered_hormones_confounds,abcd_lt01,by=c("subjectkey","src_subject_id","eventname"))
site_check$site_id_l <- as.factor(site_check$site_id_l)
ggplot(data=site_check)+geom_violin(aes(x=site_id_l,y=hours_since_wake,color=site_id_l)) # a handful of sites focused on morning data collection, the rest spread it throughout the day
ggplot(data=site_check)+geom_violin(aes(x=site_id_l,y=minutes_to_collect,color=site_id_l)) # no clear site differences
ggplot(data=site_check)+geom_violin(aes(x=site_id_l,y=as.numeric(time_to_freeze*1440),color=site_id_l))
ggplot(data=site_check,aes(x=hours_since_wake,y=dhea,color=site_id_l))+geom_point()+geom_smooth() #no sites that deviate
ggplot(data=site_check,aes(x=hours_since_wake,y=test,color=site_id_l))+geom_point()+geom_smooth() #no sites that deviate
ggplot(data=site_check,aes(x=hours_since_wake,y=estr,color=site_id_l))+geom_point()+geom_smooth() #site 22 deviates but is a small site


###############
# Calculate mentrual cycle phase based on self-reported start of last period
# using counting method in Schmalenberger et al 2021,assuming cycle of 28 days since we only have last period and (in most cases) unknown cycle length
# 
###############

header_abcd_ypdms01 <- read.table(paste0(data_dir,"abcd_ypdms01.txt", sep=''), header = F, nrows = 1, as.is = T) 
abcd_ypdms01 <- read.csv(paste0(data_dir,"abcd_ypdms01.txt"), skip = 2, header = F, sep="\t",stringsAsFactors = F) 
colnames(abcd_ypdms01)<- header_abcd_ypdms01
#find post-menarche girls
abcd_ypdms01 <- abcd_ypdms01 %>% 
  filter(eventname=="baseline_year_1_arm_1" & pds_f5_y==4|eventname=="2_year_follow_up_y_arm_1" & pds_f5_y==4) %>% 
  select(subjectkey,src_subject_id,eventname, interview_date,sex,pds_f5_y, pds_f6_y, pds_f6_y_dk,menstrualcycle1_y,menstrualcycle2_y,
         menstrualcycle2_y_dk,menstrualcycle3_y,menstrualcycle4_y,menstrualcycle5_y,menstrualcycle6_y)
#calculate number of days since last period
abcd_ypdms01 <- abcd_ypdms01 %>% 
  mutate(days_since_period=dates(interview_date,format="m/d/Y") - dates(menstrualcycle1_y,format="Y-m-d"))
#sort into menstrual cycle phases 
#Phases recommended by Schmalenberger et al. 2021 (https://doi.org/10.1016/j.psyneuen.2020.1048950) did not capture the whole cycle,
# and left 60-ish% of girls unclassifiable. Adapted this with Birdie Shirtcliff to cover all days up to 36 days since last period.
abcd_ypdms01 <- abcd_ypdms01 %>% 
  mutate(cycle_phase = as.factor(ifelse(days_since_period>=0 & days_since_period<4,"perimenstrual",
                                   ifelse(days_since_period>=4&days_since_period<9,"midfollicular",
                                      ifelse(days_since_period>=9&days_since_period<13,"latefollicular",
                                        ifelse(days_since_period>=13&days_since_period<17,"ovulatory",
                                           ifelse(days_since_period>=17&days_since_period<24,"midluteal",
                                              ifelse(days_since_period>=23&days_since_period<37,"lateluteal",
                                                     ifelse(days_since_period>=37,"over36daysago",NA)))))))))
#menstrual3=Is your menstrual cycle regular? 0=no, 1=yes,2=don't know, 3=refuse. #count 2 as irregular, 3 as missing
#Also count as irregular cycling if last period was over 36 days ago
abcd_ypdms01 <- abcd_ypdms01 %>% 
  mutate(cycle_regularity = as.factor(ifelse(menstrualcycle3_y==1&days_since_period<37,"regular",
                                        ifelse(menstrualcycle3_y==0|menstrualcycle3_y==2|days_since_period>=37,"irregular",
                                               NA))))

abcd_ypdms01 %>% group_by(cycle_phase) %>% summarize(no=n()) 
abcd_ypdms01 %>% group_by(cycle_regularity) %>% summarize(no=n()) 
# about 40% is irregular in our classification system

#Get these variables into the full dataset
postmen_hormones_confounds <-  merge(filtered_hormones_confounds,abcd_ypdms01[c("subjectkey","eventname","days_since_period","cycle_phase")],by=c("subjectkey","eventname"),all.y=T)
allmen_hormones_confounds <-  merge(filtered_hormones_confounds,
                                    abcd_ypdms01[c("subjectkey","eventname","days_since_period","cycle_phase","cycle_regularity")],
                                    by=c("subjectkey","eventname"),all=T) %>%
                                    filter(sex=="F")
#Make binary variables to use as control variables
allmen_hormones_confounds <- allmen_hormones_confounds %>% 
  mutate(irregular_cycle = as.factor(ifelse(cycle_regularity=="irregular",1,0)),
         perimentrual_phase = as.factor(ifelse(cycle_phase=="perimenstrual",1,0)),
         midfollicular_phase = as.factor(ifelse(cycle_phase=="midfollicular",1,0)),
         latefollicular_phase = as.factor(ifelse(cycle_phase=="latefollicular",1,0)),
         ovulatory_phase = as.factor(ifelse(cycle_phase=="ovulatory",1,0)),
         midluteal_phase = as.factor(ifelse(cycle_phase=="midluteal",1,0)),
         lateluteal_phase = as.factor(ifelse(cycle_phase=="lateluteal",1,0))) %>%
  replace_na(list(irregular_cycle=0,perimentrual_phase=0,midfollicular_phase=0,latefollicular_phase=0,
                  ovulatory_phase=0,midluteal_phase=0,lateluteal_phase=0))


#Plots to see distributions od cycle phase, irregularity, etc
ggplot(data=abcd_ypdms01,aes(x=cycle_phase,color=eventname))+geom_bar(fill=NA)+ 
  theme_minimal() +  labs(x="Cycle phase")
ggplot(data=abcd_ypdms01,aes(x=menstrualcycle3_y,color=eventname))+geom_bar(fill=NA)+
  theme_minimal() +  labs(x="Cycle regularity (0=no, 1=yes,2=don't know, 3=refuse)")
ggplot(data=abcd_ypdms01,aes(x=menstrualcycle2_y,color=eventname))+geom_histogram(fill=NA)+ 
  theme_minimal() +  labs(x="Average length cycle")
ggplot(data=abcd_ypdms01[which(abcd_ypdms01$days_since_period<150),],aes(x=days_since_period,color=eventname))+geom_histogram(fill=NA)+ 
  theme_minimal() 
ggplot(data=abcd_ypdms01[which(abcd_ypdms01$menstrualcycle3_y==0&abcd_ypdms01$days_since_period<150),],aes(x=days_since_period,color=eventname))+geom_histogram(fill=NA)+ 
  theme_minimal()  +  labs(title="Days since period in irregularly cycling girls only")
ggplot(data=abcd_ypdms01[which(abcd_ypdms01$menstrualcycle3_y==1&abcd_ypdms01$days_since_period<150),],aes(x=days_since_period,color=eventname))+geom_histogram(fill=NA)+ 
  theme_minimal()  +  labs(title="Days since period in regularly cycling girls only")
ggplot(data=allmen_hormones_confounds,aes(y=estr,color=as.factor(cycle_phase)))+geom_boxplot(fill=NA)+ 
  theme_minimal() + facet_grid(rows = vars(eventname))

#estradiol by cycle phase split by irregular vs regular
ggplot(data=allmen_hormones_confounds,aes(y=estr,color=as.factor(cycle_phase)))+geom_boxplot(fill=NA)+ 
  theme_minimal() + facet_grid(rows = vars(cycle_regularity),cols = vars(eventname))
estr_cyclephase <- aov(log(estr+1) ~ cycle_phase, data=postmen_hormones_confounds)
summary(estr_cyclephase) #no significant differences but starting to look like the expected pattern in 2y follow up
                       

###############
# Log-transform 
# Using natural log after adding a constant of 1 to prevent infinite values for levels of 0 pg/ml 
###############

filtered_hormones_confounds <- filtered_hormones_confounds %>%
  mutate(dhea_log=log(dhea+1),
         test_log=log(test+1))
allmen_hormones_confounds <- allmen_hormones_confounds %>%
  mutate(estr_log=log(estr+1))


###############
#Regress hormones on confounds and extract residuals
#Confounds included: hours since wake, collection duration, recent caffeine use, recent physical activity, glucocorticoid med use
#not adding time to freeze because it's short (M=2min13s, IQR=1min) and was not related to hormone levels in Herting et al.
###############


#Get clean data file and Get complete data files by hormone
clean_hormones_confounds <- filtered_hormones_confounds %>%
  select(subjectkey,eventname,sex,sex_rec,dhea,test,dhea_log,test_log,hours_since_wake, minutes_to_collect, hormone_sal_caff_y,hormone_sal_active,glucocorticoid, contraceptive)
complete_dhea_confounds <- clean_hormones_confounds %>% select(-test,test_log) %>% na.omit()
complete_test_confounds <- clean_hormones_confounds %>% select(-dhea,-dhea_log) %>% na.omit()
complete_estr_confounds <- allmen_hormones_confounds %>% 
  select(subjectkey,eventname,sex,sex_rec,estr,estr_log,hours_since_wake, minutes_to_collect, hormone_sal_caff_y,hormone_sal_active,
         glucocorticoid, contraceptive,irregular_cycle,perimentrual_phase,midfollicular_phase,latefollicular_phase,
         ovulatory_phase,midluteal_phase,lateluteal_phase) %>% 
  na.omit()

#fractional polynomial models
dhea_mfp <- mfp(dhea_log ~ fp(minutes_to_collect)+fp(hours_since_wake)+hormone_sal_caff_y+hormone_sal_active+glucocorticoid+contraceptive, 
                family = gaussian, data = complete_dhea_confounds, select=0.05)
summary(dhea_mfp)
test_mfp <- mfp(test_log ~ fp(minutes_to_collect)+fp(hours_since_wake)+hormone_sal_caff_y+hormone_sal_active+glucocorticoid+contraceptive, 
                family = gaussian, data = complete_test_confounds, select=0.05)
summary(test_mfp)
estr_mfp <- mfp(estr_log ~ fp(minutes_to_collect)+fp(hours_since_wake)+hormone_sal_caff_y+hormone_sal_active+glucocorticoid+contraceptive+
                  irregular_cycle+perimentrual_phase+midfollicular_phase+latefollicular_phase+ovulatory_phase+midluteal_phase+lateluteal_phase, 
                family = gaussian, data = complete_estr_confounds, select=0.05)
summary(estr_mfp)


#extract residuals
complete_dhea_confounds$dhea_corrected <- resid(dhea_mfp)
complete_test_confounds$test_corrected <- resid(test_mfp)
complete_estr_confounds$estr_corrected <- resid(estr_mfp)


#complete_estr_confounds2 <- merge(complete_estr_confounds,complete_postmen_estr_confounds[c("subjectkey","eventname","estr_corrected2")],
#                                  by=c("subjectkey","eventname"),all=T)
#complete_estr_confounds2 <- complete_estr_confounds2 %>%
#  mutate(estr_corrected=ifelse(is.na(estr_corrected2), estr_corrected,estr_corrected2))

###############
# Save final hormone variables 
###############

hormones_cleaned_2tps <- merge(complete_dhea_confounds[c("subjectkey","eventname","sex","sex_rec","dhea","dhea_log","dhea_corrected")],
                          complete_test_confounds[c("subjectkey","eventname","sex","sex_rec","test","test_log","test_corrected")],
                          by=c("subjectkey","eventname","sex","sex_rec"),all=T) %>%
                    merge(.,complete_estr_confounds[c("subjectkey","eventname","sex","sex_rec","estr","estr_log","estr_corrected")],
                          by=c("subjectkey","eventname","sex","sex_rec"),all=T)
write.csv(hormones_cleaned_2tps,paste0(hormones_dir,"hormones_cleaned_baseline_2year.csv"),row.names = F)

