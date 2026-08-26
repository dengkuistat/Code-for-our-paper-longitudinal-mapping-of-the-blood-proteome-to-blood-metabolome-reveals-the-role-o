#####################################################################################################################################################################################
########## Code used in the paper "Longitudinal mapping of the blood proteome to blood metabolome reveals the role of the protein-metabolite axes in metabolic health" ##############
####### Dataset construction for discovery set
#### Transform the data into wide format
## Proteomics data
library(plyr)
data_temp1 <-  data_protein_phenotype_repeated_match_discovery[,c("id","Follow_up",protein_all[1])]

data_protein_phenotype_repeated_wide_discovery <- reshape(data_temp1, idvar = "id", timevar = "Follow_up", direction = "wide",sep="_")

for (i in 2:length(protein_all)){
  cat("##########",i,"############\n")
  data_temp <-  data_protein_phenotype_repeated_matc_discovery[,c("id","Follow_up",protein_all[i])]
  
  wide_data_temp <- reshape(data_temp, idvar = "id", timevar = "Follow_up", direction = "wide",sep="_")
  
  data_protein_phenotype_repeated_wide_discovery <- merge(data_protein_phenotype_repeated_wide_discovery,wide_data_temp,by="id")
  
}

# Merge with the phenotype data
data_protein_phenotype_repeated_wide_discovery_new <- merge(data_phenotype,data_protein_phenotype_repeated_wide_discovery,by="id")

# Merge with phase infor
data_phase_discovery <- read.dta13(file="Data/protein/New/phase_id.dta",convert.factors=F)

data_phase_discovery <- data_phase_discovery[,c("id","SampleID","phase")]

data_phase_discovery_NL <- data_phase_discovery[str_sub(data_phase_discovery$SampleID,1,2)=="NL",]
data_phase_discovery_NL <- data_phase_discovery_NL[,c("id","phase")]
colnames(data_phase_discovery_NL)[2] <- "phase_NL"

data_phase_discovery_F2 <- data_phase_discovery[str_sub(data_phase_discovery$SampleID,1,2)=="F2",]
data_phase_discovery_F2 <- data_phase_discovery_F2[,c("id","phase")]
colnames(data_phase_discovery_F2)[2] <- "phase_F2"

data_phase_discovery_F3 <- data_phase_discovery[str_sub(data_phase_discovery$SampleID,1,2)=="F3",]
data_phase_discovery_F3 <- data_phase_discovery_F3[,c("id","phase")]
colnames(data_phase_discovery_F3)[2] <- "phase_F3"


data_protein_phenotype_repeated_wide_discovery_new <- merge(data_protein_phenotype_repeated_wide_discovery_new,data_phase_discovery_NL,by="id")
data_protein_phenotype_repeated_wide_discovery_new <- merge(data_protein_phenotype_repeated_wide_discovery_new,data_phase_discovery_F2,by="id")
data_protein_phenotype_repeated_wide_discovery_new <- merge(data_protein_phenotype_repeated_wide_discovery_new,data_phase_discovery_F3,by="id")


rownames(data_protein_phenotype_repeated_wide_discovery_new) <- data_protein_phenotype_repeated_wide_discovery_new$id


## Metabolomics data
Repeated_four_ID <- names(table(data_metabolomics_repeated_new$id)[table(data_metabolomics_repeated_new$id)==4])

data_metabolomics_repeated_four_new <- data_metabolomics_repeated_new[data_metabolomics_repeated_new$id %in% Repeated_four_ID,]

library(plyr)
data_temp1 <-  data_metabolomics_repeated_four_new[,c("id","Follow_up",metabolite_all[1])]

data_metabolomics_repeated_wide <- reshape(data_temp1, idvar = "id", timevar = "Follow_up", direction = "wide",sep="_")

for (i in 2:length(metabolite_all)){
  cat("##########",i,"############\n")
  data_temp <-  data_metabolomics_repeated_four_new[,c("id","Follow_up",metabolite_all[i])]
  
  wide_data_temp <- reshape(data_temp, idvar = "id", timevar = "Follow_up", direction = "wide",sep="_")
  
  data_metabolomics_repeated_wide <- merge(data_metabolomics_repeated_wide,wide_data_temp,by="id")
  
}

rownames(data_metabolomics_repeated_wide) <- data_metabolomics_repeated_wide$id

overlap_id_discovery <- intersect(data_protein_phenotype_repeated_wide_discovery_new$id,data_metabolomics_repeated_wide$id)

data_protein_phenotype_repeated_wide_discovery_new <- data_protein_phenotype_repeated_wide_discovery_new[overlap_id_discovery,]

data_metabolomics_repeated_wide_discovery <- data_metabolomics_repeated_wide[overlap_id_discovery,]

# Calculating Follow-up time
data_protein_phenotype_repeated_wide_discovery_new$visitdate_v1 <- as.Date(data_protein_phenotype_repeated_wide_discovery_new$visitdate_v1,format="%Y/%m/%d")
data_protein_phenotype_repeated_wide_discovery_new$visitdate_v2 <- as.Date(data_protein_phenotype_repeated_wide_discovery_new$visitdate_v2,format="%Y/%m/%d")
data_protein_phenotype_repeated_wide_discovery_new$visitdate_v3 <- as.Date(data_protein_phenotype_repeated_wide_discovery_new$visitdate_v3,format="%Y/%m/%d")
data_protein_phenotype_repeated_wide_discovery_new$visitdate_v4 <- as.Date(data_protein_phenotype_repeated_wide_discovery_new$visitdate_v4,format="%Y-%m-%d")


data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v2_v1 <- difftime(data_protein_phenotype_repeated_wide_discovery_new$visitdate_v2,data_protein_phenotype_repeated_wide_discovery_new$visitdate_v1,units="days")/365
data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v2_v1 <- as.numeric(data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v2_v1)

data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v3_v2 <- difftime(data_protein_phenotype_repeated_wide_discovery_new$visitdate_v3,data_protein_phenotype_repeated_wide_discovery_new$visitdate_v2,units="days")/365
data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v3_v2 <- as.numeric(data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v3_v2)

data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v4_v3 <- difftime(data_protein_phenotype_repeated_wide_discovery_new$visitdate_v4,data_protein_phenotype_repeated_wide_discovery_new$visitdate_v3,units="days")/365
data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v4_v3 <- as.numeric(data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v4_v3)


## Generating new datasets
# Baseline
library(stringr)
data_protein_discovery_temp_NL <- data_protein_phenotype_repeated_wide_discovery_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_discovery_new),"_NL")]
data_metabolite_discovery_temp_NL <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_NL")]
data_metabolite_discovery_temp_F2 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F2")]

colnames(data_protein_discovery_temp_NL) <- str_replace_all(colnames(data_protein_discovery_temp_NL),"_NL","_base")
colnames(data_metabolite_discovery_temp_NL) <- str_replace_all(colnames(data_metabolite_discovery_temp_NL),"_NL","_base")
colnames(data_metabolite_discovery_temp_F2) <- str_replace_all(colnames(data_metabolite_discovery_temp_F2),"_F2","_follow")

data_protein_metabolite_discovery_NL_F2 <- data.frame(SampleID=data_protein_phenotype_repeated_wide_discovery_new$id,id=data_protein_phenotype_repeated_wide_discovery_new$id,
                                                      age=data_protein_phenotype_repeated_wide_discovery_new$age,sex=data_protein_phenotype_repeated_wide_discovery_new$sex,
                                                      follow_up_time = data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v2_v1,
                                                      T2D_drug=data_protein_phenotype_repeated_wide_discovery_new$T2D_med_new_v1,
                                                      hypertension_drug=data_protein_phenotype_repeated_wide_discovery_new$hypertension_med_new_v1,
                                                      dyslipi_drug=data_protein_phenotype_repeated_wide_discovery_new$dyslipi_med_new_v1,
                                                      data_protein_discovery_temp_NL,data_metabolite_discovery_temp_NL,data_metabolite_discovery_temp_F2,check.names=F,stringsAsFactors = F)


# F2
data_protein_discovery_temp_F2 <- data_protein_phenotype_repeated_wide_discovery_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_discovery_new),"_F2")]
data_metabolite_discovery_temp_F2 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F2")]
data_metabolite_discovery_temp_F3 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F3")]

colnames(data_protein_discovery_temp_F2) <- str_replace_all(colnames(data_protein_discovery_temp_F2),"_F2","_base")
colnames(data_metabolite_discovery_temp_F2) <- str_replace_all(colnames(data_metabolite_discovery_temp_F2),"_F2","_base")
colnames(data_metabolite_discovery_temp_F3) <- str_replace_all(colnames(data_metabolite_discovery_temp_F3),"_F3","_follow")

data_protein_metabolite_discovery_F2_F3 <- data.frame(SampleID=paste("F2",data_protein_phenotype_repeated_wide_discovery_new$id,sep=""),id=data_protein_phenotype_repeated_wide_discovery_new$id,
                                                      age=data_protein_phenotype_repeated_wide_discovery_new$age_v2,sex=data_protein_phenotype_repeated_wide_discovery_new$sex,
                                                      follow_up_time = data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v3_v2,
                                                      T2D_drug=data_protein_phenotype_repeated_wide_discovery_new$T2D_med_new_v2,
                                                      hypertension_drug=data_protein_phenotype_repeated_wide_discovery_new$hypertension_med_new_v2,
                                                      dyslipi_drug=data_protein_phenotype_repeated_wide_discovery_new$dyslipi_med_new_v2,
                                                      data_protein_discovery_temp_F2,data_metabolite_discovery_temp_F2,data_metabolite_discovery_temp_F3,check.names=F,stringsAsFactors = F)


# F3
data_protein_discovery_temp_F3 <- data_protein_phenotype_repeated_wide_discovery_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_discovery_new),"_F3")]
data_metabolite_discovery_temp_F3 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F3")]
data_metabolite_discovery_temp_F4 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F4")]

colnames(data_protein_discovery_temp_F3) <- str_replace_all(colnames(data_protein_discovery_temp_F3),"_F3","_base")
colnames(data_metabolite_discovery_temp_F3) <- str_replace_all(colnames(data_metabolite_discovery_temp_F3),"_F3","_base")
colnames(data_metabolite_discovery_temp_F4) <- str_replace_all(colnames(data_metabolite_discovery_temp_F4),"_F4","_follow")

data_protein_metabolite_discovery_F3_F4 <- data.frame(SampleID=paste("F3",data_protein_phenotype_repeated_wide_discovery_new$id,sep=""),id=data_protein_phenotype_repeated_wide_discovery_new$id,
                                                      age=data_protein_phenotype_repeated_wide_discovery_new$age_v3,sex=data_protein_phenotype_repeated_wide_discovery_new$sex,
                                                      follow_up_time = data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v4_v3,
                                                      T2D_drug=data_protein_phenotype_repeated_wide_discovery_new$T2D_med_new_v3,
                                                      hypertension_drug=data_protein_phenotype_repeated_wide_discovery_new$hypertension_med_new_v3,
                                                      dyslipi_drug=data_protein_phenotype_repeated_wide_discovery_new$dyslipi_med_new_v3,
                                                      data_protein_discovery_temp_F3,data_metabolite_discovery_temp_F3,data_metabolite_discovery_temp_F4,check.names=F,stringsAsFactors = F)


data_discovery_new <- rbind(data_protein_metabolite_discovery_NL_F2,data_protein_metabolite_discovery_F2_F3,data_protein_metabolite_discovery_F3_F4)


####### Dataset construction for validation set
#### Transform the data into wide format
## Proteomics data
library(plyr)
data_temp1 <-  data_protein_phenotype_repeated_match_validation[,c("id","Follow_up",protein_all[1])]

data_protein_phenotype_repeated_wide_validation <- reshape(data_temp1, idvar = "id", timevar = "Follow_up", direction = "wide",sep="_")

for (i in 2:length(protein_all)){
  cat("##########",i,"############\n")
  data_temp <-  data_protein_phenotype_repeated_match_validation[,c("id","Follow_up",protein_all[i])]
  
  wide_data_temp <- reshape(data_temp, idvar = "id", timevar = "Follow_up", direction = "wide",sep="_")
  
  data_protein_phenotype_repeated_wide_validation <- merge(data_protein_phenotype_repeated_wide_validation,wide_data_temp,by="id")
  
}

## Merge with the phenotype data
data_protein_phenotype_repeated_wide_validation_new <- merge(data_phenotype,data_protein_phenotype_repeated_wide_validation,by="id")
rownames(data_protein_phenotype_repeated_wide_validation_new) <- data_protein_phenotype_repeated_wide_validation_new$id


overlap_id_validation <- intersect(data_protein_phenotype_repeated_wide_validation_new$id,data_metabolomics_repeated_wide$id)

data_protein4_phenotype_repeated_wide_new <- data_protein_phenotype_repeated_wide_validation_new[overlap_id_validation,]

data_metabolomics_repeated_wide_validation <- data_metabolomics_repeated_wide[overlap_id_validation,]


# Calculating Follow-up time
data_protein_phenotype_repeated_wide_validation_new$visitdate_v1 <- as.Date(data_protein_phenotype_repeated_wide_validation_new$visitdate_v1,format="%Y/%m/%d")
data_protein_phenotype_repeated_wide_validation_new$visitdate_v2 <- as.Date(data_protein_phenotype_repeated_wide_validation_new$visitdate_v2,format="%Y/%m/%d")
data_protein_phenotype_repeated_wide_validation_new$visitdate_v3 <- as.Date(data_protein_phenotype_repeated_wide_validation_new$visitdate_v3,format="%Y/%m/%d")
data_protein_phenotype_repeated_wide_validation_new$visitdate_v4 <- as.Date(data_protein_phenotype_repeated_wide_validation_new$visitdate_v4,format="%Y-%m-%d")


data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v2_v1 <- difftime(data_protein_phenotype_repeated_wide_validation_new$visitdate_v2,data_protein_phenotype_repeated_wide_validation_new$visitdate_v1,units="days")/365
data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v2_v1 <- as.numeric(data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v2_v1)

data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v3_v2 <- difftime(data_protein_phenotype_repeated_wide_validation_new$visitdate_v3,data_protein_phenotype_repeated_wide_validation_new$visitdate_v2,units="days")/365
data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v3_v2 <- as.numeric(data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v3_v2)

data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v4_v3 <- difftime(data_protein_phenotype_repeated_wide_validation_new$visitdate_v4,data_protein_phenotype_repeated_wide_validation_new$visitdate_v3,units="days")/365
data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v4_v3 <- as.numeric(data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v4_v3)


## Generating new datasets
# Baseline
library(stringr)
data_protein_validation_temp_NL <- data_protein_phenotype_repeated_wide_validation_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_validation_new),"_NL")]
data_metabolite_validation_temp_NL <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_NL")]
data_metabolite_validation_temp_F2 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F2")]

colnames(data_protein_validation_temp_NL) <- str_replace_all(colnames(data_protein_validation_temp_NL),"_NL","_base")
colnames(data_metabolite_validation_temp_NL) <- str_replace_all(colnames(data_metabolite_validation_temp_NL),"_NL","_base")
colnames(data_metabolite_validation_temp_F2) <- str_replace_all(colnames(data_metabolite_validation_temp_F2),"_F2","_follow")

data_protein_metabolite_validation_NL_F2 <- data.frame(SampleID=data_protein_phenotype_repeated_wide_validation_new$id,id=data_protein_phenotype_repeated_wide_validation_new$id,
                                                       age=data_protein_phenotype_repeated_wide_validation_new$age,sex=data_protein_phenotype_repeated_wide_validation_new$sex,
                                                       follow_up_time = data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v2_v1,
                                                       T2D_drug=data_protein_phenotype_repeated_wide_validation_new$T2D_med_new_v1,
                                                       hypertension_drug=data_protein_phenotype_repeated_wide_validation_new$hypertension_med_new_v1,
                                                       dyslipi_drug=data_protein_phenotype_repeated_wide_validation_new$dyslipi_med_new_v1,
                                                       data_protein_validation_temp_NL,data_metabolite_validation_temp_NL,data_metabolite_validation_temp_F2,check.names=F,stringsAsFactors = F)


# F2
data_protein_validation_temp_F2 <- data_protein_phenotype_repeated_wide_validation_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_validation_new),"_F2")]
data_metabolite_validation_temp_F2 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F2")]
data_metabolite_validation_temp_F3 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F3")]

colnames(data_protein_validation_temp_F2) <- str_replace_all(colnames(data_protein_validation_temp_F2),"_F2","_base")
colnames(data_metabolite_validation_temp_F2) <- str_replace_all(colnames(data_metabolite_validation_temp_F2),"_F2","_base")
colnames(data_metabolite_validation_temp_F3) <- str_replace_all(colnames(data_metabolite_validation_temp_F3),"_F3","_follow")

data_protein_metabolite_validation_F2_F3 <- data.frame(SampleID=paste("F2",data_protein_phenotype_repeated_wide_validation_new$id,sep=""),id=data_protein_phenotype_repeated_wide_validation_new$id,
                                                       age=data_protein_phenotype_repeated_wide_validation_new$age_v2,sex=data_protein_phenotype_repeated_wide_validation_new$sex,
                                                       follow_up_time = data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v3_v2,
                                                       T2D_drug=data_protein_phenotype_repeated_wide_validation_new$T2D_med_new_v2,
                                                       hypertension_drug=data_protein_phenotype_repeated_wide_validation_new$hypertension_med_new_v2,
                                                       dyslipi_drug=data_protein_phenotype_repeated_wide_validation_new$dyslipi_med_new_v2,
                                                       data_protein_validation_temp_F2,data_metabolite_validation_temp_F2,data_metabolite_validation_temp_F3,check.names=F,stringsAsFactors = F)


# F3
data_protein_validation_temp_F3 <- data_protein_phenotype_repeated_wide_validation_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_validation_new),"_F3")]
data_metabolite_validation_temp_F3 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F3")]
data_metabolite_validation_temp_F4 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F4")]

colnames(data_protein_validation_temp_F3) <- str_replace_all(colnames(data_protein_validation_temp_F3),"_F3","_base")
colnames(data_metabolite_validation_temp_F3) <- str_replace_all(colnames(data_metabolite_validation_temp_F3),"_F3","_base")
colnames(data_metabolite_validation_temp_F4) <- str_replace_all(colnames(data_metabolite_validation_temp_F4),"_F4","_follow")

data_protein_metabolite_validation_F3_F4 <- data.frame(SampleID=paste("F3",data_protein_phenotype_repeated_wide_validation_new$id,sep=""),id=data_protein_phenotype_repeated_wide_validation_new$id,
                                                       age=data_protein_phenotype_repeated_wide_validation_new$age_v3,sex=data_protein_phenotype_repeated_wide_validation_new$sex,
                                                       follow_up_time = data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v4_v3,
                                                       T2D_drug=data_protein_phenotype_repeated_wide_validation_new$T2D_med_new_v3,
                                                       hypertension_drug=data_protein_phenotype_repeated_wide_validation_new$hypertension_med_new_v3,
                                                       dyslipi_drug=data_protein_phenotype_repeated_wide_validation_new$dyslipi_med_new_v3,
                                                       data_protein_validation_temp_F3,data_metabolite_validation_temp_F3,data_metabolite_validation_temp_F4,check.names=F,stringsAsFactors = F)


data_validation_new <- rbind(data_protein_metabolite_validation_NL_F2,data_protein_metabolite_validation_F2_F3,data_protein_metabolite_validation_F3_F4)


###### Longitudinal prospective associations between serum proteins and metabolites
## Discovery set 
protein_metabolite_prospective_discovery <- c()
for (i in 1:length(protein_all)){
  for (j in 1:length(metabolite_all)){
    cat("###############",i,"protein","|",j,"metabolite","#####################\n")
    
    protein_temp <- protein_all[i]
    metabolite_temp <- metabolite_all[j]
    
    protein_temp_base <- paste(protein_temp,"_base",sep="")
    metabolite_temp_base <- paste(metabolite_temp,"_base",sep="")
    metabolite_temp_follow <- paste(metabolite_temp,"_follow",sep="")
    
    data_fit_temp <- data.frame(protein_base=data_discovery_new[,protein_temp_base],
                                metabolite_base=data_discovery_new[,metabolite_temp_base],
                                metabolite_follow=data_discovery_new[,metabolite_temp_follow],
                                data_discovery_new,stringsAsFactors = F)
    
    library(lme4)
    library(lmerTest)
    mix_model_temp <- lmer(metabolite_follow~protein_base+metabolite_base+age+sex+follow_up_time+factor(phase)+(1|id),data=data_fit_temp)
    mix_model_temp_summary <- summary(mix_model_temp)
    coef_mixed_model <- mix_model_temp_summary$coefficients[2,]
    var_name_temp <- c(protein=protein_temp,metabolite=metabolite_temp)
    coef_mixed_model <- c(var_name_temp,coef_mixed_model)
    
    protein_metabolite_prospective_discovery <- rbind(protein_metabolite_prospective_discovery,coef_mixed_model)
    
  }
}

protein_metabolite_prospective_discovery <- as.data.frame(protein_metabolite_prospective_discovery)
protein_metabolite_prospective_discovery[,-c(1:2)] <- apply(protein_metabolite_prospective_discovery[,-c(1:2)],2,as.numeric)
protein_metabolite_prospective_discovery$id <- paste(protein_metabolite_prospective_discovery$protein," & ",protein_metabolite_prospective_discovery$metabolite,sep="")


## validation set 
protein_metabolite_prospective_validation <- c()
for (i in 1:length(protein_all)){
  for (j in 1:length(metabolite_all)){
    cat("###############",i,"protein","|",j,"metabolite","#####################\n")
    
    protein_temp <- protein_all[i]
    metabolite_temp <- metabolite_all[j]
    
    protein_temp_base <- paste(protein_temp,"_base",sep="")
    metabolite_temp_base <- paste(metabolite_temp,"_base",sep="")
    metabolite_temp_follow <- paste(metabolite_temp,"_follow",sep="")
    
    data_fit_temp <- data.frame(protein_base=data_validation_new[,protein_temp_base],
                                metabolite_base=data_validation_new[,metabolite_temp_base],
                                metabolite_follow=data_validation_new[,metabolite_temp_follow],
                                data_validation_new,stringsAsFactors = F)
    
    library(lme4)
    library(lmerTest)
    mix_model_temp <- lmer(metabolite_follow~protein_base+metabolite_base+age+sex+follow_up_time+(1|id),data=data_fit_temp)
    mix_model_temp_summary <- summary(mix_model_temp)
    coef_mixed_model <- mix_model_temp_summary$coefficients[2,]
    var_name_temp <- c(protein=protein_temp,metabolite=metabolite_temp)
    coef_mixed_model <- c(var_name_temp,coef_mixed_model)
    
    protein_metabolite_prospective_validation <- rbind(protein_metabolite_prospective_validation,coef_mixed_model)
    
  }
}

protein_metabolite_prospective_validation <- as.data.frame(protein_metabolite_prospective_validation)
protein_metabolite_prospective_validation[,-c(1:2)] <- apply(protein_metabolite_prospective_validation[,-c(1:2)],2,as.numeric)
protein_metabolite_prospective_validation$id <- paste(protein_metabolite_prospective_validation$protein," & ",protein_metabolite_prospective_validation$metabolite,sep="")


## Meta analysis 
data_discovery <- protein_metabolite_prospective_discovery
data_validation <- protein_metabolite_prospective_validation

library(metafor)
protein_metabolite_prospective_meta <- c()
for (i in 1:nrow(data_discovery)){
  cat("#########",i,"###########\n")
  data_discovery_temp <- data_discovery[i,c("Estimate","Std. Error")]
  data_validation_temp  <- data_validation[i,c("Estimate","Std. Error")]
  
  data_test_temp <- rbind(data_discovery_temp,data_validation_temp)
  meta_result_temp <- rma(yi=Estimate,sei=`Std. Error`,data=data_test_temp,method="FE")
  result_temp <- data.frame(id=data_discovery$id[i],protein=data_discovery$protein[i],metabolite=data_discovery$metabolite[i],
                            B_discovery=data_discovery$Estimate[i],SE_discovery=data_discovery$`Std. Error`[i],P_discovery = data_discovery$`Pr(>|t|)`[i],
                            B_validation=data_validation$Estimate[i],SE_validation=data_validation$`Std. Error`[i],P_validation = data_validation$`Pr(>|t|)`[i],
                            heterogeneity=meta_result_temp$QEp,B_meta=meta_result_temp$beta,SE_meta=meta_result_temp$se,CI_lower_meta=meta_result_temp$ci.lb,CI_upper_meta=meta_result_temp$ci.ub,p_meta=meta_result_temp$pval,stringsAsFactors=F
  )
  protein_metabolite_prospective_meta <- rbind(protein_metabolite_prospective_meta,result_temp)
}

protein_metabolite_prospective_meta$FDR_meta <- p.adjust(protein_metabolite_prospective_meta$p_meta,method="fdr")

protein_metabolite_prospective_meta_select <- protein_metabolite_prospective_meta[protein_metabolite_prospective_meta$P_discovery < 0.05 &
                                                                                    protein_metabolite_prospective_meta$P_validation < 0.05 & 
                                                                                    protein_metabolite_prospective_meta$FDR_meta < 0.05 &
                                                                                    protein_metabolite_prospective_meta$heterogeneity > 0.05,]


##### Interaction analysis exploring whether the protein-metabolite associations change over time
## Discovery set 
data_discovery_new$interval_time[str_sub(data_discovery_new$SampleID,1,2)=="NL"] <- 0
data_discovery_new$interval_time[str_sub(data_discovery_new$SampleID,1,2)=="F2"] <- data_discovery_new$follow_up_time[str_sub(data_discovery_new$SampleID,1,2)=="NL"]
data_discovery_new$interval_time[str_sub(data_discovery_new$SampleID,1,2)=="F3"] <- data_discovery_new$follow_up_time[str_sub(data_discovery_new$SampleID,1,2)=="NL"] + data_discovery_new$follow_up_time[str_sub(data_discovery_new$SampleID,1,2)=="F2"]

Interaction_test_discovery <- c()
for (i in 1:nrow(protein_metabolite_prospective_meta_select)){
  cat("###############",i,"#####################\n")
  
  protein_temp <- protein_metabolite_prospective_meta_select$protein[i]
  metabolite_temp <- protein_metabolite_prospective_meta_select$metabolite[i]
  
  protein_temp_base <- paste(protein_temp,"_base",sep="")
  metabolite_temp_base <- paste(metabolite_temp,"_base",sep="")
  metabolite_temp_follow <- paste(metabolite_temp,"_follow",sep="")
  
  data_fit_temp <- data.frame(protein_base=data_discovery_new[,protein_temp_base],
                              metabolite_base=data_discovery_new[,metabolite_temp_base],
                              metabolite_follow=data_discovery_new[,metabolite_temp_follow],
                              data_discovery_new,stringsAsFactors = F)
  
  
  library(lme4)
  library(lmerTest)
  mix_model_interaction_temp <- lmer(metabolite_follow~protein_base*interval_time+metabolite_base+age+sex+follow_up_time+factor(phase)+(1|id),data=data_fit_temp)
  
  mix_model_interaction_temp_summary <- summary(mix_model_interaction_temp)
  Interaction_P <- mix_model_interaction_temp_summary$coefficients[nrow(mix_model_interaction_temp_summary$coefficients),]
  Interaction_test_discovery <- rbind(Interaction_test_discovery,Interaction_P)
  
}

Interaction_test_discovery <- as.data.frame(Interaction_test_discovery)
rownames(Interaction_test_discovery) <- protein_metabolite_prospective_meta_select$id


## Validation set 
data_validation_new$interval_time[str_sub(data_validation_new$SampleID,1,2)=="NL"] <- 0
data_validation_new$interval_time[str_sub(data_validation_new$SampleID,1,2)=="F2"] <- data_validation_new$follow_up_time[str_sub(data_validation_new$SampleID,1,2)=="NL"]
data_validation_new$interval_time[str_sub(data_validation_new$SampleID,1,2)=="F3"] <- data_validation_new$follow_up_time[str_sub(data_validation_new$SampleID,1,2)=="NL"] + data_validation_new$follow_up_time[str_sub(data_validation_new$SampleID,1,2)=="F2"]

Interaction_test_validation <- c()
for (i in 1:nrow(protein_metabolite_prospective_meta_select)){
  cat("###############",i,"#####################\n")
  
  protein_temp <- protein_metabolite_prospective_meta_select$protein[i]
  metabolite_temp <- protein_metabolite_prospective_meta_select$metabolite[i]
  
  protein_temp_base <- paste(protein_temp,"_base",sep="")
  metabolite_temp_base <- paste(metabolite_temp,"_base",sep="")
  metabolite_temp_follow <- paste(metabolite_temp,"_follow",sep="")
  
  data_fit_temp <- data.frame(protein_base=data_validation_new[,protein_temp_base],
                              metabolite_base=data_validation_new[,metabolite_temp_base],
                              metabolite_follow=data_validation_new[,metabolite_temp_follow],
                              data_validation_new,stringsAsFactors = F)
  
  
  library(lme4)
  library(lmerTest)
  mix_model_interaction_temp <- lmer(metabolite_follow~protein_base*interval_time+metabolite_base+age+sex+follow_up_time+(1|id),data=data_fit_temp)
  
  mix_model_interaction_temp_summary <- summary(mix_model_interaction_temp)
  Interaction_P <- mix_model_interaction_temp_summary$coefficients[nrow(mix_model_interaction_temp_summary$coefficients),]
  Interaction_test_validation <- rbind(Interaction_test_validation,Interaction_P)
  
}

Interaction_test_validation <- as.data.frame(Interaction_test_validation)
rownames(Interaction_test_validation) <- protein_metabolite_prospective_meta_select$id


############ Meta analysis 
data_discovery <- Interaction_test_discovery
data_validation <- Interaction_test_validation

library(metafor)
protein_metabolite_interaction_meta <- c()
for (i in 1:nrow(data_discovery)){
  cat("#########",i,"###########\n")
  data_discovery_temp <- data_discovery[i,c("Estimate","Std. Error")]
  data_validation_temp  <- data_validation[i,c("Estimate","Std. Error")]
  
  data_test_temp <- rbind(data_discovery_temp,data_validation_temp)
  meta_result_temp <- rma(yi=Estimate,sei=`Std. Error`,data=data_test_temp,method="FE")
  result_temp <- data.frame(id=rownames(data_discovery)[i],
                            B_discovery=data_discovery$Estimate[i],SE_discovery=data_discovery$`Std. Error`[i],P_discovery = data_discovery$`Pr(>|t|)`[i],
                            B_validation=data_validation$Estimate[i],SE_validation=data_validation$`Std. Error`[i],P_validation = data_validation$`Pr(>|t|)`[i],
                            heterogeneity=meta_result_temp$QEp,B_meta=meta_result_temp$beta,SE_meta=meta_result_temp$se,CI_lower_meta=meta_result_temp$ci.lb,CI_upper_meta=meta_result_temp$ci.ub,p_meta=meta_result_temp$pval,stringsAsFactors=F
  )
  protein_metabolite_interaction_meta <- rbind(protein_metabolite_interaction_meta,result_temp)
}

protein_metabolite_interaction_meta$FDR_meta <- p.adjust(protein_metabolite_interaction_meta$p_meta,method="fdr")

protein_metabolite_interaction_meta_select <- protein_metabolite_interaction_meta[protein_metabolite_interaction_meta$P_discovery < 0.05 &
                                                                                    protein_metabolite_interaction_meta$P_validation < 0.05 & 
                                                                                    protein_metabolite_interaction_meta$FDR_meta < 0.05 &
                                                                                    protein_metabolite_interaction_meta$heterogeneity > 0.05,]


##### Sensitivity analysis by further adjusting for medications for T2D, hypertension, and hyperlipidemia
## Discovery set 
protein_metabolite_prospective_discovery_sensitivity <- c()
for (i in 1:length(protein_all)){
  for (j in 1:length(metabolite_all)){
    cat("###############",i,"protein","|",j,"metabolite","#####################\n")
    
    protein_temp <- protein_all[i]
    metabolite_temp <- metabolite_all[j]
    
    protein_temp_base <- paste(protein_temp,"_base",sep="")
    metabolite_temp_base <- paste(metabolite_temp,"_base",sep="")
    metabolite_temp_follow <- paste(metabolite_temp,"_follow",sep="")
    
    data_fit_temp <- data.frame(protein_base=data_discovery_new[,protein_temp_base],
                                metabolite_base=data_discovery_new[,metabolite_temp_base],
                                metabolite_follow=data_discovery_new[,metabolite_temp_follow],
                                data_discovery_new,stringsAsFactors = F)
    
    library(lme4)
    library(lmerTest)
    mix_model_temp <- lmer(metabolite_follow~protein_base+metabolite_base+age+sex+follow_up_time+factor(phase)+T2D_drug+hypertension_drug+dyslipi_drug+(1|id),data=data_fit_temp)
    mix_model_temp_summary <- summary(mix_model_temp)
    coef_mixed_model <- mix_model_temp_summary$coefficients[2,]
    var_name_temp <- c(protein=protein_temp,metabolite=metabolite_temp)
    coef_mixed_model <- c(var_name_temp,coef_mixed_model)
    
    protein_metabolite_prospective_discovery_sensitivity <- rbind(protein_metabolite_prospective_discovery_sensitivity,coef_mixed_model)
    
  }
}

protein_metabolite_prospective_discovery_sensitivity <- as.data.frame(protein_metabolite_prospective_discovery_sensitivity)

protein_metabolite_prospective_discovery_sensitivity[,-c(1:2)] <- apply(protein_metabolite_prospective_discovery_sensitivity[,-c(1:2)],2,as.numeric)


## validation set 
protein_metabolite_prospective_validation_sensitivity <- c()
for (i in 1:length(protein_all)){
  for (j in 1:length(metabolite_all)){
    cat("###############",i,"protein","|",j,"metabolite","#####################\n")
    
    protein_temp <- protein_all[i]
    metabolite_temp <- metabolite_all[j]
    
    protein_temp_base <- paste(protein_temp,"_base",sep="")
    metabolite_temp_base <- paste(metabolite_temp,"_base",sep="")
    metabolite_temp_follow <- paste(metabolite_temp,"_follow",sep="")
    
    data_fit_temp <- data.frame(protein_base=data_validation_new[,protein_temp_base],
                                metabolite_base=data_validation_new[,metabolite_temp_base],
                                metabolite_follow=data_validation_new[,metabolite_temp_follow],
                                data_validation_new,stringsAsFactors = F)
    
    library(lme4)
    library(lmerTest)
    mix_model_temp <- lmer(metabolite_follow~protein_base+metabolite_base+age+sex+follow_up_time+T2D_drug+hypertension_drug+dyslipi_drug+(1|id),data=data_fit_temp)
    mix_model_temp_summary <- summary(mix_model_temp)
    coef_mixed_model <- mix_model_temp_summary$coefficients[2,]
    var_name_temp <- c(protein=protein_temp,metabolite=metabolite_temp)
    coef_mixed_model <- c(var_name_temp,coef_mixed_model)
    
    protein_metabolite_prospective_validation_sensitivity <- rbind(protein_metabolite_prospective_validation_sensitivity,coef_mixed_model)
    
  }
}

protein_metabolite_prospective_validation_sensitivity <- as.data.frame(protein_metabolite_prospective_validation_sensitivity)

protein_metabolite_prospective_validation_sensitivity[,-c(1:2)] <- apply(protein_metabolite_prospective_validation_sensitivity[,-c(1:2)],2,as.numeric)


###################################################################################################
###### Longitudinal cross-sectional associations between serum proteins and metabolites
### Dataset construction for discovery set
## NL
data_protein_discovery_temp_NL <- data_protein_phenotype_repeated_wide_discovery_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_discovery_new),"_NL")]
data_metabolite_discovery_temp_NL <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_NL")]

colnames(data_protein_discovery_temp_NL) <- str_replace_all(colnames(data_protein_discovery_temp_NL),"_NL","_base")
colnames(data_metabolite_discovery_temp_NL) <- str_replace_all(colnames(data_metabolite_discovery_temp_NL),"_NL","_base")


data_protein_metabolite_discovery_NL <- data.frame(SampleID=data_protein_phenotype_repeated_wide_discovery_new$id,id=data_protein_phenotype_repeated_wide_discovery_new$id,
                                                        age=data_protein_phenotype_repeated_wide_discovery_new$age,sex=data_protein_phenotype_repeated_wide_discovery_new$sex,
                                                        T2D_drug=data_protein_phenotype_repeated_wide_discovery_new$T2D_med_new_v1,
                                                        hypertension_drug=data_protein_phenotype_repeated_wide_discovery_new$hypertension_med_new_v1,
                                                        dyslipi_drug=data_protein_phenotype_repeated_wide_discovery_new$dyslipi_med_new_v1,
                                                        data_protein_discovery_temp_NL,data_metabolite_discovery_temp_NL,check.names=F,stringsAsFactors = F)


## F2
data_protein_discovery_temp_F2 <- data_protein_phenotype_repeated_wide_discovery_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_discovery_new),"_F2")]
data_metabolite_discovery_temp_F2 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F2")]

colnames(data_protein_discovery_temp_F2) <- str_replace_all(colnames(data_protein_discovery_temp_F2),"_F2","_base")
colnames(data_metabolite_discovery_temp_F2) <- str_replace_all(colnames(data_metabolite_discovery_temp_F2),"_F2","_base")

data_protein_metabolite_discovery_F2 <- data.frame(SampleID=paste("F2",data_protein_phenotype_repeated_wide_discovery_new$id,sep=""),id=data_protein_phenotype_repeated_wide_discovery_new$id,
                                                        age=data_protein_phenotype_repeated_wide_discovery_new$age_v2,sex=data_protein_phenotype_repeated_wide_discovery_new$sex,
                                                        T2D_drug=data_protein_phenotype_repeated_wide_discovery_new$T2D_med_new_v2,
                                                        hypertension_drug=data_protein_phenotype_repeated_wide_discovery_new$hypertension_med_new_v2,
                                                        dyslipi_drug=data_protein_phenotype_repeated_wide_discovery_new$dyslipi_med_new_v2,
                                                        data_protein_discovery_temp_F2,data_metabolite_discovery_temp_F2,check.names=F,stringsAsFactors = F)


## F3
data_protein_discovery_temp_F3 <- data_protein_phenotype_repeated_wide_discovery_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_discovery_new),"_F3")]
data_metabolite_discovery_temp_F3 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F3")]

colnames(data_protein_discovery_temp_F3) <- str_replace_all(colnames(data_protein_discovery_temp_F3),"_F3","_base")
colnames(data_metabolite_discovery_temp_F3) <- str_replace_all(colnames(data_metabolite_discovery_temp_F3),"_F3","_base")

data_protein_metabolite_discovery_F3 <- data.frame(SampleID=paste("F3",data_protein_phenotype_repeated_wide_discovery_new$id,sep=""),id=data_protein_phenotype_repeated_wide_discovery_new$id,
                                                        age=data_protein_phenotype_repeated_wide_discovery_new$age_v3.,sex=data_protein_phenotype_repeated_wide_discovery_new$sex,
                                                        T2D_drug=data_protein_phenotype_repeated_wide_discovery_new$T2D_med_new_v3,
                                                        hypertension_drug=data_protein_phenotype_repeated_wide_discovery_new$hypertension_med_new_v3,
                                                        dyslipi_drug=data_protein_phenotype_repeated_wide_discovery_new$dyslipi_med_new_v3,
                                                        data_protein_discovery_temp_F3,data_metabolite_discovery_temp_F3,check.names=F,stringsAsFactors = F)


data_discovery_cross_sectional_new <- rbind(data_protein_metabolite_discovery_NL,data_protein_metabolite_discovery_F2,data_protein_metabolite_discovery_F3)

### Dataset construction for validation set
## NL
data_protein_validation_temp_NL <- data_protein_phenotype_repeated_wide_validation_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_validation_new),"_NL")]
data_metabolite_validation_temp_NL <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_NL")]

colnames(data_protein_validation_temp_NL) <- str_replace_all(colnames(data_protein_validation_temp_NL),"_NL","_base")
colnames(data_metabolite_validation_temp_NL) <- str_replace_all(colnames(data_metabolite_validation_temp_NL),"_NL","_base")


data_protein_metabolite_validation_NL <- data.frame(SampleID=data_protein_phenotype_repeated_wide_validation_new$id,id=data_protein_phenotype_repeated_wide_validation_new$id,
                                                    age=data_protein_phenotype_repeated_wide_validation_new$age,sex=data_protein_phenotype_repeated_wide_validation_new$sex,
                                                    T2D_drug=data_protein_phenotype_repeated_wide_validation_new$T2D_med_new_v1,
                                                    hypertension_drug=data_protein_phenotype_repeated_wide_validation_new$hypertension_med_new_v1,
                                                    dyslipi_drug=data_protein_phenotype_repeated_wide_validation_new$dyslipi_med_new_v1,
                                                    data_protein_validation_temp_NL,data_metabolite_validation_temp_NL,check.names=F,stringsAsFactors = F)


## F2
data_protein_validation_temp_F2 <- data_protein_phenotype_repeated_wide_validation_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_validation_new),"_F2")]
data_metabolite_validation_temp_F2 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F2")]

colnames(data_protein_validation_temp_F2) <- str_replace_all(colnames(data_protein_validation_temp_F2),"_F2","_base")
colnames(data_metabolite_validation_temp_F2) <- str_replace_all(colnames(data_metabolite_validation_temp_F2),"_F2","_base")

data_protein_metabolite_validation_F2 <- data.frame(SampleID=paste("F2",data_protein_phenotype_repeated_wide_validation_new$id,sep=""),id=data_protein_phenotype_repeated_wide_validation_new$id,
                                                    age=data_protein_phenotype_repeated_wide_validation_new$age_v2,sex=data_protein_phenotype_repeated_wide_validation_new$sex,
                                                    T2D_drug=data_protein_phenotype_repeated_wide_validation_new$T2D_med_new_v2,
                                                    hypertension_drug=data_protein_phenotype_repeated_wide_validation_new$hypertension_med_new_v2,
                                                    dyslipi_drug=data_protein_phenotype_repeated_wide_validation_new$dyslipi_med_new_v2,
                                                    data_protein_validation_temp_F2,data_metabolite_validation_temp_F2,check.names=F,stringsAsFactors = F)


## F3
data_protein_validation_temp_F3 <- data_protein_phenotype_repeated_wide_validation_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_validation_new),"_F3")]
data_metabolite_validation_temp_F3 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F3")]

colnames(data_protein_validation_temp_F3) <- str_replace_all(colnames(data_protein_validation_temp_F3),"_F3","_base")
colnames(data_metabolite_validation_temp_F3) <- str_replace_all(colnames(data_metabolite_validation_temp_F3),"_F3","_base")

data_protein_metabolite_validation_F3 <- data.frame(SampleID=paste("F3",data_protein_phenotype_repeated_wide_validation_new$id,sep=""),id=data_protein_phenotype_repeated_wide_validation_new$id,
                                                    age=data_protein_phenotype_repeated_wide_validation_new$age_v3.,sex=data_protein_phenotype_repeated_wide_validation_new$sex,
                                                    T2D_drug=data_protein_phenotype_repeated_wide_validation_new$T2D_med_new_v3,
                                                    hypertension_drug=data_protein_phenotype_repeated_wide_validation_new$hypertension_med_new_v3,
                                                    dyslipi_drug=data_protein_phenotype_repeated_wide_validation_new$dyslipi_med_new_v3,
                                                    data_protein_validation_temp_F3,data_metabolite_validation_temp_F3,check.names=F,stringsAsFactors = F)


data_validation_cross_sectional_new <- rbind(data_protein_metabolite_validation_NL,data_protein_metabolite_validation_F2,data_protein_metabolite_validation_F3)


###### Longitudinal cross-sectional associations between serum proteins and metabolites
## Discovery set 
protein_metabolite_cross_sectional_discovery <- c()
for (i in 1:length(protein_all)){
  for (j in 1:length(metabolite_all)){
    cat("###############",i,"protein","|",j,"metabolite","#####################\n")
    
    protein_temp <- protein_all[i]
    metabolite_temp <- metabolite_all[j]
    
    protein_temp_base <- paste(protein_temp,"_base",sep="")
    metabolite_temp_base <- paste(metabolite_temp,"_base",sep="")
    
    data_fit_temp <- data.frame(protein_base=data_discovery_cross_sectional_new[,protein_temp_base],
                                metabolite_base=data_discovery_cross_sectional_new[,metabolite_temp_base],
                                data_discovery_cross_sectional_new,stringsAsFactors = F)
    
    lm_model_temp <- lm(metabolite_base~protein_base+age+sex+factor(phase),data=data_fit_temp)
    lm_model_temp_summary <- summary(lm_model_temp)
    coef_lm_model <- lm_model_temp_summary$coefficients[2,]
    var_name_temp <- c(protein=protein_temp,metabolite=metabolite_temp)
    coef_lm_model <- c(var_name_temp,coef_lm_model)
    
    protein_metabolite_cross_sectional_discovery <- rbind(protein_metabolite_cross_sectional_discovery,coef_lm_model)
    
  }
}

protein_metabolite_cross_sectional_discovery <- as.data.frame(protein_metabolite_cross_sectional_discovery)
protein_metabolite_cross_sectional_discovery[,-c(1:2)] <- apply(protein_metabolite_cross_sectional_discovery[,-c(1:2)],2,as.numeric)
protein_metabolite_cross_sectional_discovery$id <- paste(protein_metabolite_cross_sectional_discovery$protein," & ",protein_metabolite_cross_sectional_discovery$metabolite,sep="")


## Validation set 
protein_metabolite_cross_sectional_validation <- c()
for (i in 1:length(protein_all)){
  for (j in 1:length(metabolite_all)){
    cat("###############",i,"protein","|",j,"metabolite","#####################\n")
    
    protein_temp <- protein_all[i]
    metabolite_temp <- metabolite_all[j]
    
    protein_temp_base <- paste(protein_temp,"_base",sep="")
    metabolite_temp_base <- paste(metabolite_temp,"_base",sep="")
    
    data_fit_temp <- data.frame(protein_base=data_validation_cross_sectional_new[,protein_temp_base],
                                metabolite_base=data_validation_cross_sectional_new[,metabolite_temp_base],
                                data_validation_cross_sectional_new,stringsAsFactors = F)
    
    lm_model_temp <- lm(metabolite_base~protein_base+age+sex,data=data_fit_temp)
    lm_model_temp_summary <- summary(lm_model_temp)
    coef_lm_model <- lm_model_temp_summary$coefficients[2,]
    var_name_temp <- c(protein=protein_temp,metabolite=metabolite_temp)
    coef_lm_model <- c(var_name_temp,coef_lm_model)
    
    protein_metabolite_cross_sectional_validation <- rbind(protein_metabolite_cross_sectional_validation,coef_lm_model)
    
  }
}

protein_metabolite_cross_sectional_validation <- as.data.frame(protein_metabolite_cross_sectional_validation)
protein_metabolite_cross_sectional_validation[,-c(1:2)] <- apply(protein_metabolite_cross_sectional_validation[,-c(1:2)],2,as.numeric)
protein_metabolite_cross_sectional_validation$id <- paste(protein_metabolite_cross_sectional_validation$protein," & ",protein_metabolite_cross_sectional_validation$metabolite,sep="")


## Meta analysis 
data_discovery <- protein_metabolite_cross_sectional_discovery
data_validation <- protein_metabolite_cross_sectional_validation

library(metafor)
protein_metabolite_cross_sectional_meta <- c()
for (i in 1:nrow(data_discovery)){
  cat("#########",i,"###########\n")
  data_discovery_temp <- data_discovery[i,c("Estimate","Std. Error")]
  data_validation_temp  <- data_validation[i,c("Estimate","Std. Error")]
  
  data_test_temp <- rbind(data_discovery_temp,data_validation_temp)
  meta_result_temp <- rma(yi=Estimate,sei=`Std. Error`,data=data_test_temp,method="FE")
  result_temp <- data.frame(id=data_discovery$id[i],protein=data_discovery$protein[i],metabolite=data_discovery$metabolite[i],
                            B_discovery=data_discovery$Estimate[i],SE_discovery=data_discovery$`Std. Error`[i],P_discovery = data_discovery$`Pr(>|t|)`[i],
                            B_validation=data_validation$Estimate[i],SE_validation=data_validation$`Std. Error`[i],P_validation = data_validation$`Pr(>|t|)`[i],
                            heterogeneity=meta_result_temp$QEp,B_meta=meta_result_temp$beta,SE_meta=meta_result_temp$se,CI_lower_meta=meta_result_temp$ci.lb,CI_upper_meta=meta_result_temp$ci.ub,p_meta=meta_result_temp$pval,stringsAsFactors=F
  )
  protein_metabolite_cross_sectional_meta <- rbind(protein_metabolite_cross_sectional_meta,result_temp)
}

protein_metabolite_cross_sectional_meta$FDR_meta <- p.adjust(protein_metabolite_cross_sectional_meta$p_meta,method="fdr")

protein_metabolite_cross_sectional_meta_select <- protein_metabolite_cross_sectional_meta[protein_metabolite_cross_sectional_meta$P_discovery < 0.05 &
                                                                                            protein_metabolite_cross_sectional_meta$P_validation < 0.05 & 
                                                                                            protein_metabolite_cross_sectional_meta$FDR_meta < 0.05 &
                                                                                            protein_metabolite_cross_sectional_meta$heterogeneity > 0.05,]

#### Sensitivity analysis by further adjusting for medications for T2D, hypertension, and hyperlipidemia
## Discovery set 
protein_metabolite_cross_sectional_discovery_sensitivity <- c()
for (i in 1:length(protein_all)){
  for (j in 1:length(metabolite_all)){
    cat("###############",i,"protein","|",j,"metabolite","#####################\n")
    
    protein_temp <- protein_all[i]
    metabolite_temp <- metabolite_all[j]
    
    protein_temp_base <- paste(protein_temp,"_base",sep="")
    metabolite_temp_base <- paste(metabolite_temp,"_base",sep="")
    
    data_fit_temp <- data.frame(protein_base=data_discovery_cross_sectional_new[,protein_temp_base],
                                metabolite_base=data_discovery_cross_sectional_new[,metabolite_temp_base],
                                data_discovery_cross_sectional_new,stringsAsFactors = F)
    
    lm_model_temp <- lm(metabolite_base~protein_base+age+sex+factor(phase)+T2D_drug+hypertension_drug+dyslipi_drug,data=data_fit_temp)
    lm_model_temp_summary <- summary(lm_model_temp)
    coef_lm_model <- lm_model_temp_summary$coefficients[2,]
    var_name_temp <- c(protein=protein_temp,metabolite=metabolite_temp)
    coef_lm_model <- c(var_name_temp,coef_lm_model)
    
    protein_metabolite_cross_sectional_discovery_sensitivity <- rbind(protein_metabolite_cross_sectional_discovery_sensitivity,coef_lm_model)
    
  }
}

protein_metabolite_cross_sectional_discovery_sensitivity <- as.data.frame(protein_metabolite_cross_sectional_discovery_sensitivity)
protein_metabolite_cross_sectional_discovery_sensitivity[,-c(1:2)] <- apply(protein_metabolite_cross_sectional_discovery_sensitivity[,-c(1:2)],2,as.numeric)


## Validation set 
protein_metabolite_cross_sectional_validation_sensitivity <- c()
for (i in 1:length(protein_all)){
  for (j in 1:length(metabolite_all)){
    cat("###############",i,"protein","|",j,"metabolite","#####################\n")
    
    protein_temp <- protein_all[i]
    metabolite_temp <- metabolite_all[j]
    
    protein_temp_base <- paste(protein_temp,"_base",sep="")
    metabolite_temp_base <- paste(metabolite_temp,"_base",sep="")
    
    data_fit_temp <- data.frame(protein_base=data_validation_cross_sectional_new[,protein_temp_base],
                                metabolite_base=data_validation_cross_sectional_new[,metabolite_temp_base],
                                data_validation_cross_sectional_new,stringsAsFactors = F)
    
    lm_model_temp <- lm(metabolite_base~protein_base+age+sex+T2D_drug+hypertension_drug+dyslipi_drug,data=data_fit_temp)
    lm_model_temp_summary <- summary(lm_model_temp)
    coef_lm_model <- lm_model_temp_summary$coefficients[2,]
    var_name_temp <- c(protein=protein_temp,metabolite=metabolite_temp)
    coef_lm_model <- c(var_name_temp,coef_lm_model)
    
    protein_metabolite_cross_sectional_validation_sensitivity <- rbind(protein_metabolite_cross_sectional_validation_sensitivity,coef_lm_model)
    
  }
}

protein_metabolite_cross_sectional_validation_sensitivity <- as.data.frame(protein_metabolite_cross_sectional_validation_sensitivity)
protein_metabolite_cross_sectional_validation_sensitivity[,-c(1:2)] <- apply(protein_metabolite_cross_sectional_validation_sensitivity[,-c(1:2)],2,as.numeric)


#########################################################################################################
######### Associations of proteins in the identified protein-metabolite pairs with metabolic traits
### Discovery set
# Transform the skewed distribution to normal distribution
data_protein_metabolite_discovery$Ins_new <- log(data_protein_metabolite_discovery$Ins_new+0.01)
data_protein_metabolite_discovery$HOMA_IR_new <- log(data_protein_metabolite_discovery$HOMA_IR_new+0.001)
data_protein_metabolite_discovery$TG_new <- log(data_protein_metabolite_discovery$TG_new+0.001)

Metabolic_traits <- c("BMI_new","wc_new","WHR_new","Glu_new","Ins_new","HOMA_IR_new","HbA1c_new","TC_new","TG_new","LDL_new","HDL_new","SBP_new","DBP_new")
Protein_select <- unique(protein_metabolite_prospective_meta_select$protein)


protein_trait_association_discovery <- c()
for (i in Protein_select){
  for (j in Metabolic_traits){
    cat("########",i,"protein","|",j,"trait","##########\n")
    
    data_temp <- data.frame(trait=scale(data_protein_metabolite_discovery[,j]),
                            protein=data_protein_metabolite_discovery[,i],
                            data_protein_metabolite_discovery,stringsAsFactors = F)
    
    data_temp1 <- data_temp[,c("id","trait","protein","age","sex","phase")]
    data_temp1 <- na.omit(data_temp1)
    
    N_samples <- nrow(data_temp1)
    N_participants <- length(unique(data_temp1$id))
    
    library(lme4)
    library(lmerTest)
    mix_model_temp <- lmer(trait~protein+age+sex+factor(phase)+(1|id),data=data_temp)
    mix_model_temp_summary <- summary(mix_model_temp)
    coef_mixed_model <- mix_model_temp_summary$coefficients[2,]
    
    
    var_name_temp <- c(protein=i,trait=j)
    coef_mixed_model <- c(var_name_temp,N_participants=N_participants,N_samples=N_samples,coef_mixed_model)
    
    protein_trait_association_discovery <- rbind(protein_trait_association_discovery,coef_mixed_model)
    
  }
}
protein_trait_association_discovery <- as.data.frame(protein_trait_association_discovery)
protein_trait_association_discovery[,-c(1:2)] <- apply(protein_trait_association_discovery[,-c(1:2)],2,as.numeric)

protein_trait_association_discovery$id <- paste(protein_trait_association_discovery$protein," & ",protein_trait_association_discovery$trait,sep="")
rownames(protein_trait_association_discovery) <- protein_trait_association_discovery$id


protein_trait_association_discovery$FDR <- p.adjust(protein_trait_association_discovery$`Pr(>|t|)`,method="fdr")
protein_trait_association_discovery_select <- protein_trait_association_discovery[protein_trait_association_discovery$FDR < 0.05,]


### Validation set
# Transform the skewed distribution to normal distribution
data_protein_metabolite_validation$Ins_new <- log(data_protein_metabolite_validation$Ins_new+0.01)
data_protein_metabolite_validation$HOMA_IR_new <- log(data_protein_metabolite_validation$HOMA_IR_new+0.001)
data_protein_metabolite_validation$TG_new <- log(data_protein_metabolite_validation$TG_new+0.001)

Metabolic_traits <- c("BMI_new","wc_new","WHR_new","Glu_new","Ins_new","HOMA_IR_new","HbA1c_new","TC_new","TG_new","LDL_new","HDL_new","SBP_new","DBP_new")
Protein_select <- unique(protein_metabolite_prospective_meta_select$protein)


protein_trait_association_validation <- c()
for (i in 1:nrow(protein_trait_association_discovery_select)){
  cat("########",i,"##########\n")
  protein_temp <- protein_trait_association_discovery_select$protein[i]
  trait_temp <- protein_trait_association_discovery_select$trait[i]
  
  data_temp <- data.frame(trait=scale(data_protein_metabolite_validation[,trait_temp]),
                          protein=data_protein_metabolite_validation[,protein_temp],
                          data_protein_metabolite_validation,stringsAsFactors = F)
  
  data_temp1 <- data_temp[,c("id","trait","protein","age","sex")]
  data_temp1 <- na.omit(data_temp1)
  
  N_samples <- nrow(data_temp1)
  N_participants <- length(unique(data_temp1$id))
  
  library(lme4)
  library(lmerTest)
  mix_model_temp <- lmer(trait~protein+age+sex+(1|id),data=data_temp)
  mix_model_temp_summary <- summary(mix_model_temp)
  coef_mixed_model <- mix_model_temp_summary$coefficients[2,]
  
  var_name_temp <- c(protein=i,trait=j)
  coef_mixed_model <- c(var_name_temp,N_participants=N_participants,N_samples=N_samples,coef_mixed_model)
  
  protein_trait_association_validation <- rbind(protein_trait_association_validation,coef_mixed_model)
  
}

protein_trait_association_validation <- as.data.frame(protein_trait_association_validation)
protein_trait_association_validation[,-c(1:2)] <- apply(protein_trait_association_validation[,-c(1:2)],2,as.numeric)

protein_trait_association_validation$id <- paste(protein_trait_association_validation$protein," & ",protein_trait_association_validation$trait,sep="")
rownames(protein_trait_association_validation) <- protein_trait_association_validation$id

protein_trait_association_validation_select <- protein_trait_association_validation[protein_trait_association_validation$`Pr(>|t|)` < 0.05 &
                                                                              sign(protein_trait_association_validation$Estimate)==sign(protein_trait_association_discovery_select$Estimate),]


##### Sensitivity analysis by further adjusting for medications for T2D, hypertension, and hyperlipidemia
### Discovery set
protein_trait_association_discovery_sensitivity <- c()
for (i in Protein_select){
  for (j in Metabolic_traits){
    cat("########",i,"protein","|",j,"trait","##########\n")
    
    data_temp <- data.frame(trait=scale(data_protein_metabolite_discovery[,j]),
                            protein=data_protein_metabolite_discovery[,i],
                            data_protein_metabolite_discovery,stringsAsFactors = F)
    
    data_temp1 <- data_temp[,c("id","trait","protein","age","sex","phase")]
    data_temp1 <- na.omit(data_temp1)
    
    N_samples <- nrow(data_temp1)
    N_participants <- length(unique(data_temp1$id))
    
    library(lme4)
    library(lmerTest)
    mix_model_temp <- lmer(trait~protein+age+sex+factor(phase)+T2D_drug+hypertension_drug+dyslipi_drug+(1|id),data=data_temp)
    mix_model_temp_summary <- summary(mix_model_temp)
    coef_mixed_model <- mix_model_temp_summary$coefficients[2,]
    
    
    var_name_temp <- c(protein=i,trait=j)
    coef_mixed_model <- c(var_name_temp,N_participants=N_participants,N_samples=N_samples,coef_mixed_model)
    
    protein_trait_association_discovery_sensitivity <- rbind(protein_trait_association_discovery_sensitivity,coef_mixed_model)
    
  }
}

protein_trait_association_discovery_sensitivity <- as.data.frame(protein_trait_association_discovery_sensitivity)
protein_trait_association_discovery_sensitivity[,-c(1:2)] <- apply(protein_trait_association_discovery_sensitivity[,-c(1:2)],2,as.numeric)

protein_trait_association_discovery_sensitivity$id <- paste(protein_trait_association_discovery_sensitivity$protein," & ",protein_trait_association_discovery_sensitivity$trait,sep="")
rownames(protein_trait_association_discovery_sensitivity) <- protein_trait_association_discovery_sensitivity$id


### Validation set
protein_trait_association_validation_sensitivity <- c()
for (i in 1:nrow(protein_trait_association_discovery_select)){
  cat("########",i,"##########\n")
  protein_temp <- protein_trait_association_discovery_select$protein[i]
  trait_temp <- protein_trait_association_discovery_select$trait[i]
  
  data_temp <- data.frame(trait=scale(data_protein_metabolite_validation[,trait_temp]),
                          protein=data_protein_metabolite_validation[,protein_temp],
                          data_protein_metabolite_validation,stringsAsFactors = F)
  
  data_temp1 <- data_temp[,c("id","trait","protein","age","sex")]
  data_temp1 <- na.omit(data_temp1)
  
  N_samples <- nrow(data_temp1)
  N_participants <- length(unique(data_temp1$id))
  
  library(lme4)
  library(lmerTest)
  mix_model_temp <- lmer(trait~protein+age+sex+T2D_drug+hypertension_drug+dyslipi_drug+(1|id),data=data_temp)
  mix_model_temp_summary <- summary(mix_model_temp)
  coef_mixed_model <- mix_model_temp_summary$coefficients[2,]
  
  var_name_temp <- c(protein=i,trait=j)
  coef_mixed_model <- c(var_name_temp,N_participants=N_participants,N_samples=N_samples,coef_mixed_model)
  
  protein_trait_association_validation_sensitivity <- rbind(protein_trait_association_validation_sensitivity,coef_mixed_model)
  
}

protein_trait_association_validation_sensitivity <- as.data.frame(protein_trait_association_validation_sensitivity)
protein_trait_association_validation_sensitivity[,-c(1:2)] <- apply(protein_trait_association_validation_sensitivity[,-c(1:2)],2,as.numeric)

protein_trait_association_validation_sensitivity$id <- paste(protein_trait_association_validation_sensitivity$protein," & ",protein_trait_association_validation_sensitivity$trait,sep="")
rownames(protein_trait_association_validation_sensitivity) <- protein_trait_association_validation_sensitivity$id


#####################################################################################################################################
######### Associations of proteins in the identified protein-metabolite pairs with metabolic diseases
#### Assessed only in the discovery set, due to the small number of incident metabolic disease cases in the validation set
## Based on baseline proteome data in the discovery set
Metabolic_disease_follow <-  c("DM_diagnosis_follow","Hyper_diagnosis_follow","MetS_diagnosis_follow","Obesity_follow")
Metabolic_disease_v1 <-  c("DM_diagnosis_v1","Hyper_diagnosis_v1","MetS_diagnosis_v1","Obesity_v1")
Protein_select <- unique(protein_metabolite_prospective_meta_select$protein)


protein_disease_association_discovery <- c()
for (i in Protein_select){
  for (j in 1:length(Metabolic_disease_follow)){
    cat("########",i,"protein","|",j,"disease","##########\n")
    disease_v1_temp <- Metabolic_disease_v1[j]
    disease_follow_temp <- Metabolic_disease_follow[j]
    
    
    data_temp <- data.frame(disease_follow=data_protein_metabolite_discovery_baseline[,disease_follow_temp],
                            disease_v1 = data_protein_metabolite_discovery_baseline[,disease_v1_temp],
                            protein=data_protein_metabolite_discovery_baseline[,i],
                            data_protein_metabolite_discovery_baseline,stringsAsFactors = F)
    data_temp_new <- data_temp[data_temp$disease_v1==0,]
    
    N_case <- table(data_temp_new$disease_follow)[names(table(data_temp_new$disease_follow))==1]
    
    
    fit_temp <-glm(factor(disease_follow)~protein+age+sex+factor(phase),data=data_temp_new,family=binomial("logit"))
    fit_temp_summary <- summary(fit_temp)
    Coef_temp <- fit_temp_summary$coefficients[2,]
    
    var_name_temp <- c(protein=i,disease=disease_follow_temp,N_total=nrow(data_temp_new),N_case=N_case)
    
    fit_coef_temp <- c(var_name_temp,Coef_temp)
    
    protein_disease_association_discovery <- rbind(protein_disease_association_discovery,fit_coef_temp)
  }
}

protein_disease_association_discovery <- as.data.frame(protein_disease_association_discovery)
protein_disease_association_discovery[,-c(1:2)] <- apply(protein_disease_association_discovery[,-c(1:2)],2,as.numeric)

protein_disease_association_discovery$id <- paste(protein_disease_association_discovery$protein," & ",protein_disease_association_discovery$disease,sep="")
rownames(protein_disease_association_discovery) <- protein_disease_association_discovery$id

protein_disease_association_discovery$FDR <- p.adjust(protein_disease_association_discovery$`Pr(>|z|)`,method="fdr")
protein_disease_association_discovery_select <- protein_disease_association_discovery[protein_disease_association_discovery$FDR < 0.05,]


##### Sensitivity analysis by further adjusting for medications for T2D, hypertension, and hyperlipidemia
## Based on baseline proteome data in the discovery set
protein_disease_association_discovery_sensitivity <- c()
for (i in Protein_select){
  for (j in 1:length(Metabolic_disease_follow)){
    cat("########",i,"protein","|",j,"disease","##########\n")
    disease_v1_temp <- Metabolic_disease_v1[j]
    disease_follow_temp <- Metabolic_disease_follow[j]
    
    
    data_temp <- data.frame(disease_follow=data_protein_metabolite_discovery_baseline[,disease_follow_temp],
                            disease_v1 = data_protein_metabolite_discovery_baseline[,disease_v1_temp],
                            protein=data_protein_metabolite_discovery_baseline[,i],
                            data_protein_metabolite_discovery_baseline,stringsAsFactors = F)
    data_temp_new <- data_temp[data_temp$disease_v1==0,]
    
    N_case <- table(data_temp_new$disease_follow)[names(table(data_temp_new$disease_follow))==1]
    
    
    fit_temp <-glm(factor(disease_follow)~protein+age+sex+factor(phase)+T2D_drug+hypertension_drug+dyslipi_drug,data=data_temp_new,family=binomial("logit"))
    fit_temp_summary <- summary(fit_temp)
    Coef_temp <- fit_temp_summary$coefficients[2,]
    
    var_name_temp <- c(protein=i,disease=disease_follow_temp,N_total=nrow(data_temp_new),N_case=N_case)
    
    fit_coef_temp <- c(var_name_temp,Coef_temp)
    
    protein_disease_association_discovery_sensitivity <- rbind(protein_disease_association_discovery_sensitivity,fit_coef_temp)
  }
}



protein_disease_association_discovery_sensitivity <- as.data.frame(protein_disease_association_discovery_sensitivity)
protein_disease_association_discovery_sensitivity[,-c(1:2)] <- apply(protein_disease_association_discovery_sensitivity[,-c(1:2)],2,as.numeric)

protein_disease_association_discovery_sensitivity$id <- paste(protein_disease_association_discovery_sensitivity$protein," & ",protein_disease_association_discovery_sensitivity$disease,sep="")
rownames(protein_disease_association_discovery_sensitivity) <- protein_disease_association_discovery_sensitivity$id


#########################################################################################################
######### Associations of metabolites in the identified protein-metabolite pairs with metabolic traits
### Discovery set---based on repeated-measured metabolome data
data_metabolite_repeated_measures_discovery$Ins_new <- log(data_metabolite_repeated_measures_discovery$Ins_new+0.01)
data_metabolite_repeated_measures_discovery$HOMA_IR_new <- log(data_metabolite_repeated_measures_discovery$HOMA_IR_new+0.001)
data_metabolite_repeated_measures_discovery$TG_new <- log(data_metabolite_repeated_measures_discovery$TG_new+0.001)

Metabolic_traits <- c("BMI_new","wc_new","WHR_new","Glu_new","Ins_new","HOMA_IR_new","HbA1c_new","TC_new","TG_new","LDL_new","HDL_new","SBP_new","DBP_new")
Metabolite_select <- unique(protein_metabolite_prospective_meta_select$metabolite)


metabolite_trait_association_discovery <- c()
for (i in Metabolite_select){
  for (j in Metabolic_traits){
    cat("########",i,"metabolite","|",j,"trait","##########\n")
    
    data_temp <- data.frame(trait=scale(data_metabolite_repeated_measures_discovery[,j]),
                            metabolite=data_metabolite_repeated_measures_discovery[,i],
                            data_metabolite_repeated_measures_discovery,stringsAsFactors = F)
    
    data_temp1 <- data_temp[,c("id","trait","metabolite","age","sex")]
    data_temp1 <- na.omit(data_temp1)
    
    N_samples <- nrow(data_temp1)
    N_participants <- length(unique(data_temp1$id))
    
    library(lme4)
    library(lmerTest)
    mix_model_temp <- lmer(trait~metabolite+age+sex+(1|id),data=data_temp)
    mix_model_temp_summary <- summary(mix_model_temp)
    coef_mixed_model <- mix_model_temp_summary$coefficients[2,]
    
    
    var_name_temp <- c(metabolite=i,trait=j)
    coef_mixed_model <- c(var_name_temp,N_participants=N_participants,N_samples=N_samples,coef_mixed_model)
    
    metabolite_trait_association_discovery <- rbind(metabolite_trait_association_discovery,coef_mixed_model)
    
  }
}
metabolite_trait_association_discovery <- as.data.frame(metabolite_trait_association_discovery)
metabolite_trait_association_discovery[,-c(1:2)] <- apply(metabolite_trait_association_discovery[,-c(1:2)],2,as.numeric)

metabolite_trait_association_discovery$id <- paste(metabolite_trait_association_discovery$metabolite," & ",metabolite_trait_association_discovery$trait,sep="")
rownames(metabolite_trait_association_discovery) <- metabolite_trait_association_discovery$id

metabolite_trait_association_discovery$FDR <- p.adjust(metabolite_trait_association_discovery$`Pr(>|t|)`,method="fdr")
metabolite_trait_association_discovery_select <- metabolite_trait_association_discovery[metabolite_trait_association_discovery$FDR < 0.05,]


### Validation set---based on cross-sectional metabolome data
data_metabolite_cross_section_validation$Ins_new <- log(data_metabolite_cross_section_validation$Ins_new+0.01)
data_metabolite_cross_section_validation$HOMA_IR_new <- log(data_metabolite_cross_section_validation$HOMA_IR_new+0.001)
data_metabolite_cross_section_validation$TG_new <- log(data_metabolite_cross_section_validation$TG_new+0.001)

metabolite_trait_association_validation <- c()
for (i in 1:nrow(metabolite_trait_association_discovery_select)){
  cat("########",i,"##########\n")
  
  trait_temp <- metabolite_trait_association_discovery_select$trait[i]
  metabolite_temp <- metabolite_trait_association_discovery_select$metabolite[i]
  
  data_temp <- data.frame(trait=scale(data_metabolite_cross_section_validation[,trait_temp]),
                          metabolite=data_metabolite_cross_section_validation[,metabolite_temp],
                          data_metabolite_cross_section_validation,stringsAsFactors = F)
  
  data_temp1 <- data_temp[,c("id","trait","metabolite","age","sex","follow")]
  data_temp1 <- na.omit(data_temp1)
  
  N_participants <- nrow(data_temp1)
  
  lm_model_temp <- lm(trait~metabolite+age+sex+factor(follow),data=data_temp)
  lm_model_temp_summary <- summary(lm_model_temp)
  coef_lm_model <- lm_model_temp_summary$coefficients[2,]
  
  var_name_temp <- c(metabolite=metabolite_temp,trait=trait_temp,N_participants=N_participants)
  coef_lm_model <- c(var_name_temp,coef_lm_model)
  
  metabolite_trait_association_validation <- rbind(metabolite_trait_association_validation,coef_lm_model)
  
}

metabolite_trait_association_validation <- as.data.frame(metabolite_trait_association_validation)
metabolite_trait_association_validation[,-c(1:2)] <- apply(metabolite_trait_association_validation[,-c(1:2)],2,as.numeric)

metabolite_trait_association_validation$id <- paste(metabolite_trait_association_validation$metabolite," & ",metabolite_trait_association_validation$trait,sep="")
rownames(metabolite_trait_association_validation) <- metabolite_trait_association_validation$id

metabolite_trait_association_validation_select <- metabolite_trait_association_validation[metabolite_trait_association_validation$`Pr(>|t|)` < 0.05 & 
                                                                                  sign(metabolite_trait_association_validation$Estimate)==sign(metabolite_trait_association_discovery_select$Estimate),]


##### Sensitivity analysis by further adjusting for medications for T2D, hypertension, and hyperlipidemia
### Discovery set---based on repeated-measured metabolome data
metabolite_trait_association_discovery_sensitivity <- c()
for (i in Metabolite_select){
  for (j in Metabolic_traits){
    cat("########",i,"metabolite","|",j,"trait","##########\n")
    
    data_temp <- data.frame(trait=scale(data_metabolite_repeated_measures_discovery[,j]),
                            metabolite=data_metabolite_repeated_measures_discovery[,i],
                            data_metabolite_repeated_measures_discovery,stringsAsFactors = F)
    
    data_temp1 <- data_temp[,c("id","trait","metabolite","age","sex")]
    data_temp1 <- na.omit(data_temp1)
    
    N_samples <- nrow(data_temp1)
    N_participants <- length(unique(data_temp1$id))
    
    library(lme4)
    library(lmerTest)
    mix_model_temp <- lmer(trait~metabolite+age+sex+T2D_drug+hypertension_drug+dyslipi_drug+(1|id),data=data_temp)
    mix_model_temp_summary <- summary(mix_model_temp)
    coef_mixed_model <- mix_model_temp_summary$coefficients[2,]
    
    
    var_name_temp <- c(metabolite=i,trait=j)
    coef_mixed_model <- c(var_name_temp,N_participants=N_participants,N_samples=N_samples,coef_mixed_model)
    
    metabolite_trait_association_discovery_sensitivity <- rbind(metabolite_trait_association_discovery_sensitivity,coef_mixed_model)
    
  }
}
metabolite_trait_association_discovery_sensitivity <- as.data.frame(metabolite_trait_association_discovery_sensitivity)
metabolite_trait_association_discovery_sensitivity[,-c(1:2)] <- apply(metabolite_trait_association_discovery_sensitivity[,-c(1:2)],2,as.numeric)

metabolite_trait_association_discovery_sensitivity$id <- paste(metabolite_trait_association_discovery_sensitivity$metabolite," & ",metabolite_trait_association_discovery_sensitivity$trait,sep="")
rownames(metabolite_trait_association_discovery_sensitivity) <- metabolite_trait_association_discovery_sensitivity$id


### Validation set---based on cross-sectional metabolome data
metabolite_trait_association_validation_sensitivity <- c()
for (i in 1:nrow(metabolite_trait_association_discovery_select)){
  cat("########",i,"##########\n")
  
  trait_temp <- metabolite_trait_association_discovery_select$trait[i]
  metabolite_temp <- metabolite_trait_association_discovery_select$metabolite[i]
  
  data_temp <- data.frame(trait=scale(data_metabolite_cross_section_validation[,trait_temp]),
                          metabolite=data_metabolite_cross_section_validation[,metabolite_temp],
                          data_metabolite_cross_section_validation,stringsAsFactors = F)
  
  data_temp1 <- data_temp[,c("id","trait","metabolite","age","sex","follow")]
  data_temp1 <- na.omit(data_temp1)
  
  N_participants <- nrow(data_temp1)
  
  lm_model_temp <- lm(trait~metabolite+age+sex+factor(follow)+T2D_drug+hypertension_drug+dyslipi_drug,data=data_temp)
  lm_model_temp_summary <- summary(lm_model_temp)
  coef_lm_model <- lm_model_temp_summary$coefficients[2,]
  
  var_name_temp <- c(metabolite=metabolite_temp,trait=trait_temp,N_participants=N_participants)
  coef_lm_model <- c(var_name_temp,coef_lm_model)
  
  metabolite_trait_association_validation_sensitivity <- rbind(metabolite_trait_association_validation_sensitivity,coef_lm_model)
  
}

metabolite_trait_association_validation_sensitivity <- as.data.frame(metabolite_trait_association_validation_sensitivity)
metabolite_trait_association_validation_sensitivity[,-c(1:2)] <- apply(metabolite_trait_association_validation_sensitivity[,-c(1:2)],2,as.numeric)

metabolite_trait_association_validation_sensitivity$id <- paste(metabolite_trait_association_validation_sensitivity$metabolite," & ",metabolite_trait_association_validation_sensitivity$trait,sep="")
rownames(metabolite_trait_association_validation_sensitivity) <- metabolite_trait_association_validation_sensitivity$id


############################################################################################################
######### Associations of metabolites in the identified protein-metabolite pairs with metabolic diseases
### Discovery set---based on baseline metabolome data in the repeated-measured metabolome data
Metabolic_disease_follow <-  c("DM_diagnosis_follow","Hyper_diagnosis_follow","MetS_diagnosis_follow","Obesity_follow")
Metabolic_disease_v1 <-  c("DM_diagnosis_v1","Hyper_diagnosis_v1","MetS_diagnosis_v1","Obesity_v1")

Metabolite_select <- unique(protein_metabolite_prospective_meta_select$metabolite)

metabolite_disease_association_discovery <- c()
for (i in Metabolite_select){
  for (j in 1:length(Metabolic_disease_follow)){
    cat("########",i,"metabolite","|",j,"disease","##########\n")
    
    disease_v1_temp <- Metabolic_disease_v1[j]
    disease_follow_temp <- Metabolic_disease_follow[j]
    
    data_temp <- data.frame(disease_follow=data_metabolite_repeated_measures_discovery_baseline[,disease_follow_temp],
                            disease_v1 = data_metabolite_repeated_measures_discovery_baseline[,disease_v1_temp],
                            metabolite=data_metabolite_repeated_measures_discovery_baseline[,i],
                            data_metabolite_repeated_measures_discovery_baseline,stringsAsFactors = F)
    
    data_temp_new <- data_temp[data_temp$disease_v1==0,]
    N_case <- table(data_temp_new$disease_follow)[names(table(data_temp_new$disease_follow))==1]
    
    fit_temp <-glm(factor(disease_follow)~metabolite+age+sex,data=data_temp_new,family=binomial("logit"))
    fit_temp_summary <- summary(fit_temp)
    Coef_temp <- fit_temp_summary$coefficients[2,]
    
    var_name_temp <- c(metabolite=i,disease=disease_follow_temp,N_total=nrow(data_temp_new),N_case=N_case)
    
    fit_coef_temp <- c(var_name_temp,Coef_temp)
    
    metabolite_disease_association_discovery <- rbind(metabolite_disease_association_discovery,fit_coef_temp)
    
  }
}

metabolite_disease_association_discovery <- as.data.frame(metabolite_disease_association_discovery)
metabolite_disease_association_discovery[,-c(1:2)] <- apply(metabolite_disease_association_discovery[,-c(1:2)],2,as.numeric)

metabolite_disease_association_discovery$id <- paste(metabolite_disease_association_discovery$metabolite," & ",metabolite_disease_association_discovery$disease,sep="")
rownames(metabolite_disease_association_discovery) <- metabolite_disease_association_discovery$id

metabolite_disease_association_discovery$FDR <- p.adjust(metabolite_disease_association_discovery$`Pr(>|z|)`,method="fdr")
metabolite_disease_association_discovery_select <- metabolite_disease_association_discovery[metabolite_disease_association_discovery$FDR < 0.05,]


### Validation set---based on cross-sectional metabolome data
Metabolic_disease_follow <-  c("DM_diagnosis_follow","Hyper_diagnosis_follow","MetS_diagnosis_follow","Obesity_follow")
Metabolic_disease_new <-  c("DM_diagnosis_new","Hyper_diagnosis_new","MetS_diagnosis_new","Obesity_new")

metabolite_disease_association_validation <- c()
for (i in Metabolite_select){
  for (j in 1:length(Metabolic_disease_follow)){
    cat("########",i,"metabolite","|",j,"disease","##########\n")
    
    disease_new_temp <- Metabolic_disease_new[j]
    disease_follow_temp <- Metabolic_disease_follow[j]
    
    data_temp <- data.frame(disease_follow=data_metabolite_cross_section_validation[,disease_follow_temp],
                            disease_new = data_metabolite_cross_section_validation[,disease_new_temp],
                            metabolite=data_metabolite_cross_section_validation[,i],
                            data_metabolite_cross_section_validation,stringsAsFactors = F)
    data_temp_new <- data_temp[data_temp$disease_new==0,]
  
    N_case <- table(data_temp_new$disease_follow)[names(table(data_temp_new$disease_follow))==1]
    
    fit_temp <-glm(factor(disease_follow)~metabolite+age+sex+factor(follow),data=data_temp_new,family=binomial("logit"))
    fit_temp_summary <- summary(fit_temp)
    Coef_temp <- fit_temp_summary$coefficients[2,]
    
    var_name_temp <- c(metabolite=i,disease=disease_follow_temp,N_total=nrow(data_temp_new),N_case=N_case)
    
    fit_coef_temp <- c(var_name_temp,Coef_temp)
    
    metabolite_disease_association_validation <- rbind(metabolite_disease_association_validation,fit_coef_temp)
    
  }
}

metabolite_disease_association_validation <- as.data.frame(metabolite_disease_association_validation)
metabolite_disease_association_validation[,-c(1:2)] <- apply(metabolite_disease_association_validation[,-c(1:2)],2,as.numeric)

metabolite_disease_association_validation$id <- paste(metabolite_disease_association_validation$metabolite," & ",metabolite_disease_association_validation$disease,sep="")
rownames(metabolite_disease_association_validation) <- metabolite_disease_association_validation$id

metabolite_disease_association_validation_new <- metabolite_disease_association_validation[metabolite_disease_association_discovery_select$id,]
metabolite_disease_association_validation_new <- metabolite_disease_association_validation_new[!is.na(metabolite_disease_association_validation_new$id),]


metabolite_disease_association_discovery_select_temp <- metabolite_disease_association_discovery_select[metabolite_disease_association_validation_new$id,]
metabolite_disease_association_validation_new_select <-  metabolite_disease_association_validation_new[metabolite_disease_association_validation_new$`Pr(>|z|)` < 0.05 & 
                                                                                           sign(metabolite_disease_association_validation_new$Estimate)==sign(metabolite_disease_association_discovery_select_temp$Estimate),]


##### Sensitivity analysis by further adjusting for medications for T2D, hypertension, and hyperlipidemia
### Discovery set---based on baseline metabolome data in the repeated-measured metabolome data
Metabolic_disease_follow <-  c("DM_diagnosis_follow","Hyper_diagnosis_follow","MetS_diagnosis_follow","Obesity_follow")
Metabolic_disease_v1 <-  c("DM_diagnosis_v1","Hyper_diagnosis_v1","MetS_diagnosis_v1","Obesity_v1")

Metabolite_select <- unique(protein_metabolite_prospective_meta_select$metabolite)

metabolite_disease_association_discovery_sensitivity <- c()
for (i in Metabolite_select){
  for (j in 1:length(Metabolic_disease_follow)){
    cat("########",i,"metabolite","|",j,"disease","##########\n")
    
    disease_v1_temp <- Metabolic_disease_v1[j]
    disease_follow_temp <- Metabolic_disease_follow[j]
    
    data_temp <- data.frame(disease_follow=data_metabolite_repeated_measures_discovery_baseline[,disease_follow_temp],
                            disease_v1 = data_metabolite_repeated_measures_discovery_baseline[,disease_v1_temp],
                            metabolite=data_metabolite_repeated_measures_discovery_baseline[,i],
                            data_metabolite_repeated_measures_discovery_baseline,stringsAsFactors = F)
    
    data_temp_new <- data_temp[data_temp$disease_v1==0,]
    N_case <- table(data_temp_new$disease_follow)[names(table(data_temp_new$disease_follow))==1]
    
    fit_temp <-glm(factor(disease_follow)~metabolite+age+sex+T2D_drug+hypertension_drug+dyslipi_drug,data=data_temp_new,family=binomial("logit"))
    fit_temp_summary <- summary(fit_temp)
    Coef_temp <- fit_temp_summary$coefficients[2,]
    
    var_name_temp <- c(metabolite=i,disease=disease_follow_temp,N_total=nrow(data_temp_new),N_case=N_case)
    
    fit_coef_temp <- c(var_name_temp,Coef_temp)
    
    metabolite_disease_association_discovery_sensitivity <- rbind(metabolite_disease_association_discovery_sensitivity,fit_coef_temp)
    
  }
}

metabolite_disease_association_discovery_sensitivity <- as.data.frame(metabolite_disease_association_discovery_sensitivity)
metabolite_disease_association_discovery_sensitivity[,-c(1:2)] <- apply(metabolite_disease_association_discovery_sensitivity[,-c(1:2)],2,as.numeric)

metabolite_disease_association_discovery_sensitivity$id <- paste(metabolite_disease_association_discovery_sensitivity$metabolite," & ",metabolite_disease_association_discovery_sensitivity$disease,sep="")
rownames(metabolite_disease_association_discovery_sensitivity) <- metabolite_disease_association_discovery_sensitivity$id


### Validation set---based on cross-sectional metabolome data
Metabolic_disease_follow <-  c("DM_diagnosis_follow","Hyper_diagnosis_follow","MetS_diagnosis_follow","Obesity_follow")
Metabolic_disease_new <-  c("DM_diagnosis_new","Hyper_diagnosis_new","MetS_diagnosis_new","Obesity_new")

metabolite_disease_association_validation_sensitivity <- c()
for (i in Metabolite_select){
  for (j in 1:length(Metabolic_disease_follow)){
    cat("########",i,"metabolite","|",j,"disease","##########\n")
    
    disease_new_temp <- Metabolic_disease_new[j]
    disease_follow_temp <- Metabolic_disease_follow[j]
    
    data_temp <- data.frame(disease_follow=data_metabolite_cross_section_validation[,disease_follow_temp],
                            disease_new = data_metabolite_cross_section_validation[,disease_new_temp],
                            metabolite=data_metabolite_cross_section_validation[,i],
                            data_metabolite_cross_section_validation,stringsAsFactors = F)
    data_temp_new <- data_temp[data_temp$disease_new==0,]
    
    N_case <- table(data_temp_new$disease_follow)[names(table(data_temp_new$disease_follow))==1]
    
    fit_temp <-glm(factor(disease_follow)~metabolite+age+sex+factor(follow)+T2D_drug+hypertension_drug+dyslipi_drug,data=data_temp_new,family=binomial("logit"))
    fit_temp_summary <- summary(fit_temp)
    Coef_temp <- fit_temp_summary$coefficients[2,]
    
    var_name_temp <- c(metabolite=i,disease=disease_follow_temp,N_total=nrow(data_temp_new),N_case=N_case)
    
    fit_coef_temp <- c(var_name_temp,Coef_temp)
    
    metabolite_disease_association_validation_sensitivity <- rbind(metabolite_disease_association_validation_sensitivity,fit_coef_temp)
    
  }
}

metabolite_disease_association_validation_sensitivity <- as.data.frame(metabolite_disease_association_validation_sensitivity)
metabolite_disease_association_validation_sensitivity[,-c(1:2)] <- apply(metabolite_disease_association_validation_sensitivity[,-c(1:2)],2,as.numeric)

metabolite_disease_association_validation_sensitivity$id <- paste(metabolite_disease_association_validation_sensitivity$metabolite," & ",metabolite_disease_association_validation_sensitivity$disease,sep="")
rownames(metabolite_disease_association_validation_sensitivity) <- metabolite_disease_association_validation_sensitivity$id



####################################################################################################################
###################### Summarize the above results---potential protein-metabolite-metabolic trait/disease pathways
### For potential protein-metabolite-metabolic trait pathways
protein_metabolite_traits <- merge(protein_metabolite_prospective_meta_select[,c("id","protein","metabolite","B_meta")],protein_trait_association_validation_select[,1:3],by="protein",suffixes = c("_protien_metabolite","_protein_trait"))
protein_metabolite_traits_new <- merge(protein_metabolite_traits,metabolite_trait_association_validation_select,by="metabolite",suffixes = c("_protien_trait","_metabolite_trait"))

protein_metabolite_traits_new <- protein_metabolite_traits_new[protein_metabolite_traits_new$trait_protien_trait==protein_metabolite_traits_new$trait_metabolite_trait,]

protein_metabolite_traits_new$sign_mediate <- sign(protein_metabolite_traits_new$B_meta) * sign(protein_metabolite_traits_new$Estimate_metabolite_trait)
protein_metabolite_traits_select <- protein_metabolite_traits_new[sign(protein_metabolite_traits_new$Estimate_protien_trait)==sign(protein_metabolite_traits_new$sign_mediate),]


### For potential protein-metabolite-metabolic disease pathways
protein_metabolite_diseases <- merge(protein_metabolite_prospective_meta_select[,c("id","protein","metabolite","B_meta")],protein_disease_association_discovery_select[,c("protein","disease","Estimate")],by="protein",suffixes = c("_protien_metabolite","_protein_disease"))
protein_metabolite_diseases_new <- merge(protein_metabolite_diseases,metabolite_disease_association_validation_new_select[,c("metabolite","disease","Estimate")],by="metabolite",suffixes = c("_protien_disease","_metabolite_disease"))

protein_metabolite_diseases_new <- protein_metabolite_diseases_new[protein_metabolite_diseases_new$disease_protien_disease==protein_metabolite_diseases_new$disease_metabolite_disease,]

protein_metabolite_diseases_new$sign_mediate <- sign(protein_metabolite_diseases_new$B_meta) * sign(protein_metabolite_diseases_new$Estimate_metabolite_disease)
protein_metabolite_diseases_select <- protein_metabolite_diseases_new[sign(protein_metabolite_diseases_new$Estimate_protien_disease)==sign(protein_metabolite_diseases_new$sign_mediate),]


#################################################################################################################
######################## Mediation analysis for protein-metabolite-metabolic trait/disease pathways 
##### Dataset construction for discovery set
## NL
data_protein_discovery_temp_NL <- data_protein_phenotype_repeated_wide_discovery_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_discovery_new),"_NL")]
data_metabolite_discovery_temp_NL <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_NL")]
data_metabolite_discovery_temp_F2 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F2")]

colnames(data_protein_discovery_temp_NL) <- str_replace_all(colnames(data_protein_discovery_temp_NL),"_NL","_base")
colnames(data_metabolite_discovery_temp_NL) <- str_replace_all(colnames(data_metabolite_discovery_temp_NL),"_NL","_base")
colnames(data_metabolite_discovery_temp_F2) <- str_replace_all(colnames(data_metabolite_discovery_temp_F2),"_F2","_follow")

data_traits_discovery_temp_F2 <- data_protein_phenotype_repeated_wide_discovery_new[,c("BMI_v2","wc_v2","WHR_v2","GLU_v2","ins_v2","HOMA_IR_v2","hba1c_v2",
                                                                                       "TC_v2","TG_v2","LDL_v2","HDL_v2","SBP_v2","DBP_v2")]

colnames(data_traits_discovery_temp_F2)[colnames(data_traits_discovery_temp_F2)=="ins_v2"] <- "Insulin_v2"

colnames(data_traits_discovery_temp_F2) <- str_replace_all(colnames(data_traits_discovery_temp_F2),"_v2","_follow")

data_disease_discovery_temp_F2 <- data_protein_phenotype_repeated_wide_discovery_new[,c("DM_diagnosis_v2","Hyper_diagnosis_v2","MetS_diagnosis_v2","Obesity_v2")]
colnames(data_disease_discovery_temp_F2) <- str_replace_all(colnames(data_disease_discovery_temp_F2),"_v2","_follow")


data_protein_metabolite_NL_F2_discovery_new <- data.frame(SampleID=data_protein_phenotype_repeated_wide_discovery_new$id,id=data_protein_phenotype_repeated_wide_discovery_new$id,
                                                          age=data_protein_phenotype_repeated_wide_discovery_new$age,sex=data_protein_phenotype_repeated_wide_discovery_new$sex,
                                                          follow_up_time = data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v2_v1,
                                                          data_protein_discovery_temp_NL,data_metabolite_discovery_temp_NL,data_metabolite_discovery_temp_F2,data_traits_discovery_temp_F2,data_disease_discovery_temp_F2,check.names=F,stringsAsFactors = F)


## F2
data_protein_discovery_temp_F2 <- data_protein_phenotype_repeated_wide_discovery_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_discovery_new),"_F2")]
data_metabolite_discovery_temp_F2 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F2")]
data_metabolite_discovery_temp_F3 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F3")]

colnames(data_protein_discovery_temp_F2) <- str_replace_all(colnames(data_protein_discovery_temp_F2),"_F2","_base")
colnames(data_metabolite_discovery_temp_F2) <- str_replace_all(colnames(data_metabolite_discovery_temp_F2),"_F2","_base")
colnames(data_metabolite_discovery_temp_F3) <- str_replace_all(colnames(data_metabolite_discovery_temp_F3),"_F3","_follow")

data_traits_discovery_temp_F3 <- data_protein_phenotype_repeated_wide_discovery_new[,c("BMI_v3","wc_v3","WHR_v3","GLU_v3","Insulin_v3","HOMA_IR_v3","hba1c_v3",
                                                                                       "TC_v3","TG_v3","LDL_v3","HDL_v3","SBP_v3","DBP_v3")]

colnames(data_traits_discovery_temp_F3) <- str_replace_all(colnames(data_traits_discovery_temp_F3),"_v3","_follow")

data_disease_discovery_temp_F3 <- data_protein_phenotype_repeated_wide_discovery_new[,c("DM_diagnosis_v3","Hyper_diagnosis_v3","MetS_diagnosis_v3","Obesity_v3")]
colnames(data_disease_discovery_temp_F3) <- str_replace_all(colnames(data_disease_discovery_temp_F3),"_v3","_follow")


data_protein_metabolite_F2_F3_discovery_new <- data.frame(SampleID=paste("F2",data_protein_phenotype_repeated_wide_discovery_new$id,sep=""),id=data_protein_phenotype_repeated_wide_discovery_new$id,
                                                          age=data_protein_phenotype_repeated_wide_discovery_new$age_v2,sex=data_protein_phenotype_repeated_wide_discovery_new$sex,
                                                          follow_up_time = data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v3_v2,
                                                          data_protein_discovery_temp_F2,data_metabolite_discovery_temp_F2,data_metabolite_discovery_temp_F3,data_traits_discovery_temp_F3,data_disease_discovery_temp_F3,check.names=F,stringsAsFactors = F)


## F3
data_protein_discovery_temp_F3 <- data_protein_phenotype_repeated_wide_discovery_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_discovery_new),"_F3")]
data_metabolite_discovery_temp_F3 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F3")]
data_metabolite_discovery_temp_F4 <- data_metabolomics_repeated_wide_discovery[,str_detect(colnames(data_metabolomics_repeated_wide_discovery),"_F4")]

colnames(data_protein_discovery_temp_F3) <- str_replace_all(colnames(data_protein_discovery_temp_F3),"_F3","_base")
colnames(data_metabolite_discovery_temp_F3) <- str_replace_all(colnames(data_metabolite_discovery_temp_F3),"_F3","_base")
colnames(data_metabolite_discovery_temp_F4) <- str_replace_all(colnames(data_metabolite_discovery_temp_F4),"_F4","_follow")

data_traits_discovery_temp_F4 <- data_protein_phenotype_repeated_wide_discovery_new[,c("BMI_v4","wc_v4","WHR_v4","GLU_v4","Insulin_v4","HOMA_IR_v4","hba1c_v4",
                                                                                       "TC_v4","TG_v4","LDL_v4","HDL_v4","SBP_v4","DBP_v4")]

colnames(data_traits_discovery_temp_F4) <- str_replace_all(colnames(data_traits_discovery_temp_F4),"_v4","_follow")

data_disease_discovery_temp_F4 <- data_protein_phenotype_repeated_wide_discovery_new[,c("DM_diagnosis_v4","Hyper_diagnosis_v4","MetS_diagnosis_v4","Obesity_v4")]
colnames(data_disease_discovery_temp_F4) <- str_replace_all(colnames(data_disease_discovery_temp_F4),"_v4","_follow")


data_protein_metabolite_F3_F4_discovery_new <- data.frame(SampleID=paste("F3",data_protein_phenotype_repeated_wide_discovery_new$id,sep=""),id=data_protein_phenotype_repeated_wide_discovery_new$id,
                                                          age=data_protein_phenotype_repeated_wide_discovery_new$age_v3,sex=data_protein_phenotype_repeated_wide_discovery_new$sex,
                                                          follow_up_time = data_protein_phenotype_repeated_wide_discovery_new$follow_up_time_v4_v3,
                                                          data_protein_discovery_temp_F3,data_metabolite_discovery_temp_F3,data_metabolite_discovery_temp_F4,data_traits_discovery_temp_F4,data_disease_discovery_temp_F4,check.names=F,stringsAsFactors = F)


data_metabolite_protein_discovery_new_all <- rbind(data_protein_metabolite_NL_F2_discovery_new,data_protein_metabolite_F2_F3_discovery_new,data_protein_metabolite_F3_F4_discovery_new)

# Transform the skewed distribution to normal distribution
data_metabolite_protein_discovery_new_all$Insulin_follow <- log(data_metabolite_protein_discovery_new_all$Insulin_follow+0.01)

data_metabolite_protein_discovery_new_all$HOMA_IR_follow <- log(data_metabolite_protein_discovery_new_all$HOMA_IR_follow+0.001)

data_metabolite_protein_discovery_new_all$TG_follow <- log(data_metabolite_protein_discovery_new_all$TG_follow+0.001)


colnames(data_metabolite_protein_discovery_new_all)[colnames(data_metabolite_protein_discovery_new_all)=="Insulin_follow"] <- "Ins_follow"
colnames(data_metabolite_protein_discovery_new_all)[colnames(data_metabolite_protein_discovery_new_all)=="GLU_follow"] <- "Glu_follow"
colnames(data_metabolite_protein_discovery_new_all)[colnames(data_metabolite_protein_discovery_new_all)=="hba1c_follow"] <- "HbA1c_follow"


##### Dataset construction for validation set
## NL
data_protein_validation_temp_NL <- data_protein_phenotype_repeated_wide_validation_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_validation_new),"_NL")]
data_metabolite_validation_temp_NL <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_NL")]
data_metabolite_validation_temp_F2 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F2")]

colnames(data_protein_validation_temp_NL) <- str_replace_all(colnames(data_protein_validation_temp_NL),"_NL","_base")
colnames(data_metabolite_validation_temp_NL) <- str_replace_all(colnames(data_metabolite_validation_temp_NL),"_NL","_base")
colnames(data_metabolite_validation_temp_F2) <- str_replace_all(colnames(data_metabolite_validation_temp_F2),"_F2","_follow")

data_traits_validation_temp_F2 <- data_protein_phenotype_repeated_wide_validation_new[,c("BMI_v2","wc_v2","WHR_v2","GLU_v2","ins_v2","HOMA_IR_v2","hba1c_v2",
                                                                                         "TC_v2","TG_v2","LDL_v2","HDL_v2","SBP_v2","DBP_v2")]

colnames(data_traits_validation_temp_F2)[colnames(data_traits_validation_temp_F2)=="ins_v2"] <- "Insulin_v2"

colnames(data_traits_validation_temp_F2) <- str_replace_all(colnames(data_traits_validation_temp_F2),"_v2","_follow")

data_disease_validation_temp_F2 <- data_protein_phenotype_repeated_wide_validation_new[,c("DM_diagnosis_v2","Hyper_diagnosis_v2","MetS_diagnosis_v2","Obesity_v2")]
colnames(data_disease_validation_temp_F2) <- str_replace_all(colnames(data_disease_validation_temp_F2),"_v2","_follow")


data_protein_metabolite_NL_F2_validation_new <- data.frame(SampleID=data_protein_phenotype_repeated_wide_validation_new$id,id=data_protein_phenotype_repeated_wide_validation_new$id,
                                                           age=data_protein_phenotype_repeated_wide_validation_new$age,sex=data_protein_phenotype_repeated_wide_validation_new$sex,
                                                           follow_up_time = data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v2_v1,
                                                           data_protein_validation_temp_NL,data_metabolite_validation_temp_NL,data_metabolite_validation_temp_F2,data_traits_validation_temp_F2,data_disease_validation_temp_F2,check.names=F,stringsAsFactors = F)


## F2
data_protein_validation_temp_F2 <- data_protein_phenotype_repeated_wide_validation_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_validation_new),"_F2")]
data_metabolite_validation_temp_F2 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F2")]
data_metabolite_validation_temp_F3 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F3")]

colnames(data_protein_validation_temp_F2) <- str_replace_all(colnames(data_protein_validation_temp_F2),"_F2","_base")
colnames(data_metabolite_validation_temp_F2) <- str_replace_all(colnames(data_metabolite_validation_temp_F2),"_F2","_base")
colnames(data_metabolite_validation_temp_F3) <- str_replace_all(colnames(data_metabolite_validation_temp_F3),"_F3","_follow")

data_traits_validation_temp_F3 <- data_protein_phenotype_repeated_wide_validation_new[,c("BMI_v3","wc_v3","WHR_v3","GLU_v3","Insulin_v3","HOMA_IR_v3","hba1c_v3",
                                                                                         "TC_v3","TG_v3","LDL_v3","HDL_v3","SBP_v3","DBP_v3")]

colnames(data_traits_validation_temp_F3) <- str_replace_all(colnames(data_traits_validation_temp_F3),"_v3","_follow")

data_disease_validation_temp_F3 <- data_protein_phenotype_repeated_wide_validation_new[,c("DM_diagnosis_v3","Hyper_diagnosis_v3","MetS_diagnosis_v3","Obesity_v3")]
colnames(data_disease_validation_temp_F3) <- str_replace_all(colnames(data_disease_validation_temp_F3),"_v3","_follow")


data_protein_metabolite_F2_F3_validation_new <- data.frame(SampleID=paste("F2",data_protein_phenotype_repeated_wide_validation_new$id,sep=""),id=data_protein_phenotype_repeated_wide_validation_new$id,
                                                           age=data_protein_phenotype_repeated_wide_validation_new$age_v2,sex=data_protein_phenotype_repeated_wide_validation_new$sex,
                                                           follow_up_time = data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v3_v2,
                                                           data_protein_validation_temp_F2,data_metabolite_validation_temp_F2,data_metabolite_validation_temp_F3,data_traits_validation_temp_F3,data_disease_validation_temp_F3,check.names=F,stringsAsFactors = F)


## F3
data_protein_validation_temp_F3 <- data_protein_phenotype_repeated_wide_validation_new[,str_detect(colnames(data_protein_phenotype_repeated_wide_validation_new),"_F3")]
data_metabolite_validation_temp_F3 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F3")]
data_metabolite_validation_temp_F4 <- data_metabolomics_repeated_wide_validation[,str_detect(colnames(data_metabolomics_repeated_wide_validation),"_F4")]

colnames(data_protein_validation_temp_F3) <- str_replace_all(colnames(data_protein_validation_temp_F3),"_F3","_base")
colnames(data_metabolite_validation_temp_F3) <- str_replace_all(colnames(data_metabolite_validation_temp_F3),"_F3","_base")
colnames(data_metabolite_validation_temp_F4) <- str_replace_all(colnames(data_metabolite_validation_temp_F4),"_F4","_follow")

data_traits_validation_temp_F4 <- data_protein_phenotype_repeated_wide_validation_new[,c("BMI_v4","wc_v4","WHR_v4","GLU_v4","Insulin_v4","HOMA_IR_v4","hba1c_v4",
                                                                                         "TC_v4","TG_v4","LDL_v4","HDL_v4","SBP_v4","DBP_v4")]

colnames(data_traits_validation_temp_F4) <- str_replace_all(colnames(data_traits_validation_temp_F4),"_v4","_follow")

data_disease_validation_temp_F4 <- data_protein_phenotype_repeated_wide_validation_new[,c("DM_diagnosis_v4","Hyper_diagnosis_v4","MetS_diagnosis_v4","Obesity_v4")]
colnames(data_disease_validation_temp_F4) <- str_replace_all(colnames(data_disease_validation_temp_F4),"_v4","_follow")


data_protein_metabolite_F3_F4_validation_new <- data.frame(SampleID=paste("F3",data_protein_phenotype_repeated_wide_validation_new$id,sep=""),id=data_protein_phenotype_repeated_wide_validation_new$id,
                                                           age=data_protein_phenotype_repeated_wide_validation_new$age_v3,sex=data_protein_phenotype_repeated_wide_validation_new$sex,
                                                           follow_up_time = data_protein_phenotype_repeated_wide_validation_new$follow_up_time_v4_v3,
                                                           data_protein_validation_temp_F3,data_metabolite_validation_temp_F3,data_metabolite_validation_temp_F4,data_traits_validation_temp_F4,data_disease_validation_temp_F4,check.names=F,stringsAsFactors = F)


data_metabolite_protein_validation_new_all <- rbind(data_protein_metabolite_NL_F2_validation_new,data_protein_metabolite_F2_F3_validation_new,data_protein_metabolite_F3_F4_validation_new)

# Transform the skewed distribution to normal distribution
data_metabolite_protein_validation_new_all$Insulin_follow <- log(data_metabolite_protein_validation_new_all$Insulin_follow+0.01)

data_metabolite_protein_validation_new_all$HOMA_IR_follow <- log(data_metabolite_protein_validation_new_all$HOMA_IR_follow+0.001)

data_metabolite_protein_validation_new_all$TG_follow <- log(data_metabolite_protein_validation_new_all$TG_follow+0.001)


colnames(data_metabolite_protein_validation_new_all)[colnames(data_metabolite_protein_validation_new_all)=="Insulin_follow"] <- "Ins_follow"
colnames(data_metabolite_protein_validation_new_all)[colnames(data_metabolite_protein_validation_new_all)=="GLU_follow"] <- "Glu_follow"
colnames(data_metabolite_protein_validation_new_all)[colnames(data_metabolite_protein_validation_new_all)=="hba1c_follow"] <- "HbA1c_follow"


############### For protein-metabolite-metabolic trait pathways
## Discovery set
Mediation_protein_metabolite_trait_discovery <- c()
for (i in 1:nrow(protein_metabolite_traits_select)){
  cat("############",i,"|",nrow(protein_metabolite_traits_select),"#################\n")
  protein_temp <- protein_metabolite_traits_select$protein[i]
  metabolite_temp <- protein_metabolite_traits_select$metabolite[i]
  trait_temp <- protein_metabolite_traits_select$trait_protien_trait[i]
  
  protein_base_temp <- paste(protein_temp,"_base",sep="")
  metabolite_base_temp <- paste(metabolite_temp,"_base",sep="")
  metabolite_follow_temp <- paste(metabolite_temp,"_follow",sep="")
  trait_follow_temp <- str_replace(trait_temp,"_new","_follow")
  
  data_temp <- data.frame(protein_base = data_metabolite_protein_discovery_new_all[,protein_base_temp],
                          metabolite_base = data_metabolite_protein_discovery_new_all[,metabolite_base_temp],
                          metabolite_follow = data_metabolite_protein_discovery_new_all[,metabolite_follow_temp],
                          trait_follow = scale(data_metabolite_protein_discovery_new_all[,trait_follow_temp]),
                          data_metabolite_protein_discovery_new_all,stringsAsFactors = F)
  
  data_temp <- data_temp[!is.na(data_temp$trait_follow),]
  
  N_participants <- length(unique(data_temp$id))
  N_samples <- nrow(data_temp)
  
  library(mediation)
  library(lme4)
 
  Med_fit_temp <- lmer(metabolite_follow~protein_base+metabolite_base+age+sex+follow_up_time+factor(phase)+(1|id),data=data_temp)
  Out_fit_temp <- lmer(trait_follow~metabolite_follow+protein_base+age+sex+follow_up_time+factor(phase)+(1|id),data=data_temp)
  
  set.seed(1234)
  Med_protein_metabolite_trait_dis<- mediate(model.m=Med_fit_temp,model.y=Out_fit_temp,treat="protein_base",mediator="metabolite_follow")
  Med_protein_metabolite_trait_dis_summary <- summary(Med_protein_metabolite_trait_dis)
  
  results_temp <- c(Total_Effect=Med_protein_metabolite_trait_dis_summary$tau.coef,Total_Effect_P = Med_protein_metabolite_trait_dis_summary$tau.p,
                    ADE=Med_protein_metabolite_trait_dis_summary$z.avg,ADE_P = Med_protein_metabolite_trait_dis_summary$z.avg.p,
                    ACME=Med_protein_metabolite_trait_dis_summary$d.avg,
                    proportions_mediated=Med_protein_metabolite_trait_dis_summary$n.avg,ACME_P=Med_protein_metabolite_trait_dis_summary$d.avg.p)
  
  results_temp <- c(protein=protein_temp,metabolite=metabolite_temp,trait = trait_temp,N_participants=N_participants,N_samples=N_samples,results_temp)
  
  Mediation_protein_metabolite_trait_discovery <- rbind(Mediation_protein_metabolite_trait_discovery,results_temp)
  
}

Mediation_protein_metabolite_trait_discovery <- as.data.frame(Mediation_protein_metabolite_trait_discovery)
Mediation_protein_metabolite_trait_discovery[,-c(1:3)] <- apply(Mediation_protein_metabolite_trait_discovery[,-c(1:3)],2,as.numeric)

Mediation_protein_metabolite_trait_discovery$FDR <- p.adjust(Mediation_protein_metabolite_trait_discovery$ACME_P,method="fdr")
Mediation_protein_metabolite_trait_discovery$id <- paste(Mediation_protein_metabolite_trait_discovery$protein,"&",Mediation_protein_metabolite_trait_discovery$metabolite,"&",Mediation_protein_metabolite_trait_discovery$trait,sep=" ")

Mediation_protein_metabolite_trait_discovery_select <- Mediation_protein_metabolite_trait_discovery[Mediation_protein_metabolite_trait_discovery$FDR < 0.05,]


## Validation set
Mediation_protein_metabolite_trait_validation <- c()
for (i in 1:nrow(protein_metabolite_traits_select)){
  cat("############",i,"|",nrow(protein_metabolite_traits_select),"#################\n")
  protein_temp <- protein_metabolite_traits_select$protein[i]
  metabolite_temp <- protein_metabolite_traits_select$metabolite[i]
  trait_temp <- protein_metabolite_traits_select$trait_protien_trait[i]
  
  protein_base_temp <- paste(protein_temp,"_base",sep="")
  metabolite_base_temp <- paste(metabolite_temp,"_base",sep="")
  metabolite_follow_temp <- paste(metabolite_temp,"_follow",sep="")
  trait_follow_temp <- str_replace(trait_temp,"_new","_follow")
  
  data_temp <- data.frame(protein_base = data_metabolite_protein_validation_new_all[,protein_base_temp],
                          metabolite_base = data_metabolite_protein_validation_new_all[,metabolite_base_temp],
                          metabolite_follow = data_metabolite_protein_validation_new_all[,metabolite_follow_temp],
                          trait_follow = scale(data_metabolite_protein_validation_new_all[,trait_follow_temp]),
                          data_metabolite_protein_validation_new_all,stringsAsFactors = F)
  
  data_temp <- data_temp[!is.na(data_temp$trait_follow),]
  
  N_participants <- length(unique(data_temp$id))
  N_samples <- nrow(data_temp)
  
  library(mediation)
  library(lme4)
  
  Med_fit_temp <- lmer(metabolite_follow~protein_base+metabolite_base+age+sex+follow_up_time+(1|id),data=data_temp)
  Out_fit_temp <- lmer(trait_follow~metabolite_follow+protein_base+age+sex+follow_up_time+(1|id),data=data_temp)
  
  set.seed(1234)
  Med_protein_metabolite_trait_val<- mediate(model.m=Med_fit_temp,model.y=Out_fit_temp,treat="protein_base",mediator="metabolite_follow")
  Med_protein_metabolite_trait_val_summary <- summary(Med_protein_metabolite_trait_val)
  
  results_temp <- c(Total_Effect=Med_protein_metabolite_trait_val_summary$tau.coef,Total_Effect_P = Med_protein_metabolite_trait_val_summary$tau.p,
                    ADE=Med_protein_metabolite_trait_val_summary$z.avg,ADE_P = Med_protein_metabolite_trait_val_summary$z.avg.p,
                    ACME=Med_protein_metabolite_trait_val_summary$d.avg,
                    proportions_mediated=Med_protein_metabolite_trait_val_summary$n.avg,ACME_P=Med_protein_metabolite_trait_val_summary$d.avg.p)
  
  results_temp <- c(protein=protein_temp,metabolite=metabolite_temp,trait = trait_temp,N_participants=N_participants,N_samples=N_samples,results_temp)
  
  Mediation_protein_metabolite_trait_validation <- rbind(Mediation_protein_metabolite_trait_validation,results_temp)
  
}

Mediation_protein_metabolite_trait_validation <- as.data.frame(Mediation_protein_metabolite_trait_validation)
Mediation_protein_metabolite_trait_validation[,-c(1:3)] <- apply(Mediation_protein_metabolite_trait_validation[,-c(1:3)],2,as.numeric)

Mediation_protein_metabolite_trait_validation$FDR <- p.adjust(Mediation_protein_metabolite_trait_validation$ACME_P,method="fdr")
Mediation_protein_metabolite_trait_validation$id <- paste(Mediation_protein_metabolite_trait_validation$protein,"&",Mediation_protein_metabolite_trait_validation$metabolite,"&",Mediation_protein_metabolite_trait_validation$trait,sep=" ")

Mediation_protein_metabolite_trait_validation_temp <- Mediation_protein_metabolite_trait_validation[Mediation_protein_metabolite_trait_validation$id %in% Mediation_protein_metabolite_trait_discovery_select$id,]
Mediation_protein_metabolite_trait_validation_select <- Mediation_protein_metabolite_trait_validation_temp[Mediation_protein_metabolite_trait_validation_temp$ACME_P < 0.05,]
Mediation_protein_metabolite_trait_validation_select <- Mediation_protein_metabolite_trait_validation_select[Mediation_protein_metabolite_trait_validation_select$proportions_mediated > 0,]


############### For protein-metabolite-metabolic disease pathways
##  Assessed only in the discovery set, due to the small number of incident metabolic disease cases in the validation set.
Mediation_protein_metabolite_disease_discovery <- c()
for (i in 1:nrow(protein_metabolite_diseases_select)){
  cat("############",i,"|",nrow(protein_metabolite_diseases_select),"#################\n")
  protein_temp <- protein_metabolite_diseases_select$protein[i]
  metabolite_temp <- protein_metabolite_diseases_select$metabolite[i]
  disease_temp <- protein_metabolite_diseases_select$disease_protien_disease[i]
  
  
  protein_base_temp <- paste(protein_temp,"_base",sep="")
  metabolite_base_temp <- paste(metabolite_temp,"_base",sep="")
  metabolite_follow_temp <- paste(metabolite_temp,"_follow",sep="")
  
  data_temp <- data.frame(disease=data_metabolite_protein_discovery_new_all[,disease_temp],
                          protein_base=data_metabolite_protein_discovery_new_all[,protein_base_temp],
                          metabolite_base = data_metabolite_protein_discovery_new_all[,metabolite_base_temp],
                          metabolite_follow = data_metabolite_protein_discovery_new_all[,metabolite_follow_temp],
                          data_metabolite_protein_discovery_new_all,stringsAsFactors = F)
  
  data_temp <- data_temp[!is.na(data_temp$disease),]
  
  N_participants <- length(unique(data_temp$id))
  N_samples <- nrow(data_temp)
  
  library(mediation)
  library(lme4)
  Med_fit_temp <- lmer(metabolite_follow~protein_base+metabolite_base+age+sex+follow_up_time+factor(phase)+(1|id),data=data_temp)
  Out_fit_temp <- glmer(disease~metabolite_follow+protein_base+age+sex+follow_up_time+factor(phase)+(1|id),data=data_temp,family = binomial(link = "logit"))
  

  set.seed(1234)
  temp_model <- try(Med_protein_metabolite_disease_dis <- mediate(model.m=Med_fit_temp,model.y=Out_fit_temp,treat="protein_base",mediator="metabolite_follow"),silent = T)
  
  if (class(temp_model)=="try-error"){
    
    results_temp <- c(Total_Effect=NA,Total_Effect_P = NA,
                      ADE=NA,ADE_P = NA,
                      ACME=NA,
                      proportions_mediated=NA,ACME_P=NA)
  } else {
    Med_protein_metabolite_disease_dis_summary <- summary(Med_protein_metabolite_disease_dis)
    
    results_temp <- c(Total_Effect=Med_protein_metabolite_disease_dis_summary$tau.coef,Total_Effect_P = Med_protein_metabolite_disease_dis_summary$tau.p,
                      ADE=Med_protein_metabolite_disease_dis_summary$z.avg,ADE_P = Med_protein_metabolite_disease_dis_summary$z.avg.p,
                      ACME=Med_protein_metabolite_disease_dis_summary$d.avg,
                      proportions_mediated=Med_protein_metabolite_disease_dis_summary$n.avg,ACME_P=Med_protein_metabolite_disease_dis_summary$d.avg.p)
  }
  
  results_temp <- c(protein=protein_temp,metabolite=metabolite_temp,disease = disease_temp,N_participants=N_participants,N_samples=N_samples,results_temp)
  
  Mediation_protein_metabolite_disease_discovery <- rbind(Mediation_protein_metabolite_disease_discovery,results_temp)
  
}

Mediation_protein_metabolite_disease_discovery <- as.data.frame(Mediation_protein_metabolite_disease_discovery)
Mediation_protein_metabolite_disease_discovery[,-c(1:3)] <- apply(Mediation_protein_metabolite_disease_discovery[,-c(1:3)],2,as.numeric)
Mediation_protein_metabolite_disease_discovery <- Mediation_protein_metabolite_disease_discovery[!is.na(Mediation_protein_metabolite_disease_discovery$proportions_mediated),]
Mediation_protein_metabolite_disease_discovery$FDR <- p.adjust(Mediation_protein_metabolite_disease_discovery$ACME_P,method="fdr")

Mediation_protein_metabolite_disease_discovery_select <- Mediation_protein_metabolite_disease_discovery[Mediation_protein_metabolite_disease_discovery$FDR < 0.05,]

