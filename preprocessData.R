library(tidyverse)

# Loads raw clinical data
clin <- read.csv('./Data/clinical.csv')

# Remove missing stage information
clin <- filter(clin, ajcc_pathologic_stage != 'Stage X')

# Merges sub stages into their parent stages 
clin$stage <- mapvalues(clin$ajcc_pathologic_stage,
                        from = c('Stage I', 'Stage IA', 'Stage IB', 'Stage II',
                                 'Stage IIA', 'Stage IIB', 'Stage III', 'Stage IIIA',
                                 'Stage IIIB', 'Stage IIIC', 'Stage IV'),
                        to = c('Stage I', 'Stage I', 'Stage I', 'Stage II',
                               'Stage II', 'Stage II', 'Stage III', 'Stage III',
                               'Stage III', 'Stage III', 'Stage IV'))

# Finds and assigns molecular sub-types to normal samples
normal_patient_ids <- clin[clin$shortLetterCode == 'NT', ]
corres_tumor_ids <- clin[clin$patient %in% normal_patient_ids$patient, ]
corres_tumor_ids <- corres_tumor_ids[corres_tumor_ids$shortLetterCode == 'TP', ]
corres_tumor_ids <- corres_tumor_ids[!duplicated(corres_tumor_ids$patient), ]
normal_patient_ids <- normal_patient_ids[normal_patient_ids$patient %in% corres_tumor_ids$patient, ]
normal_patient_ids <- normal_patient_ids[order(normal_patient_ids$patient), ]
corres_tumor_ids <- corres_tumor_ids[order(corres_tumor_ids$patient), ]
normal_patient_ids$paper_BRCA_Subtype_PAM50 <- corres_tumor_ids$paper_BRCA_Subtype_PAM50

clin <- rbind(clin, normal_patient_ids)
clin <- clin %>%
  drop_na(paper_BRCA_Subtype_PAM50)

write.csv(clin, './Data/clinical_selected.csv')
