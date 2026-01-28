# 1. Data Preprocessing

library(psych)

protein_df <- read.csv("~/Desktop/SurvivalAnalysis/LGG.csv")

# remove proteins with excessive missingness
protein_data_clean <- protein_df[, !(colnames(protein_df) %in% c("ALPHACATENIN", "PARP1"))]

# locate first protein column
start_col_clean <- which(colnames(protein_data_clean) == "X1433BETA")

# isolate only the numeric protein data 
protein_data_final <- protein_data_clean[, start_col_clean:ncol(protein_data_clean)]

# Final check for any NAs. This should now be 0.
print(paste("Total missing values:", sum(is.na(protein_data_final))))


