

library(readxl)
library(dplyr)
library(openxlsx)


# PALEOZOIC DATA
###################################################################################

# Define the path to your Excel file
file_path <- "Phan/PhanData/assignLatLonSite/paleozoicAssignSites.xlsx"

# Read data from Sheet 1 and Sheet 2
sheet1 <- read_excel(file_path, sheet = "paleozoic comp")
sheet2 <- read_excel(file_path, sheet = "site data by study")

# Join sheet1 with sheet2 using source and study as keys
merged_data_paleozoic <- sheet1 %>%
  left_join(sheet2 %>% 
              select(study, `assigned location`, `assigned lat`, `assigned lon`, category),
            by = c("source" = "study", "assigned location" = "assigned location"))

write.xlsx(merged_data_paleozoic, file = "Phan/PhanData/assignLatLonSite/updated_paleozoicAssignSites.xlsx", sheetName = "paleozoic comp updated", overwrite = TRUE)



# NON-PALEOZOIC DATA
###################################################################################

file_path <- "Phan/PhanData/assignLatLonSite/nonpaleozoicAssignSites.xlsx"

# Read the two sheets
comp_data <- read_excel(file_path, sheet = "nonpaleozoic_comp")
site_data <- read_excel(file_path, sheet = "site data by study")

# Join based on matching source → study AND location → assigned location
merged_data_nonpaleozoic <- comp_data %>%
  left_join(site_data %>%
              select(study, `assigned location`, `assigned lat`, `assigned lon`, category),
            by = c("source" = "study", "assigned location" = "assigned location"))

write.xlsx(merged_data_nonpaleozoic, file = "Phan/PhanData/assignLatLonSite/updated_nonpaleozoicAssignSites.xlsx", sheetName = "nonpaleozoic_comp", overwrite = TRUE)


# COMBINE
###################################################################################

df1 <- as.data.frame(merged_data_nonpaleozoic)
df2 <- as.data.frame(merged_data_paleozoic)

# Reorder columns of df2 to match df1
df2 <- df2[, colnames(df1)]

# Combine them
combined <- rbind(df1, df2)

# Convert back to matrix if needed
PhanCompUpdated <- as.matrix(combined)

# Save the result to a new Excel file
write.csv(PhanCompUpdated, file = "Phan/PhanData/assignLatLonSite/PhanCompUpdated.csv")


