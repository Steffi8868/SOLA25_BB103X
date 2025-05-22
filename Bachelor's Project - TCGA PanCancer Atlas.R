# Load packages
library(readr) # For reading CSV and TSV-files
library(tidyverse) # For data manipulation
library(dplyr) 
library(genefu) # For subtyping of breast cancer samples
data("pam50.robust") # Data related to PAM50 from genefu
library(reshape2) # For using the "melt"-function
library(pheatmap) # For creating heatmaps
library(ggpubr) # For visualization
library(stringr) # For working with strings
library(survival) # For survival analysis
library(survminer)
library(patchwork)

### SUBTYPING -------------------------------------------------------------------------------------------------------------------------------------------------

# Read mRNA data
TCGA_mrna_data <- read_tsv("validation_data_mrna_seq_v2_rsem.txt") # Read the tab separated file
View(TCGA_mrna_data) # Show the data
head(TCGA_mrna_data) # Display the 6 first rows
dim(TCGA_mrna_data) # Check dimensions of the object: 20531 rows (n.o genes + row names) and 1084 columns (n.o patients + column names)

# Create expression matrix
TCGA_mrna_matrix <- TCGA_mrna_data[, -c(1,2)] # Remove the first two columns with Hugo symbols and Entrez ID
head(TCGA_mrna_matrix)
dim(TCGA_mrna_matrix) # Check dimensions of the object: 20531 rows (n.o genes) and 1082 columns (n.o patients)

# Normalize data
TCGA_mrna_zscore <- t(scale(t(TCGA_mrna_matrix), center = TRUE, scale = TRUE))
# Normalize data using Z-score normalization
# The dataframe is transposed bc scale works column-wise --> rows: patients, columns: genes
# center = TRUE --> data is subtracted by mean of each column
# scale = TRUE --> data is divided by stdv od each column
# Output: matrix
# Transpose back --> rows: genes, columns: patients
TCGA_mrna_zscore_tibble <- as_tibble(TCGA_mrna_zscore) # Convert matrix w. Z-scored data to a tibble to observe the results
head(TCGA_mrna_zscore_tibble)
dim(TCGA_mrna_zscore_tibble) # Check dimensions of the object: 20531 rows (n.o genes) and 1082 columns (n.o patients)

# Add the Hugo symbols to the Z-scored data
gene_symbols <- TCGA_mrna_data$Hugo_Symbol
rownames(TCGA_mrna_zscore) <- gene_symbols # The Hugo symbols are set as rownames for the Z-scored data
expression_data_df <- as.data.frame(TCGA_mrna_zscore) # Convert the tibble to a dataframe 
View(expression_data_df)
dim(expression_data_df) # Check dimensions of the object: 20531 rows (n.o genes) and 1082 columns (n.o patients)

gene_info <- data.frame(Genes = TCGA_mrna_data$Hugo_Symbol) # Create data frame with Hugo symbols with column name "Genes"

# Transpose the expression matrix (subtyping function only accepts data with rows: patients and columns: genes)
texp <- t(expression_data_df)
View(as.data.frame(texp)) 
colnames(texp) <- gene_info$Genes # Rename the columns with Hugo symbols

# Subtyping
  # sbt.model: subtyping model (PAM50)
  # data: gene expression data with rows: patients and columns: genes
  # annot: annotation data (Hugo symbols)
  # do.mapping = FALSE since gene names recognized by PAM50-model have been assigned
pam50_predictions <- molecular.subtyping(
  sbt.model = "pam50",
  data = texp,
  annot = gene_info,
  do.mapping = FALSE)

# Display subtypes
subtype <- as.data.frame(pam50_predictions$subtype)
View(subtype)
dim(subtype) # Check dimensions of the object: 1082 rows (n.o patients) and 1 column (subtypes)

# Bar chart with number of patients of each subtype according to genefu
bar_chart_subtypes <- ggplot(subtype, 
                             aes(x=pam50_predictions$subtype))+ 
  geom_bar(stat="count", fill="steelblue")+
  labs(title="Number of Patient per Subtype in TCGA PanCancer Atlas dataset", 
       x="Subtypes", y="Count")
bar_chart_subtypes

# Table with n.o. patients in each subtype and percentage of patients in each subtype
subtype_count <- table(subtype)
subtype_df <- as.data.frame(subtype_count) # Convert the table to a data frame 
colnames(subtype_df) <- c("Subtype", "Number_of_Patients") # Rename columns
total_patients <- sum(subtype_df$Number_of_Patients)
subtype_df$Percentage_of_Total <- (subtype_df$Number_of_Patients / total_patients) * 100 # Calculate the percentage of patients in each subtype
print(subtype_df) # Display the table

# Read clinical data
TCGA_clinical_data_i <- read_tsv("validation_data_clinical_patient.txt") # Read the tab separated file
View(TCGA_clinical_data_i) # Show the data
dim(TCGA_clinical_data_i)
TCGA_clinical_data <- TCGA_clinical_data_i[-c(1:4),] # Remove the first 4 rows that do not contain data
View(TCGA_clinical_data)
dim(TCGA_clinical_data) # Check dimensions of the object: 1084 rows (n.o patients) and 38 columns (each category) --> two more patients in the clinical dataset

TCGA_clinical_data <- as.data.frame(TCGA_clinical_data)

TCGA_clinical_data$'#Patient Identifier' <- paste0(TCGA_clinical_data$'#Patient Identifier', "-01") # Add -01 to the end of Sample IDs to match the mRNA dataset
head(TCGA_clinical_data$'#Patient Identifier')

# Clean patient IDs in case there are spaces or mismatches
mRNA_id_clean <- trimws(toupper(colnames(TCGA_mrna_matrix)))
clinical_ids_clean <- trimws(toupper(TCGA_clinical_data$`#Patient Identifier`))
View(clinical_ids_clean)

# Get the common patient IDs between the mRNA data and the clinical dataset
common_ids <- intersect(mRNA_id_clean, clinical_ids_clean)
View(common_ids)

# Filter the clinical dataset to only keep the patients present in the mRNA dataset
TCGA_clinical_data_filtered <- TCGA_clinical_data[clinical_ids_clean %in% common_ids, ]
View(TCGA_clinical_data_filtered)
dim(TCGA_clinical_data_filtered)

# Bar chart with number of patients of each subtype
bar_chart_clinical <- ggplot(TCGA_clinical_data_filtered, 
                             aes(x= `Subtype`))+ 
  geom_bar(stat="count", fill="steelblue")+
  labs(title="Number of Patient per Subtype According to Clinical Data in TCGA PanCancer Atlas dataset", 
       x="Subtypes", y="Count")
bar_chart_clinical

# Table with n.o. patients in each subtype and percentage of patients in each subtype
subtype_count_2 <- table(TCGA_clinical_data_filtered$Subtype, useNA = 'ifany')
subtype2_df <- as.data.frame(subtype_count_2) # Convert the table to a data frame 
colnames(subtype2_df) <- c("Subtype", "Number of Patients") # Rename columns
total_patients2 <- sum(subtype2_df$`Number of Patients`) # Calculate the percentage of patients in each subtype
subtype2_df$`Percentage of Total` <- (subtype2_df$`Number of Patients` / total_patients2) * 100
print(subtype2_df) # Calculate the percentage of patients in each subtype

# Create contingency table
contingency_table <- table(TCGA_clinical_data_filtered$Subtype, pam50_predictions$subtype, useNA ='ifany') 
contingency_table <- addmargins(contingency_table) # Add sums for rows and columns
print(contingency_table)

### IMMUNE DECONVOLUTION -------------------------------------------------------------------------------------------------------------------------------------------------

# Prepare mRNA data for CIBERSORTx

# Try to add Hugo symbols to all rows
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # Connect to the Ensembl database

conversion <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),     # Conversion from Entrez ID to Hugo Symbol
                    filters = "entrezgene_id", 
                    values = TCGA_mrna_data$Entrez_Gene_Id,  
                    mart = ensembl)
View(conversion)

missing_hugo <- is.na(TCGA_mrna_data$Hugo_Symbol) # Shows which rows have missing Hugo symbols
print(missing_hugo)
sum(missing_hugo) # 13 genes lack Hugo symbol

# Add the Hugo symbols
TCGA_mrna_data$Hugo_Symbol[missing_hugo] <- conversion$hgnc_symbol[match(TCGA_mrna_data$Entrez_Gene_Id[missing_hugo], conversion$entrezgene_id)]
View(TCGA_mrna_data)
missing_hugo2 <- is.na(TCGA_mrna_data$Hugo_Symbol) 
print(missing_hugo2)
sum(missing_hugo2) # Still 13 genes with no Hugo symbol --> these Entrez ID:s have no matching Hugo Symbol
dim(TCGA_mrna_data)

sum(is.na(TCGA_mrna_data)) # Number of NAs in the dataset --> Only NAs in the "Hugo Symbol" column
TCGA_mrna_data_cleaned <- na.omit(TCGA_mrna_data) # Removes rows with NAs
sum(is.na(TCGA_mrna_data_cleaned)) # Controls that the number of NAs in the cleaned dataset is 0
dim(TCGA_mrna_data_cleaned) # Controls the dimensions of the dataset: 20518 genes and 1082 patients --> 13 genes have been removed

# Check for duplicate genes
duplicate_genes <- duplicated(TCGA_mrna_data_cleaned$Hugo_Symbol)
View(duplicate_genes)
sum_duplicate_genes <- sum(duplicate_genes)
print(sum_duplicate_genes) # TRUE = 1 and FALSE = 0 --> sum is 7 --> there are duplicates

df_unique <- TCGA_mrna_data_cleaned %>% # Removes duplicates based on the column 'Hugo_symbols' and only keeps the first occurrence
  distinct(Hugo_Symbol, .keep_all = TRUE)
dim(df_unique) # Controls the dimensions of the dataset: 20511 x 1084

# Save the processed dataset as a TSV-file
#write_tsv(df_unique, "TCGA_mrna_data_cleaned.tsv")

# Load the CIBERSORTx results
TCGA_deconv <- read_tsv("CIBERSORTx_Results_TCGA.txt") # Read the tab separated file
View(TCGA_deconv)
dim(TCGA_deconv) # Dimensions 1083 rows (Entrez ID + patients) and 26 columns (22 celltypes + 4 additional categories)

# Prepare the CIBERSORTx results
filtered_results <- TCGA_deconv[TCGA_deconv$`P-value` < 0.05, ] # Filter out unreliable results
dim(filtered_results) # No unreliable results were found
TCGA_deconv_cleaned <- filtered_results[-1, -((ncol(filtered_results)-2):ncol(filtered_results))] # Remove the first row (Entrez ID) and the three last columns (do not contain cell data)
colnames(TCGA_deconv_cleaned)[colnames(TCGA_deconv_cleaned) == "Mixture"] <- "SampleID" # Change column name from "Mixture" to SampleID
View(TCGA_deconv_cleaned)
dim(TCGA_deconv_cleaned) # 1082 x 23

# Prepare the clinical dataset
TCGA_filtered_data <- subset(TCGA_clinical_data_filtered, !is.na(`Subtype`)) # Removes NAs in the Subtype column
View(TCGA_filtered_data)
colnames(TCGA_filtered_data)[colnames(TCGA_filtered_data) == "#Patient Identifier"] <- "SampleID" # Change column name from "#Patient Identifier" to SampleID
View(TCGA_filtered_data)
dim(TCGA_filtered_data) # 981 x 38

# Merge the CIBERSORTx results with the PAM50-classification of the clinical data, based on Sample ID 
df_combined <- merge(TCGA_deconv_cleaned, TCGA_filtered_data[, c("SampleID", "Subtype")], by = "SampleID")
View(df_combined)# Column with PAM50 subtyping for each patient added
dim(df_combined) # 981 x 24

# Present results in heatmap
selected_cells <- c(               # Abundant immune celltypes according to METABRIC heatmap
  "T cells follicular helper", 
  "Macrophages M1", 
  "T cells CD8", 
  "T cells CD4 memory resting", 
  "Mast cells resting", 
  "Macrophages M2"
)
deconv_matrix <- as.matrix(df_combined[,selected_cells]) # Extract data for selected cells
head(deconv_matrix)
normalized_matrix <- t(scale(t(deconv_matrix))) # Normalize for each cell type

normalized_matrix[normalized_matrix < 0] <- 0 ##force vals för 0

pheatmap(
  normalized_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Immune Cell Abundance (Validation Cohort - TCGA PanCancer Atlas)",
  fontsize_row = 8,
  fontsize_col = 8,
  color = colorRampPalette(c("white", "red"))(100),
  #filename = "TCGA Heatmap.jpeg"
)

# Sets the PAM50-classification to a factor, i.e. a categorical variable
df_combined$Subtype <- as.factor(df_combined$Subtype)

# Reshape the data frame so that the data values are organized in a long data frame format
df_melted <- melt(df_combined, id.vars = "Subtype", measure.vars = colnames(df_combined)[2:23])
# Arguments:
# The dataset to be transformed: df_combined
# The identifier variable with categorical data: Subtype
# The columns that should be melted into the "variable" and "value" columns of the output
View(df_melted)
dim(df_melted) # 21582 rows x 3 columns: PAM50-classification, variable (column names from df_combined), value (the data from the melted columns)

# VIOLON PLOTS - Active celltypes
df_melted_selected <- df_melted %>%      # Extract data for relevant immune cell types
  filter(variable %in% selected_cells)
View(df_melted_selected)
dim(df_melted_selected) # 5886 x 3

subtype_comparisons <- list( c("BRCA_Basal", "BRCA_Her2"), c("BRCA_Basal", "BRCA_LumA"),        # Pairwise comparisons
                             c("BRCA_Basal", "BRCA_LumB"), c("BRCA_Basal", "BRCA_Normal"), 
                             c("BRCA_Her2", "BRCA_LumA"), c("BRCA_Her2", "BRCA_LumB"), 
                             c("BRCA_Her2", "BRCA_Normal"), c("BRCA_LumA", "BRCA_LumB"), 
                             c("BRCA_LumA", "BRCA_Normal"), c("BRCA_LumB", "BRCA_Normal"))

ggplot(df_melted_selected, aes(x = Subtype, y = value, fill = Subtype)) +                              # Violin plot of immune cell fractions grouped by PAM50 subtype
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 0.5, fontface = 'bold') +  # Show p-value from Kruskal-Wallis test
  stat_compare_means(comparisons = subtype_comparisons, label = "p.signif") +                         # Shows significant and non-significant values from Dunn's test
  facet_wrap(~variable, scales = "free_y") +               # Creates separate plots for each immune cell type
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immune Infiltration by PAM50 Subtype (Validation Cohort - TCGA PanCancer Atlas)",
       y = "Immune Cell Fraction",
       x = "PAM50 Subtype"
  )

ggsave("violin_abundant_cells_TCGA.jpeg", width = 12, height = 9, dpi = 600) # Save plot as JPEG-file

# VIOLIN PLOT - Immune Activity Score
anti_tumor_cells <- c("T cells CD8", "T cells follicular helper", "Macrophages M1", "T cells CD4 memory resting") # Stimulate immune response against tumor
pro_tumor_cells  <- c("Macrophages M2", "Mast cells resting") # Suppress immune response against tumor

df_combined$`Anti-tumor cells` <- rowSums(df_combined[, anti_tumor_cells], na.rm = TRUE) # Sum of immune cell abundance for anti-tumor cells for each sample
df_combined$`Pro-tumor cells` <- rowSums(df_combined[, pro_tumor_cells], na.rm = TRUE) # Sum of immune cell abundance for pro-tumor cells for each sample
df_combined$`Immune Activity Score` <- df_combined$`Anti-tumor cells` - df_combined$`Pro-tumor cells` # Immune Activity Score calculated as difference in abundance between anti-tumor cells and pro-tumor cells
View(df_combined)

df_melted_bal <- melt(df_combined, id.vars = "Subtype", measure.vars = colnames(df_combined)[27]) # Reshapes the data to a long format
View(df_melted_bal)

ggplot(df_melted_bal, aes(x = Subtype, y = value, fill = Subtype)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  stat_compare_means(comparisons = subtype_comparisons, label = "p.signif") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 1.3, fontface = 'bold') +
  facet_wrap(~variable, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immune Activity Score by PAM50 Subtype (Validation Cohort - TCGA PanCancer Atlas)",
       y = "Immune Acitivty Score",
       x = "PAM50 Subtype"
  )

ggsave("violin__immune_activity_score_TCGA.jpeg", width = 12, height = 9, dpi = 600) # Save plot as JPEG-file

### SURVIVAL ANALYSIS: Univariate  -------------------------------------------------------------------------------------------------------------------------------------------------

# Combine immune deconvolution results with clinical data
df_merged <- merge(TCGA_deconv_cleaned, TCGA_filtered_data, by = "SampleID") 
dim(df_merged) # 981 x 60
View(df_merged)


selected_cells <- c(              # Abundant immune celltypes according to METABRIC heatmap
  "T cells follicular helper", 
  "Macrophages M1", 
  "T cells CD8", 
  "T cells CD4 memory resting", 
  "Mast cells resting", 
  "Macrophages M2"
)

# Calculate IAS
if (!"Immune Activity Score" %in% colnames(df_merged)) {
  df_merged <- df_merged %>%
    mutate(`Immune Activity Score` = 
             rowSums(across(all_of(c("T cells CD8", 
                                     "T cells follicular helper", 
                                     "Macrophages M1", 
                                     "T cells CD4 memory resting"))), na.rm = TRUE) -
             rowSums(across(all_of(c("Macrophages M2", 
                                     "Mast cells resting"))), na.rm = TRUE))
}

all_variables <- c(selected_cells, "Immune Activity Score")
df_all <- df_merged %>%
  filter(!is.na(Subtype))

# Kaplan-Meier curves for the IAS and 6 active celltypes according to the heatmap (Note: no subtype classification here)

plot_list <- list() # List to store the plots
i <- 1

# Loop over all immune variables (cells + IAS)
for (var in all_variables) {
  # Extract survival data
  df <- df_all %>%
    dplyr::select(`Overall Survival (Months)`, `Overall Survival Status`, !!sym(var)) %>%
    rename(cell_value = !!sym(var)) %>%
    mutate(
      `Overall Survival (Months)` = as.numeric(`Overall Survival (Months)`),
      `Overall Survival Status` = case_when(
        `Overall Survival Status` == "1:DECEASED" ~ 1,                # Change Overall Survival Status to numerical values 
        `Overall Survival Status` == "0:LIVING" ~ 0,
        TRUE ~ NA_real_                                               # If conditions above are not met                                             
      )
    ) %>%
    na.omit()                                   # Remove samples with NAs in survival data

  # Only if enough data  
  if (nrow(df) > 10) {
    
    # Find optimal cutoff
    cutoff <- surv_cutpoint(df,
                            time = "Overall Survival (Months)",
                            event = "Overall Survival Status",
                            variables = "cell_value")
    
    # Categorize samples as high/low based on cutoff
    df$Group <- ifelse(df$cell_value > cutoff$cutpoint$cutpoint, "High", "Low")
    
    # Fit survival model
    surv_obj <- Surv(df$`Overall Survival (Months)`, df$`Overall Survival Status`)
    fit <- survfit(surv_obj ~ Group, data = df)
    
    # Plot Kaplan-Meier curve
    kaplan_curve <- ggsurvplot(fit,
                               data = df,
                               pval = TRUE,
                               conf.int = TRUE,
                               title = paste("All Subtypes -", var),
                               legend.title = var,
                               legend.labs = c("High", "Low"),
                               risk.table = TRUE,
                               xlab = "Time (Months)")
    
    kaplan_curve$plot <- kaplan_curve$plot +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
      )
    
    print(kaplan_curve)
    
    plot_list[[i]] <- kaplan_curve$plot
    i <- i + 1
  }
}


# Kaplan-Meier curves for each subtype for the 6 active celltypes according to the METABRIC heatmap

plot_list <- list() # List to store the plots
i <- 1

for (sub in unique(df_merged$Subtype)) {
  for (cell in selected_cells) {
    
    # Extract data for current subtype
    df <- df_merged %>%
      filter(Subtype == sub) %>%
      dplyr::select(`Overall Survival (Months)`, `Overall Survival Status`, !!sym(cell)) %>%
      rename(cell_value = !!sym(cell)) %>%
      mutate(
        `Overall Survival (Months)` = as.numeric(`Overall Survival (Months)`),
        `Overall Survival Status`= case_when(
          `Overall Survival Status` == "1:DECEASED" ~ 1,                 # Change Overall Survival Status to numerical values 
          `Overall Survival Status` == "0:LIVING" ~ 0,                    
          TRUE ~ NA_real_                                                # If conditions above are not met 
        )
      ) %>% 
      na.omit()                            # Remove samples with NAs in survival data
    
    # Only if enough data
    if (nrow(df) > 10) {
      
      # Find optimal cutoff
      cutoff <- surv_cutpoint(df,
                              time = "Overall Survival (Months)",
                              event = "Overall Survival Status",
                              variables = "cell_value")
      
      # Categorize samples as high/low based on cutoff
      df$Group <- ifelse(df$cell_value > cutoff$cutpoint$cutpoint, "High", "Low")
      View(df)
      
      # Fit survival model
      surv_obj <- Surv(df$`Overall Survival (Months)`, df$`Overall Survival Status`)
      fit <- survfit(surv_obj ~ Group, data = df)
      
      # Kaplan-Meier plot
      kaplan_curve <- ggsurvplot(fit,
                                 data = df,
                                 pval = TRUE,
                                 conf.int = TRUE,
                                 title = paste(sub, "-", cell),
                                 legend.title = cell,
                                 legend.labs = c("High", "Low"),
                                 risk.table = TRUE,
                                 xlab = "Time (Months)")
      
      # Plot Kaplan-Meier curve
      kaplan_curve$plot <- kaplan_curve$plot +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12)
        )
      
      print(kaplan_curve)
      
      plot_list[[i]] <- kaplan_curve$plot
      i <- i + 1
    }
  }
}

# Kaplan-Meier curve for each subtype for the 6 active celltypes according to the heatmap and response to radiotherapy

selected_cells <- c(                       # Abundant immune celltypes according to METABRIC heatmap
  "T cells follicular helper",
  "Macrophages M1",
  "T cells CD8",
  "T cells CD4 memory resting",
  "Mast cells resting",
  "Macrophages M2"
)

anti_tumor_cells <- c("T cells follicular helper", "Macrophages M1", "T cells CD8", "T cells CD4 memory resting") 
pro_tumor_cells  <- c("Mast cells resting", "Macrophages M2") 

# Extract relevant data
df_radio_subtype <- df_merged %>%
  filter(`Radiation Therapy` == "Yes") %>%   # Only include patients that have received radiotherapy
  filter(!is.na(Subtype)) %>%                # Exclude samples with no PAM50-subtype
  mutate(
    `Overall Survival (Months)` = as.numeric(`Overall Survival (Months)`),
    `Overall Survival Status` = case_when(
      `Overall Survival Status` == "1:DECEASED" ~ 1,
      `Overall Survival Status` == "0:LIVING" ~ 0,
      TRUE ~ NA_real_
    ),
    Immune.Activity.Score = rowSums(across(all_of(anti_tumor_cells)), na.rm = TRUE) -            # Calculation of immune activity score
      rowSums(across(all_of(pro_tumor_cells)), na.rm = TRUE)
  ) %>%
  filter(!is.na(`Overall Survival (Months)`), !is.na(`Overall Survival Status`))  # Exclude samples with NAs in survival data

subtypes <- unique(df_radio_subtype$Subtype) # The PAM50-subtypes
plot_list <- list() # List to store the plots
i <- 1

for (sub in subtypes) {
  df_subtype <- df_radio_subtype %>% filter(Subtype == sub)
  
  for (cell in c(selected_cells, "Immune.Activity.Score")) {
    # Extract relevant data
    col_name <- ifelse(cell == "Immune.Activity.Score", "Immune.Activity.Score", cell)
    
    df_cell <- df_subtype %>%
      dplyr::select(`Overall Survival (Months)`, `Overall Survival Status`, !!sym(col_name)) %>%
      rename(cell_value = !!sym(col_name)) %>%
      na.omit()
    
    # Only if enough data
    if (nrow(df_cell) > 10) {
      
      # Optimal cutoff
      cutoff <- surv_cutpoint(
        df_cell,
        time = "Overall Survival (Months)",
        event = "Overall Survival Status",
        variables = "cell_value"
      )
      
      df_cell$Group <- ifelse(df_cell$cell_value > cutoff$cutpoint$cutpoint, "High", "Low")
      
      # Fit survival model
      surv_obj <- Surv(df_cell$`Overall Survival (Months)`, df_cell$`Overall Survival Status`)
      fit <- survfit(surv_obj ~ Group, data = df_cell)
      
      # Kaplan-Meier plot
      kaplan_curve <- ggsurvplot(
        fit,
        data = df_cell,
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        title = paste("Subtype:", sub, "\n", ifelse(cell == "Immune.Activity.Score", "Immune Activity Score", cell)),
        legend.title = cell,
        legend.labs = c("High", "Low"),
        xlab = "Time (Months)"
      )
     
      print(kaplan_curve)
      plot_list[[i]] <- kaplan_curve$plot
      i <- i + 1
    }
  }
}


### SURVIVAL ANALYSIS: Multivariate -------------------------------------------------------------------------------------------------------------------------------------------------

# Select relevant clinical variables
df_TCGA_selct <- TCGA_filtered_data %>% 
  dplyr::select (`SampleID`,`Subtype`,`Diagnosis Age`,
                 `Radiation Therapy`,`Overall Survival Status`,
                 `Overall Survival (Months)`,
                 `Neoplasm Disease Stage American Joint Committee on Cancer Code`,
                 `Primary Lymph Node Presentation Assessment`)
view(df_TCGA_selct)

#Select relevant immune cells
df_immune_selct <- df_combined %>% 
  dplyr::select(`Macrophages M2`,`Mast cells resting`,
                `T cells CD4 memory resting`,`Macrophages M1`,
                `Immune Activity Score`,`SampleID`)
view(df_immune_selct)

# Merge the clinical data with the immune deconvolution data
df_cox_val <- merge(df_immune_selct, df_TCGA_selct, by = "SampleID")
view(df_cox_val)

# Rename variables for convenience
df_cox_val <- df_cox_val %>% dplyr::rename(OS_months = `Overall Survival (Months)`,
                                           OS_status = `Overall Survival Status`, 
                                           `Neoplasm Disease Stage`=`Neoplasm Disease Stage American Joint Committee on Cancer Code`)
view(df_cox_val)

# Change the Time and State varibles to numeric.
 # Time varible = OS_months
 # State varible = OS_status

df_cox_val$OS_months <- as.numeric(df_cox_val$OS_months)
df_cox_val$OS_status <- as.numeric(gsub(".*1.*", 1, gsub(".*0.*", 0, df_cox_val$OS_status)))

# Check that the changes worked
str(df_cox_val$OS_months)
str(df_cox_val$OS_status)

# Regroup the Neoplasm Disease Stage American Joint Committee on Cancer Code
 # This will reduce the amount of levels and present a more reliable result
df_cox_val$`Grouped Neoplasm Disease Stage` <- with(df_cox_val,
                                                    ifelse(`Neoplasm Disease Stage` %in% c("STAGE I", "STAGE IA", "STAGE IB"),
                                                           "STAGE 1",
                                                           ifelse(`Neoplasm Disease Stage` %in% c("STAGE IIA", "STAGE IIB", "STAGE II"),
                                                                  "STAGE 2",
                                                                  ifelse(`Neoplasm Disease Stage` %in% c("STAGE III", "STAGE IIIA", "STAGE IIIB", "STAGE IIIC"),
                                                                         "STAGE 3",
                                                                         "Other"))))
view(df_cox_val)

# Change variables with multiple levels to factors.
df_cox_val$Subtype <- as.factor(df_cox_val$Subtype)
df_cox_val$`Radiation Therapy`<- as.factor(df_cox_val$`Radiation Therapy`)
df_cox_val$`Grouped Neoplasm Disease Stage` <- as.factor(df_cox_val$`Grouped Neoplasm Disease Stage`)

# Change other variables from character to numeric.
df_cox_val$`Diagnosis Age` <- as.numeric(df_cox_val$`Diagnosis Age`)

# Remove all rows that has an NA value
df_cox_val <- df_cox_val[complete.cases(df_cox_val), ]

# Establish desired reference levels for all the factor variables.
df_cox_val$Subtype <- relevel(as.factor(df_cox_val$Subtype), ref = "BRCA_LumA")
df_cox_val$`Radiation Therapy` <- relevel(as.factor(df_cox_val$`Radiation Therapy`), ref = "No")
df_cox_val$`Grouped Neoplasm Disease Stage` <- relevel(as.factor(df_cox_val$`Grouped Neoplasm Disease Stage`), ref = "STAGE 1")
df_cox_val$`Primary Lymph Node Presentation Assessment` <- relevel(as.factor(df_cox_val$`Primary Lymph Node Presentation Assessment`), ref = "No")

cox_result <- coxph(Surv(OS_months, OS_status) ~ `Macrophages M2`+
                      `Mast cells resting`+`T cells CD4 memory resting`+
                      `Macrophages M1`+`Immune Activity Score`+`Subtype`+
                      `Diagnosis Age`+`Radiation Therapy`+
                      `Primary Lymph Node Presentation Assessment`+
                      `Grouped Neoplasm Disease Stage`,
                    data = df_cox_val)
summary(cox_result)

# Generate Forestplot
ggforest(cox_result, data = df_cox_val, fontsize = 0.5,
         main = "Hazard Ratios from Cox Proportional Hazards Model (Validation)",
         refLabel = "reference")











