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

# Read mRNA-Seq data
metabric_mrna_data <- read_tsv("metabric_data_mrna_illumina_microarray.txt") # Read the tab separated file
View(metabric_mrna_data) # Show the data
head(metabric_mrna_data) # Display the 6 first rows
dim(metabric_mrna_data) # Check dimensions of the object: 20603 rows (n.o genes + row names) and 1982 columns (n.o patients + column names)

# Create expression matrix
metabric_mrna_matrix <- metabric_mrna_data[, -c(1,2)] # Remove the first two columns with Hugo symbols and Entrez ID
head(metabric_mrna_matrix)
dim(metabric_mrna_matrix) # Check dimensions of the object: 20603 rows (n.o genes) and 1980 columns (n.o patients)

# Normalize data
metabric_mrna_zscore <- t(scale(t(metabric_mrna_matrix), center = TRUE, scale = TRUE))
# Normalize data using Z-score normalization
# The dataframe is transposed bc scale works column-wise --> rows: patients, columns: genes
# center = TRUE --> data is subtracted by mean of each column
# scale = TRUE --> data is divided by stdv of each column
# Output: matrix
# Transpose back --> rows: genes, columns: patients
metabric_mrna_zscore_tibble <- as_tibble(metabric_mrna_zscore) # Convert matrix w. Z-scored data to a tibble to observe the results
head(metabric_mrna_zscore_tibble)
dim(metabric_mrna_zscore_tibble) # Check dimensions of the object: 20603 rows (n.o genes) and 1980 columns (n.o patients)

# Add the Hugo symbols to the Z-scored data
gene_symbols <- metabric_mrna_data$Hugo_Symbol
rownames(metabric_mrna_zscore) <- gene_symbols # The Hugo symbols are set as rownames for the Z-scored data
expression_data_df <- as.data.frame(metabric_mrna_zscore) # Convert the tibble to a dataframe 
View(expression_data_df)
dim(expression_data_df) # Check dimensions of the object: 20603 rows (n.o genes) and 1980 columns (n.o patients)

gene_info <- data.frame(Genes = metabric_mrna_data$Hugo_Symbol) # Create data frame with Hugo symbols with column name "Genes"

# Transpose the expression matrix (subtyping function only accepts data with rows: patients and columns: genes)
texp <- t(expression_data_df)
colnames(texp) <- gene_info$Genes # Rename the columns with Hugo symbols
View(as.data.frame(texp))

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
dim(subtype) # Check dimensions of the object: 1980 rows (n.o patients) and 1 column (subtypes)

# Bar chart with number of patients of each subtype according to genefu
bar_chart_subtypes <- ggplot(subtype, 
                             aes(x=pam50_predictions$subtype))+ 
  geom_bar(stat="count", fill="steelblue")+
  labs(title="Number of Patient per Subtype in METABRIC dataset", 
       x="Subtypes", y="Count")
bar_chart_subtypes

# Table with n.o. patients in each subtype and percentage of patients in each subtype
subtype_count <- table(subtype)
subtype_df <- as.data.frame(subtype_count) # Convert the table to a data frame 
colnames(subtype_df) <- c("Subtype", "Number of Patients") # Rename columns
total_patients <- sum(subtype_df$`Number of Patients`)
subtype_df$`Percentage of Total` <- (subtype_df$`Number of Patients` / total_patients) * 100 # Calculate the percentage of patients in each subtype
print(subtype_df) # Display the table

# Read clinical data
metabric_clinical_data_i <- read_tsv("metabric_data_clinical_patient.txt") # Read the tab separated file
View(metabric_clinical_data_i) # Show the data
metabric_clinical_data <- metabric_clinical_data_i[-c(1:4),] # Remove the first 4 rows that do not contain data
View(metabric_clinical_data)
dim(metabric_clinical_data) # Check dimensions of the object: 2509 rows (n.o patients) and 24 columns (each category) --> Note: More patients compared to mRNA dataset

metabric_clinical_data <- as.data.frame(metabric_clinical_data)
colnames(metabric_clinical_data)[colnames(metabric_clinical_data) == "Pam50 + Claudin-low subtype"] <- "PAM50" # Change column name to PAM50
View(metabric_clinical_data)

# Bar chart with number of patients of each subtype according to clinical data
bar_chart_clinical <- ggplot(metabric_clinical_data, 
                             aes(x= `PAM50`))+ 
  geom_bar(stat="count", fill="steelblue")+
  labs(title="Number of Patient per Subtype in METABRIC Clinical Dataset", 
       x="Subtypes", y="Count")
bar_chart_clinical

# Table with n.o patients in each subtype and percentage of patients in each subtype
subtype_count_2 <- table(metabric_clinical_data$PAM50, useNA = 'ifany') # N.o patients per subtype
subtype2_df <- as.data.frame(subtype_count_2) # Concert to dataframe
colnames(subtype2_df) <- c("Subtype", "Number of Patients") # Rename columns
total_patients2 <- sum(subtype2_df$`Number of Patients`)
subtype2_df$`Percentage of Total` <- (subtype2_df$`Number of Patients` / total_patients2) * 100 # Percentage of patients in each subtype
print(subtype2_df) # Display the table

# Bar chart with number of patients of each subtype according to clinical data (NAs excluded)
metabric_clinical_data_cleaned <- metabric_clinical_data %>% filter(!is.na(`PAM50`)) # Exclude NAs in PAM50 column
bar_chart_clinical <- ggplot(metabric_clinical_data_cleaned, 
                             aes(x= `PAM50`))+ 
  geom_bar(stat="count", fill="steelblue")+
  labs(title="Number of Patient per Subtype in METABRIC Clinical Dataset", 
       x="Subtypes", y="Count")
bar_chart_clinical

# Table with n.o. patients in each subtype and percentage of patients in each subtype (NAs excluded)
subtype_count_3 <- table(metabric_clinical_data$PAM50) # N.o. patients
subtype3_df <- as.data.frame(subtype_count_3) # Convert the table to a data frame
colnames(subtype3_df) <- c("Subtype", "Number of Patients") # Rename columns
total_patients3 <- sum(subtype3_df$`Number of Patients`)
subtype3_df$`Percentage of Total` <- (subtype3_df$`Number of Patients` / total_patients3) * 100 # Calculate the percentage of patients in each subtype
print(subtype3_df) # Display the table

# Contingency table to compare the classifications
metabric_clinical_data_no_NA <- na.omit(metabric_clinical_data$PAM50) # Remove NAs since both vectors need to be of equal length
metabric_clinical_data_no_NA_df <- as.data.frame(metabric_clinical_data_no_NA) # Convert dataframe to check dimension --> 1980 patients
dim(metabric_clinical_data_no_NA_df)
View(metabric_clinical_data_no_NA_df)

contingency_table <- table(metabric_clinical_data_no_NA, pam50_predictions$subtype) # Create contingency table¨
contingency_table <- addmargins(contingency_table) # Add sums for rows and columns
print(contingency_table)

### IMMUNE DECONVOLUTION -------------------------------------------------------------------------------------------------------------------------------------------------

# Prepare mRNA data for CIBERSORTx
sum(is.na(metabric_mrna_data)) # Number of NAs in the dataset
metabric_mrna_data_cleaned <- metabric_mrna_data[, colSums(is.na(metabric_mrna_data)) == 0] # Removes columns with NAs
sum(is.na(metabric_mrna_data_cleaned)) # Controls that the number of NAs in the cleaned dataset is 0
dim(metabric_mrna_data_cleaned) # Dimensions of the dataset: 20603 rows (genes) and 1966 columns (patients)

df_unique <- metabric_mrna_data_cleaned %>% # Removes duplicates based on the column 'Hugo_symbols' and only keeps the first occurrence
  distinct(Hugo_Symbol, .keep_all = TRUE)
View(df_unique)
dim(df_unique) # Controls the dimensions of the dataset: 20385 genes and 1966 patients

# Check for duplicate genes
duplicate_genes <- duplicated(df_unique$Hugo_Symbol)
sum_duplicate_genes <- sum(duplicate_genes)
print(sum_duplicate_genes) # TRUE = 1 and FALSE = 0 --> Sum is 0 --> no duplicates

# Save the processed dataset as a TSV-file
#write_tsv(df_unique, "metabric_mRNA_cleaned.tsv")

# Load the CIBERSORTx results
metabric_deconv <- read_tsv("CIBERSORTx_Results_METABRIC.txt") # Read the tab separated file
View(metabric_deconv)
dim(metabric_deconv) # Dimensions 1965 rows (Entrez ID + patients) and 26 columns (22 celltypes + 4 additional categories)
head(metabric_deconv)

# Prepare the CIBERSORTx results
filtered_results <- metabric_deconv[metabric_deconv$`P-value` < 0.05, ] # Filter out unreliable results
dim(filtered_results) # Dimension: 1965 rows and 26 columns --> no unreliable results were found

metabric_deconv_cleaned <- filtered_results[-1, -((ncol(filtered_results)-2):ncol(filtered_results))] # Remove the first row (Entrez ID) and the last three columns (do not contain cell data)
View(metabric_deconv_cleaned)
dim(metabric_deconv_cleaned) # 1964 x 23
colnames(metabric_deconv_cleaned)[colnames(metabric_deconv_cleaned) == "Mixture"] <- "SampleID" # Change column name from "Mixture" to SampleID
View(metabric_deconv_cleaned)

# Prepare the clinical dataset
filtered_data <- subset(metabric_clinical_data, !PAM50 %in% c("claudin-low", "NC") & !is.na(PAM50)) # Removes claudin-low and NC from dataframe with clinical data as well as NAs in the PAM50 column
View(filtered_data)
colnames(filtered_data)[colnames(filtered_data) == "#Patient Identifier"] <- "SampleID" # Change column name from "#Patient Identifier" to SampleID
View(filtered_data)

# Merge the CIBERSORTx results with the PAM50-classification of the clinical data, based on Sample ID 
df_combined <- merge(metabric_deconv_cleaned, filtered_data[, c("SampleID", "PAM50")], by = "SampleID")
View(df_combined)# Column with PAM50 subtyping for each patient added
dim(df_combined) # 1744 x 24

# Heatmap of immune cell abundance
deconv_matrix <- as.matrix(df_combined[2:23]) # Exclude first column since it contains SampleID
head(deconv_matrix)
normalized_matrix <- t(scale(t(deconv_matrix))) # Normalize for each cell type

breaks <- seq(min(normalized_matrix), max(normalized_matrix), length.out = 100) # Define breaks for color scale
midpoint <- which.min(abs(breaks)) # Value closest to zero
colors <- c(
  colorRampPalette(c("white", "white"))(midpoint),  # Values below zero are white
  colorRampPalette(c("white", "red"))(100 - midpoint)    # Values above zero are red 
)

# Create heatmap
pheatmap(
  normalized_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colors,
  breaks = breaks,
  main = "Immune Cell Abundance Heatmap (METABRIC)",
  fontsize_row = 6,
  fontsize_col = 6,
  #filename = "Heatmap METABRIC.jpeg"                   # Remove comment to save as JPEG
)
 
# Sets the PAM50-classification to a factor, i.e. a categorical variable
df_combined$PAM50 <- as.factor(df_combined$PAM50)

# Reshape the data frame so that the data values are organized in a long data frame format
df_melted <- melt(df_combined, id.vars = "PAM50", measure.vars = colnames(df_combined)[2:23])
# Arguments:
# The dataset to be transformed: df_combined
# The identifier variable with categorical data: PAM50
# The columns that should be melted into the "variable" and "value" columns of the output
View(df_melted)
dim(df_melted) # 38368 rows x 3 columns: PAM50-classification, variable (column names from df_combined), value (the data from the melted columns)

# VIOLIN PLOT - Abundant celltypes
selected_cells <- c("Macrophages M2", "T cells CD8", "T cells follicular helper", 
                    "Macrophages M1", "T cells CD4 memory resting", "Mast cells resting") # Abundant immune cells acc. to heatmap

df_melted_selected <- df_melted %>%      # Extract data for relevant immune cell types
  filter(variable %in% selected_cells)
View(df_melted_selected)
dim(df_melted_selected) # 10464 x 3

subtype_comparisons <- list( c("Basal", "Her2"), c("Basal", "LumA"), c("Basal", "LumB"), # Pairwise comparisons
                             c("Basal", "Normal"), c("Her2", "LumA"), 
                             c("Her2", "LumB"), c("Her2", "Normal"), c("LumA", "LumB"), 
                             c("LumA", "Normal"), c("LumB", "Normal"))

ggplot(df_melted_selected, aes(x = PAM50, y = value, fill = PAM50)) +   # Violin plot of immune cell fractions grouped by PAM50 subtype
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 0.5, fontface = 'bold') + # Show p-value from Kruskal-Wallis test
  stat_compare_means(comparisons = subtype_comparisons, method = "wilcox.test", label = "p.signif") +  # Shows significant and non-significant values from pairwise comparison
  facet_wrap(~variable, scales = "free_y") +         # Creates separate plots for each immune cell type
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immune Infiltration by PAM50 Subtype (METABRIC)",
       y = "Immune Cell Abundance",
       x = "PAM50 Subtype"
  )

ggsave("violin_abundant_cells.jpeg", width = 12, height = 9, dpi = 600) # Save plot as JPEG-file

# VIOLIN PLOT - Immune Activity Score
anti_tumor_cells <- c("T cells CD8", "T cells follicular helper", "Macrophages M1", "T cells CD4 memory resting") # Stimulate immune response against tumor
pro_tumor_cells  <- c("Macrophages M2", "Mast cells resting") # Suppress immune response against tumor

df_combined$`Anti-tumor cells` <- rowSums(df_combined[, anti_tumor_cells], na.rm = TRUE) # Sum of immune cell abundance for anti-tumor cells for each sample
df_combined$`Pro-tumor cells` <- rowSums(df_combined[, pro_tumor_cells], na.rm = TRUE) # Sum of immune cell abundance for pro-tumor cells for each sample
df_combined$`Immune Activity Score` <- df_combined$`Anti-tumor cells` - df_combined$`Pro-tumor cells` # Immune Activity Score calculated as difference in abundance between anti-tumor cells and pro-tumor cells
View(df_combined)

df_melted_bal <- melt(df_combined, id.vars = "PAM50", measure.vars = colnames(df_combined)[27]) # Reshapes the data to a long format
View(df_melted_bal)

ggplot(df_melted_bal, aes(x = PAM50, y = value, fill = PAM50)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  stat_compare_means(comparisons = subtype_comparisons, label = "p.signif") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 1.3, fontface = 'bold') +
  facet_wrap(~variable, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immune Activity Score by PAM50 Subtype (METABRIC)",
       y = "Immune Acitivty Score",
       x = "PAM50 Subtype"
  )

ggsave("violin__immune_activity_score.jpeg", width = 12, height = 9, dpi = 600) # Save plot as JPEG-file

### SURVIVAL ANALYSIS: Univariate  -------------------------------------------------------------------------------------------------------------------------------------------------

# Combine immune deconvolution results with clinical data
df_merged <- merge(metabric_deconv_cleaned, filtered_data, by = "SampleID") 
dim(df_merged) # 1744 x 46
View(df_merged)

# Create Immune Activity Score
df_merged$`Anti-tumor cells` <- rowSums(df_merged[, anti_tumor_cells], na.rm = TRUE) # Sum of immune cell abundance for anti-tumor cells for each sample
df_merged$`Pro-tumor cells` <- rowSums(df_merged[, pro_tumor_cells], na.rm = TRUE) # Sum of immune cell abundance for pro-tumor cells for each sample
df_merged$`Immune Activity Score` <- df_merged$`Anti-tumor cells` - df_merged$`Pro-tumor cells` # Immune Activity Score calculated as difference in abundance between anti-tumor cells and pro-tumor cells
View(df_merged)

# Kaplan-Meier curve for IAS (Note: no subtype classification here)

# Extract survival data
df <- df_merged %>%
  dplyr::select("Overall Survival (Months)", "Overall Survival Status", "Immune Activity Score") %>%
  mutate(`Overall Survival (Months)` = as.numeric(`Overall Survival (Months)`)) %>%
  mutate(`Overall Survival Status` = case_when(
    `Overall Survival Status` %in% c("1:DECEASED") ~ 1,        # Change Overall Survival Status to numerical values    
    `Overall Survival Status` %in% c("0:LIVING") ~ 0)) %>% 
  na.omit()                                                    # Remove samples with NAs in survival data
View(df)
dim(df)

# Find optimal cutoff
cutoff <- surv_cutpoint(df,
                        time = "Overall Survival (Months)",
                        event = "Overall Survival Status",
                        variables = "Immune Activity Score")

# Categorize samples as high/low based on cutoff
df$Group <- ifelse(df$`Immune Activity Score` > cutoff$cutpoint$cutpoint, "High", "Low")

# Fit survival model
surv_obj <- Surv(df$`Overall Survival (Months)`, df$`Overall Survival Status`)
fit <- survfit(surv_obj ~ Group, data = df)

# Plot Kaplan-Meier curve
kaplan_curve <- ggsurvplot(fit, 
                           data = df, 
                           pval = TRUE,
                           conf.int = TRUE,
                           title = "Kaplan-Meier: Immune Activity Score (METABRIC)", 
                           legend.title = "Immune Activity Score",
                           legend.labs = c("High", "Low"), 
                           risk.table = TRUE, 
                           xlab = "Time (Months)")
kaplan_curve$plot <- kaplan_curve$plot + 
  theme(
    plot.title = element_text(size = 14),
    axis.title.x = element_text(size = 11), 
    axis.title.y = element_text(size = 11)   
  )

print(kaplan_curve)
ggsave("KM_Immune_Activity_Score_METABRIC.png", plot = kaplan_curve$plot, width = 8, height = 6)

# Kaplan-Meier curves for the 6 active celltypes according to the heatmap (Note: no subtype classification here)

plot_list <- list() # List to store the plots
i <- 1

for (cell in selected_cells) {
  df <- df_merged %>%
    dplyr::select("Overall Survival (Months)", "Overall Survival Status", !!sym(cell)) %>%
    rename(cell_value = !!sym(cell)) %>%
    mutate(`Overall Survival (Months)` = as.numeric(`Overall Survival (Months)`)) %>%
    mutate(`Overall Survival Status` = case_when(
      `Overall Survival Status` %in% c("1:DECEASED") ~ 1,             # Change Overall Survival Status to numerical values  
      `Overall Survival Status` %in% c("0:LIVING") ~ 0)) %>%          
    na.omit()                                                         # Remove samples with NAs in survival data
    
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
                              title = paste(cell), 
                              legend.title = cell,
                              legend.labs = c("High", "Low"), 
                              risk.table = TRUE, 
                              xlab = "Time (months)")
  kaplan_curve$plot <- kaplan_curve$plot + 
    theme(
      plot.title = element_text(size = 14),
      axis.title.x = element_text(size = 11),  
      axis.title.y = element_text(size = 11)   
    )
    
  print(kaplan_curve)
    
  plot_list[[i]] <- kaplan_curve$plot
  i <- i + 1
}

wrap_plots(plotlist = plot_list, ncol = length(selected_cells)) # All plots in the same figure

# Kaplan-Meier curves for each subtype for the 6 active celltypes according to the heatmap

plot_list <- list() # List to store the plots
i <- 1

for (sub in unique(df_merged$PAM50)) {
  for (cell in selected_cells) {
    
    # Extract data for current subtype
    df <- df_merged %>%
      filter(PAM50 == sub) %>%
      dplyr::select("Overall Survival (Months)", "Overall Survival Status", !!sym(cell)) %>%
      rename(cell_value = !!sym(cell)) %>%
      mutate(`Overall Survival (Months)` = as.numeric(`Overall Survival (Months)`)) %>%
      mutate(`Overall Survival Status` = case_when(
        `Overall Survival Status` %in% c("1:DECEASED") ~ 1,           # Change Overall Survival Status to numerical values 
        `Overall Survival Status` %in% c("0:LIVING") ~ 0)) %>% 
      na.omit()                                                       # Remove samples with NAs in survival data
    
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
                               title = paste(sub, "-", cell), 
                               legend.title = cell,
                               legend.labs = c("High", "Low"), 
                               risk.table = TRUE, 
                               xlab = "Time (months)")
    kaplan_curve$plot <- kaplan_curve$plot + 
      theme(
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 11),  
        axis.title.y = element_text(size = 11)   
      )
    
    print(kaplan_curve)
    
    plot_list[[i]] <- kaplan_curve$plot
    i <- i + 1
    
  }
}

# Kaplan-Meier curve for each subtype for the 6 active celltypes according to the heatmap and response to radiotherapy

treatment_list = treatment_list <- c("Radio Therapy")  # List with selected therapies

plot_list <- list()
i <- 1

for (sub in unique(df_merged$PAM50)) {
  for (cell in selected_cells) {
    for (treatment in treatment_list) {
      # Extract data for current subtype
      df <- df_merged %>%
        filter(PAM50 == sub) %>%
        filter(.data[[treatment]] == "YES") %>%
        dplyr::select("Overall Survival (Months)", "Overall Survival Status", !!sym(cell)) %>%
        rename(cell_value = !!sym(cell)) %>%
        mutate(`Overall Survival (Months)` = as.numeric(`Overall Survival (Months)`)) %>%
        mutate(`Overall Survival Status` = case_when(
          `Overall Survival Status` %in% c("1:DECEASED") ~ 1,
          `Overall Survival Status` %in% c("0:LIVING") ~ 0)) %>% 
        na.omit()
      
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
                                 title = paste(sub, "-", cell, "\n-", treatment), 
                                 legend.title = cell, 
                                 legend.labs = c("High", "Low"), 
                                 risk.table = TRUE, xlab = "Time (months)")
      
      kaplan_curve$plot <- kaplan_curve$plot + 
        theme(
          plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 11),  # for x-axis title
          axis.title.y = element_text(size = 11)   # for y-axis title
        )
      print(kaplan_curve)
      plot_list[[i]] <- kaplan_curve$plot
      i <- i + 1
      
    }
  }
  plots <- wrap_plots(plotlist = plot_list, ncol = length(selected_cells))
  print(plots)
  plot_list <- list()
  i <- 1
}

# Kaplan-Meier curve for each subtype for the IAS

score_list <- c("Immune Activity Score") # List with selected scores

plot_list <- list()
i <- 1

for (sub in unique(df_merged$PAM50)) {
  for (score in score_list) {
    for (treatment in treatment_list) {
      # Extract data for current subtype
      df <- df_merged %>%
        filter(PAM50 == sub) %>%
        filter(.data[[treatment]] == "YES") %>%
        dplyr::select("Overall Survival (Months)", "Overall Survival Status", !!sym(score)) %>%
        rename(score_value = !!sym(score)) %>%
        mutate(`Overall Survival (Months)` = as.numeric(`Overall Survival (Months)`)) %>%
        mutate(`Overall Survival Status` = case_when(
          `Overall Survival Status` %in% c("1:DECEASED") ~ 1,
          `Overall Survival Status` %in% c("0:LIVING") ~ 0)) %>% 
        na.omit()
      
      # Find optimal cutoff
      cutoff <- surv_cutpoint(df,
                              time = "Overall Survival (Months)",
                              event = "Overall Survival Status",
                              variables = "score_value")
      
      # Categorize samples as high/low based on cutoff
      df$Group <- ifelse(df$score_value > cutoff$cutpoint$cutpoint, "High", "Low")
      
      # Fit survival model
      surv_obj <- Surv(df$`Overall Survival (Months)`, df$`Overall Survival Status`)
      fit <- survfit(surv_obj ~ Group, data = df)
      
      # Plot Kaplan-Meier curve
      kaplan_curve <- ggsurvplot(fit, 
                                 data = df, 
                                 pval = TRUE,
                                 conf.int = TRUE,
                                 title = paste(sub, "-", score, "\n-", treatment), 
                                 legend.title = score, 
                                 legend.labs = c("High", "Low"), 
                                 risk.table = TRUE, xlab = "Time (months)")
      
      kaplan_curve$plot <- kaplan_curve$plot + 
        theme(
          plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 11),  # for x-axis title
          axis.title.y = element_text(size = 11)   # for y-axis title
        )
      
      plot_list[[i]] <- kaplan_curve$plot
      i <- i + 1
    }
  }
  plots <- wrap_plots(plotlist = plot_list, ncol = length(score_list))
  print(plots)
  plot_list <- list()
  i <- 1
}

### SURVIVAL ANALYSIS: Multivariate -------------------------------------------------------------------------------------------------------------------------------------------------

# Extract relevant data from df_merged dataframe: 
    # Deconvolution data for 6 most abundant immune cells
    # Clinical characteristics
    # Survival data
    # Subtyping
    # Response to therapy
    # IAS
df_cox <- df_merged %>% dplyr::select(`Immune Activity Score`,
                                      `Lymph nodes examined positive`,
                                      `Nottingham prognostic index`,`Cellularity`,
                                      `ER status measured by IHC`,
                                      `HER2 status measured by SNP6`,
                                      `Hormone Therapy`,`Age at Diagnosis`,
                                      `Overall Survival (Months)`,
                                      `Overall Survival Status`,`PAM50`,
                                      `Radio Therapy`)

# Rename variables for convenience
df_cox <- df_cox %>%
  dplyr::rename(OS_months = `Overall Survival (Months)`,
                OS_status = `Overall Survival Status`)

# Change the Time and State variables to numeric.
    # Time varible = OS_months
    # State varible = OS_status
df_cox$OS_months <- as.numeric(df_cox$OS_months)
df_cox$OS_status <- as.numeric(gsub(".*1.*", 1, gsub(".*0.*", 0, df_cox$OS_status)))

# Check that the changes worked
str(df_cox$OS_months)
str(df_cox$OS_status)

# Change variables with multiple levels to factors.
df_cox$Cellularity <- as.factor(df_cox$Cellularity)
df_cox$`ER status measured by IHC` <- as.factor(df_cox$`ER status measured by IHC`)
df_cox$`HER2 status measured by SNP6` <- as.factor(df_cox$`HER2 status measured by SNP6`)
df_cox$`Hormone Therapy` <- as.factor(df_cox$`Hormone Therapy`)
df_cox$PAM50 <- as.factor(df_cox$PAM50)
df_cox$`Radio Therapy` <- as.factor(df_cox$`Radio Therapy`)

# Change other variables from character to numeric.
df_cox$`Lymph nodes examined positive` <- as.numeric(df_cox$`Lymph nodes examined positive`)
df_cox$`Nottingham prognostic index` <- as.numeric(df_cox$`Nottingham prognostic index`)
df_cox$`Age at Diagnosis` <- as.numeric(df_cox$`Age at Diagnosis`)

# Remove all rows that has an NA value
df_cox <- df_cox[complete.cases(df_cox), ]

# Establish desired reference levels for all the factor variables.
df_cox$PAM50 <- relevel(as.factor(df_cox$PAM50), ref = "LumA")
df_cox$Cellularity <- relevel(as.factor(df_cox$Cellularity), ref = "Low")
df_cox$`ER status measured by IHC` <- relevel(as.factor(df_cox$`ER status measured by IHC`), ref = "Negative")
df_cox$`HER2 status measured by SNP6` <- relevel(as.factor(df_cox$`HER2 status measured by SNP6`), ref = "NEUTRAL")
df_cox$`Hormone Therapy` <- relevel(as.factor(df_cox$`Hormone Therapy`), ref = "NO")
df_cox$`Radio Therapy` <- relevel(as.factor(df_cox$`Radio Therapy`), ref = "NO")

# Run the Cox regression model
cox_result <- coxph(Surv(OS_months, OS_status) ~ ., data = df_cox)
summary(cox_result)

# Generate Forestplot
ggforest(cox_result, data = df_cox, fontsize = 0.5,
         main = "Hazard Ratios from Cox Proportional Hazards Model",
         refLabel = "reference",  cpositions = c(0.22, 0.35, 0.43))





