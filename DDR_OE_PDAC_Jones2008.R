#load necessary packages
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(stringr)

#import data
rawdata <- read_excel("C://Users/jayra/Desktop/practice_data/hpne_heatmap/pdac_heatmap/table_s9/table_S9_1164368tableS9.xls") 


# understand the structure of the data
str(rawdata)
head(rawdata)
tail(rawdata)
summary(rawdata)

# identify and separate variables and observations
#check for column names
print(names(rawdata))
head(rawdata, 10)
#extract desires column names (alternatively assign desired column names)
column_name <- rawdata[7, 1:34]

#assign correct or desired column names
names(rawdata) <- column_name
#remove symbols from column names and assign desired name
names(rawdata)[33] <- "overexpressed"
names(rawdata)[34] <- "extracellular"

# tidy data

#remove unnecessary column using slice() and select() functions of dplyr for rows and column respectively
tidy_rawdata <- rawdata %>%
  slice(8:nrow(rawdata)) %>%
  select(-34)
 
# print(head(tidy_rawdata), width = Inf) to confirm removal

# Remove rows with any missing values )
tidy_rawdata_clean <- na.omit(tidy_rawdata)
#check if the NA rows were removed successfully
#print(which(!complete.cases(tidy_rawdata_clean))) 

# Convert "+" and "-" to "yes" and "no"
data_mutated <- tidy_rawdata_clean %>%
  mutate(overexpressed = if_else(overexpressed == "+", "yes", "no"))

#filter rows with genes overexpressed as "yes"
genes_overexpressed <- data_mutated |> filter(overexpressed == "yes")
#convert all gene names to UPPERCASE
genes_overexpressed$Gene <- toupper(genes_overexpressed$Gene)

str(genes_overexpressed)
#541 genes are overexpressed
#this matches with the article
#at least 10-fold overexpressed in >90% pf the 24 cancers (compared to normal pancreatic duct cells or HPDE cells)

# Remove duplicate rows based on the "Gene" column
go_df_gene <- genes_overexpressed %>%
  distinct(Gene, .keep_all = TRUE)

#reshape the data from wide to long
gene_overexpressed_long <- pivot_longer(go_df_gene,
                                        cols = c(2:32),
                                        names_to = "Sample_ID",
                                        values_to = "fold_change")
gene_overexpressed_long <- gene_overexpressed_long |> select(-2)


#create metadata 
metadata <- as.data.frame(t(rawdata |> slice(1, 2))) 
names(metadata)[1] <- "Sample_ID"
names(metadata)[2] <- "Sample_type"
str(metadata)
#remove unnecessary rows
metadata <- metadata |> slice(2:32)

#join metadata and gene_overexpressed_long data 
genes_oe_df <- left_join(gene_overexpressed_long, metadata, by = "Sample_ID")


# download HGNC data
#get the gene symbol and alternative names from HGNC
#https://www.genenames.org/download/archive/#!/#tocAnchor-1-3
# Specify the URL to download the file from
hgnc_url <- "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"

# Specify the destination file path to save the downloaded file
gene_set_2024apr16 <- "C://Users/jayra/Desktop/practice_data/hpne_heatmap/pdac_heatmap/hgnc_complete_set.txt"
download.file(hgnc_url, gene_set_2024apr16)

#open gene set detab file
hgnc_genelist <- read.delim(gene_set_2024apr16, sep = "\t")

#subset desired columns
hgnc_genelist1 <- hgnc_genelist |> select (2, 9, 11)
names(hgnc_genelist1)[1] <- "Gene"


hgnc_genelist2 <- hgnc_genelist1 %>%
  mutate(alias_pieces = str_count(alias_symbol, "\\|") + 1) |>
  mutate(prev_pieces = str_count(prev_symbol, "\\|") + 1)
print(max(hgnc_genelist2$alias_pieces))
print(max(hgnc_genelist2$prev_pieces))

# Create a vector of column names dynamically
alias_symbol_cols <- paste0("as", seq_len(22))
prev_symbol_cols <- paste0("ps", seq_len(18))


hgnc_genelist3 <- hgnc_genelist1 %>%
  separate(
  col = alias_symbol,
  into = alias_symbol_cols,
  sep = "\\|"
)

hgnc_genelist4 <- hgnc_genelist1 %>%
  separate(
    col = prev_symbol,
    into = prev_symbol_cols,
    sep = "\\|"
  )

hgnc_genelist_5 <- left_join(hgnc_genelist3, hgnc_genelist4, by = "Gene")
hgnc_genelist_5 <- hgnc_genelist_5 |> select(-24, -25)
hgnc_long <- pivot_longer(hgnc_genelist_5,
                          cols = c(2:41),
                          names_to = "asps_col",
                          values_to = "as_or_ps")
hgnc_long <- select(hgnc_long, 1, 3)
hgnc_long <- hgnc_long |> distinct()
hgnc_long_clean <- filter(hgnc_long, !is.na(as_or_ps))
#convert all the gene names in as_or_ps column to upper case (merge function is case sensitive)
hgnc_long_clean$Gene <- toupper(hgnc_long_clean$Gene)
hgnc_long_clean$as_or_ps <- toupper(hgnc_long_clean$as_or_ps)
--------
#join over-expressed and HGNC list#####

# Perform joins for each column: Gene, alias_symbol, and prev_symbol

#find matches between OE gene list and HGNC list

# find matches between gene column and create new dataframe of matches by using merge() function
joined_gene_oe <- merge(hgnc_long_clean, genes_oe_df, by = "Gene")
joined_gene_oe <- select(joined_gene_oe, 1)
joined_gene_oe <- joined_gene_oe |> distinct() #number of genes is 446 (remaining unmatched genes should be 95 (541-446))

# find matches between Alis/Prev name column of HGNC and Gene column of OE genes
# first identify OE genes for which matches have not been found
column_gene_oe <- genes_overexpressed$Gene
column_joined_oe <- joined_gene_oe$Gene
#these are the oe genes for which match has not been  found
unique_gene_oe <- as.data.frame(setdiff(column_gene_oe, column_joined_oe))
unique_gene_oe <- unique_gene_oe |> distinct()
names(unique_gene_oe)[1] <- "Gene"
#sum of unique gene not matched plus matched genes should be equal to total genes over-expressed.
#in this case it is same so all is well!

#find hgnc genes for which match with oe gene has not been found
column_hgnc_long_clean <- hgnc_long_clean$Gene
unique_hgnc_long_clean_oe <- as.data.frame(setdiff(column_hgnc_long_clean, column_joined_oe))
names(unique_hgnc_long_clean_oe)[1] <- "Gene"
#gene the hgnc genes for which match with ddr genes has not been found with the ALias names column for hgnc
merge_unique_hgnc_oe <- merge(unique_hgnc_long_clean_oe, hgnc_long_clean, by = "Gene")

#find matches between unique hgnc and unique oe for which matches have not been found
joined_alias_oe <- merge(merge_unique_hgnc_oe, unique_gene_oe, by.x = "as_or_ps", by.y = "Gene")

#gene names for ddr genes matching with the alias
joined_alias_oe1 <- select(joined_alias_oe, 1) #83
names(joined_alias_oe1)[1] <- "Gene"
joined_alias_oe1 <- joined_alias_oe1 |> distinct() #79 
#gene names for ddr genes matching with the alias
joined_gene_oe1 <- joined_gene_oe |> select(1)
joined_df_oe1 <- data.frame()
joined_df_oe1 <- rbind(joined_gene_oe1, joined_alias_oe1)
joined_df_oe1 <- joined_df_oe1 |> distinct() #this number should be 541 but it is 525, which means 16 values are missing
column_joined_df_oe1 <- joined_df_oe1$Gene
#missing should be zero and if not then it may be because it is a transcript or old name not in HGNC database
missing_oe_gene1 <- as.data.frame(setdiff(column_gene_oe, column_joined_df_oe1)) # unmatched genes are 16.

--------------------
#gene names for ddr genes matching with the HGNC gene symbol of the alias
joined_alias_oe2 <- select(joined_alias_oe, 2)
joined_alias_oe2 <- joined_alias_oe2 |> distinct() #this number should be 95, but it is 83, 12 unmatched genes.
#unmatched genes are either transcript or unknown 

-----------------
# Combine all the results into one dataframe
#gene names for oe genes matching with the HGNC gene symbol of the alias
joined_df_oe <- data.frame()
joined_df_oe <- rbind(joined_gene_oe, joined_alias_oe2)
joined_df_oe <- joined_df_oe |> distinct()
#joined_df_ddr <- joined_df_ddr |> filter(Gene != "BHLHE40" & Gene != "IRF4")

----------------
#Final dataframe of DDR genes matching HGNC data
oe_final <- as.data.frame(joined_df_oe$Gene)
names(oe_final)[1] <- "Gene"
-------------------------
----------


#join DDR and HGNC list#####
#Identify over expressed DDR genes
#import DDR genes list
ddr_genes <- as.data.frame(read_excel("C://Users/jayra/Desktop/practice_data/hpne_heatmap/pdac_heatmap/DDRlib_guide_RNA.xlsx"))
ddr_gene_list <- as.data.frame(read_excel("C://Users/jayra/Desktop/practice_data/hpne_heatmap/pdac_heatmap/DDRlib_guide_RNA.xlsx"))

-------------------
#use the str_replace() function of stringr to remove the suffixes from the gene_ID column:
ddr_gene_list$gene_ID <- str_replace(ddr_gene_list$gene_ID, "_\\d+$", "")
#Alternatively, you can use the gsub() function from base R to achieve the same result:
#ddr_gene_list$gene_ID <- gsub("_\\d+$", "", ddr_gene_list$gene_ID)

----------------
# Remove duplicate rows based on the specified columns
ddr_gene_list_distinct <- ddr_gene_list %>%
  distinct(gene_ID, .keep_all = TRUE)
names(ddr_gene_list_distinct)[2] <- "Gene"
#remove non-target row
ddr_gene_list_distinct <- ddr_gene_list_distinct |> slice(1:504)
#convert values in Gene column to UPPERCASE
ddr_gene_list_distinct$Gene <- toupper(ddr_gene_list_distinct$Gene)
write_xlsx(ddr_gene_list_distinct, path = "C:/Users/jayra/Desktop/practice_data/hpne_heatmap/pdac_heatmap/ddr_gene_list_distinct.xlsx", col_names = TRUE)
--------------------
#find matches between DDR gene list and HGNC list

# find matches between gene column and create new dataframe of matches by using merge() function
joined_gene_ddr <- merge(hgnc_long_clean, ddr_gene_list_distinct, by = "Gene")
joined_gene_ddr <- select(joined_gene_ddr, 1)
joined_gene_ddr <- joined_gene_ddr |> distinct() #number of genes is 473 (remaining unmatched genes should be 31 (504-473))
------------------------
# find matches between Alis/Prev name column of HGNC and Gene column of DDR
# first identify ddr genes for which matches have not been found
column_ddr_distinct <- ddr_gene_list_distinct$Gene
column_joined <- joined_gene_ddr$Gene
#these are the ddr genes for which match has not been  found
unique_ddr_distinct <- as.data.frame(setdiff(column_ddr_distinct, column_joined))
names(unique_ddr_distinct)[1] <- "Gene" #matches correct number 31 unmatched genes

#find hgnc genes for which match with ddr gene has not been found
column_hgnc_long_clean <- hgnc_long_clean$Gene
unique_hgnc_long_clean <- as.data.frame(setdiff(column_hgnc_long_clean, column_joined))
names(unique_hgnc_long_clean)[1] <- "Gene"
#gene the hgnc genes for which match with ddr genes has not been found with the ALias names column for hgnc
merge_unique_hgnc <- merge(unique_hgnc_long_clean, hgnc_long_clean, by = "Gene")
--------------------
#find matches between unique hgnc and unique ddr for which matches have not been found
joined_alias_ddr <- merge(merge_unique_hgnc, unique_ddr_distinct, by.x = "as_or_ps", by.y = "Gene")

#gene names for ddr genes matching with the alias
joined_alias_ddr1 <- select(joined_alias_ddr, 1)
names(joined_alias_ddr1)[1] <- "Gene"

#gene names for ddr genes matching with the alias
joined_df_ddr1 <- data.frame()
joined_gene_ddr1 <- joined_gene_ddr |> select(1) 
joined_alias_ddr1 <- joined_alias_ddr1 |> distinct() #this number should be 31, which it is!
joined_df_ddr1 <- rbind(joined_gene_ddr1, joined_alias_ddr1)
joined_df_ddr1 <- joined_df_ddr1 |> distinct() #this should be 504 genes, which it is!
column_joined_df_ddr1 <- joined_df_ddr1$Gene
#missing should be zero
missing_ddr_gene1 <- as.data.frame(setdiff(column_ddr_distinct, column_joined_df_ddr1)) #missing should be zero and it is!

---------------
#gene names for ddr genes matching with the HGNC gene symbol of the alias
joined_alias_ddr2 <- select(joined_alias_ddr, 2)
joined_alias_ddr2 <- joined_alias_ddr2 |> distinct() #this number should be 31, which it is not! 2 additional gene.
#two genes may have common alias or previous name
#CENPX (previous name STRA13) (is a DDR gene)
#BHLHE40 (previous name STRA13) (not a DDR gene)
#PWWP3A (previous name MUM1) (is a DDR gene)
#IRF4 (previous name MUM1) (not a DDR gene)


-----------------
# Combine all the results into one dataframe
#gene names for ddr genes matching with the HGNC gene symbol of the alias
joined_df_ddr <- data.frame()
joined_df_ddr <- rbind(joined_gene_ddr, joined_alias_ddr2)
joined_df_ddr <- joined_df_ddr |> distinct()
joined_df_ddr <- joined_df_ddr |> filter(Gene != "BHLHE40" & Gene != "IRF4")

----------------
#Final fataframe of DDR genes matching HGNC data
ddr_final <- as.data.frame(joined_df_ddr$Gene)
write_xlsx(ddr_final, path = "C:/Users/jayra/Desktop/practice_data/hpne_heatmap/pdac_heatmap/ddr_final.xlsx", col_names = TRUE)
-------------------------
----------

# Combine all the results into one dataframe######


#find overexpressed ddr genes
oe_match_ddr <- merge(joined_df_oe, joined_df_ddr, "Gene")

#find expression data of overexpressed ddr genes
oe_ddr_expression <- merge(oe_match_ddr, gene_overexpressed_long, "Gene")

#plot overexpressed ddr gene data########
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#headmap
Heatmap(oe_ddr_expression,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_labels = colnames(oe_ddr_expression),``
        row_labels = rownames(oe_ddr_expression),
        name = "fold_change")
plot <- ggplot(data= oe_ddr_expression, aes(x = Sample_ID, y = fold_change)) + geom_point()

------------------

-----------
#scRNAseq##########

#open data
data_sc <- read_excel("C://Users/jayra/Desktop/practice_data/hpne_heatmap/pdac_heatmap/peng_2019_pdac_scrnaseq/table_s3.xlsx") 

#match and merge ddr genes with scRNAseq data
data_sc_ddr <- merge(data_sc, joined_df_ddr, "Gene")