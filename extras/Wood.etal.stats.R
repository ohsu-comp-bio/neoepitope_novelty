### Load data ###

setwd("/Users/wooma/Documents/TCGA_data_analysis/")

# Allele overlap for novel binding peptides
immunogenic_overlap <- read.table("newer.immunogenic.overlap.data.tsv", header=F, sep="\t")
colnames(immunogenic_overlap)  <- c("Allele_set", "One", "Two", "Three", "Four", "Five", "Six")

# HLA allele frequencies
allele_freqs <- read.table("allele_freqs.tsv", header=FALSE, sep="\t")
colnames(allele_freqs) <- c("Allele", "Frequency")
allele_freqs <- allele_freqs[order(allele_freqs$Allele),]

# TCGA data that did not undergo blast
BRCA <- read.table("BRCA.noblast.all.tsv", header=TRUE, sep="\t", quote = "\"")
BRCA <- BRCA[!duplicated(BRCA[,c(2,3,27)]),]

CESC <- read.table("CESC.noblast.all.tsv", header=TRUE, sep="\t", quote = "\"")
CESC <- CESC[!duplicated(CESC[,c(2,3,27)]),]

PRAD <- read.table("PRAD.noblast.all.tsv", header=TRUE, sep="\t", quote = "\"")
PRAD <- PRAD[!duplicated(PRAD[,c(2,3,27)]),]

SARC <- read.table("SARC.noblast.all.tsv", header=TRUE, sep="\t", quote = "\"")
SARC <- SARC[!duplicated(SARC[,c(2,3,27)]),]

# TCGA data that did undergo blast
ACC <- read.table("ACC.all.tsv", header=TRUE, sep="\t", quote = "\"", fill = T)
ACC <- ACC[!duplicated(ACC[,c(2,3,27)]),]
ACC <- ACC[-which(is.na(ACC$Adjusted_PS)),]

CHOL <- read.table("CHOL.all.tsv", header=TRUE, sep="\t", quote = "\"")
CHOL <- CHOL[!duplicated(CHOL[,c(2,3,27)]),]

DLBC <- read.table("DLBC.all.tsv", header=TRUE, sep="\t", quote = "\"", fill = T)
DLBC <- DLBC[!duplicated(DLBC[,c(2,3,27)]),]
DLBC <- DLBC[-which(is.na(DLBC$Adjusted_PS)),]

KICH <- read.table("KICH.all.tsv", header=TRUE, sep="\t", quote = "\"")
KICH <- KICH[!duplicated(KICH[,c(2,3,27)]),]

LAML <- read.table("LAML.all.tsv", header=TRUE, sep="\t", quote = "\"", fill=T)
LAML <- LAML[!duplicated(LAML[,c(2,3,27)]),]
LAML <- LAML[-which(is.na(LAML$Adjusted_PS)),]

MESO <- read.table("MESO.all.tsv", header=TRUE, sep="\t", quote = "\"")
MESO <- MESO[!duplicated(MESO[,c(2,3,27)]),]

PCPG <- read.table("PCPG.all.tsv", header=TRUE, sep="\t", quote = "\"")
PCPG <- PCPG[!duplicated(PCPG[,c(2,3,27)]),]

TGCT <- read.table("TGCT.all.tsv", header=TRUE, sep="\t", quote = "\"")
TGCT <- TGCT[!duplicated(TGCT[,c(2,3,27)]),]

THYM <- read.table("THYM.all.tsv", header=TRUE, sep="\t", quote = "\"")
THYM <- THYM[!duplicated(THYM[,c(2,3,27)]),]

UCS <- read.table("UCS.all.tsv", header=TRUE, sep="\t", quote = "\"", fill = T)
UCS <- UCS[!duplicated(UCS[,c(2,3,27)]),]
UCS <- UCS[-which(is.na(UCS$Adjusted_PS)),]

UVM <- read.table("UVM.all.tsv", header=TRUE, sep="\t", quote = "\"")
UVM <- UVM[!duplicated(UVM[,c(2,3,27)]),]

TCGA_data <- rbind(ACC, CHOL, DLBC, KICH, LAML, MESO, PCPG, TGCT, THYM, UCS, UVM)


### Adjust novel/non-novel binding status for new criteria ###

adj_stat <- rep("NA", length(TCGA_data$Tumor_affinity))
for (i in 1:length(TCGA_data$Tumor_affinity)){
  tumor_aff <- TCGA_data$Tumor_affinity[i]
  norm_aff <- TCGA_data$Normal_affinity[i]
  fold_change <- norm_aff/tumor_aff
  if (tumor_aff < 500 & norm_aff > 500 & fold_change >= 5){
    stat <- "immunogenic"
    adj_stat[i] <- stat
  } else {
    stat <- "nonimmunogenic"
    adj_stat[i] <- stat
  }
}
TCGA_data$Adj_stat <- adj_stat
TCGA_data$Colour <- "gray"
TCGA_data$Colour[TCGA_data$Adj_stat == "immunogenic"] <- "red"

adj_statBRCA <- rep("NA", length(BRCA$Tumor_affinity))
for (i in 1:length(BRCA$Tumor_affinity)){
  tumor_aff <- BRCA$Tumor_affinity[i]
  norm_aff <- BRCA$Normal_affinity[i]
  fold_change <- norm_aff/tumor_aff
  if (tumor_aff < 500 & norm_aff > 500 & fold_change >= 5){
    stat <- "immunogenic"
    adj_statBRCA[i] <- stat
  } else {
    stat <- "nonimmunogenic"
    adj_statBRCA[i] <- stat
  }
}
BRCA$Adj_stat <- adj_statBRCA
BRCA$Colour <- "gray"
BRCA$Colour[BRCA$Adj_stat == "immunogenic"] <- "red"

adj_statCESC <- rep("NA", length(CESC$Tumor_affinity))
for (i in 1:length(CESC$Tumor_affinity)){
  tumor_aff <- CESC$Tumor_affinity[i]
  norm_aff <- CESC$Normal_affinity[i]
  fold_change <- norm_aff/tumor_aff
  if (tumor_aff < 500 & norm_aff > 500 & fold_change >= 5){
    stat <- "immunogenic"
    adj_statCESC[i] <- stat
  } else {
    stat <- "nonimmunogenic"
    adj_statCESC[i] <- stat
  }
}
CESC$Adj_stat <- adj_statCESC
CESC$Colour <- "gray"
CESC$Colour[CESC$Adj_stat == "immunogenic"] <- "red"

adj_statPRAD <- rep("NA", length(PRAD$Tumor_affinity))
for (i in 1:length(PRAD$Tumor_affinity)){
  tumor_aff <- PRAD$Tumor_affinity[i]
  norm_aff <- PRAD$Normal_affinity[i]
  fold_change <- norm_aff/tumor_aff
  if (tumor_aff < 500 & norm_aff > 500 & fold_change >= 5){
    stat <- "immunogenic"
    adj_statPRAD[i] <- stat
  } else {
    stat <- "nonimmunogenic"
    adj_statPRAD[i] <- stat
  }
}
PRAD$Adj_stat <- adj_statPRAD
PRAD$Colour <- "gray"
PRAD$Colour[PRAD$Adj_stat == "immunogenic"] <- "red"

adj_statSARC <- rep("NA", length(SARC$Tumor_affinity))
for (i in 1:length(SARC$Tumor_affinity)){
  tumor_aff <- SARC$Tumor_affinity[i]
  norm_aff <- SARC$Normal_affinity[i]
  fold_change <- norm_aff/tumor_aff
  if (tumor_aff < 500 & norm_aff > 500 & fold_change >= 5){
    stat <- "immunogenic"
    adj_statSARC[i] <- stat
  } else {
    stat <- "nonimmunogenic"
    adj_statSARC[i] <- stat
  }
}
SARC$Adj_stat <- adj_statSARC
SARC$Colour <- "gray"
SARC$Colour[SARC$Adj_stat == "immunogenic"] <- "red"

# Data frame for all diseases
TCGA_alldis <- rbind(TCGA_data, BRCA, CESC, PRAD, SARC)

# Data frame for only those that went through blast
TCGA_9merscors <- read.table("updated.bac.vir.9mers.tsv", header=TRUE, sep="\t")
TCGA_data <- cbind(TCGA_data, TCGA_9merscors)


### Neoepitope burden by disease ###

# All neoepitopes
ACC_pat <- as.data.frame(table(ACC$Tumor_sample))
ACC_pat$Var2 <- "ACC"
BRCA_pat <- as.data.frame(table(BRCA$Tumor_sample))
BRCA_pat$Var2 <- "BRCA"
CESC_pat <- as.data.frame(table(CESC$Tumor_sample))
CESC_pat$Var2 <- "CESC"
CHOL_pat <- as.data.frame(table(CHOL$Tumor_sample))
CHOL_pat$Var2 <- "CHOL"
DLBC_pat <- as.data.frame(table(DLBC$Tumor_sample))
DLBC_pat$Var2 <- "DLBC"
KICH_pat <- as.data.frame(table(KICH$Tumor_sample))
KICH_pat$Var2 <- "KICH"
LAML_pat <- as.data.frame(table(LAML$Tumor_sample))
LAML_pat$Var2 <- "LAML"
MESO_pat <- as.data.frame(table(MESO$Tumor_sample))
MESO_pat$Var2 <- "MESO"
PCPG_pat <- as.data.frame(table(PCPG$Tumor_sample))
PCPG_pat$Var2 <- "PCPG"
PRAD_pat <- as.data.frame(table(PRAD$Tumor_sample))
PRAD_pat$Var2 <- "PRAD"
SARC_pat <- as.data.frame(table(SARC$Tumor_sample))
SARC_pat$Var2 <- "SARC"
TGCT_pat <- as.data.frame(table(TGCT$Tumor_sample))
TGCT_pat$Var2 <- "TGCT"
THYM_pat <- as.data.frame(table(THYM$Tumor_sample))
THYM_pat$Var2 <- "THYM"
UCS_pat <- as.data.frame(table(UCS$Tumor_sample))
UCS_pat$Var2 <- "UCS"
UVM_pat <- as.data.frame(table(UVM$Tumor_sample))
UVM_pat$Var2 <- "UVM"
all_pat <- rbind(ACC_pat, BRCA_pat, CESC_pat, CHOL_pat, DLBC_pat, KICH_pat, LAML_pat, MESO_pat, PCPG_pat, PRAD_pat, SARC_pat, TGCT_pat, THYM_pat, UCS_pat, UVM_pat)

# Novel-binding neoepitopes
ACC_immpat <- as.data.frame(table(ACC$Tumor_sample[ACC$Tumor_affinity < 500 & ACC$Normal_affinity > 500 & ACC$Normal_affinity >= 5*ACC$Tumor_affinity]))
ACC_immpat$Var2 <- "ACC"
BRCA_immpat <- as.data.frame(table(BRCA$Tumor_sample[BRCA$Tumor_affinity < 500 & BRCA$Normal_affinity > 500 & BRCA$Normal_affinity >= 5*BRCA$Tumor_affinity]))
BRCA_immpat$Var2 <- "BRCA"
CESC_immpat <- as.data.frame(table(CESC$Tumor_sample[CESC$Tumor_affinity < 500 & CESC$Normal_affinity > 500 & CESC$Normal_affinity >= 5*CESC$Tumor_affinity]))
CESC_immpat$Var2 <- "CESC"
CHOL_immpat <- as.data.frame(table(CHOL$Tumor_sample[CHOL$Tumor_affinity < 500 & CHOL$Normal_affinity > 500 & CHOL$Normal_affinity >= 5*CHOL$Tumor_affinity]))
CHOL_immpat$Var2 <- "CHOL"
DLBC_immpat <- as.data.frame(table(DLBC$Tumor_sample[DLBC$Tumor_affinity < 500 & DLBC$Normal_affinity > 500 & DLBC$Normal_affinity >= 5*DLBC$Tumor_affinity]))
DLBC_immpat$Var2 <- "DLBC"
KICH_immpat <- as.data.frame(table(KICH$Tumor_sample[KICH$Tumor_affinity < 500 & KICH$Normal_affinity > 500 & KICH$Normal_affinity >= 5*KICH$Tumor_affinity]))
KICH_immpat$Var2 <- "KICH"
LAML_immpat <- as.data.frame(table(LAML$Tumor_sample[LAML$Tumor_affinity < 500 & LAML$Normal_affinity > 500 & LAML$Normal_affinity >= 5*LAML$Tumor_affinity]))
LAML_immpat$Var2 <- "LAML"
MESO_immpat <- as.data.frame(table(MESO$Tumor_sample[MESO$Tumor_affinity < 500 & MESO$Normal_affinity > 500 & MESO$Normal_affinity >= 5*MESO$Tumor_affinity]))
MESO_immpat$Var2 <- "MESO"
PCPG_immpat <- as.data.frame(table(PCPG$Tumor_sample[PCPG$Tumor_affinity < 500 & PCPG$Normal_affinity > 500 & PCPG$Normal_affinity >= 5*PCPG$Tumor_affinity]))
PCPG_immpat$Var2 <- "PCPG"
PRAD_immpat <- as.data.frame(table(PRAD$Tumor_sample[PRAD$Tumor_affinity < 500 & PRAD$Normal_affinity > 500 & PRAD$Normal_affinity >= 5*PRAD$Tumor_affinity]))
PRAD_immpat$Var2 <- "PRAD"
SARC_immpat <- as.data.frame(table(SARC$Tumor_sample[SARC$Tumor_affinity < 500 & SARC$Normal_affinity > 500 & SARC$Normal_affinity >= 5*SARC$Tumor_affinity]))
SARC_immpat$Var2 <- "SARC"
TGCT_immpat <- as.data.frame(table(TGCT$Tumor_sample[TGCT$Tumor_affinity < 500 & TGCT$Normal_affinity > 500 & TGCT$Normal_affinity >= 5*TGCT$Tumor_affinity]))
TGCT_immpat$Var2 <- "TGCT"
THYM_immpat <- as.data.frame(table(THYM$Tumor_sample[THYM$Tumor_affinity < 500 & THYM$Normal_affinity > 500 & THYM$Normal_affinity >= 5*THYM$Tumor_affinity]))
THYM_immpat$Var2 <- "THYM"
UCS_immpat <- as.data.frame(table(UCS$Tumor_sample[UCS$Tumor_affinity < 500 & UCS$Normal_affinity > 500 & UCS$Normal_affinity >= 5*UCS$Tumor_affinity]))
UCS_immpat$Var2 <- "UCS"
UVM_immpat <- as.data.frame(table(UVM$Tumor_sample[UVM$Tumor_affinity < 500 & UVM$Normal_affinity > 500 & UVM$Normal_affinity >= 5*UVM$Tumor_affinity]))
UVM_immpat$Var2 <- "UVM"
all_immpat <- rbind(ACC_immpat, BRCA_immpat, CESC_immpat, CHOL_immpat, DLBC_immpat, KICH_immpat, LAML_immpat, MESO_immpat, PCPG_immpat, PRAD_immpat, SARC_immpat, TGCT_immpat, THYM_immpat, UCS_immpat, UVM_immpat)

# Number of patients per diseases
disease_counts <- as.data.frame(rev(sort(table(all_pat$Var2))))

# Median burden per disease
med_burdens <- c(median(ACC_pat$Freq), median(BRCA_pat$Freq), median(CESC_pat$Freq), median(CHOL_pat$Freq), median(DLBC_pat$Freq), median(KICH_pat$Freq), median(LAML_pat$Freq), median(MESO_pat$Freq), median(PCPG_pat$Freq), median(PRAD_pat$Freq), median(SARC_pat$Freq), median(TGCT_pat$Freq), median(THYM_pat$Freq), median(UCS_pat$Freq), median(UVM_pat$Freq))
dis_names <- c("ACC", "BRCA", "CESC", "CHOL", "DLBC", "KICH", "LAML", "MESO", "PCPG", "PRAD", "SARC", "TGCT", "THYM", "UCS", "UVM")
median_burdens <- as.data.frame(cbind(dis_names, as.numeric(med_burdens)))

# Proportion of novel binding neoepitopes per disease
proportions <- rep(0, length(all_pat$Var1))
for (i in 1:length(all_pat$Var1)){
  this_prop <- all_immpat$Freq[i]/all_pat$Freq[i]
  proportions[i] <- this_prop
}
all_pat$Proportions <- proportions
mean(all_pat$Proportions)
max(aggregate(all_pat$Proportions, list(all_pat$Var2), mean)$x)
min(aggregate(all_pat$Proportions, list(all_pat$Var2), mean)$x)


### Neoepitope burden by HLA allele ###

TCGA_filtered <- subset(TCGA_alldis, TCGA_data$Adj_stat == "immunogenic")
TCGA_alleles_filtered <- as.data.frame(rev(sort(table(TCGA_filtered$Allele))))
colors <- NULL
A_counts <- NULL
B_counts <- NULL
C_counts <- NULL
for (i in 1:length(TCGA_alleles_filtered$Var1)){
  if (grepl("A", TCGA_alleles_filtered$Var1[i]) == TRUE){
    colors <- append(colors, "green")
    A_counts <- append(A_counts, TCGA_alleles_filtered$Freq[i])
  } else if (grepl("B", TCGA_alleles_filtered$Var1[i]) == TRUE) {
    colors <- append(colors, "blue")
    B_counts <- append(B_counts, TCGA_alleles_filtered$Freq[i])
  } else if (grepl("C", TCGA_alleles_filtered$Var1[i]) == TRUE){
    colors <- append(colors, "yellow")
    C_counts <- append(C_counts, TCGA_alleles_filtered$Freq[i])
  }
}
TCGA_alleles_filtered$Allele_colors <- colors
TCGA_alleles_filtered <- TCGA_alleles_filtered[order(TCGA_alleles_filtered$Var1),]
frequencies <- NULL
for (i in 1:length(TCGA_alleles_filtered$Var1)){
  allele_f <- allele_freqs$Frequency[as.character(allele_freqs$Allele) == TCGA_alleles_filtered$Var1[i]]
  frequencies <- append(frequencies, allele_f)
}
TCGA_alleles_filtered$Freq2 <- frequencies

# Range
TCGA_alleles_disease <- as.data.frame(table(TCGA_filtered$Allele, TCGA_filtered$Disease))
quant1 <- c(quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "ACC"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "BRCA"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "CESC"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "CHOL"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "DLBC"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "KICH"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "LAML"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "MESO"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "PCPG"],  0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "PRAD"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "SARC"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "TGCT"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "THYM"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "UCS"], 0.25), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "UVM"], 0.25))
quant2 <- c(quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "ACC"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "BRCA"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "CESC"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "CHOL"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "DLBC"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "KICH"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "LAML"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "MESO"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "PCPG"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "PRAD"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "SARC"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "TGCT"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "THYM"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "UCS"], 0.75), quantile(TCGA_alleles_disease$Freq[TCGA_alleles_disease$Var2 == "UVM"], 0.75))
inter_quant <- quant2/quant1
median(inter_quant)

# Correlation of allele frequency and novel binding neoepitope burden
cor.test(TCGA_alleles_filtered$Freq, TCGA_alleles_filtered$Freq2, alternative = "two.sided", method="pearson")


### Proportion of neoepitopes belonging to HLA-A vs HLA-B vs HLA-C alleles ###

ACC_A <- 0
ACC_B <- 0
ACC_C <- 0
ACC_total <- 0
for (i in 1:length(ACC$Allele)){
  ACC_total <- ACC_total + 1
  allele <- ACC$Allele[i]
  if (grepl("A", allele) == T) {
    ACC_A <- ACC_A + 1
  } else if (grepl("B", allele) == T){
    ACC_B <- ACC_B + 1
  } else if (grepl("C", allele) == T) {
    ACC_C <- ACC_C + 1
  }
}

BRCA_A <- 0
BRCA_B <- 0
BRCA_C <- 0
BRCA_total <- 0
for (i in 1:length(BRCA$Allele)){
  BRCA_total <- BRCA_total + 1
  allele <- BRCA$Allele[i]
  if (grepl("A", allele) == T) {
    BRCA_A <- BRCA_A + 1
  } else if (grepl("B", allele) == T){
    BRCA_B <- BRCA_B + 1
  } else if (grepl("C", allele) == T) {
    BRCA_C <- BRCA_C + 1
  }
}

CESC_A <- 0
CESC_B <- 0
CESC_C <- 0
CESC_total <- 0
for (i in 1:length(CESC$Allele)){
  CESC_total <- CESC_total + 1
  allele <- CESC$Allele[i]
  if (grepl("A", allele) == T) {
    CESC_A <- CESC_A + 1
  } else if (grepl("B", allele) == T){
    CESC_B <- CESC_B + 1
  } else if (grepl("C", allele) == T) {
    CESC_C <- CESC_C + 1
  }
}

CHOL_A <- 0
CHOL_B <- 0
CHOL_C <- 0
CHOL_total <- 0
for (i in 1:length(CHOL$Allele)){
  CHOL_total <- CHOL_total + 1
  allele <- CHOL$Allele[i]
  if (grepl("A", allele) == T) {
    CHOL_A <- CHOL_A + 1
  } else if (grepl("B", allele) == T){
    CHOL_B <- CHOL_B + 1
  } else if (grepl("C", allele) == T) {
    CHOL_C <- CHOL_C + 1
  }
}

DLBC_A <- 0
DLBC_B <- 0
DLBC_C <- 0
DLBC_total <- 0
for (i in 1:length(DLBC$Allele)){
  DLBC_total <- DLBC_total + 1
  allele <- DLBC$Allele[i]
  if (grepl("A", allele) == T) {
    DLBC_A <- DLBC_A + 1
  } else if (grepl("B", allele) == T){
    DLBC_B <- DLBC_B + 1
  } else if (grepl("C", allele) == T) {
    DLBC_C <- DLBC_C + 1
  }
}

KICH_A <- 0
KICH_B <- 0
KICH_C <- 0
KICH_total <- 0
for (i in 1:length(KICH$Allele)){
  KICH_total <- KICH_total + 1
  allele <- KICH$Allele[i]
  if (grepl("A", allele) == T) {
    KICH_A <- KICH_A + 1
  } else if (grepl("B", allele) == T){
    KICH_B <- KICH_B + 1
  } else if (grepl("C", allele) == T) {
    KICH_C <- KICH_C + 1
  }
}

LAML_A <- 0
LAML_B <- 0
LAML_C <- 0
LAML_total <- 0
for (i in 1:length(LAML$Allele)){
  LAML_total <- LAML_total + 1
  allele <- LAML$Allele[i]
  if (grepl("A", allele) == T) {
    LAML_A <- LAML_A + 1
  } else if (grepl("B", allele) == T){
    LAML_B <- LAML_B + 1
  } else if (grepl("C", allele) == T) {
    LAML_C <- LAML_C + 1
  }
}

MESO_A <- 0
MESO_B <- 0
MESO_C <- 0
MESO_total <- 0
for (i in 1:length(MESO$Allele)){
  MESO_total <- MESO_total + 1
  allele <- MESO$Allele[i]
  if (grepl("A", allele) == T) {
    MESO_A <- MESO_A + 1
  } else if (grepl("B", allele) == T){
    MESO_B <- MESO_B + 1
  } else if (grepl("C", allele) == T) {
    MESO_C <- MESO_C + 1
  }
}

PCPG_A <- 0
PCPG_B <- 0
PCPG_C <- 0
PCPG_total <- 0
for (i in 1:length(PCPG$Allele)){
  PCPG_total <- PCPG_total + 1
  allele <- PCPG$Allele[i]
  if (grepl("A", allele) == T) {
    PCPG_A <- PCPG_A + 1
  } else if (grepl("B", allele) == T){
    PCPG_B <- PCPG_B + 1
  } else if (grepl("C", allele) == T) {
    PCPG_C <- PCPG_C + 1
  }
}

PRAD_A <- 0
PRAD_B <- 0
PRAD_C <- 0
PRAD_total <- 0
for (i in 1:length(PRAD$Allele)){
  PRAD_total <- PRAD_total + 1
  allele <- PRAD$Allele[i]
  if (grepl("A", allele) == T) {
    PRAD_A <- PRAD_A + 1
  } else if (grepl("B", allele) == T){
    PRAD_B <- PRAD_B + 1
  } else if (grepl("C", allele) == T) {
    PRAD_C <- PRAD_C + 1
  }
}

SARC_A <- 0
SARC_B <- 0
SARC_C <- 0
SARC_total <- 0
for (i in 1:length(SARC$Allele)){
  SARC_total <- SARC_total + 1
  allele <- SARC$Allele[i]
  if (grepl("A", allele) == T) {
    SARC_A <- SARC_A + 1
  } else if (grepl("B", allele) == T){
    SARC_B <- SARC_B + 1
  } else if (grepl("C", allele) == T) {
    SARC_C <- SARC_C + 1
  }
}

TGCT_A <- 0
TGCT_B <- 0
TGCT_C <- 0
TGCT_total <- 0
for (i in 1:length(TGCT$Allele)){
  TGCT_total <- TGCT_total + 1
  allele <- TGCT$Allele[i]
  if (grepl("A", allele) == T) {
    TGCT_A <- TGCT_A + 1
  } else if (grepl("B", allele) == T){
    TGCT_B <- TGCT_B + 1
  } else if (grepl("C", allele) == T) {
    TGCT_C <- TGCT_C + 1
  }
}

THYM_A <- 0
THYM_B <- 0
THYM_C <- 0
THYM_total <- 0
for (i in 1:length(THYM$Allele)){
  THYM_total <- THYM_total + 1
  allele <- THYM$Allele[i]
  if (grepl("A", allele) == T) {
    THYM_A <- THYM_A + 1
  } else if (grepl("B", allele) == T){
    THYM_B <- THYM_B + 1
  } else if (grepl("C", allele) == T) {
    THYM_C <- THYM_C + 1
  }
}

UCS_A <- 0
UCS_B <- 0
UCS_C <- 0
UCS_total <- 0
for (i in 1:length(UCS$Allele)){
  UCS_total <- UCS_total + 1
  allele <- UCS$Allele[i]
  if (grepl("A", allele) == T) {
    UCS_A <- UCS_A + 1
  } else if (grepl("B", allele) == T){
    UCSC_B <- UCS_B + 1
  } else if (grepl("C", allele) == T) {
    UCS_C <- UCS_C + 1
  }
}

UVM_A <- 0
UVM_B <- 0
UVM_C <- 0
UVM_total <- 0
for (i in 1:length(UVM$Allele)){
  UVM_total <- UVM_total + 1
  allele <- UVM$Allele[i]
  if (grepl("A", allele) == T) {
    UVM_A <- UVM_A + 1
  } else if (grepl("B", allele) == T){
    UVM_B <- UVM_B + 1
  } else if (grepl("C", allele) == T) {
    UVM_C <- UVM_C + 1
  }
}

A_proportions <- c((ACC_A/ACC_total), (BRCA_A/BRCA_total), (CESC_A/CESC_total), (CHOL_A/CHOL_total), (DLBC_A/DLBC_total), (KICH_A/KICH_total), (LAML_A/LAML_total), (MESO_A/MESO_total), (PCPG_A/PCPG_total), (PRAD_A/PRAD_total), (SARC_A/SARC_total), (TGCT_A/TGCT_total), (THYM_A/THYM_total), (UCS_A/UCS_total), (UVM_A/UVM_total))
mean(A_proportions)
sd(A_proportions)

B_proportions <- c((ACC_B/ACC_total), (BRCA_B/BRCA_total), (CESC_B/CESC_total), (CHOL_B/CHOL_total), (DLBC_B/DLBC_total), (KICH_B/KICH_total), (LAML_B/LAML_total), (MESO_B/MESO_total), (PCPG_B/PCPG_total), (PRAD_B/PRAD_total), (SARC_B/SARC_total), (TGCT_B/TGCT_total), (THYM_B/THYM_total), (UCS_B/UCS_total), (UVM_B/UVM_total))
mean(B_proportions)
sd(B_proportions)

C_proportions <- c((ACC_C/ACC_total), (BRCA_C/BRCA_total), (CESC_C/CESC_total), (CHOL_C/CHOL_total), (DLBC_C/DLBC_total), (KICH_C/KICH_total), (LAML_C/LAML_total), (MESO_C/MESO_total), (PCPG_C/PCPG_total), (PRAD_C/PRAD_total), (SARC_C/SARC_total), (TGCT_C/TGCT_total), (THYM_C/THYM_total), (UCS_C/UCS_total), (UVM_C/UVM_total))
mean(C_proportions)
sd(C_proportions)


### Top alleles ###

head(TCGA_alleles_filtered, n=10)

B1503_eps <- subset(TCGA_alldis, TCGA_alldis$Allele=="B1503")
B1503_eps$Tumor_peptide <- factor(B1503_eps$Tumor_peptide)
ep_counts <- as.data.frame(table(B1503_eps$Tumor_peptide))

A0206_eps <- subset(TCGA_alldis, TCGA_alldis$Allele=="A0206")
A0206_eps$Tumor_peptide <- factor(A0206_eps$Tumor_peptide)
ep_counts2 <- as.data.frame(table(A0206_eps$Tumor_peptide))

A3001_eps <- subset(TCGA_alldis, TCGA_alldis$Allele=="A3001")
A3001_eps$Tumor_peptide <- factor(A3001_eps$Tumor_peptide)
ep_counts5 <- as.data.frame(table(A3001_eps$Tumor_peptide))

B1517_eps <- subset(TCGA_alldis, TCGA_alldis$Allele=="B1517")
B1517_eps$Tumor_peptide <- factor(B1517_eps$Tumor_peptide)
ep_counts6 <- as.data.frame(table(B1517_eps$Tumor_peptide))

A0202_eps <- subset(TCGA_alldis, TCGA_alldis$Allele=="A0202")
A0202_eps$Tumor_peptide <- factor(A0202_eps$Tumor_peptide)
ep_counts3 <- as.data.frame(table(A0202_eps$Tumor_peptide))

A0203_eps <- subset(TCGA_alldis, TCGA_alldis$Allele=="A0203")
A0203_eps$Tumor_peptide <- factor(A0203_eps$Tumor_peptide)
ep_counts4 <- as.data.frame(table(A0203_eps$Tumor_peptide))

A0205_eps <- subset(TCGA_alldis, TCGA_alldis$Allele=="A0205")
A0205_eps$Tumor_peptide <- factor(A0205_eps$Tumor_peptide)
ep_counts8 <- as.data.frame(table(A0205_eps$Tumor_peptide))

C1203_eps <- subset(TCGA_alldis, TCGA_alldis$Allele=="C1203")
C1203_eps$Tumor_peptide <- factor(C1203_eps$Tumor_peptide)
ep_counts9 <- as.data.frame(table(C1203_eps$Tumor_peptide))

C1402_eps <- subset(TCGA_alldis, TCGA_alldis$Allele=="C1402")
C1402_eps$Tumor_peptide <- factor(C1402_eps$Tumor_peptide)
ep_counts10 <- as.data.frame(table(C1402_eps$Tumor_peptide))

C1403_eps <- subset(TCGA_alldis, TCGA_alldis$Allele=="C1403")
C1403_eps$Tumor_peptide <- factor(C1403_eps$Tumor_peptide)
ep_counts7 <- as.data.frame(table(C1403_eps$Tumor_peptide))

ep_means <- c(mean(ep_counts$Freq), mean(ep_counts2$Freq), mean(ep_counts3$Freq), mean(ep_counts4$Freq), mean(ep_counts5$Freq), mean(ep_counts6$Freq), mean(ep_counts7$Freq), mean(ep_counts8$Freq), mean(ep_counts9$Freq), mean(ep_counts10$Freq))
ep_max <- c(max(ep_counts$Freq), max(ep_counts2$Freq), max(ep_counts3$Freq), max(ep_counts4$Freq), max(ep_counts5$Freq), max(ep_counts6$Freq), max(ep_counts7$Freq), max(ep_counts8$Freq), max(ep_counts9$Freq), max(ep_counts10$Freq))

range(ep_max)
mean(ep_means)


### Allelic verlap of novel binding epitopes ###

immunogenic_overlap$Total <- immunogenic_overlap$One + immunogenic_overlap$Two + immunogenic_overlap$Three + immunogenic_overlap$Four + immunogenic_overlap$Five + immunogenic_overlap$Six

mean(immunogenic_overlap$Total)
mean(immunogenic_overlap$One/immunogenic_overlap$Total)
mean(immunogenic_overlap$Two/immunogenic_overlap$Total)
mean(immunogenic_overlap$Three/immunogenic_overlap$Total)
mean(immunogenic_overlap$Four/immunogenic_overlap$Total)
mean(immunogenic_overlap$Five/immunogenic_overlap$Total)
mean(immunogenic_overlap$Six/immunogenic_overlap$Total)
max(immunogenic_overlap$Six/immunogenic_overlap$Total)

immunogenic_overlap$Prp_one <- 100*immunogenic_overlap$One/immunogenic_overlap$Total
immunogenic_overlap$Prp_two <- 100*immunogenic_overlap$Two/immunogenic_overlap$Total
immunogenic_overlap$Prp_three <- 100*immunogenic_overlap$Three/immunogenic_overlap$Total
immunogenic_overlap$Prp_four <- 100*immunogenic_overlap$Four/immunogenic_overlap$Total
immunogenic_overlap$Prp_five <- 100*immunogenic_overlap$Five/immunogenic_overlap$Total
immunogenic_overlap$Prp_six <- 100*immunogenic_overlap$Six/immunogenic_overlap$Total


### Binding difference ###

mutant_pos <- rep("Nonanchor", length(TCGA_alldis$Tumor_peptide))
for (i in 1:length(TCGA_alldis$Tumor_peptide)){
  tum_pep <- strsplit(as.character(TCGA_alldis$Tumor_peptide[i]),"")[[1]]
  norm_pep <- strsplit(as.character(TCGA_alldis$Normal_peptide[i]),"")[[1]]
  if (tum_pep[2] != norm_pep[2]) {
    mutant_pos[i] <- "Anchor"
  } else if (tum_pep[9] != norm_pep[9]) {
    mutant_pos[i] <- "Anchor"
  }
}
TCGA_alldis$Mutant_pos <- mutant_pos

wilcox.test(TCGA_alldis$Binding_difference~TCGA_alldis$Mutant_pos)
median(TCGA_alldis$Binding_difference[TCGA_alldis$Mutant_pos == "Anchor"])
median(TCGA_alldis$Binding_difference[TCGA_alldis$Mutant_pos == "Nonanchor"])


### Protein similarity to paired normal epitope ###

mean(TCGA_alldis$Adjusted_PS[TCGA_alldis$Mutant_pos == "Nonanchor"])
sd(TCGA_alldis$Adjusted_PS[TCGA_alldis$Mutant_pos == "Nonanchor"])
max(TCGA_alldis$Adjusted_PS[TCGA_alldis$Mutant_pos == "Nonanchor"])
min(TCGA_alldis$Adjusted_PS[TCGA_alldis$Mutant_pos == "Nonanchor"])

wilcox.test(TCGA_alldis$Adjusted_PS[TCGA_alldis$Mutant_pos == "Nonanchor"]~TCGA_alldis$Adj_stat[TCGA_alldis$Mutant_pos == "Nonanchor"])

TCGA_nonanch <- subset(TCGA_alldis, TCGA_alldis$Mutant_pos == "Nonanchor")
TCGA_nonanch$Mutant_pos <- factor(TCGA_nonanch$Mutant_pos)

wilcox.test(TCGA_nonanch$Adjusted_PS[TCGA_nonanch$Tumor_affinity < 100]~TCGA_nonanch$Adj_stat[TCGA_nonanch$Tumor_affinity < 100])
wilcox.test(TCGA_nonanch$Adjusted_PS[TCGA_nonanch$Tumor_affinity < 200 & TCGA_nonanch$Tumor_affinity > 99]~TCGA_nonanch$Adj_stat[TCGA_nonanch$Tumor_affinity < 200 & TCGA_nonanch$Tumor_affinity > 99])
wilcox.test(TCGA_nonanch$Adjusted_PS[TCGA_nonanch$Tumor_affinity < 300 & TCGA_nonanch$Tumor_affinity > 199]~TCGA_nonanch$Adj_stat[TCGA_nonanch$Tumor_affinity < 300 & TCGA_nonanch$Tumor_affinity > 199])
wilcox.test(TCGA_nonanch$Adjusted_PS[TCGA_nonanch$Tumor_affinity < 400 & TCGA_nonanch$Tumor_affinity > 299]~TCGA_nonanch$Adj_stat[TCGA_nonanch$Tumor_affinity < 400 & TCGA_nonanch$Tumor_affinity > 299])
wilcox.test(TCGA_nonanch$Adjusted_PS[TCGA_nonanch$Tumor_affinity < 500 & TCGA_nonanch$Tumor_affinity > 399]~TCGA_nonanch$Adj_stat[TCGA_nonanch$Tumor_affinity < 500 & TCGA_nonanch$Tumor_affinity > 399])

mean(TCGA_alldis$Adjusted_PS[TCGA_alldis$Mutant_pos == "Nonanchor" & TCGA_alldis$Adj_stat == "immunogenic"])
mean(TCGA_alldis$Adjusted_PS[TCGA_alldis$Mutant_pos == "Nonanchor" & TCGA_alldis$Adj_stat == "nonimmunogenic"])


### Protein similarity to all human peptides ###

TCGA_blast <- TCGA_data[complete.cases(TCGA_data$Blast_score),]

length(TCGA_blast$Blast_score)/length(TCGA_data$Blast_score)
mean(TCGA_blast$Adjusted_blast_PS)
sd(TCGA_blast$Adjusted_blast_PS)

# Find matching vs nonmatching hits
TCGA_blast_complete <- TCGA_blast[TCGA_blast$Orig_gene != "missing",]
gene_match <- rep("NA", length(TCGA_blast_complete$Orig_gene))
for (i in 1:length(TCGA_blast_complete$Orig_gene)){
  if (grepl(TCGA_blast_complete$Orig_gene[i], TCGA_blast_complete$Match_gene[i]) == T){
    gene_match[i] <- "Matching"
  } else {
    gene_match[i] <- "Nonmatching"
  }
}
TCGA_blast_complete$Blast_match <- gene_match
num_gene_matches <- table(TCGA_blast_complete$Blast_match, TCGA_blast_complete$Disease)

# Find exact vs inexact matches
TCGA_blast_complete$Match_exactness <- rep("inexact", length(TCGA_blast_complete$Tumor_peptide))
for (i in 1:length(TCGA_blast_complete$Match_exactness)){
  if (as.character(TCGA_blast_complete$Tumor_peptide[i]) == as.character(TCGA_blast_complete$Match_seq[i])){
    TCGA_blast_complete$Match_exactness[i] <- "exact"
  }
}

# Per patient by disease
ACC_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "ACC")
ACC_blast$Tumor_sample <- factor(ACC_blast$Tumor_sample)
ACC_match <- as.data.frame(table(ACC_blast$Tumor_sample[ACC_blast$Match_Stat == "Matching"]))
ACC_match$Var2 <- "ACC"
ACC_nonmatch <- as.data.frame(table(ACC_blast$Tumor_sample[ACC_blast$Match_Stat == "Nonmatching"]))
ACC_nonmatch$Var2 <- "ACC"
CHOL_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "CHOL")
CHOL_blast$Tumor_sample <- factor(CHOL_blast$Tumor_sample)
CHOL_match <- as.data.frame(table(CHOL_blast$Tumor_sample[CHOL_blast$Match_Stat == "Matching"]))
CHOL_match$Var2 <- "CHOL"
CHOL_nonmatch <- as.data.frame(table(CHOL_blast$Tumor_sample[CHOL_blast$Match_Stat == "Nonmatching"]))
CHOL_nonmatch$Var2 <- "CHOL"
DLBC_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "DLBC")
DLBC_blast$Tumor_sample <- factor(DLBC_blast$Tumor_sample)
DLBC_match <- as.data.frame(table(DLBC_blast$Tumor_sample[DLBC_blast$Match_Stat == "Matching"]))
DLBC_match$Var2 <- "DLBC"
DLBC_nonmatch <- as.data.frame(table(DLBC_blast$Tumor_sample[DLBC_blast$Match_Stat == "Nonmatching"]))
DLBC_nonmatch$Var2 <- "DLBC"
KICH_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "KICH")
KICH_blast$Tumor_sample <- factor(KICH_blast$Tumor_sample)
KICH_match <- as.data.frame(table(KICH_blast$Tumor_sample[KICH_blast$Match_Stat == "Matching"]))
KICH_match$Var2 <- "KICH"
KICH_nonmatch <- as.data.frame(table(KICH_blast$Tumor_sample[KICH_blast$Match_Stat == "Nonmatching"]))
KICH_nonmatch$Var2 <- "KICH"
LAML_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "LAML")
LAML_blast$Tumor_sample <- factor(LAML_blast$Tumor_sample)
LAML_match <- as.data.frame(table(LAML_blast$Tumor_sample[LAML_blast$Match_Stat == "Matching"]))
LAML_match$Var2 <- "LAML"
LAML_nonmatch <- as.data.frame(table(LAML_blast$Tumor_sample[LAML_blast$Match_Stat == "Nonmatching"]))
LAML_nonmatch$Var2 <- "LAML"
MESO_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "MESO")
MESO_blast$Tumor_sample <- factor(MESO_blast$Tumor_sample)
MESO_match <- as.data.frame(table(MESO_blast$Tumor_sample[MESO_blast$Match_Stat == "Matching"]))
MESO_match$Var2 <- "MESO"
MESO_nonmatch <- as.data.frame(table(MESO_blast$Tumor_sample[MESO_blast$Match_Stat == "Nonmatching"]))
MESO_nonmatch$Var2 <- "MESO"
PCPG_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "PCPG")
PCPG_blast$Tumor_sample <- factor(PCPG_blast$Tumor_sample)
PCPG_match <- as.data.frame(table(PCPG_blast$Tumor_sample[PCPG_blast$Match_Stat == "Matching"]))
PCPG_match$Var2 <- "PCPG"
PCPG_nonmatch <- as.data.frame(table(PCPG_blast$Tumor_sample[PCPG_blast$Match_Stat == "Nonmatching"]))
PCPG_nonmatch$Var2 <- "PCPG"
TGCT_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "TGCT")
TGCT_blast$Tumor_sample <- factor(TGCT_blast$Tumor_sample)
TGCT_match <- as.data.frame(table(TGCT_blast$Tumor_sample[TGCT_blast$Match_Stat == "Matching"]))
TGCT_match$Var2 <- "TGCT"
TGCT_nonmatch <- as.data.frame(table(TGCT_blast$Tumor_sample[TGCT_blast$Match_Stat == "Nonmatching"]))
TGCT_nonmatch$Var2 <- "TGCT"
THYM_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "THYM")
THYM_blast$Tumor_sample <- factor(THYM_blast$Tumor_sample)
THYM_match <- as.data.frame(table(THYM_blast$Tumor_sample[THYM_blast$Match_Stat == "Matching"]))
THYM_match$Var2 <- "THYM"
THYM_nonmatch <- as.data.frame(table(CHOL_blast$Tumor_sample[THYM_blast$Match_Stat == "Nonmatching"]))
THYM_nonmatch$Var2 <- "THYM"
UCS_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "UCS")
UCS_blast$Tumor_sample <- factor(UCS_blast$Tumor_sample)
UCS_match <- as.data.frame(table(UCS_blast$Tumor_sample[UCS_blast$Match_Stat == "Matching"]))
UCS_match$Var2 <- "UCS"
UCS_nonmatch <- as.data.frame(table(UCS_blast$Tumor_sample[UCS_blast$Match_Stat == "Nonmatching"]))
UCS_nonmatch$Var2 <- "UCS"
UVM_blast <- subset(TCGA_blast_complete, TCGA_blast_complete$Disease == "UVM")
UVM_blast$Tumor_sample <- factor(UVM_blast$Tumor_sample)
UVM_match <- as.data.frame(table(UVM_blast$Tumor_sample[UVM_blast$Match_Stat == "Matching"]))
UVM_match$Var2 <- "UVM"
UVM_nonmatch <- as.data.frame(table(UVM_blast$Tumor_sample[UVM_blast$Match_Stat == "Nonmatching"]))
UVM_nonmatch$Var2 <- "UVM"
all_match <- rbind(ACC_match, CHOL_match, DLBC_match, KICH_match, LAML_match, MESO_match, PCPG_match, TGCT_match, THYM_match, UCS_match, UVM_match)
all_nonmatch <- rbind(ACC_nonmatch, CHOL_nonmatch, DLBC_nonmatch, KICH_nonmatch, LAML_nonmatch, MESO_nonmatch, PCPG_nonmatch, TGCT_nonmatch, THYM_nonmatch, UCS_nonmatch, UVM_nonmatch)

mean(all_match$Freq/(all_match$Freq+all_nonmatch$Freq))
mean(all_nonmatch$Freq/(all_match$Freq+all_nonmatch$Freq))

TCGA_nonmatch <- subset(TCGA_blast_complete, TCGA_blast_complete$Match_Stat == "Nonmatching")
exact_matches <- 0
for (i in 1:length(TCGA_nonmatch$Tumor_peptide)){
  if (as.character(TCGA_nonmatch$Tumor_peptide[i]) == as.character(TCGA_nonmatch$Match_seq[i])){
    exact_matches <- exact_matches + 1
  }
}
exact_matches/length(TCGA_nonmatch$Tumor_peptide)

more_sim <- 0
for (i in 1:length(TCGA_nonmatch$Adjusted_PS)){
  if (TCGA_nonmatch$Adjusted_PS[i] < TCGA_nonmatch$Adjusted_blast_PS[i]){
    more_sim <- more_sim + 1
  }
}
more_sim/length(TCGA_nonmatch$Tumor_peptide)


### Protein similarity to bacterial/viral peptides ###

# Bacterial
TCGA_bac <- TCGA_data[complete.cases(TCGA_data$Bac_score),]
length(TCGA_bac$Bac_score)/length(TCGA_data$Bac_score)
mean(TCGA_bac$Adjusted_Bac_PS)
sd(TCGA_bac$Adjusted_Bac_PS)

# Viral
TCGA_vir <- TCGA_data[complete.cases(TCGA_data$Vir_score),]
length(TCGA_vir$Vir_score)/length(TCGA_data$Vir_score)
mean(TCGA_vir$Adjusted_Vir_PS)
sd(TCGA_vir$Adjusted_Vir_PS)

t.test(TCGA_vir$Adjusted_Vir_PS, TCGA_bac$Adjusted_Bac_PS)

# More similar bacterial peptides
bac_better <- subset(TCGA_bac, TCGA_bac$Adjusted_Bac_PS > TCGA_bac$Adjusted_blast_PS & TCGA_bac$Adjusted_Bac_PS > TCGA_bac$Adjusted_PS)
length(bac_better$Tumor_peptide)/length(TCGA_bac$Tumor_peptide)
bac_better9mer <- subset(TCGA_bac, TCGA_bac$Bac_9mer_PS > TCGA_bac$Blast_9mer & TCGA_bac$Bac_9mer_PS > TCGA_bac$Norm_9mer)
length(bac_better9mer$Tumor_peptide)/length(TCGA_bac$Tumor_peptide)

# More similar viral peptides
vir_better <- subset(TCGA_vir, TCGA_vir$Adjusted_Vir_PS > TCGA_vir$Adjusted_blast_PS & TCGA_vir$Adjusted_Vir_PS > TCGA_vir$Adjusted_PS)
length(vir_better$Tumor_peptide)/length(TCGA_vir$Tumor_peptide)
vir_better9mer <- subset(TCGA_vir, TCGA_vir$Vir_9mer_PS. > TCGA_vir$Blast_9mer & TCGA_vir$Vir_9mer_PS. > TCGA_vir$Norm_9mer)
length(vir_better9mer$Tumor_peptide)/length(TCGA_vir$Tumor_peptide)

# Bacterial mismatch
bac_mismatches <- rep(0, length(TCGA_bac))
for(i in 1:length(TCGA_bac$Tumor_peptide)){
  a <- as.character(TCGA_bac$Tumor_peptide[i])
  b <- as.character(TCGA_bac$Bac_seq[i])
  mismatches <- length(which(unlist(strsplit(a, "")) != unlist(strsplit(b, ""))))
  bac_mismatches[i] <- mismatches
}
norm_mismatches <- rep(0, length(TCGA_bac))
for(i in 1:length(TCGA_bac$Tumor_peptide)){
  a <- as.character(TCGA_bac$Tumor_peptide[i])
  b <- as.character(TCGA_bac$Normal_peptide[i])
  mismatches <- length(which(unlist(strsplit(a, "")) != unlist(strsplit(b, ""))))
  norm_mismatches[i] <- mismatches
}
blast_mismatches <- rep(0, length(TCGA_bac))
for(i in 1:length(TCGA_bac$Tumor_peptide)){
  a <- as.character(TCGA_bac$Tumor_peptide[i])
  b <- as.character(TCGA_bac$Match_seq[i])
  mismatches <- length(which(unlist(strsplit(a, "")) != unlist(strsplit(b, ""))))
  blast_mismatches[i] <- mismatches
}
mean(bac_mismatches)
mean(norm_mismatches)
mean(blast_mismatches)
t.test(bac_mismatches, norm_mismatches)
t.test(bac_mismatches, blast_mismatches)

better_bac_species <- rev(sort(table(bac_better$Bac_match)))