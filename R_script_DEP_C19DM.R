# https://gtk-teaching.github.io/Microarrays-R/06-DifferentialGeneExpression/index.html
# scritp version was built on limma R package
# R 4.0.5

library(limma)
library(reshape2)
library(ggplot2)
library(dplyr)

# intensity normalized dataset of all peptide spectrum matches was used as input (omnidelta_new_IRS_normalized_kk.csv), but only the experiments at the time 24h were analyzed, 
to maximize the signals of infection. 
# https://www.ebi.ac.uk/pride/archive/projects/PXD037265 

setwd("./scripting_DEp_BYCOVID/")
data_row <- read.csv("omnidelta_new_IRS_normalized_kk_24h.csv", sep=";", row.names = 1)


# experiental design

conditions <- read.table("exp_desi.txt")
colnames(conditions)<-conditions[1,]
conditions<-conditions[-c(1),]
str(conditions)

# only unique proteins
data_row$Protein_Id<-row.names(data_row)
data_unique <- make_unique(data_row, "Gene_Symbol.1", "Protein_Id", delim = ";")
str(data_unique)
data_unique$name %>% duplicated() %>% any()

# logaritmic tranformation (log2)
log_data <- log2(data_unique[,2:14])


# Create design matrix for variant infection condition (B.1, Delta, Omicron, Mock)
conditions_list<-factor(conditions$condition)
design <- model.matrix(~ 0 + conditions_list)
colnames(design) <- levels(conditions_list)

# linear model fitting
fit <- lmFit(log_data, design)

# Contrast function to compare SARS-CoV-2-infected cell coltures and uninfected cell coltures

contrast.matrix_B1_mock <- makeContrasts(B.1 - Mock, levels = design)
# contrast.matrix_Delta_mock <- makeContrasts(Delta - Mock, levels = design)
# contrast.matrix_Omicron_mock <- makeContrasts(Omicron - Mock, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix_B1_mock)
# fit2 <- contrasts.fit(fit, contrast.matrix_Delta_mock)
# fit2 <- contrasts.fit(fit, contrast.matrix_Omicron_mock)


# Differential analysis of protein expression
fit2 <- eBayes(fit2)
topTable(fit2, number = 10)
res <- topTable(fit2, number = Inf)

# save all results
write.csv(res, "allproteins_B.1_Mock.csv")
#write.csv(res, "allproteins_Delta_Mock.csv")
#write.csv(res, "allproteins_Omicron_Mock.csv")


# Treasholds imposed FDR < 0.05 e logFC > 0.5, as reported in Mezek eta al., 2023.

sign_res <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 0.5,]

write.csv(sign_res, "DEP_Omicron_Mock_5%_0.5.csv")
#write.csv(sign_res, "DEP_Delta_Mock_5%_0.5.csv")
#write.csv(sign_res, "DEP_B.1_Mock_5%_0.5.csv")

# volcano plot

res$Significant <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1, "Significant", "Not Significant")
res$Protein <- rownames(res)

# print vulcano plot
ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot") +
  theme(legend.position = "right") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_text(data = subset(results, Significant == "Significant"),
      aes(label = Protein), vjust = -0.5, size = 3, check_overlap = TRUE)


# Transfomation for overlay of C19DM

overlayer = sign_res   

# Making a proportion to trasform the significant LogFC values, reporting all values in range [-1, 1].
# In overlayer_top you can find the higest value, in absolute value, of LogFC list to calculate new proportional values.

overlayer_top <- ifelse(max(abs(C19DM_overlayer$logFC)) > min(abs(C19DM_overlayer$logFC)), max(C19DM_overlayer$logFC), min(C19DM_overlayer$logFC))

overlayer$value <- if (C19DM_overlayer_top > 0) {
		(1*C19DM_overlayer$logFC)/C19DM_overlayer_top
	} else if (C19DM_overlayer_top < 0){
		-((1*C19DM_overlayer$logFC)/C19DM_overlayer_top)
	}

overlayer$name<-rownames(C19DM_overlayer)
overlayer$AveExpr<-NULL
overlayer$t<-NULL
overlayer$P.Value<-NULL
overlayer$adj.P.Val<-NULL
overlayer$B<-NULL
rownames(overlayer)<-NULL
write.table(overlayer, "C19DM_overlayer.txt", sep = "\t", quote=FALSE, row.names=FALSE)

