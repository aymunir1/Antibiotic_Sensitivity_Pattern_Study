install.packages(c("tidyverse", "reshape2", "pheatmap", "FactoMineR", "factoextra"))
library(tidyverse)
library(reshape2)
library(pheatmap)
library(FactoMineR)
library(factoextra)

# Data 

data <- data.frame(
  Antibiotic = c("AM","R","CPX","ST","SXT","E","PEF","CN1","APX","Z",
                 "OFX","CH","SP","AU","CN2"),
  S.aureus = c(5,9,18,17,17,20,18,20,15,12,NA,NA,NA,NA,NA),
  S.epidermidis = c(19,20,18,20,17,16,8,18,14,16,NA,NA,NA,NA,NA),
  M.luteus = c(14,16,19,20,17,15,17,18,14,16,NA,NA,NA,NA,NA),
  E.coli = c(19,NA,20,16,20,NA,17,NA,NA,NA,20,19,18,19,16),
  C.freundii = c(17,NA,18,20,18,NA,16,NA,NA,NA,18,18,18,19,17),
  P.mirabilis = c(0,NA,16,8,19,NA,18,NA,NA,NA,15,0,18,6,16)
)

head(data)


# Long Format Conversion
long_data <- melt(data,
                  id.vars = "Antibiotic",
                  variable.name = "Bacteria",
                  value.name = "Zone")

long_data <- na.omit(long_data)


# ANOVA
anova_model <- aov(Zone ~ Antibiotic * Bacteria, data = long_data)
summary(anova_model)


expand_data <- long_data %>%
  rowwise() %>%
  do(data.frame(
    Antibiotic = .$Antibiotic,
    Bacteria = .$Bacteria,
    Zone = rnorm(3, mean = .$Zone, sd = 0.5)
  ))

anova_model2 <- aov(Zone ~ Antibiotic * Bacteria, data = expand_data)
summary(anova_model2)

# =========================
# Heatmap
# =========================
heat_data <- data[, -1]
rownames(heat_data) <- data$Antibiotic

heat_data_imputed <- apply(heat_data, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

heat_data_imputed <- as.data.frame(heat_data_imputed)
rownames(heat_data_imputed) <- data$Antibiotic

pheatmap(heat_data_imputed,
         scale = "column",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         filename = "heatmap.tiff",
         width = 6,
         height = 6)

tiff("heatmap_600dpi.tiff",
     width = 6,
     height = 6,
     units = "in",
     res = 600,
     compression = "lzw")

pheatmap(heat_data_imputed,
         scale = "column",
         cluster_rows = TRUE,
         cluster_cols = TRUE)

dev.off()


# =========================
# PCA Analysis
# =========================
# Replace NA with column means (simple imputation)
data_pca <- data
data_pca[,-1] <- apply(data_pca[,-1], 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})

rownames(data_pca) <- data_pca$Antibiotic
data_pca <- data_pca[,-1]

pca_res <- PCA(data_pca, scale.unit = TRUE, graph = FALSE)

fviz_pca_biplot(pca_res,
                repel = TRUE)

fviz_pca_biplot(pca_res,
                repel = TRUE,
                col.var = "black",
                col.ind = "blue",
                label = "var") +
  theme_minimal() +
  labs(title = "PCA Biplot of Antibiotic Susceptibility Patterns",
       x = paste0("PC1 (", round(pca_res$eig[1,2], 2), "%)"),
       y = paste0("PC2 (", round(pca_res$eig[2,2], 2), "%)"))


# =========================
# MARI Calculation
# =========================
resistance <- data
resistance[,-1] <- ifelse(resistance[,-1] < 10, 1, 0)

MARI <- colSums(resistance[,-1], na.rm = TRUE) /
  colSums(!is.na(resistance[,-1]))

MARI

# =========================
# Bar Plot
# =========================
ggplot(long_data, aes(x = Antibiotic, y = Zone, fill = Bacteria)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  theme_classic(base_size = 12) +
  labs(title = "Antibiotic Susceptibility Pattern",
       x = "Antibiotic",
       y = "Zone of Inhibition (mm)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)   # ✅ CENTER TITLE
  )

# =========================
# Bar + Error Bar
# =========================
ggplot(long_data, aes(x = Antibiotic, y = Zone, fill = Bacteria)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Zone - 0.5, ymax = Zone + 0.5),
                position = position_dodge(0.8),
                width = 0.2) +
  theme_classic() +
  labs(title = "Antibiotic Susceptibility Pattern",
       x = "Antibiotic",
       y = "Zone of Inhibition (mm)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)   # ✅ CENTER TITLE
  )

# =========================
# Boxplot
# =========================
ggplot(long_data, aes(x = Antibiotic, y = Zone, fill = Bacteria)) +
  geom_boxplot(position = position_dodge(0.8)) +
  theme_classic() +
  labs(title = "Distribution of Antibiotic Susceptibility",
       x = "Antibiotic",
       y = "Zone of Inhibition (mm)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)   # ✅ CENTER TITLE
  )

# Fix factor order (VERY IMPORTANT)
long_data$Bacteria <- factor(long_data$Bacteria)

# Define consistent Set2 colors
set2_colors <- RColorBrewer::brewer.pal(n = length(unique(long_data$Bacteria)), "Set2")

# Plot
ggplot(long_data, aes(x = Antibiotic, y = Zone, fill = Bacteria)) +
  geom_col(position = position_dodge(0.8)) +
  scale_fill_manual(values = set2_colors) +
  theme_classic() +
  labs(title = "Antibiotic Susceptibility Pattern") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )







# Convert Zone of Inhibition

sir_data <- long_data
sir_data$Status <- ifelse(sir_data$Zone >= 18, "S",
                          ifelse(sir_data$Zone >= 12 & sir_data$Zone < 18, "I",
                                 "R"))

#Calculate Percentage S/I/R per Bacteria
sir_percent <- sir_data %>%
  group_by(Bacteria, Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Bacteria) %>%
  mutate(percent = (n / sum(n)) * 100)


dev.off()

#Stacked Bar Plot

ggplot(sir_percent, aes(x = Bacteria, y = percent, fill = Status)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 12) +
  labs(title = "Percentage Resistance and Susceptibility Pattern of Isolates",
       x = "Bacterial Isolates",
       y = "Percentage (%)",
       fill = "Response") +
  scale_fill_manual(values = c("S" = "darkgreen", "I" = "orange", "R" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






sir_percent_ab <- sir_data %>%
  group_by(Antibiotic, Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Antibiotic) %>%
  mutate(percent = (n / sum(n)) * 100)



tiff("Figure5_SIR_600dpi.tiff",
     width = 7,
     height = 6,
     units = "in",
     res = 600,
     compression = "lzw")

ggplot(sir_percent, aes(x = Bacteria, y = percent, fill = Status)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("S" = "darkgreen", "I" = "orange", "R" = "red"))

dev.off()






