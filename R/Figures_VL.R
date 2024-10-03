


# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Create the data frame
data <- data.frame(
  AnimalID = c("A1=37374", "A2=37636", "A3=35644", "A4=36497", "A5=36135", 
               "A6=36622", "B1=37163", "B2=38180", "B3=37166", "B4=37167", 
               "B5=36590", "B6=37448", "C1=35112", "C2=36131", "C3=37246", 
               "C4=37826", "C5=36686", "C6=36767", "D1=38089", "D2=38090", 
               "D3=35547", "D4=36562", "D5=36683", "D6=37012"),
  Testis_vRNA = c(0.01, 0.01, 0.01, 5.32, 0.01, 0.41, 0.01, 1.36, 0.07, 1.23, 
                  0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.86, 0.01, 0.01, 
                  0.16, 0.01, 0.21, 0.01),
  Testis_vDNA = c(0.032186667, 0.032186667, 0.032186667, 5.323, 0.032186667, 
                  0.032186667, 0.032186667, 0.032186667, 0.032186667, 24.42, 
                  3.037, 0.032186667, 0.032186667, 0.032186667, 0.032186667, 
                  0.032186667, 0.032186667, 0.032186667, 0.032186667, 
                  0.032186667, 0.032186667, 0.032186667, 1.51, 0.032186667),
  Testis_RNA_DNA_ratio = c(0.310687652, 0.310687652, 0.310687652, 0.999436408, 
                           0.310687652, 12.73819374, 0.310687652, 42.25352069, 
                           2.174813565, 0.05036855, 0.003292723, 0.310687652, 
                           0.310687652, 0.310687652, 0.310687652, 0.310687652, 
                           0.310687652, 26.71913808, 0.310687652, 0.310687652, 
                           4.971002434, 0.310687652, 0.139072848, 0.310687652),
  PBMC_vRNA = c(0.01, 0.61, 41.68, 31.14, 15.28, 3.63, 1.21, 1.1, 2.11, 37.98, 
                52.38, 2.11, 11.09, 0.7, 5.04, 1.3, 0.01, 41.55, 7.6, 1.53, 
                0.01, 0.01, 680.3, 0),
  PBMC_vDNA = c(3.475, 50.48, 88.32, 274.5, 209.8, 84.4, 29.84, 26.96, 82.99, 
                72.21, 48.14, 41.57, 239.7, 27.62, 49.02, 21.71, 15.03, 47.8, 
                70.63, 18.75, 127.2, 13.67, 120.4, 0),
  PBMC_RNA_DNA_ratio = c(0.002877698, 0.012083994, 0.47192029, 0.113442623, 
                         0.072831268, 0.043009479, 0.040549598, 0.040801187, 
                         0.02542475, 0.525965933, 1.088076444, 0.050757758, 
                         0.046266166, 0.025343954, 0.102815177, 0.05988024, 
                         0.000665336, 0.869246862, 0.107603002, 0.0816, 
                         7.86164E-05, 0.000731529, 5.650332226, 0)
)

# Transform the data into long format
data_long <- data %>%
  select(AnimalID, Testis_RNA_DNA_ratio, PBMC_RNA_DNA_ratio) %>%
  gather(key = "Tissue", value = "RNA_DNA_ratio", -AnimalID)

# Rename levels for clarity
data_long$Tissue <- factor(data_long$Tissue, levels = c("Testis_RNA_DNA_ratio", "PBMC_RNA_DNA_ratio"),
                           labels = c("Testis", "PBMC"))

# Plot the data
ggplot(data_long, aes(x = AnimalID, y = (RNA_DNA_ratio), fill = Tissue)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_fill_manual(values = c("Testis" = "#1f78b4", "PBMC" = "#33a02c")) +
  theme_bw(base_size = 14)  + 
  labs(title = "RNA/DNA in Testis vs. PBMC",
       x = "Available Chronic SIV+ animals",
       y = "RNA/DNA Ratio",
       fill = "Tissue") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) 




# Transform the data into long format foasinh()# Transform the data into long format for RNA/DNA ratio
data_long <- data %>%
  select(AnimalID, Testis_RNA_DNA_ratio, PBMC_RNA_DNA_ratio) %>%
  gather(key = "Tissue", value = "RNA_DNA_ratio", -AnimalID)

# Rename levels for clarity
data_long$Tissue <- factor(data_long$Tissue, levels = c("Testis_RNA_DNA_ratio", "PBMC_RNA_DNA_ratio"),
                           labels = c("Testis", "PBMC"))

# Plot the data
ggplot() +
  geom_bar(data = data_long, aes(x = AnimalID, y = RNA_DNA_ratio, fill = Tissue), 
           stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  scale_fill_manual(values = c("Testis" = "#1f78b4", "PBMC" = "#33a02c")) +
  # geom_line(data = data, aes(x = AnimalID, y = (PBMC_RNA_DNA_ratio), group = 1), color = "red", size = 1) +
  # geom_line(data = data, aes(x = AnimalID, y = (Testis_vDNA), group = 1), color = "black", size = 1) +
  # geom_line(data = data, aes(x = AnimalID, y = (Testis_vRNA), group = 1), color = "dodgerblue", size = 1) +
  # geom_point(data = data, aes(x = AnimalID, y = sqrt(PBMC_vDNA)), color = "red", size = 2) +
  scale_y_continuous(
    name = "RNA/DNA Ratio",
    # sec.axis = sec_axis(~ .^2, name = "PBMC vDNA") 
    sec.axis = sec_axis(~ ., name = "PBMC vDNA")
  ) +
  theme_minimal(base_size = 14) +
  labs(title = "Comparison of RNA/DNA Ratios in Testis and PBMC",
       x = "Animal ID",
       fill = "Tissue") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.y.right = element_text(margin = margin(l = 10), color = "red"),
    axis.text.y.right = element_text(color = "red"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )
# 


