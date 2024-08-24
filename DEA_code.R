#Load Required Libraries:
library(GEOquery)
library(limma)
library(cluster)
library(factoextra)
library(ggplot2) 
library(writexl)
library(ggfortify)
library(dplyr)
############################################################################### 
#Download and Load GEO Data:
geo_data <- "GSE225292" #Retrieves the GEO dataset with accession number "GSE225292".
geodf <- getGEO(geo_data, GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]] #GSEMatrix = TRUE loads expression data, and AnnotGPL = TRUE loads annotation information.The result is stored in geodf.
################################################################################################################################################
#Explore and View the Data:
View(geodf) #Inspects and views the dataset using View.
dim(geodf) #dim displays the dimensions of the data (number of rows and columns).
colnames(geodf) #colnames shows the column names.
rownames(geodf)
###############################################################################################################################################
# Extract Feature Data and Phenotype Data:
f_data <- fData(geodf) #Extracts feature data (f_data) from the GEO dataset.
View(f_data)
p_data <- pData(geodf) # Extracts phenotype data (p_data) from the GEO dataset.
View(p_data)
################################################################################################################################################
#Extract Expression Data and Perform Log Transformation:
exp <- exprs(geodf) #Extracts expression data (exp)
View(exp)
################################################################################################################################################
exp_log2 <- log2(exp+1)
nan_indices <- which(is.nan(exp_log2))
exp_log2[nan_indices] <- 0
View(exp_log2)
#############################################################################################################
# Create a Density Plot for the Expression Data:
par(mar = c(5, 5, 2, 2)) # Adjust margins
plot(density(exp_log2), main = "Density Plot of exp_log2", xlab = "Value", ylab = "Density")
########################################################################################################################################################
#Normalization using Quantile Normalization:
exp_norm <- normalizeQuantiles(exp_log2) #Performs quantile normalization on the log-transformed data.
####################################################################################################################
plot(density(exp_norm), main = "Density Plot of Log-Normalised Expression Data", xlab = "Log-Normalised", ylab = "Density")

################################################################################################################################################
#Perform K-Means Clustering:
K_means_clustering <- exp_norm
K_means_clustering <- na.omit(exp_norm)
K_means_clustering<-t(K_means_clustering)
K_means_clustering2<-data.frame(K_means_clustering)
View(K_means_clustering2)

################################################################################################################################################
# Visualize Elbow Plot:
fviz_nbclust(K_means_clustering2, kmeans, method = "wss") #Visualizes the elbow plot to help determine an appropriate number of clusters.
km <- kmeans(K_means_clustering2, centers = 2, nstart = 25) #performs k-means clustering with 2 centers and 25 different initial configurations.
#View and Assign Clusters:
View(km) #Views the k-means clustering results

################################################################################################################################################
#add the cluster assignment to the data
cluster_data <- cbind(K_means_clustering2, cluster = km$cluster) 
table(km$cluster)
# Combine the expression data with the cluster assignments
cluster_data <- data.frame(K_means_clustering2, cluster = km$cluster)
# View the final data frame
View(cluster_data)
################################################################################################################################################
gene_data_pca <- cluster_data[, -c(1, ncol(cluster_data))]
# Perform PCA
pca_result <- prcomp(gene_data_pca)
# Print summary of PCA
summary(pca_result)
# Access PC scores
pc_scores <- pca_result$x
# Add cluster information to PC scores
pc_scores_with_cluster <- cbind(pc_scores, Cluster = cluster_data$cluster)
# Create a data frame for plotting
plot_data <- as.data.frame(pc_scores_with_cluster)
# Define colors for each cluster
cluster_colors <- c("deeppink", "darkslateblue")
# Plot PCA results using ggplot2 with different colors for each cluster
ggplot(plot_data, aes(x = PC1, y = PC2, color = as.factor(Cluster))) +
  geom_point(size = 2) +
  labs(title = "Principal component analysis (PCA) of Gene Expression Data",
       x = "PC1",
       y = "PC2") +
  scale_color_manual(values = cluster_colors, 
                     labels = c("Cluster 1", "Cluster 2")) +  # Set custom colors
  labs(color = "Clusters")
theme_minimal()
################################################################################################################################################
# Specify design matrix for differential expression analysis
design_matrix <- model.matrix(~ 0 + factor(cluster_data$cluster)) #design matrix
colnames(design_matrix)=c("CL1","CL2")
cluster_data_dea<-cluster_data[,(1:ncol(cluster_data)-1)]
dea_data<-data.frame(t(cluster_data_dea)) #convert final_data to dataframe
write.csv(dea_data,"DEA_data.csv")
fit <- lmFit(dea_data, design_matrix) #Fits a linear model (lmFit) to the gene expression data (exp) using the previously defined design matrix
con<-makeContrasts(CL1-CL2, levels = design_matrix)
fit_contrasts <- contrasts.fit(fit, con)
fit_ebayes <- eBayes(fit_contrasts) #Empirical Bayes moderation is applied using eBayes
#Extracts the differential expression results for each cluster using topTable.
#coef specifies the contrast of interest
DEA_results <- topTable(fit_ebayes, number = Inf)

#######################################################################################################################################################
# Volcano Plot function
ggplot(data = DEA_results, aes(x = DEA_results[, "logFC"], y = -log10(DEA_results[, "P.Value"])))+
  geom_vline(xintercept = c(-0.6, 0.6), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + 
  geom_point() + 
  theme_set(theme_classic(base_size = 12) +
              theme(
                axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1), color = 'black'),
                plot.title = element_text(hjust = 0.5)
              ))

# Adding 'diffexpressed' column to DEA_results
DEA_results$diffexpressed <- "NO"
DEA_results$diffexpressed[DEA_results$logFC > 0.6 & DEA_results$P.Value < 0.05] <- "UP"
DEA_results$diffexpressed[DEA_results$logFC < -0.6 & DEA_results$P.Value < 0.05] <- "DOWN"

# Updated ggplot with the correct dataframe (DEA_results) and labels
ggplot(data = DEA_results, aes(x = DEA_results[, "logFC"], y = -log10(DEA_results[, "P.Value"]), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("purple", "grey", "orange"),
                     labels = c("Downregulated", "Non-significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) +
  labs(color = 'Expression ', x = expression("Log Fold Change"), y = expression("-log"[10]*"P-value")) +
  scale_x_continuous(breaks = seq(-10, 10, 2))
##########################################################################################################################################################################
# Filter DEA_results based on logFC and adjusted p-value
filtered_results <- DEA_results[abs(DEA_results$logFC) > 1.5 & DEA_results$adj.P.Val < 0.05, ]
# Print the filtered results
View(filtered_results)
#######################################################################################################################################################
write.csv(geodf,"GEO_data.csv")
write.csv(exp_norm,"Normalized Expression Data.csv")
write.csv(K_means_clustering2,"K_means Clustering Data.csv")
write.csv(cluster_data,"Cluster Data.csv")
write.csv(plot_data,"PCA Plot data.csv")
write.csv(DEA_results,"DEA_Results.csv")
write.csv(filtered_results, "filtered_results3.csv")
####################################################################################################################################################################