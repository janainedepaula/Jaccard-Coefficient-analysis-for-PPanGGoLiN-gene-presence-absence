######################################################################################################
# APPLYING JACCARD COEFFICIENT TO BUILD A DENDOGRAM BASED ON PPANGGOLIN GENE PRESENCE/ABSENCE OUTPUT #
######################################################################################################

#df = data frame. Your data frame must be gene_presence_absence.csv table derived from PPanGGoLiN. 
df <- read.table("gene_presence_absence_sem.csv", header = TRUE, row.names = 1, check.names = FALSE)

#This function computes Jaccard coefficient between two lineages and repeats the calculation for all pairs
jaccard <- function(lineage1, lineage2, df) {
  shared <- sum(df[[lineage1]] == 1 & df[[lineage2]] == 1) #Compute CDS shared by both lineages (intersection)
  union <- sum(df[[lineage1]] == 1 | df[[lineage2]] == 1)  #Count CDS present in at least one of the two lineages (union)
  
  if (union == 0) {
    return(0)  #If union value is zero, the two lineages share no CDS.
  }
  
  return(shared / union)  #Apply Jaccard similarity formula
}

lineages <- colnames(df) #Retrieve lineage names from column names

#Prepare an empty matrix to compute the Jaccard similarity between all pairs of lineages
matrix_jaccard <- matrix(0, nrow = length(lineages), ncol = length(lineages), 
                         dimnames = list(lineages, lineages))

#Compute the Jaccard similarity matrix for each pair of lineages
for (i in 1:(length(lineages)-1)) {
  for (j in (i+1):length(lineages)) {
    jaccard <- jaccard(lineages[i], lineages[j], df)
    matrix_jaccard[i, j] <- jaccard
    matrix_jaccard[j, i] <- jaccard  #Symmetric matrix
  }
}

#Convert similarity matrix to distance matrix (1 - Jaccard) and display it 
matrix_distance <- 1 - matrix_jaccard
print(matrix_distance)

#install.packages("ape")
library(ape) #Load 'ape' library to generate the dendrogram based on the distance matrix 
dist_matrix <- as.dist(matrix_distance) #Convert the distance matrix to a distance object
dendrograma <- hclust(dist_matrix) #Perform hierarchical clustering to build the dendrogram 
plot(dendrograma, main = "Dendrogram based on CDS sharing matrix using Jaccard Coefficient", xlab = "Lineages", sub = "") #Plot the dendrogram

dist_matrix_df <- as.data.frame(matrix_distance) #convert the distance matrix to dataframe
print(dist_matrix_df)