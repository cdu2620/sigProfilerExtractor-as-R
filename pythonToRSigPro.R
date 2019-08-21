setwd('"/Users/cdu2620"')

library(tidyverse)
library(hrbrthemes)
library(scales)
library(ggpubr)
library(ggsci)
library(janitor)
library(factoextra)
library(ggplot2)
library(dplyr)
library(rlist)
library(entropy)
cos_sim <- function(a, b) {
  if (sum(a) == 0 | sum(b) == 0) {
    value = 0
  }
  dot_product = a %*% b
  norm_a <- norm(as.matrix(a), "2")
  norm_b <- norm(as.matrix(b), "2")
  value = dot_product / (norm_a * norm_b)
  return(value)
}

calculate_similarities <- function(genomes, est_genomes, sample_names) {
  cosine_sim_list = c()
  kl_divergence_list = c()
  l1_norm_list = c()
  l2_norm_list = c()
  total_mutations = c()
  relative_l1_list = c()
  relative_l2_list = c()
  for (i in 1:ncol(genomes)) {
    p_i <- as.numeric(genomes[, i])
    q_i = (est_genomes[, i])
    cosine_sim_list = append(cosine_sim_list, round(cos_sim(p_i, q_i), digits=3))
    kl_divergence_list = append(kl_divergence_list, round(KL.empirical(p_i, q_i), digits=4))
    l1_norm_list = append(l1_norm_list, round(norm(as.matrix(p_i-q_i), "1"), digits=3))
    relative_l1_list = append(relative_l1_list, round((dplyr::last(l1_norm_list)/norm(as.matrix(p_i), "1"))*100, digits=3))
    l2_norm_list = append(l2_norm_list, round(norm(as.matrix(p_i-q_i), "2"), digits=3))
    relative_l2_list = append(relative_l2_list, round((dplyr::last(l2_norm_list)/norm(as.matrix(p_i), "2"))*100, digits=3))
    total_mutations = append(total_mutations, sum(p_i))
  }
  kl_divergence_list[!is.finite(kl_divergence_list)] <- 1000
  similarities_dataframe = data.frame("Sample Names"=sample_names,
                                      "Total Mutations"=total_mutations,
                                      "Cosine similarity"=cosine_sim_list,
                                      "L1 Norm"=l1_norm_list,
                                      "L1_Norm_%"=relative_l1_list,
                                      "L2 Norm"=l2_norm_list,
                                      "L2_Norm_%"=relative_l2_list,
                                      "KL Divergence"= kl_divergence_list)
  #write.csv(similarities_dataframe, file="MyData.csv")
  write.table(similarities_dataframe, file="MyData.txt", sep="\t", row.names = FALSE)
}

calculateSims <- function(text_file, text_file2, text_file3) {
  data2 <- read.delim(text_file)
  data2 <- data.frame(data2)
  data2 <- data2[,!is.na(colSums(data2 != 0)) & colSums(data2 != 0) > 0]
  genomes <- data2[, 2:length(data2)]
  allcolnames <- colnames(data2[,2:ncol(data2)])
  data3 <- read.delim(text_file2)
  data3 <- data3[,2:ncol(data3)]
  data4 <- read.delim(text_file3)
  data4 <- data4[,2:ncol(data4)]
  est_genomes <- as.data.frame(as.matrix(data3) %*% as.matrix(t(data4)))
  calculate_similarities(genomes, est_genomes, allcolnames)
}

calculateSims("orignal_genomes.txt", "Wsignature_sigs.txt", "Wsignature_activaties.txt")
