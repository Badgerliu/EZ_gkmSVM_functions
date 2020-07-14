## before start the program, resize the training set and remove the repeat-masked region (>70%) ###


setwd("/Users/huanliu/OneDrive/HPC/oral_biology_paper")

##########Start running gkmsvm#########
library(gkmSVM) 
set.seed(12)
library(BSgenome.Hsapiens.UCSC.hg19.masked)

source("easy_gkmSVM_vector_hg19.R")
easy_gkmSVM_vector_hg19("HIOEC_intersected_trainingset.bed",20,1) #example usage

