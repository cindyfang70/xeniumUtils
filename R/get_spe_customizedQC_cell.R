
####################################
library(data.table)
library(Biostrings)
library(ggplot2)
library(gridExtra)
library(spatialLIBD)
require(colorout)
library(default)
library(MoleculeExperiment)
library(magrittr)
library(scater)
library(scran)
library(BiocParallel)
library(job)
library(SpatialExperiment)

# need to set the default path with enough space (~1G per sample)
###################################

get_spe_customizedQC_cell <- function(dir_rawData, cutoff_QC = 20){
  
  # create tmp folder
  dir_tmp = paste0("tmp",runif(1, 1, 99),"_cellLevel")
  dir.create(dir_tmp)
  
  ####################################
  samples_filderName = list.files(paste0(dir_rawData))[1]
  
  ####################################
  for(sample_tmp in samples_filderName){
    print(sample_tmp)
    ## read transcripts.csv 
    # gunzip -k transcripts.csv.gz
    df_transcripts <- fread(paste0(dir_rawData,sample_tmp,
                                   # "/output-XETG00089__00011586__Region_3__20231019__161214/",
                                   "/transcripts.csv.gz"), header = T)
    ###### filter the transcripts
    df_transcripts_QC20 = df_transcripts[df_transcripts$qv >= cutoff_QC,]
    
    ###### create a csv file with only filtered transcripts
    dir.create(paste0(dir_tmp,"/",sample_tmp,"_cellLevel_qv"))
    
    print("writing csv file")
    
    fwrite(df_transcripts_QC20, file= paste0(dir_tmp,"/",sample_tmp,"_cellLevel_qv/", "transcripts.csv"))
    
    ###### keep with-cell data
    dir.create(paste0(dir_tmp,"/",sample_tmp,"_cellLevel_qv"))
    
    df_transcripts_QC20_cell = df_transcripts_QC20[df_transcripts_QC20$overlaps_cell ==1  ,]
    fwrite(df_transcripts_QC20_cell, file= paste0(dir_tmp,"/",sample_tmp,"_cellLevel_qv/", "transcripts.csv"))
    
  }
  
  
  ##################################################################################################
  ##### copy and paste the files to a organized folder set, grouped by sample
  ##################################################################################################
  # unchanged data
  for(sample_tmp in samples_filderName){
    print(sample_tmp)
    
    files_cell = c(paste0(dir_rawData,sample_tmp, "/cell_boundaries.csv.gz"),
                      paste0(dir_rawData,sample_tmp, "/cell_boundaries.parquet") )
    
    ###### cell-level data
    file.copy(files_cell, paste0(dir_tmp,"/",sample_tmp,"_cellLevel_qv/"))
    
  }
  
  ##################################################################################################
  ### cell level processing
  ##################################################################################################
  ######## load raw data
  xens <- MoleculeExperiment::readXenium(paste0(dir_tmp,"/"),
                                         keepCols = "essential",
                                         addBoundaries = "cell")  
  
  
  ######## 
  for (i in c(1:length(samples_filderName))){
    names(xens@molecules$detected)[i] <- samples_filderName[i]
    names(xens@boundaries$cell)[i] <- samples_filderName[i]
  }
  
  ######## this step takes ~30min for running 13 samples
  xens <- countMolecules(xens, moleculesAssay = "detected",boundariesAssay = "cell")
  
  ######## remove temp folder
  unlink(dir_tmp, recursive=TRUE) 
  ######## return spe
  return(xens)
  
}
