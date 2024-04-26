
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

get_spe_customizedQC_nucleus <- function(dir_rawData, cutoff_QC = 20){
  
  # create tmp folder
  dir_tmp = paste0("tmp",runif(1, 1, 99),"_nucleusLevel")
  dir.create(dir_tmp)

  ####################################
  samples_filderName = list.files(paste0(dir_rawData))
  ####################################
  for(sample_tmp in samples_filderName){
    print(sample_tmp)
    ## read transcripts.csv 
    # gunzip -k transcripts.csv.gz
    df_transcripts <- fread(paste0(dir_rawData, sample_tmp,
                                   # "/output-XETG00089__00011586__Region_3__20231019__161214/",
                                   "/transcripts.csv.gz"), header = T)
    ###### filter the transcripts
    df_transcripts_QC20 = df_transcripts[df_transcripts$qv >= cutoff_QC,]
    
    ###### create a csv file with only filtered transcripts
    dir.create(paste0(dir_tmp,"/",sample_tmp,"_nucleusLevel_qv"))
    
    # print("writing csv file")
    
    fwrite(df_transcripts_QC20, file= paste0(dir_tmp,"/",sample_tmp,"_nucleusLevel_qv/", "transcripts.csv"))
    
    ###### keep with-nuclei data
    dir.create(paste0(dir_tmp,"/",sample_tmp,"_nucleusLevel_qv"))
    
    df_transcripts_QC20_nucleus = df_transcripts_QC20[df_transcripts_QC20$overlaps_nucleus ==1  ,]
    fwrite(df_transcripts_QC20_nucleus, file= paste0(dir_tmp,"/",sample_tmp,"_nucleusLevel_qv/", "transcripts.csv"))
    
  }
  
  
  ##################################################################################################
  ##### copy and paste the files to a organized folder set, grouped by sample
  ##################################################################################################
  # unchanged data
  for(sample_tmp in samples_filderName){
    print(sample_tmp)
    
    files_nucleus = c(paste0(dir_rawData,sample_tmp, "/nucleus_boundaries.csv.gz"),
                     paste0(dir_rawData,sample_tmp, "/nucleus_boundaries.parquet") )
   
    ###### nucleus-level data
    file.copy(files_nucleus, paste0(dir_tmp,"/",sample_tmp,"_nucleusLevel_qv/"))
    
  }
  
  ##################################################################################################
  ### nucleus level processing
  ##################################################################################################
  ######## load raw data
  xens <- MoleculeExperiment::readXenium(paste0(dir_tmp,"/"),
                                         keepCols = "essential",
                                         addBoundaries = "nucleus")  
  
  
  ######## 
  for (i in c(1:length(samples_filderName))){
    names(xens@molecules$detected)[i] <- samples_filderName[i]
    names(xens@boundaries$nucleus)[i] <- samples_filderName[i]
  }
  
  ######## this step takes ~30min for running 13 samples
  xens <- countMolecules(xens,moleculesAssay = "detected",boundariesAssay = "nucleus")
  
  ######## remove temp folder
  unlink(dir_tmp, recursive=TRUE) 
  ######## return spe
  return(xens)
  
}
