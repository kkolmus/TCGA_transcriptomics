rm(list = ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAbiolinks")

suppressPackageStartupMessages({
  library(futile.logger)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
  library(rlist)
})


flog.threshold(DEBUG)


flog.debug("Set project directory")

proj.dir = "~/Desktop/HT projects/TCGA_transcriptomics"
ifelse(test = !dir.exists(file.path(proj.dir)), 
       yes = c(dir.create(file.path(proj.dir)), 
               setwd(file.path(proj.dir))), 
       no = "Folder exists")
getwd()


flog.debug("List available TCGA projects")

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl = TRUE)]
projects

flog.debug("Download and analyze RNASeq data for TCGA project with paired normal and tumour samples")

for(proj in projects) {
  flog.debug(paste0("Begin processing the ", proj, " project")) # print project 
  flog.debug("Query GDC data")
  query <- GDCquery(project = proj, 
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq",
                    file.type = "results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)
  
  flog.debug("Query for a type of sample: normal vs. tumour")
  sampleType <- query$results[[1]]$cases
  MatchedCoupledSampleTypes <- TCGAquery_MatchedCoupledSampleTypes(sampleType, c("NT","TP"))
  SampleTP <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "TP")
  SampleNT <- TCGAquery_SampleTypes(barcode = MatchedCoupledSampleTypes, typesample = "NT")
  
  # skip the project if there's less then 3 normal tissue samples 
  if (length(SampleNT) < 3) { 
    flog.debug("Less then 3 healthy tissue samples, skipping project")
    next 
  } 
  else {
    # extract info on project, tissue source site and participant from barcode ID 
    SampleTP_short <- substr(x = SampleTP, start = 1, stop = 12)
    SampleNT_short <- substr(x = SampleNT, start = 1, stop = 12)
    
    # get clinical data
    flog.debug("Get clinical data")
    dataClin <- assign(paste0("dataClin_", proj), 
                       GDCquery_clinic(project = proj, type = "clinical"))
    odd_indexes <- seq(1, nrow(dataClin), 2)
    dataClin <- dataClin[odd_indexes, ]
    assign(paste0("dataClin_", proj), dataClin) 
    
    # match clinical information on samples with RNASeq data
    flog.debug("Query GDC data for matched clinical samples")
    
    # primary tumour with matched normal tissue
    dataClin_TP <- dplyr::filter(dataClin, dataClin$bcr_patient_barcode %in% SampleTP_short)
    SampleTP_matched <- pmatch(dataClin_TP$bcr_patient_barcode, SampleTP) # SampleTP_short
    SampleTP <- SampleTP[SampleTP_matched]
    
    # normal tissue
    dataClin_NT <- dplyr::filter(dataClin, dataClin$bcr_patient_barcode %in% SampleNT_short)
    
    # prepare clinical information
    dataClin <- rbind(dataClin_NT, dataClin_TP)
    saveRDS(dataClin, file.path(proj.dir, "data", paste0("dataClin_", proj, ".RDS")))
    
    # if there's less then 3 normal tissue samples with associated clinical data skip the project
    if (nrow(dataClin_NT) < 3) {
      next
    } 
    else {
      SampleNT_matched <- pmatch(dataClin_NT$bcr_patient_barcode, SampleNT) # SampleNT_short
      SampleNT <- SampleNT[SampleNT_matched]
      query_down <- GDCquery(project = proj,
                             data.category = "Gene expression",
                             data.type = "Gene expression quantification",
                             platform = "Illumina HiSeq",
                             file.type = "results",
                             barcode = c(SampleTP, SampleNT),
                             experimental.strategy = "RNA-Seq",
                             legacy = TRUE)
      saveRDS(query_down, file.path(proj.dir, "data", paste0("GDCquery_", proj, ".RDS")))
      
      # download data
      tryCatch(GDCdownload(query_down, 
                           method = "api", 
                           files.per.chunk = 20,
                           directory = file.path(proj.dir, "data", "GDCdata")),
               error = function(e) GDCdownload(query_down, 
                                               method = "client", 
                                               files.per.chunk = 20,
                                               directory = file.path(proj.dir, "data", "GDCdata")))
      
      # read the data downloaded and save it as an R object
      flog.debug("Prepare GDC data")
      prep <- GDCprepare(query_down, 
                         directory = file.path(proj.dir, "data", "GDCdata"))
      saveRDS(prep, file.path(proj.dir, "data", paste0("prep_", proj, ".RDS")))
      
      # Array Array Intensity Correlation (AAIC) 
      # - identification of abnormal samples
      flog.debug("Perform intensity correlation")
      dataPreproc <- TCGAanalyze_Preprocessing(object = prep, cor.cut = 0.6)
      saveRDS(prep, file.path(proj.dir, "data", paste0("prep_", proj, ".RDS")))
      
      # Normalization for GC content
      flog.debug("Perform normalization based on GC content")
      dataNorm_GC <- TCGAanalyze_Normalization(tabDF = dataPreproc, 
                                               geneInfo = geneInfo, method = "gcContent")
      saveRDS(dataNorm_GC, file.path(proj.dir, "data", paste0("normGC_", proj, ".RDS")))
      flog.debug("Perform normalization based on gene length")
      dataNorm_length <- TCGAanalyze_Normalization(tabDF = dataPreproc, 
                                                   geneInfo = geneInfo, method = "geneLength")
      saveRDS(dataNorm_length, file.path(proj.dir, "data", paste0("normGeneLength_", proj, ".RDS")))
      
      # prepare data frame with reads for cancer and normal samples 
      # - identification of abnormally expressed genes
      flog.debug("Perform data filtering based on threshold defined quantile mean across all samples")
      
      flog.debug("Filtering for samples normalized based on GC content")
      dataFilt_GC <- assign(paste0("dataFilt_GC_", proj), 
                            TCGAanalyze_Filtering(tabDF = dataNorm_GC, 
                                                  method = "quantile", qnt.cut =  0.25))
      saveRDS(dataFilt_GC, file.path(proj.dir, "data", paste0("dataFilt_GC_", proj, ".RDS")))
      
      flog.debug("Filtering for samples normalized based on gene length")
      dataFilt_length <- assign(paste0("dataFilt_length_", proj), 
                                TCGAanalyze_Filtering(tabDF = dataNorm_length,
                                                      method = "quantile", qnt.cut =  0.25))
      saveRDS(dataFilt_length, file.path(proj.dir, "data", paste0("dataFilt_length_", proj, ".RDS")))
      
      
      # create a representative patient for the project with a mean number of reads
      # for the given gene in normal and primary tumor samples 
      flog.debug("Averaged patient created from data normalized based on GC content")
      
      dataFilt_GC_transposed <- t(dataFilt_GC)
      dataFilt_GC_transposed <- as.data.frame(dataFilt_GC_transposed)
      dataFilt_GC_transposed <- rownames_to_column(dataFilt_GC_transposed, "barcode_id")
      
      # normal tissue
      dataFilt_GC_transposed_NT <- filter(dataFilt_GC_transposed, 
                                       dataFilt_GC_transposed$barcode_id %in% SampleNT)
      dataFilt_GC_transposed_NT <- column_to_rownames(dataFilt_GC_transposed_NT, "barcode_id")
      dataFilt_GC_transposed_NT <- t(dataFilt_GC_transposed_NT)
      dataFilt_GC_transposed_NT <- as.data.frame(dataFilt_GC_transposed_NT)
      dataFilt_GC_transposed_NT <- rownames_to_column(dataFilt_GC_transposed_NT, "Symbol")
      dataFilt_GC_transposed_NT <- dplyr::mutate(dataFilt_GC_transposed_NT, 
                                              Av_NT = rowMeans(
                                                dataFilt_GC_transposed_NT[, c(2:ncol(dataFilt_GC_transposed_NT))]))
      dataFilt_GC_transposed_NT <- dataFilt_GC_transposed_NT[, c(1, ncol(dataFilt_GC_transposed_NT))]
      
      # tumor tissue
      dataFilt_GC_transposed_TP <- filter(dataFilt_GC_transposed, 
                                       dataFilt_GC_transposed$barcode_id %in% SampleTP)
      dataFilt_GC_transposed_TP <- column_to_rownames(dataFilt_GC_transposed_TP, "barcode_id")
      dataFilt_GC_transposed_TP <- t(dataFilt_GC_transposed_TP)
      dataFilt_GC_transposed_TP <- as.data.frame(dataFilt_GC_transposed_TP)
      dataFilt_GC_transposed_TP <- rownames_to_column(dataFilt_GC_transposed_TP, "Symbol")
      dataFilt_GC_transposed_TP <- dplyr::mutate(dataFilt_GC_transposed_TP, 
                                              Av_TP = rowMeans(
                                                dataFilt_GC_transposed_TP[, c(2:ncol(dataFilt_GC_transposed_TP))]))
      dataFilt_GC_transposed_TP <- dataFilt_GC_transposed_TP[, c(1, ncol(dataFilt_GC_transposed_TP))]
      
      # combine data frames from normal and tumor tissue
      dataFilt_GC_rep_patient <- assign(paste0("dataFilt_GC_rep_patient_", proj), 
                                        merge(dataFilt_GC_transposed_NT, 
                                              dataFilt_GC_transposed_TP, 
                                              by = "Symbol"))
      saveRDS(dataFilt_GC_rep_patient, 
              file.path(proj.dir, "data", paste0("dataFilt_GC_rep_patient_", proj, ".RDS")))
      
      # query for normal and tumor samples in pre-processed, normalized and filtrated data
      SampleNT_GC <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_GC), 
                                           typesample = "NT")
      SampleTP_GC <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_GC), 
                                           typesample = "TP")

      # create a representative patient for the project with a mean number of reads
      # for the given gene in normal and primary tumor samples 
      flog.debug("Averaged patient created from data normalized based on gene length")
      
      dataFilt_length_transposed <- t(dataFilt_length)
      dataFilt_length_transposed <- as.data.frame(dataFilt_length_transposed)
      dataFilt_length_transposed <- rownames_to_column(dataFilt_length_transposed, "barcode_id")
      
      # normal tissue
      dataFilt_length_transposed_NT <- filter(dataFilt_length_transposed, 
                                          dataFilt_length_transposed$barcode_id %in% SampleNT)
      dataFilt_length_transposed_NT <- column_to_rownames(dataFilt_length_transposed_NT, "barcode_id")
      dataFilt_length_transposed_NT <- t(dataFilt_length_transposed_NT)
      dataFilt_length_transposed_NT <- as.data.frame(dataFilt_length_transposed_NT)
      dataFilt_length_transposed_NT <- rownames_to_column(dataFilt_length_transposed_NT, "Symbol")
      dataFilt_length_transposed_NT <- dplyr::mutate(dataFilt_length_transposed_NT, 
                                                 Av_NT = rowMeans(
                                                   dataFilt_length_transposed_NT[, c(2:ncol(dataFilt_length_transposed_NT))]))
      dataFilt_length_transposed_NT <- dataFilt_length_transposed_NT[, c(1, ncol(dataFilt_length_transposed_NT))]
      
      # tumor tissue
      dataFilt_length_transposed_TP <- filter(dataFilt_length_transposed, 
                                          dataFilt_length_transposed$barcode_id %in% SampleTP)
      dataFilt_length_transposed_TP <- column_to_rownames(dataFilt_length_transposed_TP, "barcode_id")
      dataFilt_length_transposed_TP <- t(dataFilt_length_transposed_TP)
      dataFilt_length_transposed_TP <- as.data.frame(dataFilt_length_transposed_TP)
      dataFilt_length_transposed_TP <- rownames_to_column(dataFilt_length_transposed_TP, "Symbol")
      dataFilt_length_transposed_TP <- dplyr::mutate(dataFilt_length_transposed_TP, 
                                                 Av_TP = rowMeans(
                                                   dataFilt_length_transposed_TP[, c(2:ncol(dataFilt_length_transposed_TP))]))
      dataFilt_length_transposed_TP <- dataFilt_length_transposed_TP[, c(1, ncol(dataFilt_length_transposed_TP))]
      
      # combine data frames from normal and tumor tissue
      dataFilt_length_rep_patient <- assign(paste0("dataFilt_length_rep_patient_", proj), 
                                            merge(dataFilt_length_transposed_NT, 
                                                  dataFilt_length_transposed_TP, 
                                                  by = "Symbol"))
      saveRDS(dataFilt_length_rep_patient, 
              file.path(proj.dir, "data", paste0("dataFilt_length_rep_patient_", proj, ".RDS")))
      
      # query for normal and tumor samples in pre-processed, normalized and filtrated data
      SampleNT_length <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_length), 
                                               typesample = "NT")
      SampleTP_length <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_length), 
                                               typesample = "TP")
      
      # perform differential expression analysis
      # with contrasts between tumor and normal samples and batch normalization for tissue source site 
      flog.debug("Perform differential gene expression for GC-normalized data")
      dataDEG_GC <- TCGAanalyze_DEA(mat1 = dataFilt_GC[,SampleNT_GC], 
                                    mat2 = dataFilt_GC[,SampleTP_GC],
                                    Cond1type = "Normal", 
                                    Cond2type = "Tumor",
                                    batch.factors = "TSS", 
                                    fdr.cut = 1, logFC.cut = 0, 
                                    method = "glmLRT")
      dfDEG_GC <- assign(paste0("dfDEG_", proj), 
                         tibble::rownames_to_column(dataDEG_GC, var = "Symbol"))
      saveRDS(dfDEG_GC, 
              file.path(proj.dir, "data", paste0("dfDEG_GC_", proj, ".RDS")))
      
      flog.debug("Perform differential gene expression for gene length-normalized data")
      dataDEG_length <- TCGAanalyze_DEA(mat1 = dataFilt_length[,SampleNT_length], 
                                        mat2 = dataFilt_length[,SampleTP_length],
                                        Cond1type = "Normal", 
                                        Cond2type = "Tumor",
                                        batch.factors = "TSS", 
                                        fdr.cut = 1, logFC.cut = 0,
                                        method = "glmLRT")
      dfDEG_length <- assign(paste0("dfDEG_length_", proj), 
                      tibble::rownames_to_column(dataDEG_length, var = "Symbol"))
      saveRDS(dfDEG_length, 
              file.path(proj.dir, "data", paste0("dfDEG_length_", proj, ".RDS")))
    }
  }
}