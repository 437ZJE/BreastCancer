library(TCGAbiolinks)
library(dplyr)
library(ggplot2)

CancerProject <- "TCGA-BRCA"
DataDirectory <- paste0("./GDC/",gsub("-","_",CancerProject))
FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")

query <- GDCquery(project = CancerProject, 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

dataAssy.sub <- TCGAquery_subtype(tumor = gsub("TCGA-","",CancerProject))

samplesDown <- getResults(query,cols=c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")
dataSmTP <- dataSmTP[substr(dataSmTP,1,12) %in% dataAssy.sub$patient]

#dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
#                                  typesample = "NT")


queryDown <- GDCquery(project = CancerProject, 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = dataSmTP)
GDCdownload(query = queryDown,
            directory = DataDirectory)
dataPrep <- GDCprepare(query = queryDown, 
                       save = TRUE, 
                       directory =  DataDirectory,
                       save.filename = FileNameData)

dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")                      

#get normalized counts
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent") 

#covert EnsembleID to gene symbol
gene.annot <- read.csv("/home/aa/Data/annotation/Homo_sapiens.GRCh38/GRCh38.94.tx2gene.csv",header = F)
names(gene.annot) <- c("Tx","GeneEsID","GeneSymbol")
gene.annot <- unique(gene.annot[,2:3])
dataNorm2 <- as.data.frame(dataNorm)
dataNorm2$ID <- rownames(dataNorm2)
dataNorm2 <- left_join(dataNorm2, gene.annot, by=c("ID" = "GeneEsID"))
dataNorm2$GeneSymbol[dataNorm2$GeneSymbol == "LINC-PINT"] <- c("AC058791.1", "LINC-PINT") #The annotation file has duplicates for "LINC-PINT"
rownames(dataNorm2) <- dataNorm2$GeneSymbol
dataNorm2 <- select(dataNorm2,-c(ID,GeneSymbol))
saveRDS(dataNorm2, file = "TCGA_BRCA_HTSeq_Normalized_Counts.rds")

#get normalized counts to geneLength
dataNorm3 <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "geneLength") 
dataNorm3 <- as.data.frame(dataNorm3)
dataNorm3$ID <- rownames(dataNorm3)
dataNorm3 <- left_join(dataNorm3, gene.annot, by=c("ID" = "GeneEsID"))
dataNorm3$GeneSymbol[dataNorm3$GeneSymbol == "LINC-PINT"] <- c("AC058791.1", "LINC-PINT")
rownames(dataNorm3) <- dataNorm3$GeneSymbol
dataNorm3 <- select(dataNorm3,-c(ID,GeneSymbol))
saveRDS(dataNorm3, file = "TCGA_BRCA_HTSeq_TPM.rds")


#generate clin data
dataNorm2 <- readRDS("TCGA_BRCA_HTSeq_Normalized_Counts.rds")
clin <- data.frame(sample = colnames(dataNorm2), patient =substr(colnames(dataNorm2), 1, 12) )
clin <- left_join(clin, dataAssy.sub)
clin <- select(clin, sample, patient, BRCA_Subtype_PAM50,vital_status, days_to_death, days_to_last_followup )
clin$vital_status <- ifelse(clin$vital_status=="Alive", 0 ,1)
times <- c()
for (i in 1:nrow(clin)){
  times[i] <- ifelse(clin$days_to_death[i]=="NA", as.numeric(clin$days_to_last_followup[i]), as.numeric(clin$days_to_death[i]))
 
}
clin$times <- times
saveRDS(clin,"BRCA_clinical_data.rds")

