# setwd('R/GDCData/')

#### Logging ####
### Create a log file to print in a file instead of the R's console
logFile = file('LOG.Rout', open = 'a')
sink(logFile)

source('Functions.R')

#### Initialization ####
### Download & load required packages
Initialize()

#### ID Preparation ####
### Import ENSEMBL RNAs ID, The IDs are retrieved from ENSEMBL Biomart
ensemblIDs=read.table('Files/AllRNAsID_GDC.txt',header = T)

### Import mirBase miRNAs ID, The IDs are retrieved from mirBase
mirbaseNameAndIDs = read.table('Files/AllMiRNAsNameAndID_miRBase_v21.txt', header = T)

#### Get & Set Projects Info & Data Type ####
GDCProjects = as.data.frame(getGDCprojects())
GDCProjectsID = GDCProjects$project_id

### Specify different data types to download
dataTypes = c('Counts', 'FPKM', 'miRNA')

#### Find, Query, Download, Import & Save GDC Data for each Project ####
for(selectedProject in GDCProjectsID)
{
  ### Create an empty data structure appropriate for GDC data
  GDCData = GetAnEmptyGDCDataStructure()
  
  ### Find desired cases for the project
  ### Desired cases are the ones whit both RNA-Seq & miRNA-Seq data
  selectedCases = FindDesiredCases(selectedProject)
  
  
  ### Query, download & import data for each data type separately
  for(dataType in dataTypes)
  {
    ### Query & download the data of the desired cases and selected data types
    QueryAndDownloadData(selectedProject,
                         selectedCases,
                         dataType)
    
    ### Import (read) the raw text files
    ### Each sample is a text (document) file contain the data
    ImportData(selectedProject, dataType)
  }
  
  ### Query, download & impor the clinical data of the desired cases
  queryOfClinical = QueryAndDownloadClinicalData(selectedProject, selectedCases)
  ImportClinicalData(selectedProject, queryOfClinical)
  
  ### Save the final filled GDC data structure of the project as an RDS file.
  SaveGDCData(selectedProject)
}

#### Arrange Data ####
### Arrangment means align the case IDs between RNA-Seq, miRNA-seq & Clinical data.
GDCDataRepository = ArrangeData()

saveRDS(GDCData, file = 'GDCDataRepositpory/GDCDataRepository.rds')

#### Merge GDC Data Repository ####
### Merge the GDC data for all projects
GDCData = MergeAllGDCData(GDCDataRepository)

saveRDS(GDCData, file = 'GDCDataRepositpory/GDCData.rds')

#### Generate & Get GDC Data Summary ####
GDCDataSummary = GetGDCDataSummary()

GDCData = readRDS('GDCDataRepositpory/GDCData.rds')

#### Fix A Clinical Label Issue ####
GDCData[['Clinical']][GDCData[['Clinical']]$type=='Normal',]$disease = 'Normal'

#### Get & Add GDC RNAs Info ####
AddGenesInfoToGDCData()
saveRDS(GDCData, file = 'GDCDataRepositpory/GDCData_v2.rds')

#### Extract Specific RNA Type ####
selectedRNATypeData = ExtractSpecificRNATypeData(selectedRNAType = 'protein_coding')

GDCProteinCodingData = GDCData
for(dataType in names(selectedRNATypeData))
  GDCProteinCodingData[[dataType]] = selectedRNATypeData[[dataType]]

saveRDS(GDCProteinCodingData, file = 'GDCDataRepositpory/GDCData_ProteinCoding.rds')

#### Harmonize Clinical Data ####
GDCData[['Clinical']] = HarmonizeClinicalData(GDCData[['Clinical']])

saveRDS(GDCData, file = 'GDCDataRepositpory/GDCData_ProteinCoding_v2.rds')


firstNetworkClinicalFields = c('id', 'tissue', 'type', 'disease')
firstNetworkClinicals = ExtractSpecificClinicalFields(GDCData[['Clinical']],
                                                      firstNetworkClinicalFields,
                                                      tumorsOnly = F,
                                                      removeNAs = T)

secondNetworkClinicalFields = c('id', 'disease', 'stage', 'vitalStatus', 'overallSurvivalTimeInYears',
                                'isRelapsed', 'relapseTimeInYears', 'yearsBeingAliveAfterRelapse')
secondNetworkClinicals = ExtractSpecificClinicalFields(GDCData[['Clinical']],
                                                       secondNetworkClinicalFields,
                                                       tumorsOnly = T,
                                                       removeNAs = T)

GDCDataForFirstNetwork = GDCData
GDCDataForFirstNetwork[['Clinical']] = firstNetworkClinicals
GDCDataForFirstNetwork[['AllClinical']] = NULL


GDCDataForSecondNetwork = GDCData
GDCDataForSecondNetwork[['Clinical']] = secondNetworkClinicals
dataTypes = c('Counts','FPKM','miRNA')
for(dataType in dataTypes)
  GDCDataForSecondNetwork[[dataType]] = GDCDataForSecondNetwork[[dataType]][, GDCDataForSecondNetwork[['Clinical']]$id]
GDCDataForSecondNetwork[['AllClinical']] = NULL

saveRDS(GDCDataForFirstNetwork, file = 'GDCDataRepositpory/GDCData_FirstNetwork.rds')
saveRDS(GDCDataForSecondNetwork, file = 'GDCDataRepositpory/GDCData_SecondNetwork.rds')


######################## #
#### Single Cell Data ####
######################## #

GDCData = readRDS(file = 'GDCDataRepositpory/GDCData_SecondNetwork.rds')

ensemblInfo = read.delim('Files/Ensembl_RNAs_Info.txt')

### Sample 1
s1Counts = GDCData$Counts$`TCGA-OR-A5J1-01A`
s1Fpkm = GDCData$FPKM$`TCGA-OR-A5J1-01A`
s1Lengths = (s1Counts*1e9)/(s1Fpkm*sum(s1Counts))
# s1Info = data.frame(Counts = s1Counts, FPKM = s1Fpkm, CountsB = (s1Counts*1e9), FPKMSum = (s1Fpkm*sum(s1Counts)), Length = s1Lengths)

### Sample 2
s2Counts = GDCData$Counts$`TCGA-DK-AA76-01A`
s2Fpkm = GDCData$FPKM$`TCGA-DK-AA76-01A`
s2Lengths = (s2Counts*1e9)/(s2Fpkm*sum(s2Counts))

allSInfo = data.frame(s1 = s1Lengths, s2 = s2Lengths)

rownames(allSInfo) = rownames(GDCData$Counts)
allSInfo$id = rownames(allSInfo)

ensemblExonInfo = read.delim(file = 'Files/Gene_Exon_Start_End.gz')
names(ensemblExonInfo) = c('id', 'exonID', 'start', 'end', 'rank')
ensemblExonInfo$rank = NULL
ensemblExonInfo = unique(ensemblExonInfo)

exonRanges = IRanges(ensemblExonInfo$start, ensemblExonInfo$end)

limitOfRanges = lapply(unique(ensemblExonInfo$id), function(gene) 
  range(which(ensemblExonInfo$id == gene)))

limitOfRanges = do.call(rbind, limitOfRanges)
names(limitOfRanges) = c('sIndex', 'eIndex')

genesLength = mapply(function(sIndex, eIndex){
  sum(width(reduce(exonRanges[sIndex:eIndex])))
}, limitOfRanges[,1], limitOfRanges[,2])

genesLengthDF = data.frame(id = unique(ensemblExonInfo$id), length = genesLength)


# ensemblExonInfo$ranges = IRanges(ensemblExonInfo$Exon.region.start..bp., ensemblExonInfo$Exon.region.end..bp.)
# ensemblExonInfo = aggregate(ensemblExonInfo, list('ranges'), 'reduce')
# 
# ensemblExonInfo = ensemblExonInfo[, c(1, 5)]
# ensemblExonInfo = as.data.table(ensemblExonInfo)
# setnames(ensemblExonInfo, 'Gene.stable.ID', 'id')
# ensemblGenesLength = aggregate(length ~ id, ensemblExonInfo, sum)

genesLengthDF = merge(allSInfo, genesLengthDF, all.x =T)
genesLengthDF$s1 = NULL
names(genesLengthDF) = c('id', 'byGDCData', 'byEnsemblBiomart')


# temp = ensemblExonInfo[1:100,]
# temp$ranges = IRanges(temp$Exon.region.start..bp., temp$Exon.region.end..bp.)
# 
# finalLengths = data.frame(gene=character(), length = integer())
# for(gene in unique(temp$Gene.stable.ID))
# {
#   rangesUnion = temp$ranges[temp$Gene.stable.ID==gene][1]
#   for(i in seq(2, length(temp$ranges[temp$Gene.stable.ID==gene])))
#   {
#     rangesUnion = union(rangesUnion, temp$ranges[temp$Gene.stable.ID==gene][i])
#   }
#    
#   finalLengths = rbind(finalLengths, sum(width(rangesUnion)))
# }
# 
# 
# tempUnion = tempRanges[1]
# for(i in seq(2, length(tempRanges)))
# {
#   tempUnion = union(tempUnion, tempRanges[i])
# }
#   
# sum(width(tempUnion))
# 
# temp$ranges = IRanges(temp$Exon.region.start..bp., temp$Exon.region.end..bp.)
# temp$ranges
# 
# head(temp[, c(1,6)])
# 
# tempUnion = union(tempRanges[1],tempRanges[2])
# tempUnion
# 
# reduce(tempRanges)

scData = read.delim('Files/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt.gz')
scData$gene_id = substring(scData$gene_id, 1, 15)
scData = scData[scData$gene_type == 'protein_coding', ]

scData = scData[, -c(2,3)]
setnames(scData, 'gene_id', 'id')

head(scData)

gdcGenesID = rownames(GDCData$Counts)
gdcGenesID = data.frame(id = gdcGenesID)
scData = merge(gdcGenesID, scData, all.x=T)
scData[is.na(scData)] = 0

rownames(scData) = scData$id
scData$id = NULL

lengthData = merge(gdcGenesID, genesLengthDF, all.x =T)

lBar = sweep(scData, 1, lengthData$byEnsemblBiomart, '*')
lBar = colSums(lBar)/1e6
fpkmValue = sweep(scData, 2, lBar, '/')
fpkmValue = fpkmValue*1e3

GSE75688Data = list()
GSE75688Data[['FPKM']] = fpkmValue
GSE75688Data[['TPM']] = scData
GSE75688Data[['LengthInfo']] = lengthData

saveRDS(GSE75688Data, 'GEO/GSE75688Data.rds')

########################################## #
#### Separating data for second network ####
########################################## #

GDCData = readRDS(file = 'GDCDataRepositpory/GDCData_SecondNetwork_v2.rds')
GDCData_FN = readRDS('GDCDataRepositpory/GDCData_FirstNetwork.rds')

fnCounts = GDCData_FN$Counts
fnCounts = rbind(fnCounts, seq(ncol(fnCounts)))
snCounts = fnCounts[, GDCData$Clinical$id]
selectedSamples = snCounts[nrow(snCounts), ]
selectedSamples = as.integer(selectedSamples)

firstNet = lapply(list.files('Files/FirstNetwork_Output/', full.names = T), read.csv, header = F)
names(firstNet) = list.files('Files/FirstNetwork_Output/')

firstNet = lapply(firstNet, function(x) x[selectedSamples, ])

for(aCSV in names(firstNet))
  write.table(firstNet[[aCSV]], row.names = F, col.names = F, quote = F, sep = ',', file = aCSV)

saveRDS(firstNet, 'SecondNetwork_Input.rds')

summary(firstNet)

# x = colnames(GDCData_FN$Counts)[selectedSamples]
# y = GDCData$Clinical$id
# identical(x, y)
# identical(colnames(GDCData$Counts)[selectedSamplesIndex], GDCData$Clinical$id)
# 
# head(x)
# head(y)
# 
# 
# 
# x= colnames(GDCData_FN$Counts)[selectedSamplesIndex]
# y= GDCData$Clinical$id
# setdiff(y, x)
# 
# identical(colnames(GDCData_FN$Counts)[selectedSamplesIndex], GDCData$Clinical$id)
# head(colnames(GDCData_FN$Counts)[selectedSamplesIndex])
# head(GDCData$Clinical$id)
# 
# allSamples = colnames(GDCData_FN$Counts)
# temp = strsplit2(allSamples, split = '')
# 
# snSamples = colnames(GDCData$Counts)
# 
# length(which(allSamples%in%snSamples))
# 
# commonSamples = allSamples[allSamples%in%snSamples]
# setdiff(snSamples, commonSamples)
# 
# identical(GDCData$Clinical$id, GDCData_FN$Clinical$id[rownames(GDCData$Clinical)])
# head(GDCData$NameMapping)



############################################################# #
#### Extract All Clinical Labels to Find New Useful Labels ####
############################################################# #

GDCData = readRDS(file = 'GDCDataRepositpory/GDCDataRepository.rds')

allClinicalLabels = character()
for(i in names(GDCData))
{
  clinicalLabelsOfAProject = colnames(GDCData[[i]]$Clinical)
  diseaseOfTheProject = unique(GDCData[[i]]$Clinical$disease)
  clinicalLabelsOfAProject = c(diseaseOfTheProject, clinicalLabelsOfAProject)
  
  clinicalLabelsOfAProject = paste0(clinicalLabelsOfAProject, collapse = '\t')
  
  allClinicalLabels = c(allClinicalLabels, clinicalLabelsOfAProject)
}

write.table(allClinicalLabels, file = 'AllClinicalLabels.txt',
            quote = F, sep = '\n', row.names = F, col.names = F)


tcgaClinicalLabels = character()
targetClinicalLabels = character()

tcgaCommonClinicalLabels = character()
targetCommonClinicalLabels = character()

tcgaIsFirst = T
targetIsFirst = T

for(i in names(GDCData))
{
  clinicalLabelsOfAProject = colnames(GDCData[[i]]$Clinical)
  
  if(grepl('TCGA', i))
  {
    tcgaClinicalLabels = c(tcgaClinicalLabels, clinicalLabelsOfAProject)
    
    if(tcgaIsFirst)
    {
      tcgaCommonClinicalLabels = clinicalLabelsOfAProject
      tcgaIsFirst = F
    } else {
      tcgaCommonClinicalLabels = intersect(tcgaCommonClinicalLabels, clinicalLabelsOfAProject)
    }
    
  } else {
    targetClinicalLabels = c(targetClinicalLabels, clinicalLabelsOfAProject)
    
    if(targetIsFirst)
    {
      targetCommonClinicalLabels = clinicalLabelsOfAProject
      targetIsFirst = F
    } else {
      targetCommonClinicalLabels = intersect(targetCommonClinicalLabels, clinicalLabelsOfAProject)
    }
  }
  
}
tcgaFreqs = as.data.frame(table(tcgaClinicalLabels))
tcgaFreqs = tcgaFreqs[order(tcgaFreqs$Freq, decreasing = T), ]

targetFreqs = as.data.frame(table(targetClinicalLabels))
targetFreqs = targetFreqs[order(targetFreqs$Freq, decreasing = T), ]


tcgaLabelSamplesNumber = integer()
# label = tcgaFreqs$tcgaClinicalLabels[1]
for(label in tcgaFreqs$tcgaClinicalLabels)
{
  numOfAllSamplesOfALabel = 0
  for(i in names(GDCData))
    if(grepl('TCGA', i))
    {
      result = GDCData[[i]]$Clinical[[label]]
      result = result[!is.na(result)]
      numOfAllSamplesOfALabel = numOfAllSamplesOfALabel + length(result)
    }
  tcgaLabelSamplesNumber = c(tcgaLabelSamplesNumber, numOfAllSamplesOfALabel)  
}


targetLabelSamplesNumber = integer()
# label = targetFreqs$targetClinicalLabels[1]
for(label in targetFreqs$targetClinicalLabels)
{
  numOfAllSamplesOfALabel = 0
  for(i in names(GDCData))
    if(grepl('TARGET', i))
    {
      result = GDCData[[i]]$Clinical[[label]]
      result = result[!is.na(result)]
      numOfAllSamplesOfALabel = numOfAllSamplesOfALabel + length(result)
    }
  targetLabelSamplesNumber = c(targetLabelSamplesNumber, numOfAllSamplesOfALabel) 
}

tcgaFreqs$NumOfSamples = tcgaLabelSamplesNumber
targetFreqs$NumOfSamples = targetLabelSamplesNumber


write.table(tcgaCommonClinicalLabels, file = 'ClinicalLabelsStatistics/tcgaCommonClinicalLabels.txt',
            quote = F, sep = '\n', row.names = F, col.names = F)
write.table(targetCommonClinicalLabels, file = 'ClinicalLabelsStatistics/targetCommonClinicalLabels.txt',
            quote = F, sep = '\n', row.names = F, col.names = F)

write.table(tcgaFreqs, file = 'ClinicalLabelsStatistics/tcgaFreqs.txt',
            quote = F, sep = '\t', row.names = F, col.names = F)
write.table(targetFreqs, file = 'ClinicalLabelsStatistics/targetFreqs.txt',
            quote = F, sep = '\t', row.names = F, col.names = F)


#### Number of TCGA samples ####
count = 0
for(i in names(GDCData))
{
  if(grepl('TCGA', i))
  {
    count = count + nrow(GDCData[[i]]$Clinical)
  }
}


############################## #
#### Cancer Staging Systems ####
############################## #
### The American Joint Committee on Cancer: Cancer Staging
### https://cancerstaging.org

GDCData = readRDS(file = 'GDCDataRepositpory/GDCDataRepository.rds')

### Tumor–Node–Metastasis (TNM) Cancer Staging System
### https://www.cancer.gov/about-cancer/diagnosis-staging/staging
GDCData[[i]]$Clinical$stage_event_tnm_categories

### https://www.cancer.gov/about-cancer/diagnosis-staging/prognosis/tumor-grade-fact-sheet
GDCData[[i]]$Clinical$neoplasm_histologic_grade

### Porstate Cancer : https://www.cancer.gov/about-cancer/diagnosis-staging/prognosis/tumor-grade-fact-sheet
GDCData[['TCGA-PRAD']]$Clinical$stage_event_gleason_grading

### http://oncologypro.esmo.org/Oncology-in-Practice/Practice-Tools/Performance-Scales
GDCData[[i]]$Clinical$karnofsky_performance_score


#### Generating a Table of Data Summary ####
GDCProjects = as.data.frame(getGDCprojects())
GDCProjectsID = GDCProjects$project_id

GDCData = readRDS('GDCDataRepositpory/GDCData.rds')
GDCDataSummary = GetGDCDataSummary()

gdcprojs = GDCProjects[, c('primary_site', 'disease_type', 'project_id')]
names(gdcprojs) = c('tissue', 'disease', 'project')
dataSummary = merge(GDCDataSummary, gdcprojs, all.x = T, sort = F)

dataSummary = dataSummary[, c(2,5,1,3,4)]
dataSummary$miRNA.Tumor = dataSummary$tumor
dataSummary$miRNA.Normal = dataSummary$normal

names(dataSummary) = c('Disease', 'Project', 'Tissue',
                       'mRNA.Tumor', 'mRNA.Normal',
                       'miRNA.Tumor', 'miRNA.Normal')
DataSummary = dataSummary
DataSummary[, -seq(3)] = apply(DataSummary[, -seq(3)], 2, as.integer)
DataSummary = DataSummary[DataSummary$Disease!='All', ]

saveRDS(DataSummary, 'GDCDataRepositpory/GDCDataSummry.rds')
write.xlsx(DataSummary, file = 'GDCDataRepositpory/GDCDataSummary.xlsx', row.names = F)

