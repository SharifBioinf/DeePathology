Initialize = function()
{
  # source('https://bioconductor.org/biocLite.R')
  # if(!('TCGAbiolinks' %in% rownames(installed.packages())))
  #   biocLite('TCGAbiolinks', suppressUpdates = T)
  # if(!('limma' %in% rownames(installed.packages())))
  #   biocLite('limma', suppressUpdates = T)
  # if(!('xlsx' %in% rownames(installed.packages())))
  #   biocLite('xlsx', suppressUpdates = T)
  
  ### 'TCGAbiolinks' is an R interface for GDC's API
  library(TCGAbiolinks)
  
  ### 'limma' is loaded to use it's 'strsplit2' function
  library(limma)
  
  ### 'parallel' is used to import (read) raw text files in parallel
  library(parallel)
  
  ### 'xlsx' is used to import (read) TARGET clinical Excel files.
  ### TCGA clinical files are XML, but TARGET clinical files are Excel.
  library(xlsx)
  
  options(stringsAsFactors = F)
}


GetProjectInfo = function(selectedProject, getFullNameOnly=F)
{
  GDCProjects=as.data.frame(getGDCprojects())
  projectInfo = GDCProjects[GDCProjects$project_id==selectedProject, c('id',
                                                                       'primary_site.0',
                                                                       'disease_type.0')]
  names(projectInfo) = c('project', 'tissue', 'disease')
  projectFullName=paste0('[',projectInfo$tissue,'] ', projectInfo$disease)
  
  projectInfo = if(getFullNameOnly) projectFullName else projectInfo
  projectInfo
}



FindDesiredCases = function(selectedProject)
{
  print(selectedProject)
  
  tumorSampleTypes = 'Primary solid Tumor'
  normalSampleTypes = 'Solid Tissue Normal'
  
  ### The tumor and normal types of projects which thier primary site is Blood-based, are different.
  if(GDCProjects[GDCProjects$id==selectedProject, 'primary_site.0'] %in% c('Blood','Bone Marrow'))
  {
    tumorSampleTypes = c(tumorSampleTypes,
                         'Primary Blood Derived Cancer - Peripheral Blood',
                         'Primary Blood Derived Cancer - Bone Marrow')
    normalSampleTypes = c(normalSampleTypes,
                          'Blood Derived Normal',
                          'Bone Marrow Normal')
  }
  
  ### Query the RNA-Seq data
  queryOfGeneExpressionCounts = GDCquery(project = selectedProject,
                                         data.category = 'Transcriptome Profiling',
                                         data.type = 'Gene Expression Quantification',
                                         workflow.type = 'HTSeq - Counts',
                                         sample.type = c(tumorSampleTypes, normalSampleTypes))
  
  resultOfGeneExpressionCounts = getResults(queryOfGeneExpressionCounts)
  
  ### Retrieve cases which have RNA-Seq data
  countsCases = resultOfGeneExpressionCounts[, 'cases']
  if(grepl('TCGA', selectedProject))
    countsCases = substr(countsCases, 1, 16)
  
  ### Query the miRNA-Seq data
  queryOfMiRNAExpressionIsoform = GDCquery(project = selectedProject,
                                           data.category = 'Transcriptome Profiling',
                                           data.type = 'Isoform Expression Quantification',
                                           sample.type = c(tumorSampleTypes, normalSampleTypes))
  
  resultOfMiRNAExpressionIsoform = getResults(queryOfMiRNAExpressionIsoform)
  
  ### Retrieve cases which have miRNA-Seq data
  mirCases = resultOfMiRNAExpressionIsoform[, 'cases']
  if(grepl('TCGA', selectedProject))
    mirCases = substr(mirCases, 1, 16)
  
  selectedCases = intersect(countsCases, mirCases)
  
  selectedCases
}



QueryAndDownloadData = function(selectedProject, selectedCases, dataType)
{
  print(dataType)
  
  projectInfo = GetProjectInfo(selectedProject)
  
  if(dataType=='Counts')
  {
    data.type = 'Gene Expression Quantification'
    workflow.type = 'HTSeq - Counts'
  }
  else if(dataType=='FPKM')
  {
    data.type = 'Gene Expression Quantification'
    workflow.type = 'HTSeq - FPKM'
  }
  else if(dataType=='miRNA')
  {
    data.type = 'Isoform Expression Quantification'
    workflow.type = 'BCGSC miRNA Profiling'
  }
  
  queryOfData = GDCquery(project = selectedProject,
                         data.category = 'Transcriptome Profiling',
                         data.type = data.type,
                         workflow.type = workflow.type,
                         barcode = selectedCases)
  
  
  GDCdownload(queryOfData)
  resultOfData = getResults(queryOfData)
  
  importantCols = c('cases', 'id', 'submitter_id', 'tissue.definition', 'project')
  
  dataInfo = resultOfData[, importantCols]
  
  if(grepl('TCGA', selectedProject))
    dataInfo$cases = substr(dataInfo$cases, 1, 20)
  dataInfo$submitter_id = substr(dataInfo$submitter_id, 1, 36)
  
  dataInfoAndMapping = merge(dataInfo, projectInfo, by = 'project', sort = F)
  dataInfo = dataInfoAndMapping[,c('cases', 'tissue', 'tissue.definition', 'disease')]
  names(dataInfo) = c('id', 'tissue', 'type', 'disease')
  dataInfo$type = as.character(dataInfo$type)
  
  if(length(which(dataInfo$type %in% tumorSampleTypes))!=0)
  {
    dataInfo[dataInfo$type %in% tumorSampleTypes,]$type = 'Tumor'
  }
  
  if(length(which(dataInfo$type %in% normalSampleTypes))!=0)
  {
    dataInfo[dataInfo$type %in% normalSampleTypes,]$type = 'Normal'
    dataInfo[dataInfo$type=='Normal',]$disease = NA
  }
  
  dataNameMapping = dataInfoAndMapping[,c('project', 'cases', 'id', 'submitter_id')]
  names(dataNameMapping) = c('project_id', 'case_id', 'uuid', 'submitter_id')
  
  GDCData[['NameMapping']][[dataType]] <<- dataNameMapping 
  
  if(length(GDCData[['SampleFilesInfo']])==0)
    GDCData[['SampleFilesInfo']] <<- dataInfo
}



ImportData = function(selectedProject, dataType)
{
  if(dataType %in% c('Counts','FPKM'))
  {
    
    dataFormat = if(dataType == 'Counts') 'counts' else 'FPKM.txt'
    pathOfSamplesAll=list.files('GDCdata', recursive = T, pattern = paste0("*.",dataFormat,".gz$"), full.names = T)
    
    ### Import data in serial 
    
    # expressionMatrix = data.frame(ensemblIDs)
    # names(expressionMatrix)='ENSEMBL_ID'
    
    # for(j in pathOfSamplesAll)
    # {
    #   ### Extract Submitter ID
    #   temp = strsplit2(j, '/')
    #   fileName = temp[length(temp)]
    #   #submitterID=paste0('S__',gsub('-','_',strsplit2(fileName,'[.]')[1]))
    #   submitterID=strsplit2(fileName,'[.]')[1]
    #   print(submitterID)
    #   ### Read the expression file 
    #   sampleExpression=read.delim(j, header = F, col.names = c('ENSEMBL_ID',submitterID), comment.char = '_', check.names = F)
    #   sampleExpression$ENSEMBL_ID=substr(sampleExpression$ENSEMBL_ID,1,15)
    #   sampleExpression=sampleExpression[order(sampleExpression$ENSEMBL_ID),]
    #   expressionMatrix=cbind(expressionMatrix,sampleExpression[submitterID])
    #   
    # }
    
    
    ### Import data in parallel using 'mclapply'
    
    pathsSplitted = split(pathOfSamplesAll, cut(seq_along(pathOfSamplesAll), detectCores()-1))
    names(pathsSplitted) = seq(detectCores()-1)
    expressionMatrix =
      mclapply(pathsSplitted,
               function(paths)
               {
                 tempExpressionMatrix = data.frame(ensemblIDs)
                 for(j in paths)
                 {
                   ### Extract Submitter ID
                   temp = strsplit2(j, '/')
                   fileName = temp[length(temp)]
                   submitterID=strsplit2(fileName,'[.]')[1]
                   print(submitterID)
                   
                   ### Read the expression file
                   sampleExpression=read.delim(j, header = F, col.names = c('ENSEMBL_ID',submitterID), comment.char = '_', check.names = F)
                   sampleExpression$ENSEMBL_ID=substr(sampleExpression$ENSEMBL_ID,1,15)
                   sampleExpression=sampleExpression[order(sampleExpression$ENSEMBL_ID),]
                   tempExpressionMatrix=cbind(tempExpressionMatrix,sampleExpression[submitterID])
                 }
                 tempExpressionMatrix$ENSEMBL_ID=NULL
                 tempExpressionMatrix
               },
               mc.cores = detectCores()-1)
    
    
    expressionMatrix = Reduce(cbind, expressionMatrix)
    rownames(expressionMatrix)=ensemblIDs$ENSEMBL_ID
    
  } else if(dataType=='miRNA')
  {
    pathOfSamplesAll=list.files('GDCdata', recursive = T, pattern = paste0("*.",'.isoforms.quantification.txt'), full.names = T)
    
    ### Import data in serial
    
    # expressionMatrix = mirbaseNameAndIDs
    
    # for(j in pathOfSamplesAll)
    # {
    #   ### Extract Submitter ID
    #   temp = strsplit2(j, '/')
    #   fileName = temp[length(temp)]
    #   #submitterID=paste0('S__',gsub('-','_',strsplit2(fileName,'[.]')[1]))
    #   submitterID=strsplit2(fileName,'[.]')[1]
    #   print(submitterID)
    #   ### Read the expression file 
    #   sampleExpression=read.delim(j, header = T, check.names = F)
    #   # colnames(isoformSampleExpression)
    #   sampleExpression = sampleExpression[,c('miRNA_region','read_count')]
    #   names(sampleExpression) = c('id', 'count')
    #   sampleExpression = sampleExpression[grepl('mature', sampleExpression$id),]
    #   sampleExpression$id = gsub(pattern = 'mature,', '', sampleExpression$id)
    #   sampleExpression = aggregate(count ~ id, sampleExpression, sum)
    #   sampleExpression = merge(mirbaseNameAndIDs, sampleExpression, by = 'id', all.x=T)
    #   sampleExpression$count[is.na(sampleExpression$count)]=0
    #   
    #   names(sampleExpression) = c('id', 'name', submitterID)
    #   expressionMatrix=cbind(expressionMatrix, sampleExpression[submitterID]) 
    # }
    # 
    # rownames(expressionMatrix)=expressionMatrix$id
    # expressionMatrix[,c('id','name')]=NULL
    
    ### Import data in parallel using 'mclapply'
    
    pathsSplitted = split(pathOfSamplesAll, cut(seq_along(pathOfSamplesAll), detectCores()-1))
    names(pathsSplitted) = seq(detectCores()-1)
    
    expressionMatrix =
      mclapply(pathsSplitted,
               function(paths)
               {
                 tempExpressionMatrix = mirbaseNameAndIDs
                 for(j in paths)
                 {
                   ### Extract Submitter ID
                   temp = strsplit2(j, '/')
                   fileName = temp[length(temp)]
                   submitterID=strsplit2(fileName,'[.]')[1]
                   print(submitterID)
                   
                   ### Read the expression file 
                   sampleExpression=read.delim(j, header = T, check.names = F)
                   sampleExpression = sampleExpression[,c('miRNA_region','read_count')]
                   names(sampleExpression) = c('id', 'count')
                   sampleExpression = sampleExpression[grepl('mature', sampleExpression$id),]
                   sampleExpression$id = gsub(pattern = 'mature,', '', sampleExpression$id)
                   sampleExpression = aggregate(count ~ id, sampleExpression, sum)
                   sampleExpression = merge(mirbaseNameAndIDs, sampleExpression, by = 'id', all.x=T)
                   sampleExpression$count[is.na(sampleExpression$count)]=0
                   
                   names(sampleExpression) = c('id', 'name', submitterID)
                   tempExpressionMatrix=cbind(tempExpressionMatrix, sampleExpression[submitterID])
                 }
                 tempExpressionMatrix[,c('id','name')]=NULL
                 tempExpressionMatrix
               },
               mc.cores = detectCores()-1)
    
    expressionMatrix = Reduce(cbind, expressionMatrix)
    rownames(expressionMatrix)=mirbaseNameAndIDs$id
  }
  
  submitterIDs = data.frame(submitter_id = colnames(expressionMatrix))
  fullCaseIDs = merge(submitterIDs, GDCData[['NameMapping']][[dataType]], by='submitter_id', sort = F)$case_id
  colnames(expressionMatrix) = fullCaseIDs
  
  GDCData[[dataType]] <<- expressionMatrix
  
  projectFullName = GetProjectInfo(selectedProject, getFullNameOnly = T)
  pathToData = paste(getwd(), 'GDCDataRepositpory', projectFullName, dataType, sep = '/')
  dir.create(pathToData, recursive = T)
  system(paste0('mv GDCdata/* ', paste0("'",pathToData,"'")))
  system(paste0('mv MANIFEST.txt ', paste0("'",pathToData,"'")))
}


QueryAndDownloadClinicalData = function(selectedProject, selectedCases)
{
  print('Clinical')
  
  selectedParticipants = if(grepl('TCGA', selectedProject)) unique(substr(selectedCases, 1, 12)) else unique(substr(selectedCases, 1, 16))
  
  if(grepl('TCGA', selectedProject))
  {
    queryOfClinical = GDCquery(project = selectedProject,
                               data.category = 'Clinical',
                               barcode = selectedParticipants)
  } else
  {
    queryOfClinical = GDCquery(project = selectedProject,
                               data.category = 'Clinical')
  }
  
  GDCdownload(queryOfClinical)
  resultOfClinical = getResults(queryOfClinical)
  
  importantCols = c('project', 'cases', 'file_id', 'file_name')
  clinicalNameMapping = resultOfClinical[,c('project', 'cases', 'file_id', 'file_name')]
  names(clinicalNameMapping) = c('project_id', 'case_id', 'uuid', 'file_name')
  
  GDCData[['NameMapping']][['Clinical']] <<- clinicalNameMapping
  
  queryOfClinical
  
}



ImportClinicalData=function(selectedProject, queryOfClinical)
{
  if(grepl('TCGA', selectedProject))
  {
    clinicalInfo = GDCprepare_clinic(queryOfClinical, clinical.info = "patient")
    caseIDLength  = 12
  } else
  {
    clinicalInfoFile = list.files(path = 'GDCdata', pattern = '[.]xlsx', recursive = T, full.names = T)[1]
    if(length(clinicalInfoFile)>1)
      clinicalInfoFile = list.files(path = 'GDCdata', pattern = '.*Discovery.*[.]xlsx', recursive = T, full.names = T)[1]
    clinicalInfo = read.xlsx(file = clinicalInfoFile, sheetIndex = 1, as.data.frame = T, header = T)
    names(clinicalInfo) = c('bcr_patient_barcode', names(clinicalInfo)[-1])
    caseIDLength  = 16
  }
  
  sampleFilesInfo = GDCData[['SampleFilesInfo']]
  sampleFilesInfo$cases = substr(sampleFilesInfo$id, 1, caseIDLength)
  
  clinicalInfo = merge(sampleFilesInfo, clinicalInfo, by.x = 'cases', by.y = 'bcr_patient_barcode')
  
  GDCData[['Clinical']] <<- clinicalInfo
  
  projectFullName = GetProjectInfo(selectedProject, getFullNameOnly = T)
  dataType = 'Clinical'
  
  pathToData = paste(getwd(), 'GDCDataRepositpory', projectFullName, dataType, sep = '/')
  dir.create(pathToData, recursive = T)
  system(paste0('mv GDCdata/* ', paste0("'",pathToData,"'")))
  system(paste0('mv MANIFEST.txt ', paste0("'",pathToData,"'")))
  
}



SaveGDCData = function(selectedProject)
{
  print('Saving Data')
  
  projectFullName = GetProjectInfo(selectedProject, getFullNameOnly = T)
  pathToData = paste(getwd(), 'GDCDataRepositpory', projectFullName, sep = '/')
  
  DBName = paste0('GDCData__', sub('-','_',selectedProject))
  assign(DBName, GDCData)
  
  RDataFileName = paste0(pathToData, '/', DBName, '.rds')
  
  saveRDS(object = get(DBName), file = RDataFileName)
  
  print('Done')
}


ArrangeData = function(dataPath = 'GDCDataRepositpory')
{
  GDCDataPaths = list.files(path = dataPath, pattern = 'GDCData__.*[.]rds', recursive = T, full.names = T)
  GDCDataRepository = list()
  # fieldsOfClinicalInfo = list()
  
  for(GDCDataPath in GDCDataPaths)
  {
    selectedProject = sub('_','-',strsplit2(strsplit2(GDCDataPath,'__')[,2],'[.]')[,1])
    print(selectedProject)
    
    GDCData = readRDS(GDCDataPath)
    
    dataTypes = c('Counts', 'FPKM', 'miRNA')
    
    for(dataType in dataTypes)
      colnames(GDCData[[dataType]])=substr(colnames(GDCData[[dataType]]), 1, 16)
    
    
    ### Counts-miRNA index matching
    if(length(colnames(GDCData[['Counts']])) <= length(colnames(GDCData[['miRNA']])))
    {
      GDCData[['Counts']] = GDCData[['Counts']][, order(colnames(GDCData[['Counts']]))]
      colnames(GDCData[['Counts']])=substr(colnames(GDCData[['Counts']]), 1, 16)
      
      colIndex = match(colnames(GDCData[['Counts']]), colnames(GDCData[['miRNA']]))
      colIndexRev = match(colnames(GDCData[['Counts']]), rev(colnames(GDCData[['miRNA']])))
      colIndex[duplicated(colIndex)]=length(colnames(GDCData[['miRNA']])) - colIndexRev[duplicated(colIndexRev)] +1
      
      GDCData[['miRNA']] = GDCData[['miRNA']][, colIndex]
      colnames(GDCData[['miRNA']])=substr(colnames(GDCData[['miRNA']]), 1, 16)
      
      clinicalsToRemoveCauseOfSampleDeletion = setdiff(seq(ncol(GDCData[['miRNA']])), colIndex)
      clinicalsToRemoveCauseOfSampleDeletion = colnames((GDCData[['miRNA']]))[clinicalsToRemoveCauseOfSampleDeletion]
      
    } else # if(colnames(GDCData[['Counts']]) > colnames(GDCData[['miRNA']]))
    {
      GDCData[['miRNA']] = GDCData[['miRNA']][, order(colnames(GDCData[['miRNA']]))]
      colnames(GDCData[['miRNA']])=substr(colnames(GDCData[['miRNA']]), 1, 16)
      
      colIndex = match(colnames(GDCData[['miRNA']]), colnames(GDCData[['Counts']]))
      colIndexRev = match(colnames(GDCData[['miRNA']]), rev(colnames(GDCData[['Counts']])))
      colIndex[duplicated(colIndex)]=length(colnames(GDCData[['Counts']])) - colIndexRev[duplicated(colIndexRev)] +1
      
      GDCData[['Counts']] = GDCData[['Counts']][, colIndex]
      colnames(GDCData[['Counts']])=substr(colnames(GDCData[['Counts']]), 1, 16)
      
      clinicalsToRemoveCauseOfSampleDeletion = setdiff(seq(ncol(GDCData[['Counts']])), colIndex)
      clinicalsToRemoveCauseOfSampleDeletion = colnames((GDCData[['Counts']]))[clinicalsToRemoveCauseOfSampleDeletion]
    }
    
    ### FPKM index matching
    colIndex = match(colnames(GDCData[['Counts']]), colnames(GDCData[['FPKM']]))
    GDCData[['FPKM']] = GDCData[['FPKM']][, colIndex]
    colnames(GDCData[['FPKM']])=substr(colnames(GDCData[['FPKM']]), 1, 16)
    
    ### Clinical index matching
    ### Fix samples which do NOT have clinical info
    GDCData[['Clinical']]$id = substr(GDCData[['Clinical']]$id, 1, 16)
    
    samplesID = colnames(GDCData[['Counts']])
    clinicalsID = GDCData[['Clinical']]$id
    samplesWithoutClinical = samplesID[!(samplesID%in%clinicalsID)]
    
    for(sampleID in samplesWithoutClinical)
    {
      cases = substr(sampleID, 1, 12)
      id = sampleID
      # tissue = GDCProjects[GDCProjects$project_id==selectedProject,'primary_site.0']
      tissue = GDCProjects[GDCProjects$project_id==selectedProject,'primary_site'][[1]]
      type = ifelse(substr(sampleID, 14, 15)=='11', 'Normal', 'Tumor')
      # disease = ifelse(type=='Normal', NA, GDCProjects[GDCProjects$project_id==selectedProject,'disease_type.0'])
      # disease = ifelse(type=='Normal', NA, GDCProjects[GDCProjects$project_id==selectedProject,'disease_type'][[1]])
      if(type=='Normal') disease = NA else disease = GDCProjects[GDCProjects$project_id==selectedProject,'disease_type'][[1]]
      other = rep(NA, ncol(GDCData[['Clinical']])-5)
      
      sampleClinicalInfo = c(cases, id, tissue, type, disease, other)
      GDCData[['Clinical']] = rbind(GDCData[['Clinical']], sampleClinicalInfo)
    }
    
    ### Fix clinical info which do NOT have sample
    GDCData[['Clinical']] = GDCData[['Clinical']][order(GDCData[['Clinical']]$id), ]
    clinicalsToRemoveCauseOfSampleDeletion = which(GDCData[['Clinical']]$id %in% clinicalsToRemoveCauseOfSampleDeletion)
    if(length(clinicalsToRemoveCauseOfSampleDeletion)>0 && nrow(GDCData[['Clinical']])>ncol(GDCData[['Counts']]))
      GDCData[['Clinical']]=GDCData[['Clinical']][-clinicalsToRemoveCauseOfSampleDeletion, ]
    
    ### matching index
    colIndex = match(GDCData[['Clinical']]$id, colnames(GDCData[['Counts']]))
    GDCData[['Clinical']] = GDCData[['Clinical']][colIndex, ]
    
    
    #### Fix Clinical data issues ####
    GDCData[['Clinical']] = FixClinicalDataLabelsIssues(selectedProject, GDCData[['Clinical']])
  
    
    # fieldsOfClinicalInfo[[selectedProject]]=colnames(GDCData[['Clinical']])
    GDCDataRepository[[selectedProject]] = GDCData
  }
  
  system('rm -d GDCdata')
  GDCDataRepository
}



FixClinicalDataLabelsIssues = function(selectedProject, clinicals)
{
  ### Correct clinical 'type' lables, identify metastatic & new tumors as 'Tumor'
  if(length(which(!(clinicals$type%in%c('Tumor','Normal'))))>0)
    clinicals[!(clinicals$type%in%c('Tumor','Normal')),]$type = 'Tumor'
  
  ### Correct clinical 'disease' labels, previously the 'disease' of 'Normal' cases was set to 'NA'
  diseaseAndOrNA = unique(clinicals$disease)
  if(selectedProject == 'TARGET-RT')
    diseaseAndOrNA = 'Rhabdoid Tumor'
  disease = diseaseAndOrNA[!is.na(diseaseAndOrNA)]

  clinicals$disease = disease
  
  clinicals
}



MergeAllGDCData = function(GDCDataRepository)
{
  GDCData = GetAnEmptyGDCDataStructure()
  
  clinicalInfo = list()
  
  for(i in names(GDCDataRepository))
  {
    print(i)
    GDCProjectID = i
    GDCDataOfAProject = GDCDataRepository[[i]]
    
    if(length(GDCData[['Counts']])==0)
    {
      isItFirstProject = T
      GDCData = GDCDataOfAProject
    } else
    {
      GDCData[c('Counts','FPKM','miRNA')] = 
        mapply(cbind, GDCData[c('Counts','FPKM','miRNA')], GDCDataOfAProject[c('Counts','FPKM','miRNA')], SIMPLIFY = F)
      
      GDCData[['NameMapping']] =
        mapply(rbind, GDCData[['NameMapping']], GDCDataOfAProject[['NameMapping']], SIMPLIFY = F)
      
      GDCData['SampleFilesInfo'] =
        mapply(rbind, GDCData['SampleFilesInfo'], GDCDataOfAProject['SampleFilesInfo'], SIMPLIFY = F)
    }
    
    
    clinicalInfo[[GDCProjectID]] = GDCDataOfAProject[['Clinical']]
    
    projectType = ifelse(grepl('TCGA', GDCProjectID), 'TCGA', 'TARGET')
    includedDesiredCols = colnames(GDCDataOfAProject[['Clinical']]) %in% GetDesiredFieldsOfClinicalInfo(projectType)
    notIncludedDesiredCols = setdiff(GetDesiredFieldsOfClinicalInfo(projectType), colnames(GDCDataOfAProject[['Clinical']])[includedDesiredCols])
    GDCDataOfAProject[['Clinical']][, notIncludedDesiredCols] = NA
    # desiredCols = colnames(GDCDataOfAProject[['Clinical']]) %in% GetDesiredFieldsOfClinicalInfo(projectType)
    
    clinicalInfoOfAProject = GDCDataOfAProject[['Clinical']][, GetDesiredFieldsOfClinicalInfo(projectType)]
    colnames(clinicalInfoOfAProject) = GetDesiredFieldsOfClinicalInfo('ProperNamesOnly')
    
    clinicalInfoOfAProject[,'ageAtDiagnosisInDays'] = as.integer(clinicalInfoOfAProject[,'ageAtDiagnosisInDays'])
    if(projectType=='TCGA')
      clinicalInfoOfAProject[,'ageAtDiagnosisInDays'] = -clinicalInfoOfAProject[,'ageAtDiagnosisInDays']
    
    if(isItFirstProject)
    {
      GDCData[['Clinical']] = clinicalInfoOfAProject
      isItFirstProject = F
    } else
    {
      GDCData[['Clinical']] = rbind(GDCData[['Clinical']], clinicalInfoOfAProject)
    }
    
    # fieldsOfClinicalInfo[[GDCProjectID]]=colnames(GDCDataOfAProject[['Clinical']])
  }
  
  GDCData[['AllClinical']] = clinicalInfo
  
  GDCData
}


GetDesiredFieldsOfClinicalInfo = function(projectType, method='Manual')
{
  if(method=='Manual')
  {
    if(projectType=='TCGA')
    {
      desiredFieldsOfClinicalInfo = c('cases', 'id',	'tissue',	'type',	'disease',
                                      'gender', 'race_list',	'ethnicity', 'days_to_birth',
                                      'vital_status', 'primary_pathology_year_of_initial_pathologic_diagnosis',
                                      'stage_event_pathologic_stage')
    }
    else if(projectType=='TARGET')
    {
      desiredFieldsOfClinicalInfo = c('cases', 'id',	'tissue',	'type',	'disease',
                                      'Gender', 'Race', 'Ethnicity', 'Age.at.Diagnosis.in.Days',
                                      'Vital.Status',' Year.of.Diagnosis', 'Stage')
      
    } else
    {
      desiredFieldsOfClinicalInfo = c('cases', 'id',	'tissue',	'type',	'disease',
                                      'gender', 'race', 'ethnicity', 'ageAtDiagnosisInDays',
                                      'vitalStatus', 'yearOfDiagnosis', 'Stage')
    }
  }
  
  desiredFieldsOfClinicalInfo
}



GetGDCDataSummary=function()
{
  GDCDataSummary = data.frame(tissue = character(0),
                              disease = character(0),
                              total = integer(0),
                              tumor = integer(0),
                              normal = integer(0))
  clinicals = GDCData[['Clinical']]
  tissues = unique(clinicals$tissue)
  tissues = tissues[!is.na(tissues)]
  # tissue = tissues[1]
  for(tissue in tissues)
  {
    clinicalsOfATissue = clinicals[clinicals$tissue==tissue, c('disease', 'type')]
    diseasesOfATissue = unique(clinicalsOfATissue$disease)
    diseasesOfATissue = diseasesOfATissue[!is.na(diseasesOfATissue)]
    # disease = diseasesOfATissue[1]
    for(disease in diseasesOfATissue)
    {
      clinicalsOfADisease = clinicalsOfATissue[clinicalsOfATissue$disease==disease, ]
      total = nrow(clinicalsOfADisease)
      tumor = nrow(clinicalsOfADisease[clinicalsOfADisease$type=='Tumor',])
      normal = nrow(clinicalsOfADisease[clinicalsOfADisease$type=='Normal',])
      
      GDCDataSummary = rbind(GDCDataSummary, c(tissue, disease, total, tumor, normal))
    }
  }
  
  names(GDCDataSummary) = c('tissue', 'disease', 'total', 'tumor', 'normal')
  class(GDCDataSummary[, 'total'])  = class(GDCDataSummary[, 'tumor']) = class(GDCDataSummary[, 'normal']) = 'integer'

  sumsOfAll = colSums(GDCDataSummary[, c('total', 'tumor', 'normal')])
  GDCDataSummary = rbind(GDCDataSummary, c('All', 'All', sumsOfAll))
  
  GDCDataSummary
}



AddGenesInfoToGDCData = function()
{
  ensemblRNAsInfo = read.delim('Files/Ensembl_RNAs_Info.txt', header = T, sep = '\t')
  names(ensemblRNAsInfo) = c('id', 'symbol', 'type', 'description', 'entrez')
  
  ### Some ENSEMBL IDs have more than one Entrez IDs
  ensemblRNAsInfoWithoutEntrez = ensemblRNAsInfo[!duplicated(ensemblRNAsInfo$id),-5]
  
  GDCGenesID = data.frame(id=rownames(GDCData[['Counts']]))
  
  GDCGenesInfo = merge(GDCGenesID, ensemblRNAsInfoWithoutEntrez, by='id', all.x = T)
  GDCGenesEntrezIDs = merge(GDCGenesID, ensemblRNAsInfo[,c('id','entrez')], by='id', all.x = T)
  GDCData[['GenesInfo']] <<- list(All=GDCGenesInfo, EntrezIDs = GDCGenesEntrezIDs)
}



ExtractSpecificRNATypeData = function(selectedRNAType)
{
  GDCGenesInfo = GDCData[['GenesInfo']][['All']]
  selectedRNATypeIDs = GDCGenesInfo[GDCGenesInfo$type==selectedRNAType,'id']
  selectedRNATypeIDs = selectedRNATypeIDs[!is.na(selectedRNATypeIDs)]
  
  selectedRNATypeData = list()
  selectedRNATypeData[['Counts']] = GDCData[['Counts']][rownames(GDCData[['Counts']]) %in% selectedRNATypeIDs,]
  selectedRNATypeData[['FPKM']] = GDCData[['FPKM']][rownames(GDCData[['FPKM']]) %in% selectedRNATypeIDs,]
  
  genesInfoAll = GDCData[['GenesInfo']][['All']][GDCData[['GenesInfo']][['All']]$id %in% selectedRNATypeIDs, ]
  genesEntrezIDs = GDCData[['GenesInfo']][['EntrezIDs']][GDCData[['GenesInfo']][['EntrezIDs']]$id %in% selectedRNATypeIDs, ]
  selectedRNATypeData[['GenesInfo']] = list(All = genesInfoAll, EntrezIDs = genesEntrezIDs)
  
  selectedRNATypeData
}



HarmonizeClinicalData = function(clinicals)
{
  tcgaClinicals = clinicals[grepl('TCGA', clinicals$cases),]
  targetClinicals = clinicals[grepl('TARGET', clinicals$cases),]
  
  xenaPhenotypes = read.delim('Files/TCGA_TARGET_phenotype.txt', header = T, sep = '\t', check.names = F)
  tcgaSurvivalData = xenaPhenotypes[,c('_PATIENT','_OS','_RFS','_RFS_IND')]
  names(tcgaSurvivalData) = substr(names(tcgaSurvivalData), 2, nchar(names(tcgaSurvivalData)))
  tcgaSurvivalData = tcgaSurvivalData[!duplicated(tcgaSurvivalData$PATIENT),]
  
  newTCGAClinicals = merge(tcgaClinicals, tcgaSurvivalData, by.x = 'cases', by.y = 'PATIENT', all.x=T, sort = F)
  
  
  targetRawClinicals = GDCData[['AllClinical']][grepl('TARGET', names(GDCData$AllClinical))]
  targetSurvivalData = lapply(targetRawClinicals, function(clinical) clinical[,c('cases', 'id',
                                                                                 'First.Event',
                                                                                 'Event.Free.Survival.Time.in.Days',
                                                                                 'Overall.Survival.Time.in.Days')])
  targetSurvivalData = Reduce(rbind, targetSurvivalData)
  names(targetSurvivalData) = c('cases','id', 'RFS_IND', 'RFS', 'OS')
  targetSurvivalData$RFS_IND[targetSurvivalData$RFS_IND=='Relapse'] = 1
  targetSurvivalData$RFS_IND[targetSurvivalData$RFS_IND=='Censored'] = NA
  targetSurvivalData$RFS_IND[!is.na(targetSurvivalData$RFS_IND) & targetSurvivalData$RFS_IND!=1]=0
  targetSurvivalData = targetSurvivalData[!duplicated(targetSurvivalData$cases),]
  targetSurvivalData$id = NULL
  
  newTARGETClinicals = merge(targetClinicals, targetSurvivalData, by.x = 'cases', all.x=T, sort = F)
  
  newTCGAClinicals$index = grep('TCGA', clinicals$cases)
  newTARGETClinicals$index = grep('TARGET', clinicals$cases)
  newClinicals = rbind(newTCGAClinicals, newTARGETClinicals)
  newClinicals = newClinicals[order(newClinicals$index),]
  newClinicals$index = NULL
  
  colnames(newClinicals)[colnames(newClinicals)=='RFS_IND']='isRelapsed'
  newClinicals$isRelapsed[newClinicals$isRelapsed==1]=T
  newClinicals$isRelapsed[newClinicals$isRelapsed==0]=F
  newClinicals$isRelapsed = as.factor(newClinicals$isRelapsed)
  
  colnames(newClinicals)[colnames(newClinicals)=='RFS']='relapseTimeInDays'
  newClinicals$relapseTimeInDays = as.integer(newClinicals$relapseTimeInDays)
  colnames(newClinicals)[colnames(newClinicals)=='OS']='overallSurvivalTimeInDays'
  newClinicals$overallSurvivalTimeInDays = as.integer(newClinicals$overallSurvivalTimeInDays)
  
  newClinicals$gender[newClinicals$gender=='MALE'] = 'Male'
  newClinicals$gender[newClinicals$gender=='FEMALE'] = 'Female'
  newClinicals$gender = as.factor(newClinicals$gender)
  
  newClinicals$vitalStatus[newClinicals$vitalStatus=='ALIVE'] = 'Alive'
  newClinicals$vitalStatus[newClinicals$vitalStatus=='DEAD'] = 'Dead'
  newClinicals$vitalStatus[newClinicals$vitalStatus==''] = NA
  newClinicals$vitalStatus = as.factor(newClinicals$vitalStatus)
  
  newClinicals$ethnicity = NULL
  newClinicals$race = NULL
  
  newClinicals$ageAtDiagnosisInYears = round(newClinicals$ageAtDiagnosisInDays/365,1)
  newClinicals$relapseTimeInYears = round(newClinicals$relapseTimeInDays/365,1)
  newClinicals$overallSurvivalTimeInYears = round(newClinicals$overallSurvivalTimeInDays/365,1)
  newClinicals$yearsBeingAliveAfterRelapse = newClinicals$overallSurvivalTimeInYears - newClinicals$relapseTimeInYears
  
  colnames(newClinicals)[colnames(newClinicals)=='Stage']='stage'
  newClinicals$stage = gsub('Stage| |/.*|A|B|C|S|NOS', '', newClinicals$stage)
  newClinicals$stage[newClinicals$stage=='' | newClinicals$stage=='X']=NA
  newClinicals$stage = factor(newClinicals$stage, labels = 0:4)
  
  newClinicals
}



ExtractSpecificClinicalFields = function(clinicals, selectedFields, tumorsOnly, removeNAs)
{
  
  if(tumorsOnly)  
    clinicals = clinicals[clinicals$type=='Tumor',]
  
  extractedClinicals = clinicals[,selectedFields]
  if(removeNAs)
    for(field in selectedFields)
      extractedClinicals = extractedClinicals[!is.na(extractedClinicals[[field]]),]
  
  extractedClinicals
}

