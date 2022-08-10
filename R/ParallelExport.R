library(data.table)
library(RSQLite)
DeletePrevious <- T
ExportPeaks <- T
ExportRatios <- T

if(!file.exists("Export_Vali.rda")){
  stop("No Export_Vali.rda found. Aborting Export")
}

try({
  ls <- load("Export_Vali.rda")
  (source(paste(SystemPath,"/R/Vali_Functions.R",sep = "")))
  
  library(parallel)
  
  # stopCluster(cl)
  # wd <- getwd()
  
  if(DeletePrevious){
    try({
      
      wd <- getwd()
      ls <- c(ls)
      ls <- unique(ls)
      li <- list.files("./ParallelExport/",full.names = T)
      if(length(li)>0){
        sapply(li,unlink)
      }
      li <- list.files("./export/",full.names = T)
      if(length(li)>0){
        sapply(li,unlink)
      }
    })
    
  }
  dir.create("ParallelExport",showWarnings = F)
  unlink("export.txt",force=T)
  clu <- makeCluster(threads,outfile="export.txt")
  parallel::clusterExport(clu,ls)
  
  if(ExportPeaks){
    try({parSapply(clu,1:length(anaexport),function(i){
      # try({sapply(7:length(anaexport),function(i){
      icurrent <<- i
      # options(warn=1)
      # save.image(paste("./ParallelExport/ExportParallel_",Sys.getpid(),".rda",sep = ""))
      sessionId <- paste("./ParallelExport/CurrentExportProcess",Sys.getpid(),i,sep="#")
      sessionId_grep <- paste("CurrentExportProcess",Sys.getpid(),sep="#")
      sessionId_DB <- paste("./ParallelExport/CurrentExportDB",Sys.getpid(),sep="#")
      write("","./ParallelExport/TableCompiled")
      li <- list.files("./ParallelExport",pattern=sessionId_grep,full.names = T)
      if(length(li)==0){
        write(i,gsub("Current","Initial",sessionId))
        
      }else{
        unlink(li,force = T)
        
      }
      write(i,sessionId)
      # names 
      library(RSQLite)
      library(data.table)
      library(pracma)
      
      (source(paste(SystemPath,"/R/Vali_Functions.R",sep = "")))
      
      PeaksTableName <- paste("./ParallelExport/",Sys.getpid(),"Peaks.txt",sep="")
      
      fipath <- paste("./ParallelExport/TransitionList",Sys.getpid(),".txt",sep = "_")
      pw <- 1
      
      fifu <- anaexport[i]
      cat("\rworking on ",fifu)
      session <- NULL
      # setwd("settings")
      # finame <- gsub(".txt$","",fifu)
      # # print("Loading Stuff")
      # if(file.exists(finame)){
      #   load(finame)
      # }else{
      # }
      
      # print("Export: Creating DP Table")
      # Creating DP Table
      # inputList$FDR <- FDRpep
      if(length(session) > 0){
        progress$set(message = 'Export, working on:',
                     detail = finame,value = i)
      }
      
      # setwd("../")
      # if(i==1){
      #   unlink("./export/Peaks.txt")
      #   
      # }
      
      PeaksName <- dbtaNameExport(dbp,ana = fifu)
      PeaksNameRatin <- paste("RATING_",PeaksName,sep = "")
      
      
      if(any(PeaksNameRatin== dblistT(dbp))){
        # print("check Rating")
        PrecursorRating <- dbread(PeaksNameRatin,db = dbp)$rating
        # dbread(PeaksName,dbpath)
        # print("Found Rating")
        if(length(PrecursorRating) == 0){
          PrecursorRating <- "not set"
        }
      }else{
        # print("not set Rating")
        PrecursorRating <- "not set"
      }
      
      
      # check if Peptide (including or excluding PW) was already ana in the database  
      TYPIES <- which(gsub(".$","",PeaksName)==gsub(".$","",dblistT(dbp)))
      Lo <- dblistT(dbp)[TYPIES]
      
      LoopCount <- 1
      if(length(TYPIES) >1){
        LoopCount <- length(TYPIES)
      }
      # Loop is to account several Peak tables caused by different checked PeakWidths in the database, need to manually be selected afterwards
      # print("Export: PeakCount")
      if(file.exists("tempa")){
        try(rm("tempa"),silent = T)
        
      }
      # dbp <- dbpath()
      
      for(PeakCount in 1:LoopCount){
        tempa <<- dbread(x = fifu,dbp)
        if(length(tempa$RF_Scores)==0){
          tempa$RF_Scores <- 0
        }
        if(length(tempa$DL_Scores)==0){
          tempa$DL_Scores <- 0
        }
        
        if(dim(tempa)[1]>=minthresh){
          
          if(length(TYPIES) > 0&F){
            
            cl <- class(try(DP <- dbread(Lo[PeakCount],dbp)))
            PeakType <- "Manual"
            # pw <- as.numeric(strsplitslot(Lo[PeakCount],8,"_"))
            # tempa <- dbread(x = fifu,dbpath())#
            
          }else{
            
            tempa$RT_Used <- tempa$RT_min
            tempa$FDR <- tempa$FDR_RF_all
            FDRcut <- FDRpep
            
            CANDIDATE_RT <- RTcandidatesFun(list(tempa),FDRcut,Requantify_Priority=Requantify_Priority)
            LI <- list(tempa)
            names(LI) <- fifu
            
            tempList <- list(tempa)
            names(tempList) <- fifu
            tempList <- tempList
            DPlist_XIC <- DetectPeakWrapper(ana = tempList,CANDIDATE_RT = CANDIDATE_RT,dbp = dbp,
                                            RetentionTimeWindow = RTwin,QType = "Intensities",
                                            Reanalysis = F,
                                            RT_BASED_onbestScore = T,
                                            MinPeakWidth=MinPeakWidth_vec,
                                            MaxPeakWidth=MaxPeakWidth_vec,
                                            supersmooth_I_set = F,#input$supersmooth_I_set,
                                            supersmooth_bw_set = supersmooth_bw_set_vec,
                                            ApplyMaximumWidth =  ApplyMaximumWidth_vec,
                                            Requantify_Priority= Requantify_Priority,
                                            parallel_db = paste(sessionId_DB,".sqlite",sep = "")
                                            
            )
            
            # DP <- rbind(DPlist_I[[1]],DPlist_XIC[[1]])
            DP <- data.table(DPlist_XIC[[1]],stringsAsFactors = F)
            
            tempa$PEP[is.na(tempa$PEP)] <- 1
            # cl <- class(try(DP <- dbread(Lo[PeakCount],dbp)))
            cl <- "everythingok"
            # cl <- class(try(DP <- PeakFun(tempa,RetentionTimeWindow = inputList$RetentionTimeWindow,alpha = input$FDRpep,PeakWidth = inputList$PeakWidth,SimilarRT = F,session = session)))
            PeakType <- "Automated"
          }
          # try({
          # tlnx <<- tlnx
          try(transname <- gsub("^mz","TRANS",fifu))
          # print(transname)
          dbl <- dblistT(dbp)
          # print(dbl)
          if(any(dbl==transname)){
            # print("FoundTrans")
            dbp <- dbpath()
            TransSaved <- NULL
            try(selec <- unlist(dbread(gsub("^PEAKS","TRANS",transname),db =dbp)))
            transis <- colnames(tempa)[1:(which(colnames(tempa) == "charge")-1)]
            notused <- which(is.na(match(transis,selec)))
            if(length(notused)> 0){
              showNotification("Excluding Fragments", duration = 1, closeButton =  TRUE, type = "warning")
              tempa <- tempa[,-notused]
            }
            
          }
          # })
          # Writing Out FragmentTable
          ({
            tempa_dt <- data.table(tempa)
            DPdt  <- data.table(DP)
            # tempa <- cbind(tempa,DPdt[match(tempa$rawfile,DP$rawfile),.(Q1,Q2,Q1align,Q2align)])
            tempa_dt[,Peak :=RT_min%between%DPdt[unlist(.BY)==rawfile,.(Q1,Q2)][1,],rawfile]
            tempa_dt <- tempa_dt[Peak == TRUE,]
            tempa_dt[,precursor_ID :=fifu]
            tempa_dt$Precursor_Rating <- PrecursorRating
            rf <- tempa_dt$rawfile
            tempa_dt$rawfile <- NULL
            tempa_dt$rawfile <- rf
            id <- 1:(grep("^charge$",names(tempa_dt))-1)
            
            # check order of COlumns:
            # tempa_dt <- tempa_dt
            # id <- id
            NAMES <- c(colnames(tempa_dt)[id],"charge","mz","index","RT_min","MatchCount","R","R_p","IntensitySum","SCAall","SCAcut","ZeroScore","DiffSum","FDR","PEP","SCA","CorrelationScoreMedian","CorrelationScoreMean","CorrelationScoreSD","CorrelationScoreSum","CorrelationScoreSumSpecial","CorrelationScoreSumSpecial2","Score","Peaks1","Peaks2","ScaBefore","ScaBefore2","ScaAfter","ScaAfter2","RF_Scores","DL_Scores","GBM_Scores","GLM_Scores","FDR_RFspl_peptide","FDR_RF_peptide","FDR_RF_all","FDR_DL_all","FDR_GLM_all","FDR_GBM_all","SCA_mz_fdr","SCAfdr","Peak","precursor_ID","Precursor_Rating","rawfile")
            NAMES <- unique(NAMES)
            MORDER <- match(NAMES,colnames(tempa_dt))
            if(any(is.na(MORDER))){
              try({showNotification(session,paste("missing columns in table",fifu))},silent = T)
              for(i in which((is.na(MORDER)))){
                tempa_dt <- cbind(tempa_dt,as.double(NA))
              }
              colnames(tempa_dt)[(dim(tempa_dt)[2]-(sum(is.na(MORDER))-1)):dim(tempa_dt)[2]] <- NAMES[is.na(MORDER)]
            }
            tempa_dt <- tempa_dt[,.SD,.SDcols=NAMES]
            tempa_dt <- setcolorder(tempa_dt,NAMES)
            
            HUI <- melt(tempa_dt,measure.vars = id,id.vars = setdiff(1:dim(tempa_dt)[2],id))
            HUI <- HUI[value >-1,]
            # if(dim(HUI)[2]==46){
            #   NAMES46 <<- colnames(HUI)
            # }
            # if(dim(HUI)[2]==47){
            #   NAMES47 <<- colnames(HUI)
            # }
            
            
            
            if(dim(HUI)[1]>0){
              if(file.exists(PeaksTableName)){
                header_b <- F
              }else{
                header_b <- T
                
              }
              write.table(HUI,PeaksTableName,sep = "\t",row.names = F,quote = F,append = T,col.names  = header_b)
              
            }
          })
          
          
          
          
          # print(head(DP))
          if(cl != "try-error"){
            DPA <<- DP
            DP <- DP
            if(!all(is.na(DP$DL_Scores))){
              
              DP$QuantitationType[DP$QuantitationType=="intensity"] <-  "Intensities"
              fifu <- fifu
              fifu <- gsub("__","-",fifu)
              DP <- data.frame(DP)
              # trans <- DP[,1:(which(colnames(DP) == "Q1")-1)]
              # info <- DP[,(which(colnames(DP) == "Q1")):dim(DP)[2]]
              
              
              DPmelt <- melt(data.table(DP)[!is.na(DP$DL_Scores),],id.vars=(which(colnames(DP) == "Q1")):dim(DP)[2])
              data.table::setnames(DPmelt,"value","Intensity")
              data.table::setnames(DPmelt,"variable","Fragment")
              
              AddInf <- unlist(strsplit(fifu,"_"))
              AddInf[1] <- gsub("^mz","",AddInf[1])
              AddInf[length(AddInf)] <- gsub(".txt","",AddInf[length(AddInf)])
              if(length(AddInf)< 6){
                AddInf <- c(AddInf,rep(NA,6-length(AddInf)))
                
              }
              if(length(AddInf)> 6){
                AddInf[6] <- paste(AddInf[6:length(AddInf)],collapse = "_")
                AddInf <- AddInf[1:6]
              }
              
              longtab <- DPmelt
              MaIn <- matrix(AddInf,dim(longtab)[1],length(AddInf),byrow = T)
              MaIn <- cbind(MaIn,PeakType,pw,PrecursorRating)
              longtab$PeakType <- PeakType
              longtab$PrecursorRating <- PrecursorRating
              
              InfInsert <- data.table(t(matrix((AddInf))))
              longtab[,c("m/z","z","Sequence","Accession","Gene Symbol","Library_ScanID"):={
                InfInsert
              },]
              if(length(MaIn)>0){
                
                # colnames(MaIn) <- c("m/z","z","Sequence","Accession","Gene Symbol","Library_ScanID","PeakType","Peak Detection Width","Rating")
                # longtab <- cbind(longtab,MaIn)
                # longtab <- unique(longtab)
                # print("Writing")
                # AddInf <- unlist(strsplit(fifu,"_"))
                if(!file.exists(fipath)){
                  try(write.table(longtab,file = fipath,quote = F,row.names = F,append = F,sep = "\t"))
                }else{
                  try(write.table(longtab,file = fipath,quote = F,row.names = F,col.names = F,append = T,sep = "\t"))
                }
              }
            }
            
            
          }
        }
      }
    })})
    
  }
  # Ratio Estimation
  li <- list.files("./ParallelExport",pattern=".Peaks.txt$",full.names = T)
  if(length(li)>0&ExportRatios){
    clusterExport(clu,c("PeaksProcessingLM","Fun.Fit.lm"))
    hu <-  parSapply(clu,li,function(peakspath){
      # hu <-  sapply(li,function(peakspath){
      print(peakspath)
      library(data.table)
      library(deming)
      if(file.exists(peakspath)){
        peaks <- fread(peakspath,sep ="\t",stringsAsFactors = F)
        peaks[,Precursor := gsub("_.*$","",precursor_ID)]
        peaks[,PrecursorInfo := substr(precursor_ID,gregexpr("_",precursor_ID)[[1]][1],nchar(precursor_ID))]
        peaks[,Sequence:= strsplit(precursor_ID,"_")[[1]][3],precursor_ID]
        peaks[,Label:= strsplit(Sequence,"#")[[1]][2],precursor_ID]
        peaks[,Sequence:= strsplit(Sequence,"#")[[1]][1],precursor_ID]
        
        peaks[,Charge:= strsplit(precursor_ID,"_")[[1]][[2]],precursor_ID]
        peaks[,Gene:= strsplit(precursor_ID,"_")[[1]][[5]],precursor_ID]
        
        if(any(is.na(peaks$mz))){
          peaks$mz <- as.numeric(gsub("^mz","",peaks$Precursor))
        }
        peaks <- peaks[!is.na(value),]
        Ratios <- PeaksProcessingLM(peaks,progressObj = NULL,outpath="./ParallelExport/")
        fwrite(Ratios,gsub("Peaks.txt$","Ratios.txt",peakspath),sep = "\t")
      }
      NULL
    })
  }
  
  stopCluster(clu)
  # system("open export.txt")
  # consolidating
  # combining Peaks
  write("","./ParallelExport/CombiningPeaksTable")
  li <- list.files("./ParallelExport",pattern=".Peaks.txt$",full.names = T)
  unlink("./export/Peaks.txt")
  l1 <- lapply(li,function(x){
    cat("\r",x)
    fwrite(fread(x,sep = "\t"),"./export/Peaks.txt",sep = "\t",append = T)
    NULL
  })
  # combining TransitionLists
  write("","./ParallelExport/CombiningTransitionLists")
  li <- list.files("./ParallelExport",pattern="^TransitionList.*.txt$$",full.names = T)
  unlink("export/TransitionList.txt")
  l2 <- lapply(li,function(x){
    fwrite(fread(x,sep = "\t"),"./export/TransitionList.txt",sep = "\t",append = T)
    NULL
  })
  # combining DB PeaksTables
  write("","./ParallelExport/Combining_DB_PeaksTables")
  
  sql <- list.files("./ParallelExport/",pattern="^CurrentExportDB#.*.sqlite$",full.names = T)
  db <- dbConnect(SQLite(),dbp)
  hum <- lapply(sql,function(xhu){
    library(RSQLite)
    xhu <<- xhu
    cat("\r",xhu)
    gbtemp <-dbConnect(SQLite(),xhu)
    try({
      
      
      Ta <- dbListTables(gbtemp)
      Ta <- grep("^PEAKS",Ta,value = T)
      if(length(Ta)>0){
        tabi <- sapply(Ta,function(fil){
          fil <<- fil
          try({
            dbr <- dbReadTable(gbtemp,fil)
            dbWriteTable(db,fil,dbr,overwrite=T)
          })
        })
        print(table(tabi))
        
      }
    })
    dbDisconnect(gbtemp)
    NULL
  })
  
})

# combining Ratios
write("","./ParallelExport/CombiningRatioTables")
li <- list.files("./ParallelExport",pattern="Ratios.txt$",full.names = T)
unlink("./export/Ratios.txt")
l2 <- lapply(li,function(x){
  fwrite(fread(x,sep = "\t"),"./export/Ratios.txt",sep = "\t",append = T)
  NULL
})

write("","./ParallelExport/Finished_Parallel_Export")
dbDisconnect(db)

# 
# 
# peakspath <- "./ParallelExport/Peaks.txt"
# if(file.exists(peakspath)){
#   peaks <- fread(peakspath,sep ="\t",stringsAsFactors = F)
#   peaks[,Precursor := gsub("_.*$","",precursor_ID)]
#   peaks[,PrecursorInfo := substr(precursor_ID,gregexpr("_",precursor_ID)[[1]][1],nchar(precursor_ID))]
#   peaks[,Sequence:= strsplit(precursor_ID,"_")[[1]][3],precursor_ID]
#   peaks[,Label:= strsplit(Sequence,"#")[[1]][2],precursor_ID]
#   peaks[,Sequence:= strsplit(Sequence,"#")[[1]][1],precursor_ID]
#   
#   peaks[,Charge:= strsplit(precursor_ID,"_")[[1]][[2]],precursor_ID]
#   peaks[,Gene:= strsplit(precursor_ID,"_")[[1]][[5]],precursor_ID]
#   
#   if(any(is.na(peaks$mz))){
#     peaks$mz <- as.numeric(gsub("^mz","",peaks$Precursor))
#   }
#   peaks <- peaks[!is.na(value),]
#   Ratios <- PeaksProcessingLM(peaks,progressObj = progress_Export)
#   fwrite(Ratios,"./export/Ratios.txt",sep = "\t")
# }
