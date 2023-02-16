ExpandFragmentlist_Bruker <- function(PEPTIDESexport,precursorInfo,Tselect,Reverse = F,rname="NoName"){
  # Wrapper for  Vali compatible output
  vec <- c("charge","rawfile","R","R_p","IntensitySum","SCAall","SCAcut","MatchCount","ZeroScore","DiffSum","inv_ion_mobility","RT","scan")
  
  PEPTIDESexport$charge <- precursorInfo$Charge[1] # not relevant
  PEPTIDESexport$rawfile <- rname # done
  if(length(PEPTIDESexport$retention_time)>0){
    PEPTIDESexport$RT <- PEPTIDESexport$retention_time
  }else{
    # PEPTIDESexport$RT <- PEPTIDESexport$retention_time
    
  }
  
  # Expand Rawfile with Faims information if FaimsInfo is available:
  if(length(PEPTIDESexport$FAIMS_cv)>0){
    if(any(PEPTIDESexport$FAIMS_cv!="")){
      PEPTIDESexport$rawfile <- paste(basename(r),PEPTIDESexport$FAIMS_cv,sep="#CV")
    }
  }
  PEPTIDESexport$FAIMS_cv <- NULL
  
  Fragments <- setdiff(colnames(PEPTIDESexport),vec)  
  
  
  PEPTIDESexport$R <- 0 # DONE
  PEPTIDESexport$R_p <- 1 # DONE
  PEPTIDESexport$SCAall <- 0 #DONE
  PEPTIDESexport$SCAcut <- 0 #DONE
  PEPTIDESexport$IntensitySum <- 0 # DONE
  PEPTIDESexport$ZeroScore <- 0 # UseLess
  PEPTIDESexport$DiffSum <- 0 # DONE
  PEPTIDESexport$RT_min <- PEPTIDESexport$RT#/60
  PEPTIDESexport[,SCAall:={
    temp <- .SD
    temp <- unlist(temp)
    temp[temp == -1] <- 0
    NSA <- 0
    try({NSA <- NormSpecAngle(temp,Tselect$Intensities)},silent = T)
    
  },.SDcols=Fragments,by=scan]
  
  PEPTIDESexport[,SCAcut:={
    temp <- .SD
    temp <- unlist(temp)
    temp[temp == -1] <- 0
    NormSpecAngle(temp[1:6],Tselect$Intensities[1:6])
    
  },.SDcols=Fragments,by=scan]
  
  PEPTIDESexport[,IntensitySum:={
    temp <- .SD
    apply(temp,1,function(x){sum(x[x!=-1])})
    
    
  },.SDcols=Fragments,by=scan]
  
  
  COR <- PEPTIDESexport[,{
    temp <- .SD
    temp <- unlist(temp)
    temp[temp == -1] <- 0
    LI <- list(R=0,R_p=1)
    try({
      ct <- cor.test(Tselect$Intensities,temp)
      LI <- list(R= ct$estimate,R_p= ct$p.value)
    },silent = T)
    LI
    
  },.SDcols=Fragments,by=scan]
  COR[!is.na(COR$R)]
  
  COR <- COR[match(PEPTIDESexport$scan,scan),,]
  COR$R[is.na(COR$R)] <- 0
  COR$R_p[is.na(COR$R)] <- 1
  PEPTIDESexport$R <- COR$R
  PEPTIDESexport$R_p <- COR$R_p
  
  
  ColOrder <- c(Fragments,vec)
  PEPTIDESexport <-  setcolorder(PEPTIDESexport,ColOrder)
  setdiff(ColOrder,names(PEPTIDESexport))
  
  mz <- precursorInfo$mz
  z = precursorInfo$Charge
  Use_ModifiedSequence <- T
  if(Use_ModifiedSequence){
    Modfun <- grep("modified.sequence",names(precursorInfo),ignore.case = T)
    if(length(Modfun)==1){
      se <- as.data.frame(precursorInfo)[,Modfun]
      
    }else{
      print("Cannot find corred modified sequence column. Switching to Sequence")
      se <- precursorInfo$Sequence
      
    }
    se <- gsub("^_","",se)
    se <- gsub(" ","",se)
    se <- gsub(" ","",se)
    se <- gsub(".","p",se,fixed = T)
    se <- gsub("_",".",se,fixed = T)
    se <- gsub("..",".",se,fixed = T)
    se <- gsub(".","",se,fixed = T)
    
  }else{
    se = precursorInfo$Sequence
    
  }
  id= precursorInfo$SpecID
  OutnameInfo <- paste(mz,z,se,precursorInfo$Proteins,precursorInfo$GeneSymbol,precursorInfo$SpecID,sep = "_")
  if(Reverse){
    name_pep <- paste(paste(paste("mz",mz,sep = ""),z,se,id,sep = "_"),"Rev+.txt",sep = "")
    
  }else{
    name_pep <- paste(paste(paste("mz",mz,sep = ""),z,se,id,sep = "_"),".txt",sep = "")
    
  }
  infoname <- paste(paste(paste("",mz,sep = ""),z,se,id,sep = "_"),".info",sep = "")
  write(OutnameInfo,file=infoname)
  if(file.exists(name_pep)){
    temp <- fread(name_pep,sep = "\t",nrows = 1)
    if(length(colnames(temp))!=length(colnames(PEPTIDESexport))){
      # PEPTIDESexport <<- PEPTIDESexport
      # precursorInfo <<- precursorInfo
      # Tselect <<- Tselect
      fwrite(PEPTIDESexport,paste("Wrong",name_pep,sep="_"),sep = "\t",append = F)
      
      stop("Length of tables differs.")
      
    }else{
      if(any(colnames(temp)!=colnames(PEPTIDESexport))){
        # PEPTIDESexport <<- PEPTIDESexport
        # precursorInfo <<- precursorInfo
        # Tselect <<- Tselect
        fwrite(PEPTIDESexport,paste("Wrong",name_pep,sep="_"),sep = "\t",append = F)
        
        stop("Different Column Names in tables.")
      }
    }
  }
  fwrite(PEPTIDESexport,name_pep,sep = "\t",append = T)
  
}
# ExpandFragmentlist_Bruker(PEPTIDES,precursorInfo,Tselecttemp[Decoy==FALSE],Reverse = F,rname=basename(r))
# ExpandFragmentlist_Bruker(DECOYPEPTIDES,precursorInfo,Tselecttemp[Decoy==FALSE],Reverse=T,rname=basename(r))
# stop()

Rawrr_Vali_Extracter_DIA <- function(r,dirout,outpath,ST,TT,TTdec,ppm1=10,ppm2=10,scans_atOnce=2000,ShowInfo=F,extraSink=T,SystemPath=NULL,session=NULL,DIA=T,ExtractMS1 = F,SaveScans=T,  descriptionfilter =NULL,InvertSearch=T){
  # if(extraSink){
  #   sink("SessionSink.txt",split = T,type="message")
  # }
  #
  
  setwd("E:/MaxQuant/HZ/TimsTofTest/")
  # setwd("/Users/henno/Documents/Skripte/R-Functions/Selbach-Functions/ToolTestArea/Vali_TestExperiments/Bruker")
  wd <- load("Temp_ExtractorData.rda")
  l <- load("Bruker_EVENTS_COMPLETE_Selected.rda")
  descriptionfilter =ScanStringFilter
  InvertSearch=ScanStringFilterInvert
  r<- raws[1]
  library(rawrr)
  require(compiler)
  require(data.table)
  require(Peptides)
  library(RSQLite)
  library(opentimsr)
  
  try({
    # this loop goes through the raw files
    setwd(outpath)
    # outpath <<- outpath
    # r <<- r
    source(paste(SystemPath,"R/Vali_Functions.R",sep = "/"))
    # colnames(RAW <- read.raw(file = r,rawDiag = F))
    # system.time(RAW <- read.raw(file = r,rawDiag = T))
    # rawrr::installRawFileReaderDLLs()
    # rawrr::installRawrrExe()
    # system.time(RAW <- rawrr::readIndex(r)) # rawrr replacement
    brukerFile <- OpenTIMS(r)
    RAW <-as.data.table(brukerFile@frames)
    dir.create(dirout)
    setwd(dirout)
    dir.create(rpath <- gsub(".raw$","",basename(r)))
    setwd(rpath)
    cat(paste("\rReading RawFile",r))
    if(extraSink){
      msgcon <- file("SessionSink.txt", open = "a")
      sink(msgcon,append = T,split = F,type = "message")
      sink("SessionSink:output.txt", type = "output", append = TRUE, split = TRUE)
      
    }
    if(length(session) != 0){
      Max <- TTdec[,1,.(SpecID,Matches)]
      progress_internal <- Progress$new(session,min=0, 0)
      on.exit(progress_internal$close())
      progress_internal$set(message = 'Collecting MS1 Scans',
                            detail = "",value = 0)
    }
    
    ms1scans <- RAW[RAW$MsMsType==0,] # extracts scan numbers for ms1scans # rawrr replacement
    
    if(length(session) != 0){
      Max <- TTdec[,1,.(SpecID,Matches)]
      progress_internal2 <- Progress$new(session,min=0, 0)
      on.exit(progress_internal2$close())
      progress_internal2$set(message = 'Collecting MS2 Scans',
                             detail = "",value = 0)
    }
    
    ms2scans <- RAW[RAW$MsMsType!=0,] # actual
    
    
    # if(length(descriptionfilter)!=0){
    #   if(is.character(descriptionfilter)){
    #     foundi <- grepl(descriptionfilter,RAW$ScanDescription)
    #     if(any(foundi)){
    #       if(InvertSearch){
    #         ms2scans <- RAW[RAW$MSOrder == "Ms2"&!foundi,] # Hack Surequant
    #         
    #       }else{
    #         ms2scans <- RAW[RAW$MSOrder == "Ms2"&foundi,] # Hack Surequant
    #         
    #       }
    #     }
    #     
    #   }
    # }
    
    
    
    
    
    ##
    #
    Precursors <- ST$mz
    PrecursorsWindows <- sapply(Precursors,ppmWindow,ppm = ppm1)
    # limitToMatchingMasses <- T
    # if(limitToMatchingMasses&!DIA){
    #   WorthToBeScanned <- sapply( ms2scans$precursorMass,function(x){
    #     any(PrecursorsWindows[1,1]<=x&PrecursorsWindows[2,1]<=x)
    #   })
    #   ms2scans <- ms2scans[WorthToBeScanned,]
    #   
    # }
    
    
    
    DiaFrameMsMsInfo <- opentimsr::table2df(brukerFile,"DiaFrameMsMsInfo")[[1]]
    DiaFrameMsMsWindowGroups <- opentimsr::table2df(brukerFile,"DiaFrameMsMsWindowGroups")[[1]]
    DiaFrameMsMsWindows <- data.table(opentimsr::table2df(brukerFile,"DiaFrameMsMsWindows")[[1]])
    
    
    sn <- ms2scans$Id
    sn <- sort(sn)
    AllScansCount  <- length(sn)
    scansToScan <- T
    startID <- 1
    ITXOVERchecker <<- 0
    while(scansToScan){
      cat("\r",length(sn),"Scans ToGo")
      if(length(sn)<scans_atOnce){
        scans_atOnce <- length(sn)
        scansToScan <- F
      }
      temp_sn <- sn[sel <- startID:(startID+scans_atOnce)]
      temp_sn <- temp_sn[!is.na(temp_sn)]
      sn <- sn[-sel]
      Spectra <- query(brukerFile,temp_sn)
      Spectra <- data.table(Spectra,integer64="double")
      EVENTS_COMPLETE<- Spectra[,{
        temp <- .SD
        dfm <- DiaFrameMsMsInfo[DiaFrameMsMsInfo$Frame==.BY$frame,]
        if(length(dfm$WindowGroup)>0){
          DIAWINDOWS <- data.table(DiaFrameMsMsWindows[DiaFrameMsMsWindows$WindowGroup==dfm$WindowGroup])
          li <- DIAWINDOWS[,{ScanNumBegin:ScanNumEnd},.(IsolationMz,IsolationWidth,CollisionEnergy)]
          li2 <- li[match(temp$scan,li$V1)]
          ex <- cbind(temp,li2)
        }else{
          ex <- NULL
        }
        ex
      },frame]
      EVENTS_COMPLETE <- EVENTS_COMPLETE[!is.na(IsolationMz)]
      
      # Generating DIA Windows:
      DIAwindows_subset <- data.frame(lo = EVENTS_COMPLETE$IsolationMz -EVENTS_COMPLETE$IsolationWidth/2,
                                      up= EVENTS_COMPLETE$IsolationMz +EVENTS_COMPLETE$IsolationWidth/2)
      
      
      cat("\rRead",max(temp_sn),"from",max(sn))
      sapply(1:dim(PrecursorsWindows)[2],function(it){
        
        # t
        if(length(session)!=0){
          progress_internal$set(message = r,
                                detail = "Scanning Fragment Masses",value = it)
        }
        it <<- it
        lib_precursor <- PrecursorsWindows[,it]
        precursorInfo <- ST[it,]
        if(ShowInfo){
          cat("\nWorking on window:",lib_precursor)
        }else{
          cat("\r",it,"from",dim(PrecursorsWindows)[2])
          
          if(it%%1==0){
            LI <- list.files(pattern="^SessionStatusMS2",)
            if(length(LI)>0){
              sapply(LI,unlink)
            }
            # write(it,paste("SessionStatusMS2",it,dim(PrecursorsWindows)[2],sep = "_"))
            write(it,paste("SessionStatusMS2",AllScansCount-length(sn),AllScansCount,sep="_"))
            
          }
        }
        mz <- EVENTS_COMPLETE$IsolationMz
        if(DIA){
          TimeSelect <- EVENTS_COMPLETE$retention_time >= precursorInfo$Start&EVENTS_COMPLETE$retention_time<=precursorInfo$End
          table(EVENTS_COMPLETE$IsolationMz[TimeSelect])
          if(length(TimeSelect)==0){
            TimeSelect <- rep(T,length(ms2scans_subset$startTime))
          }
          DIAselect  <- (EVENTS_COMPLETE$IsolationMz-EVENTS_COMPLETE$IsolationWidth/2)<= ST$mz[it] & 
            (EVENTS_COMPLETE$IsolationMz+EVENTS_COMPLETE$IsolationWidth/2)>= ST$mz[it] 
          
          EVENTS_COMPLETE_Selected <- EVENTS_COMPLETE[DIAselect&TimeSelect]
          
        }else{
          EVENTS_COMPLETE_Selected <- EVENTS_COMPLETE$scanNumber[which(IsolationMz>=lib_precursor[1]&IsolationMz<=lib_precursor[2])]
          # if(length(EVENTS)>0){
          #   stop()
          # }
        }
        
        data.table::setnames(EVENTS_COMPLETE_Selected,c("mz","IsolationMz"),c("masses","mz"))
        if(length(EVENTS_COMPLETE_Selected)>0){
          # new readscans section 20200810 #########
          # ScansMissing <- setdiff(paste("scans",EVENTS,sep = ""),dbListTables(dbscans))
          # ScansAvailable <- intersect(paste("scans",EVENTS,sep = ""),list.files("scans"))
          # ScansMissingNum <- gsub("scans","",ScansMissing)
          # ScansAvailableNum <- gsub("scans","",ScansAvailable)
          # cat(paste("\rScans to read",length(ScansMissing)))
          # cat(paste("\rLoadable Scans",length(ScansAvailableNum)))
          
          # if(length(ScansMissing)>0){
          #   
          #   Mass <- Spectra[match(ScansMissingNum,names(Spectra))]#readScans(r,sort(ScansMissingNum),)
          #   Worked <- lapply(Mass,function(x){
          #     if(ShowInfo){
          #       
          #       cat("\r",x$title)
          #     }
          #     # x <<- x
          #     tempx <- x
          #     masses <- (tempx$mZ)
          #     Int <- (tempx$intensity)
          #     RT <- tempx$rtinseconds
          #     mz <- tempx$mZ
          #     # if(length(mz)==0){
          #     #   Xit <- grep("mZ",names(tempx))
          #     #   mz <- tempx[,.SD,.SDcols=xit]
          #     # }
          #     st <- tempx$scanType
          #     FI <- data.table(masses=masses,Int=Int,RT=RT,mz=mz,st=st,scan = tempx$scan)
          #   })  
          #   
          #   names(Worked) <- ScansMissing
          #   if(SaveScans&F){
          #     print("Saving Tables")
          #     lapply(1:length(Worked),function(x){
          #       # cat("\rSaving Tables",names(Worked)[x])
          #       # dbWriteTable(dbscans,names(Worked)[x],Worked[[x]],overwrite=T  )
          #       workedfile=Worked[[x]]
          #       save(workedfile,file=paste("scans",names(Worked)[x],sep = "/"))
          #       
          #       return(NULL)
          #     })
          #   }
          #   
          #   Workeddt <- rbindlist(Worked) 
          #   
          # }else{
          #   cat("no missing scans")
          # }
          if(ShowInfo){
            
            # cat("\rWorking on window:",lib_precursor,"Found",length(EVENTS),"Spectra candidates")
          }
          
          
          
          # if(length(ScansAvailable)>0){
          #   print(paste("Reloading", length(ScansAvailable),"Scans"))
          #   WorkedAvailable <- lapply(paste("scans",ScansAvailable,sep="/"),function(x){
          #     # cat("\r reloading",x)
          #     load(x)
          #     workedfile
          #   })
          #   
          #   WorkedAvailable <- rbindlist(WorkedAvailable)
          #   if(length(ScansMissing)==0){
          #     Workeddt <- WorkedAvailable
          #   }
          # }
          # if(length(ScansAvailable)>0&length(ScansMissing)>0){
          #   Workeddt <- rbind(Workeddt,WorkedAvailable)
          # }
          #########
          # if(0){
          #   
          #   ##### OLD scanning section
          #   Mass <- readScans(r,sort(EVENTS))
          #   if(ShowInfo){
          #     
          #     cat("\rWorking on window:",lib_precursor,"Found",length(EVENTS),"Spectra candidates")
          #   }
          #   
          #   Worked <- lapply(Mass,function(x){
          #     if(ShowInfo){
          #       
          #       cat("\r",x$title)
          #     }
          #     # x <<- x
          #     tempx <- x
          #     masses <- (tempx$mZ)
          #     Int <- (tempx$intensity)
          #     RT <- tempx$rtinseconds
          #     mz <- tempx$mZ
          #     # if(length(mz)==0){
          #     #   Xit <- grep("mZ",names(tempx))
          #     #   mz <- tempx[,.SD,.SDcols=xit]
          #     # }
          #     st <- tempx$scanType
          #     FI <- data.table(masses=masses,Int=Int,RT=RT,mz=mz,st=st,scan = tempx$scan)
          #   })  
          #   Workeddt <- rbindlist(Worked) 
          #   
          # }
          
          
          
          
          
          # data.table with fragments in the list
          # Adding FaimsInfo 
          # Workeddt[,FAIMS_cv := gsub(" .*","",gsub(".*cv=","",st)),st]
          # Workeddt$FAIMS_cv[is.na(as.numeric(Workeddt$FAIMS_cv))] <- ""
          # if(ShowInfo){
          #   
          #   cat("\rWorking on window:",lib_precursor,"Extracted Information from",dim(EVENTS_COMPLETE_Selected)[2],"Spectra")
          # }
          
          # Mapping Matching Spectra and writing to a file
          # Check if Faims cv has been applied:
          
          Matched <- ST[,{
            cat("\r",.GRP)
            gr <-.BY
            temp <- .SD
            gr2 <<- gr
            tempi <<- temp
            wi <- ppmWindow(gr$mz,ppm1)
            
            if((wi[1]<= mean(lib_precursor)&wi[2]>= mean(lib_precursor))|DIA){
              # Extract Masses
              if(ShowInfo){
                print(gr$SpecID)
                
              }
              Tselect <- TT[SpecID== gr$SpecID,] # library
              Tselect <- Tselect[order(as.numeric(Intensities),decreasing = T),]
              Tselect$Decoy <-F
              TselectDEC <- TTdec[SpecID== gr$SpecID,]
              TselectDEC <- TselectDEC[order(as.numeric(Intensities),decreasing = T),]
              TselectDEC$Decoy <- T
              TselectDEC$Matches <- paste(TselectDEC$Matches,"_DECOY",sep = "")
              
              Tselect <- rbind(Tselect,TselectDEC)
              
              # windows for Precursor in Library
              WiMapps <- Tselect[,{
                wi <- ppmWindow(Masses,ppm2)
                list(lo = wi[1],up = wi[2])
              },.(Masses,Matches)]
              WiMapps <- WiMapps[match(Tselect$Matches,WiMapps$Matches)]
              
              
              enableJIT(3)
              LE <- length(unique(EVENTS_COMPLETE_Selected$scan))
              
              system.time( PEPTIDES_init <- EVENTS_COMPLETE_Selected[,{
                if(ShowInfo){
                  
                  cat("\r Extracting Fragments:",.GRP,"from",LE)
                }
                temp3 <- .SD
                # stop()
                # temp2 <<- temp3
                gr2 <- .BY
                if(length(temp3$Int)==0){
                  try(setnames(temp3,"intensity","Int"),silent=F)
                }
                MAPPS <- WIMAPPS(WiMapps,temp3)
                MAPPS_compiler(MAPPS)
              },.(retention_time,scan,inv_ion_mobility)])
              
              table(PEPTIDES_init$MatchCount)
              enableJIT(3)
              
              # charge <- which(colnames(PEPTIDES)=="charge")
              PEPTIDES_init$inv_ion_mobility_DECOY <- PEPTIDES_init$inv_ion_mobility
              PEPTIDES_init$RT <- PEPTIDES_init$retention_time
              PEPTIDES_init$RT_min <- PEPTIDES_init$retention_time
              
              decoys <- grepl("_DECOY$",colnames(PEPTIDES_init),perl = T)
              PEPTIDES_init <- PEPTIDES_init
              Tselecttemp <- Tselect
              DECOYPEPTIDES <- PEPTIDES_init[,.SD,.SDcols = colnames(PEPTIDES_init)[decoys]]
              PEPTIDES <- PEPTIDES_init[,.SD,.SDcols = colnames(PEPTIDES_init)[!decoys]]
              
              DECOYPEPTIDES$scan <- PEPTIDES_init$scan
              DECOYPEPTIDES$RT <- 0
              colnames(DECOYPEPTIDES) <- gsub("_DECOY$","",colnames(DECOYPEPTIDES))
              
              precursorInfo$FAIMS_CV <- NULL
              # PEPTIDES <<- PEPTIDES
              ExpandFragmentlist_Bruker(PEPTIDES,precursorInfo,Tselecttemp[Decoy==FALSE],Reverse = F,rname=basename(r))
              ExpandFragmentlist_Bruker(DECOYPEPTIDES,precursorInfo,Tselecttemp[Decoy==FALSE],Reverse=T,rname=basename(r))
              
              
            }else{
              PEPTIDES <- NULL
            }
            # as.list(PEPTIDES)
            NULL
            
          },.(mz,SpecID)]
          
        }
        
      })
      
      
    }
    
    
    
    dbDisconnect(dbscans)
    
    if(ExtractMS1&F){
      EVENTS <- ms1scans$scan#[which(mz>=lib_precursor[1]&mz<=lib_precursor[2])]
      EVENTS <- sort(EVENTS)
      
      if(length(session) != 0){
        Max <- TTdec[,1,.(SpecID,Matches)]
        progress_internal2 <- Progress$new(session,min=0, length(ST$mz))
        on.exit(progress$close())
        progress_internal2$set(message = r,
                               detail = "Scanning Precursor Masses",value = 0)
      }
      
      lapply(1:length(ST$mz),function(itx){
        mz <- ST$mz[itx]
        
        if(length(session)!=0){
          progress_internal2$set(message = r,
                                 detail = "Scanning Fragment Masses",value = itx)
        }
        
        
        cat("\rExtracting MS1 Intensities for",mz,"#",itx,"from",length(ST$mz))
        if(itx%%10==0){
          LI <- list.files(pattern="^SessionStatusMS1",)
          if(length(LI)>0){
            sapply(LI,unlink)
          }
          write(itx,paste("SessionStatusMS1",itx,length(ST$mz),sep = "_"))
        }
        
        ch <- ST$Charge[itx]
        MASSES <- mz*ch-ch*1.00784+c(0:5)
        mzs <- (MASSES+ch*1.00784)/ch
        
        system.time(XICs <- readChromatogram(r,mass = mzs,tol = ppm1,type="bpc"))# unclear how this function works
        
        Ranges <- lapply(XICs,function(x){x$times})
        TRange <- sort(unique(unlist(Ranges)))
        ms1Int <- sapply(XICs,function(x){
          x <-x
          M <- match(TRange,x$times)
          Is <- x$intensities[M]
          if(length(Is) == 0){
            Is <- rep(NA,length(TRange))
          }
          Is
        })
        
        # Preparing Data.frame for export
        colnames(ms1Int) <- paste("mz",0:5,"+",sep="")
        ms1Int <- data.frame(ms1Int,check.names = F)
        ms1Int$Charge <- ch
        ms1Int$RT <- TRange
        ms1Int$mz <- mz
        ms1Int$Modified_Sequence <- ST$Modified.sequence[itx]
        ms1Int$Sequence <- ST$Sequence[itx]
        ms1Int$ScanNumber <- NA
        ms1Int$rawfile <- r
        se <- ST$Sequence[itx]
        id <- ST$SpecID[itx]
        name_pep <- paste(paste(paste("ms1scans",mz,sep = ""),ch,se,id,sep = "_"),".txt",sep = "")
        
        fwrite(ms1Int,name_pep,sep = "\t",append = T)
        
        
      })
    }
    
    if(length(session)!=0){
      progress_internal2$close()
      progress_internal$close()
      
    }
    
  })
  write("DONE","RawFile_Finished")
  if(extraSink){
    sink(NULL, type = "message")
    sink(NULL, type = "output")
  }
  
}
