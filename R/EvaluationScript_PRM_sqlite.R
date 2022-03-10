
################# General Functions #########
NumericMaker <- function(x){
  x <- matrix(x,byrow = T,ncol = 2)
  as.numeric(paste(x[,1],x[,2],sep = "."))
}
unfactor <- function(x){as.numeric(as.character(x))}
ppmWindow <- function(x,ppm=10){
  wi <- (x/10^6)*ppm
  c(x-wi/2,x+wi/2)
}
strsplitslot <- function(x,k=1,str = ";",strreplace = ";",...){
  unlist(sapply(strsplit(as.character(x),split = str,...),function(x){paste(x[k],collapse = strreplace)}))
}
PeptideHydrophobicity <- function(Seq){
  AA <- c("W","F","L","I","M","V","Y","A","T","P","E","D","C","S","Q","G","N","R","H","K")
  Rc <- c(11,10.5,9.6,8.4,5.8,5,4,0.8,0.4,0.2,0,-0.5,-0.8,-0.8,-0.9,-0.9,-1.2,-1.3,-1.3,-1.9)
  RcNT <- c(-4,-7,-9,-8,-5.5,-5.5,-3,-1.5,5,4,7,9,4,5,1,5,5,8,4,4.6)
  rc <- cbind(AA,Rc,RcNT)
  colnames(rc) <- c("AA","Rc","RcNt")
  try(rc <- data.frame(rc,stringsAsFactors = F),silent = T)
  rc$Rc <- as.numeric(rc$Rc)
  rc$RcNt <- as.numeric(rc$RcNt)
  
  if(!exists("rc")){
    "RcConstantsTable is Missing!"
  }
  SeqL <- unlist(strsplit(Seq,""))
  rcm <- rc[match(SeqL,rc$AA),]
  
  N <- nchar(Seq)
  
  K_L = 1
  if(N< 10){K_L = 1-0.027*(10-N)}
  if(N> 20){K_L = 1-0.014*(N-20)}
  
  
  
  if(N< 38){
    H = K_L*(sum(rcm$Rc)+0.42*rcm$RcNt[1]+0.22*rcm$RcNt[2]+0.05*rcm$RcNt[3])
  }
  if(N>= 38){
    H = K_L*(sum(rcm$Rc)+0.42*rcm$RcNt[1]+0.22*rcm$RcNt[2]+0.05*rcm$RcNt[3])-0.3*(K_L*(sum(rcm$Rc)+0.42*rcm$RcNt[1]+0.22*rcm$RcNt[2]+0.05*rcm$RcNt[3])-38)
    
  }
  
  return(H)
}
SquareNorm <- function(x){
  x^0.5/sum(x,na.rm = T)^0.5
}
cdotp <- function(x,y){
  HUI <- sapply(1:length(x),function(i){
    x[i]*y[i]
    
  })
  sum(HUI,na.rm =T)
}
NormSpecAngle <- function(x,y){
  1-2*acos(cdotp(SquareNorm(x),SquareNorm(y)))/pi
}
# Reading Lines Fast 
my.read.lines2=function(fname) {
  s = file.info(fname)$size 
  buf = readChar( fname, s, useBytes=T)
  strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]
}
dbtaName <- function(ana,dbpath){
  
  dbtaName <- paste(gsub("^mz","PEAKS",names(ana)),sep = "_")
  
  return(dbtaName)
}
dbtaNameExport <- function(ana,dbpath,pw){
  
  dbtaName <- paste(gsub("^mz","PEAKS",ana),sep = "_")
  return(dbtaName)
}

# Amino Acid Mass Information (used for MS1 Isotopic Cluster Prediction)
AAtab <- function(){
  #CH
  AAtab <- 
    "Gly	G	C2H3NO	57.021464	57.05	2	3	1	1	0
  Ala	A	C3H5NO	71.037114	71.08	3	5	1	1	0
  Ser	S	C3H5NO2	87.032029	87.08	3	5	1	2	0
  Pro	P	C5H7NO	97.052764	97.12	5	7	1	1	0
  Val	V	C5H9NO	99.068414	99.07	5	9	1	1	0
  Thr	T	C4H7NO2	101.04768	101.1	4	7	1	2	0
  Cys	C	C3H5NOS	103.00919	103.1	3	5	1	1	1
  Leu	L	C6H11NO	113.08406	113.2	6	11	1	1	0
  Ile	I	C6H11NO	113.08406	113.2	6	11	1	1	0
  Asn	N	C4H6N2O2	114.04293	114	4	6	2	2	0
  Asp	D	C4H5NO3	115.02694	115.1	4	5	1	3	0
Gln	Q	C5H8N2O2	128.05858	128.1	5	8	2	2	0
Lys	K	C6H12N2O	128.09496	128.2	6	12	2	1	0
  Glu	E	C5H7NO3	129.04259	129.1	5	7	1	3	0
  Met	M	C5H9NOS	131.04048	131.2	5	9	1	1	1
  Phe	F	C9H9NO	147.06841	147.2	9	9	1	1	0
Arg	R	C6H12N4O	156.10111	156.2	6	12	4	1	0
  Tyr	Y	C9H9NO2	163.06333	163.2	6	12	4	1	0
Trp	W	C11H10N2O	186.07931	186.2	11	10	2	1	0
  His	H	C6H7N3O	137.05891	137.1	6	7	3	1	0"
  AAtab <- strsplit(AAtab,"\n")
  AAtab <- t(as.data.frame(lapply(AAtab,strsplit,split = "\t")))
  rownames(AAtab) <- NULL
  colnames(AAtab) <-c("AA","AA_L","AA_Comp","Mono","Ave","C_count","H_count","N_count","O_count","S_count")
  AAtab <- data.frame(AAtab,stringsAsFactors = F)
  return(AAtab)
}
mzcalc <- function(x,charge,proton = 1.007276){
  outvec <- (x+charge*proton)/charge
  return(outvec)
}


################# Rawfile Scanner #############

## RawFile Scanner Starter, Main Function
PrepareTransitionList <- function(mainPath,maxquant,inclusionList,ppm = 10,ppm2 = 10,threads = 2,session= NULL,pythonpath = NULL,test = F,useDIA=F,ScanStringFilter=NULL,ScanStringFilterInvert=T){
  #mainPath="E:/Projects/Paolo_ACM2/20181206_ArgC_TryP_ACM/Tryp/"
  #inclusionList = "E:/Projects/Paolo_ACM2/20181206_ArgC_TryP_ACM/Tryp/IL/"
  # maxquant <<- maxquant
  # inclusionList <<- inclusionList
  # humpf <<- "fwef"
  # mainPath <<- mainPath
  # Py_MSFR_Module <<- Py_MSFR_Module
  # pythonpath <<- pythonpath[1]
  # Py_MSFR_Module <<- Py_MSFR_Module[1]
  
  inclusionList = gsub("\\\\","/",inclusionList)
  # Py_MSFR_Module = gsub("\\\\","/",Py_MSFR_Module)
  maxquant = gsub("\\\\","/",maxquant)
  mainPath = gsub("\\\\","/",mainPath)
  
  setwd(mainPath)
  
  
  Extract_Intensities_Raw(RawFilePath = mainPath,InclusionlistPath = inclusionList,ppm1 = ppm,
                          ppm2 =ppm2,useDIA = useDIA,session=session,threads = threads,ScanStringFilter = ScanStringFilter,ScanStringFilterInvert = ScanStringFilterInvert)
  # print("Start Compiling Database")
  maxquant <- NULL
  # session <- NULL
  wd <- getwd()
  try(RETI <- CompilingDatabase(mainPath,maxquant,session,test = F,SystemPath=SystemPath,threads = threads))
  setwd(wd)
  if(ExpandingSpecificity){
    try(ExpandSpectrumSpecificity(InclusionlistPath))
  }
  # try({rm(db)})
  
  # if(length(session) != 0){
  #   progress <- Progress$new(session,min=0, max=8)
  #   # on.exit(progress$close())
  #   progress$set(message = 'compiling MS1 Scans',
  #                detail = "This may take a while",value = 0)
  # }
  
  try(CompileMS(wd = "./",dbname = "./PRM_Analyzer_Matches/PickyAnalyzer.sqlite",session = NULL))
  return(NULL)
}
## loops through list of mass windows (WiMapps) and returns intensities of matching masses in temp
WIMAPPS <- function(WiMapps,temp){
  # loops through list of mass windows (WiMapps) and returns intensities of matching masses in temp
  MASSES <- as.numeric(temp$masses)
  Intensities <- apply(WiMapps,1,function(x){
    # x <<- x
    sel <- MASSES>=as.numeric(x[3])&MASSES<=as.numeric(x[4])
    if(any(sel)){
      SELE <- max(temp$Int[sel],na.rm = T)
      # Diff
      Diff <- max(as.numeric(x[1])-MASSES[sel],na.rm = T)
      
      return(c(SELE,Diff))
    }else{return(as.double(c(-1,NA)))}
  })
  colnames(Intensities)<- WiMapps$Matches
  Intensities
}
## results of WIMAPPS to proper list, which can be converted to a table using rbindlist
MAPPS_compiler <- function(MAPPS){
  #converts results of WIMAPPS to proper list, which can be converted to a table using rbindlist
  Intensities <- MAPPS[1,]
  # names(Intensities) <- attributes(MAPPS)$names[1:length(Intensities)]
  decoy <- grepl("_DECOY$",names(Intensities),perl = T)
  DiffSum  <- MAPPS[2,]
  DiffSumFUN<- function(DiffSum){
    DiffSum <- DiffSum[DiffSum!=-1]
    DiffSum <- sum(DiffSum)/length(DiffSum)
    DiffSum
  }
  
  Intensities <- as.list(Intensities)
  
  Intensities$MatchCount= sum(Intensities[!decoy]!=-1)
  Intensities$MatchCount_DECOY= sum(Intensities[decoy]!=-1)
  
  
  Intensities$DiffSum_DECOY=             DiffSumFUN(DiffSum[decoy])
  Intensities$DiffSum=             DiffSumFUN(DiffSum[!decoy])
  Intensities
}

## Wrapper for  Vali compatible output, called in RawDiag_Vali_Extracter
ExpandFragmentlist <- function(PEPTIDESexport,precursorInfo,Tselect,Reverse = F,rname="NoName"){
  # Wrapper for  Vali compatible output
  vec <- c("charge","rawfile","R","R_p","IntensitySum","SCAall","SCAcut","MatchCount","ZeroScore","DiffSum","RT","scan")
  
  PEPTIDESexport$charge <- precursorInfo$Charge[1] # not relevant
  PEPTIDESexport$rawfile <- rname # done
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
  PEPTIDESexport$RT_min <- PEPTIDESexport$RT/60
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
  write(OutnameInfo,infoname)
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

RawDiag_Vali_Extracter_DIA <- function(r,dirout,outpath,ST,TT,TTdec,ppm1=10,ppm2=10,ShowInfo=F,extraSink=T,SystemPath=NULL,session=NULL,DIA=T,ExtractMS1 = F,SaveScans=T,descriptionfilter =NULL,InvertSearch=T){
  # if(extraSink){
  #   sink("SessionSink.txt",split = T,type="message")
  # }
  #
  require(rawDiag)
  library(rawrr)
  require(compiler)
  require(data.table)
  require(Peptides)
  library(RSQLite)
  
  try({
    # this loop goes through the raw files
    setwd(outpath)
    # outpath <<- outpath
    # r <<- r
    source(paste(SystemPath,"R/EvaluationScript_PRM_sqlite.R",sep = "/"))
    system.time(RAW <- read.raw(file = r,rawDiag = F))
    # system.time(RAW <- rawrr::readIndex( r)) # rawrr replacement
    
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
    
    ms1scans <- RAW[RAW$MSOrder == "Ms",] # extracts scan numbers for ms1scans # rawrr replacement
    
    if(length(session) != 0){
      Max <- TTdec[,1,.(SpecID,Matches)]
      progress_internal2 <- Progress$new(session,min=0, 0)
      on.exit(progress_internal2$close())
      progress_internal2$set(message = 'Collecting MS2 Scans',
                             detail = "",value = 0)
    }
    
    ms2scans <- RAW[RAW$MSOrder == "Ms2",] # actual
    

    if(length(descriptionfilter)>0){
      if(is.character(descriptionfilter)){
        foundi <- grepl(descriptionfilter,RAW$ScanDescription)
        if(any(foundi)){
          if(InvertSearch){
            ms2scans <- RAW[RAW$MSOrder == "Ms2"&!foundi,] # Hack Surequant
            
          }else{
            ms2scans <- RAW[RAW$MSOrder == "Ms2"&foundi,] # Hack Surequant
            
          }
        }
        
      }
    }
    
    
    
    removeDataDependentScans <- F
    if(removeDataDependentScans){
      MS2_undependent <- grep("NSI Full ms2",ms2scans$ScanType,fixed = T)# thermo specific, data dependent will not be searched...
      if(length(MS2_undependent)>0){
        ms2scans <- ms2scans[MS2_undependent,]
      }
    }
    # plot(density(ms2scans$StartTime))
    # extracts scan numbers for ms2scans
    # plot(as.numeric(gsub(",",".",ms2scans$StartTime)),as.numeric(gsub(",",".",ms2scans$BasePeakIntensity)),type = "l")
    # plot(as.numeric(gsub(",",".",ms1scans$StartTime)),as.numeric(gsub(",",".",ms1scans$BasePeakIntensity)),type = "l")
    
    ###
    ##
    #
    Precursors <- ST$mz
    PrecursorsWindows <- sapply(Precursors,ppmWindow,ppm = ppm1)
    
    mz <- as.numeric(gsub(",",".",ms2scans$PrecursorMass))
    mzuni <- unique(mz)
    # for DIA
    IsolationWidth_spacer <- 2*0.5
    # currently switched OFF
    IsolationWidth_spacer <- 0
    
    IsolationWindow <- ms2scans$IsolationWidth - IsolationWidth_spacer
    if(IsolationWindow==0){
      IsolationWindow <- ms2scans$IsolationWidth #- IsolationWidth_spacer
      
    }
    
    if(length(session) != 0){
      Max <- TTdec[,1,.(SpecID,Matches)]
      progress_internal2 <- Progress$new(session,min=0, dim(PrecursorsWindows)[2])
      on.exit(progress_internal2$close())
      progress_internal2$set(message = r,
                             detail = "Scanning Fragment Masses",value = 0)
    }
    DIAwindows <- data.frame(lo = ms2scans$PrecursorMass -(IsolationWindow/2),up= ms2scans$PrecursorMass +(IsolationWindow/2))
    dbscans <- dbConnect(SQLite(),scans_backup <- paste(basename(r),".sqlite",sep = ""))
    SelectedScans <- apply(DIAwindows,1,function(x){
      any(x[1]<=Precursors&x[2]>=Precursors)
    })
    
    
    sn <- ms2scans$scanNumber[SelectedScans]
    sn <- sort(sn)
    AllScansCount  <- length(sn)
    scans_atOnce <- 2000
    scansToScan <- T
    startID <- 1
    ITXOVERchecker <<- 0
    while(scansToScan){
      if(length(sn)<scans_atOnce){
        scans_atOnce <- length(sn)
        scansToScan <- F
      }
      temp_sn <- sn[sel <- startID:(startID+scans_atOnce)]
      temp_sn <- temp_sn[!is.na(temp_sn)]
      sn <- sn[-sel]
      Spectra <- readScans(r,temp_sn)
      names(Spectra)<- sapply(Spectra,function(x){x$scan})
      ms2scans_subset <- ms2scans[MFUN <- match(temp_sn,ms2scans$scanNumber),]
      DIAwindows_subset <- DIAwindows[MFUN,]
      
      cat("\rRead",max(temp_sn),"from",max(sn))
      sapply(1:dim(PrecursorsWindows)[2],function(it){
        
        # t
        if(length(session)!=0){
          progress_internal$set(message = r,
                                detail = "Scanning Fragment Masses",value = it)
        }
        # it <<- it
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
        mz <- ms2scans_subset$PrecursorMass
        if(DIA){
          TimeSelect <- ms2scans_subset$StartTime >= precursorInfo$Start&ms2scans_subset$StartTime<=precursorInfo$End
          DIAselect  <- DIAwindows_subset$lo <= ST$mz[it] & DIAwindows_subset$up >= ST$mz[it]
          EVENTS <- ms2scans_subset$scanNumber[DIAselect&TimeSelect]
          table(ms2scans_subset$PrecursorMass[DIAselect&TimeSelect])
          range(ms2scans_subset$StartTime[DIAselect&TimeSelect])
        }else{
          EVENTS <- ms2scans_subset$scanNumber[which(mz>=lib_precursor[1]&mz<=lib_precursor[2])]
          
        }
        EVENTS <- EVENTS[!is.na(EVENTS)]
        if(length(EVENTS)>0){
          # new readscans section 20200810 #########
          ScansMissing <- setdiff(paste("scans",EVENTS,sep = ""),dbListTables(dbscans))
          ScansAvailable <- intersect(paste("scans",EVENTS,sep = ""),list.files("scans"))
          ScansMissingNum <- gsub("scans","",ScansMissing)
          ScansAvailableNum <- gsub("scans","",ScansAvailable)
          cat(paste("\rScans to read",length(ScansMissing)))
          cat(paste("\rLoadable Scans",length(ScansAvailableNum)))
          
          if(length(ScansMissing)>0){
            
            Mass <- Spectra[match(ScansMissingNum,names(Spectra))]#readScans(r,sort(ScansMissingNum),)
            Worked <- lapply(Mass,function(x){
              if(ShowInfo){
                
                cat("\r",x$title)
              }
              # x <<- x
              tempx <- x
              masses <- (tempx$mZ)
              Int <- (tempx$intensity)
              RT <- tempx$rtinseconds
              mz <- tempx$mZ
              # if(length(mz)==0){
              #   Xit <- grep("mZ",names(tempx))
              #   mz <- tempx[,.SD,.SDcols=xit]
              # }
              st <- tempx$scanType
              FI <- data.table(masses=masses,Int=Int,RT=RT,mz=mz,st=st,scan = tempx$scan)
            })  
            
            names(Worked) <- ScansMissing
            if(SaveScans&F){
              print("Saving Tables")
              lapply(1:length(Worked),function(x){
                # cat("\rSaving Tables",names(Worked)[x])
                # dbWriteTable(dbscans,names(Worked)[x],Worked[[x]],overwrite=T  )
                workedfile=Worked[[x]]
                save(workedfile,file=paste("scans",names(Worked)[x],sep = "/"))
                
                return(NULL)
              })
            }
            
            Workeddt <- rbindlist(Worked) 
            
          }else{
            cat("no missing scans")
          }
          if(ShowInfo){
            
            cat("\rWorking on window:",lib_precursor,"Found",length(EVENTS),"Spectra candidates")
          }
          
          
          
          if(length(ScansAvailable)>0){
            print(paste("Reloading", length(ScansAvailable),"Scans"))
            WorkedAvailable <- lapply(paste("scans",ScansAvailable,sep="/"),function(x){
              # cat("\r reloading",x)
              load(x)
              workedfile
            })
            
            WorkedAvailable <- rbindlist(WorkedAvailable)
            if(length(ScansMissing)==0){
              Workeddt <- WorkedAvailable
            }
          }
          if(length(ScansAvailable)>0&length(ScansMissing)>0){
            Workeddt <- rbind(Workeddt,WorkedAvailable)
          }
          #########
          if(0){
            
            ##### OLD scanning section
            Mass <- readScans(r,sort(EVENTS))
            if(ShowInfo){
              
              cat("\rWorking on window:",lib_precursor,"Found",length(EVENTS),"Spectra candidates")
            }
            
            Worked <- lapply(Mass,function(x){
              if(ShowInfo){
                
                cat("\r",x$title)
              }
              # x <<- x
              tempx <- x
              masses <- (tempx$mZ)
              Int <- (tempx$intensity)
              RT <- tempx$rtinseconds
              mz <- tempx$mZ
              # if(length(mz)==0){
              #   Xit <- grep("mZ",names(tempx))
              #   mz <- tempx[,.SD,.SDcols=xit]
              # }
              st <- tempx$scanType
              FI <- data.table(masses=masses,Int=Int,RT=RT,mz=mz,st=st,scan = tempx$scan)
            })  
            Workeddt <- rbindlist(Worked) 
            
          }
          
          
          
          
          
          
          
          # data.table with fragments in the list
          # Adding FaimsInfo 
          Workeddt[,FAIMS_cv := gsub(" .*","",gsub(".*cv=","",st)),st]
          Workeddt$FAIMS_cv[is.na(as.numeric(Workeddt$FAIMS_cv))] <- ""
          if(ShowInfo){
            
            cat("\rWorking on window:",lib_precursor,"Extracted Information from",dim(Workeddt)[2],"Spectra")
          }
          
          # Mapping Matching Spectra and writing to a file
          # Check if Faims cv has been applied:
          
          
          Matched <- ST[SpecID==precursorInfo$SpecID,{
            gr <-.BY
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
              LE <- length(unique(Workeddt$scan))
              
              system.time( PEPTIDES_init <- Workeddt[,{
                if(ShowInfo){
                  
                  cat("\r Extracting Fragments:",.GRP,"from",LE)
                }
                temp <- .SD
                MAPPS <- WIMAPPS(WiMapps,temp)
                
                # sort(WiMapps$Masses)
                # MAPPS$
                # temp$masses
                # Intensities <- apply(WiMapps,1,function(x){
                #   sel <- temp$masses>=x[3]&temp$masses<=x[4]
                #   if(any(sel)){
                #     SELE <- max(temp$Int[sel])
                #     return(SELE)
                #   }else{return(as.double(-1))}
                # })
                MAPPS_compiler(MAPPS)
              },.(RT,scan,FAIMS_cv)])
              table(PEPTIDES_init$MatchCount)
              enableJIT(3)
              
              
              
              
              # charge <- which(colnames(PEPTIDES)=="charge")
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
              ExpandFragmentlist(PEPTIDES,precursorInfo,Tselecttemp[Decoy==FALSE],Reverse = F,rname=basename(r))
              ExpandFragmentlist(DECOYPEPTIDES,precursorInfo,Tselecttemp[Decoy==FALSE],Reverse=T,rname=basename(r))
              
              
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
    
    if(ExtractMS1){
      EVENTS <- ms1scans$scanNumber#[which(mz>=lib_precursor[1]&mz<=lib_precursor[2])]
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
        
        XICs <- readXICs(r,masses = mzs,tol = ppm1)# unclear how this function works
        
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


Rawrr_Vali_Extracter_DIA <- function(r,dirout,outpath,ST,TT,TTdec,ppm1=10,ppm2=10,scans_atOnce=2000,ShowInfo=F,extraSink=T,SystemPath=NULL,session=NULL,DIA=T,ExtractMS1 = F,SaveScans=T,  descriptionfilter =NULL,InvertSearch=T){
  # if(extraSink){
  #   sink("SessionSink.txt",split = T,type="message")
  # }
  #
  require(rawDiag)
  library(rawrr)
  require(compiler)
  require(data.table)
  require(Peptides)
  library(RSQLite)
  
  try({
    # this loop goes through the raw files
    setwd(outpath)
    # outpath <<- outpath
    # r <<- r
    source(paste(SystemPath,"R/EvaluationScript_PRM_sqlite.R",sep = "/"))
    # colnames(RAW <- read.raw(file = r,rawDiag = F))
    # system.time(RAW <- read.raw(file = r,rawDiag = T))
    # rawrr::installRawFileReaderDLLs()
    # rawrr::installRawrrExe()
    system.time(RAW <- rawrr::readIndex(r)) # rawrr replacement
    
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
    
    ms1scans <- RAW[RAW$MSOrder == "Ms",] # extracts scan numbers for ms1scans # rawrr replacement
    
    if(length(session) != 0){
      Max <- TTdec[,1,.(SpecID,Matches)]
      progress_internal2 <- Progress$new(session,min=0, 0)
      on.exit(progress_internal2$close())
      progress_internal2$set(message = 'Collecting MS2 Scans',
                             detail = "",value = 0)
    }
    
    ms2scans <- RAW[RAW$MSOrder == "Ms2",] # actual
    
    
    if(length(descriptionfilter)!=0){
      if(is.character(descriptionfilter)){
        foundi <- grepl(descriptionfilter,RAW$ScanDescription)
        if(any(foundi)){
          if(InvertSearch){
            ms2scans <- RAW[RAW$MSOrder == "Ms2"&!foundi,] # Hack Surequant
            
          }else{
            ms2scans <- RAW[RAW$MSOrder == "Ms2"&foundi,] # Hack Surequant
            
          }
        }
        
      }
    }
    
    
    
    removeDataDependentScans <- F
    if(removeDataDependentScans){
      MS2_undependent <- grep("NSI Full ms2",ms2scans$ScanType,fixed = T)# thermo specific, data dependent will not be searched...
      if(length(MS2_undependent)>0){
        ms2scans <- ms2scans[MS2_undependent,]
      }
    }
    
    ##
    #
    Precursors <- ST$mz
    PrecursorsWindows <- sapply(Precursors,ppmWindow,ppm = ppm1)
    limitToMatchingMasses <- T
    if(limitToMatchingMasses&!DIA){
      WorthToBeScanned <- sapply( ms2scans$precursorMass,function(x){
        any(PrecursorsWindows[1,1]<=x&PrecursorsWindows[2,1]<=x)
      })
      ms2scans <- ms2scans[WorthToBeScanned,]
      
    }
    
    
    
    
    mz <- as.numeric(gsub(",",".",ms2scans$PrecursorMass))
    mzuni <- unique(mz)
    # for DIA
    # IsolationWidth_spacer <- 2*0.5
    # # currently switched OFF
    # IsolationWidth_spacer <- 0
    # 
    # IsolationWindow <- ms2scans$IsolationWidth - IsolationWidth_spacer
    # if(IsolationWindow==0){
    #   IsolationWindow <- ms2scans$IsolationWidth #- IsolationWidth_spacer
    #   
    # }
    # 
    # if(length(session) != 0){
    #   Max <- TTdec[,1,.(SpecID,Matches)]
    #   progress_internal2 <- Progress$new(session,min=0, dim(PrecursorsWindows)[2])
    #   on.exit(progress_internal2$close())
    #   progress_internal2$set(message = r,
    #                          detail = "Scanning Fragment Masses",value = 0)
    # }
    # DIAwindows <- data.frame(lo = ms2scans$PrecursorMass -(IsolationWindow/2),up= ms2scans$PrecursorMass +(IsolationWindow/2))
    dbscans <- dbConnect(SQLite(),scans_backup <- paste(basename(r),".sqlite",sep = ""))
    # SelectedScans <- apply(DIAwindows,1,function(x){
    #   any(x[1]<=Precursors&x[2]>=Precursors)
    # })
    # 
    # 
    sn <- ms2scans$scan#[SelectedScans]
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
      Spectra <- readSpectrum(r,temp_sn)
      # cat("\r",startID,"Scanning ROund Finished.")
      
      names(Spectra) <- sapply(Spectra,function(x){x$scan})
      # Filter Spectra
      spiInfo <- lapply(Spectra,function(spi ){
        spi <- spi
        list(PrecursorMass= spi$pepmass,ScanDescription=spi$`Scan Description:`,
             Iso_Width = as.numeric(spi$`MS2 Isolation Width:`),
             Offset= as.numeric(spi$`MS2 Isolation Offset:`),
             scanNumber=spi$scan,
             startTime=spi$rtinseconds)
      })
      ms2scans_subset <- rbindlist(spiInfo)
      table(ms2scans_subset$ScanDescription)
      # Applying Filter: 
      if(length(descriptionfilter)!=0){
        if(is.character(descriptionfilter)){
          foundi <- grepl(descriptionfilter,ms2scans_subset$ScanDescription)
          if(any(foundi)){
            if(InvertSearch){
              Spectra2 <- Spectra[wh <- which(!foundi)] # Hack Surequant
              ms2scans_subset <- ms2scans_subset[wh,]
            }else{
              Spectra2 <- Spectra[wh <- which(foundi)] # Hack Surequant
              ms2scans_subset <- ms2scans_subset[wh,]
              
            }
          }else{
            next()
          }
          
        }
      }
      
      # Generating DIA Windows:
      DIAwindows_subset <- data.frame(lo = ms2scans_subset$PepMass -(ms2scans_subset$Offset+ms2scans_subset$Iso_Width/2),
                                      up= ms2scans_subset$PepMass +(ms2scans_subset$Offset+ms2scans_subset$Iso_Width/2))
      
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
        mz <- ms2scans_subset$PrecursorMass
        if(DIA){
          stop("NotImplemented")
          TimeSelect <- ms2scans_subset$StartTime >= precursorInfo$Start&ms2scans_subset$StartTime<=precursorInfo$End
          DIAselect  <- DIAwindows_subset$lo <= ST$mz[it] & DIAwindows_subset$up >= ST$mz[it]
          EVENTS <- ms2scans_subset$scanNumber[DIAselect&TimeSelect]
          table(ms2scans_subset$PrecursorMass[DIAselect&TimeSelect])
          range(ms2scans_subset$StartTime[DIAselect&TimeSelect])
        }else{
          EVENTS <- ms2scans_subset$scanNumber[which(mz>=lib_precursor[1]&mz<=lib_precursor[2])]
          # if(length(EVENTS)>0){
          #   stop()
          # }
        }
        EVENTS <- EVENTS[!is.na(EVENTS)]
        if(length(EVENTS)>0){
          # new readscans section 20200810 #########
          ScansMissing <- setdiff(paste("scans",EVENTS,sep = ""),dbListTables(dbscans))
          ScansAvailable <- intersect(paste("scans",EVENTS,sep = ""),list.files("scans"))
          ScansMissingNum <- gsub("scans","",ScansMissing)
          ScansAvailableNum <- gsub("scans","",ScansAvailable)
          cat(paste("\rScans to read",length(ScansMissing)))
          cat(paste("\rLoadable Scans",length(ScansAvailableNum)))
          
          if(length(ScansMissing)>0){
            
            Mass <- Spectra[match(ScansMissingNum,names(Spectra))]#readScans(r,sort(ScansMissingNum),)
            Worked <- lapply(Mass,function(x){
              if(ShowInfo){
                
                cat("\r",x$title)
              }
              # x <<- x
              tempx <- x
              masses <- (tempx$mZ)
              Int <- (tempx$intensity)
              RT <- tempx$rtinseconds
              mz <- tempx$mZ
              # if(length(mz)==0){
              #   Xit <- grep("mZ",names(tempx))
              #   mz <- tempx[,.SD,.SDcols=xit]
              # }
              st <- tempx$scanType
              FI <- data.table(masses=masses,Int=Int,RT=RT,mz=mz,st=st,scan = tempx$scan)
            })  
            
            names(Worked) <- ScansMissing
            if(SaveScans&F){
              print("Saving Tables")
              lapply(1:length(Worked),function(x){
                # cat("\rSaving Tables",names(Worked)[x])
                # dbWriteTable(dbscans,names(Worked)[x],Worked[[x]],overwrite=T  )
                workedfile=Worked[[x]]
                save(workedfile,file=paste("scans",names(Worked)[x],sep = "/"))
                
                return(NULL)
              })
            }
            
            Workeddt <- rbindlist(Worked) 
            
          }else{
            cat("no missing scans")
          }
          if(ShowInfo){
            
            cat("\rWorking on window:",lib_precursor,"Found",length(EVENTS),"Spectra candidates")
          }
          
          
          
          if(length(ScansAvailable)>0){
            print(paste("Reloading", length(ScansAvailable),"Scans"))
            WorkedAvailable <- lapply(paste("scans",ScansAvailable,sep="/"),function(x){
              # cat("\r reloading",x)
              load(x)
              workedfile
            })
            
            WorkedAvailable <- rbindlist(WorkedAvailable)
            if(length(ScansMissing)==0){
              Workeddt <- WorkedAvailable
            }
          }
          if(length(ScansAvailable)>0&length(ScansMissing)>0){
            Workeddt <- rbind(Workeddt,WorkedAvailable)
          }
          #########
          if(0){
            
            ##### OLD scanning section
            Mass <- readScans(r,sort(EVENTS))
            if(ShowInfo){
              
              cat("\rWorking on window:",lib_precursor,"Found",length(EVENTS),"Spectra candidates")
            }
            
            Worked <- lapply(Mass,function(x){
              if(ShowInfo){
                
                cat("\r",x$title)
              }
              # x <<- x
              tempx <- x
              masses <- (tempx$mZ)
              Int <- (tempx$intensity)
              RT <- tempx$rtinseconds
              mz <- tempx$mZ
              # if(length(mz)==0){
              #   Xit <- grep("mZ",names(tempx))
              #   mz <- tempx[,.SD,.SDcols=xit]
              # }
              st <- tempx$scanType
              FI <- data.table(masses=masses,Int=Int,RT=RT,mz=mz,st=st,scan = tempx$scan)
            })  
            Workeddt <- rbindlist(Worked) 
            
          }
          
          
          
          
          
          
          
          # data.table with fragments in the list
          # Adding FaimsInfo 
          Workeddt[,FAIMS_cv := gsub(" .*","",gsub(".*cv=","",st)),st]
          Workeddt$FAIMS_cv[is.na(as.numeric(Workeddt$FAIMS_cv))] <- ""
          if(ShowInfo){
            
            cat("\rWorking on window:",lib_precursor,"Extracted Information from",dim(Workeddt)[2],"Spectra")
          }
          
          # Mapping Matching Spectra and writing to a file
          # Check if Faims cv has been applied:
          
          
          Matched <- ST[SpecID==precursorInfo$SpecID,{
            gr <-.BY
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
              LE <- length(unique(Workeddt$scan))
              
              system.time( PEPTIDES_init <- Workeddt[,{
                if(ShowInfo){
                  
                  cat("\r Extracting Fragments:",.GRP,"from",LE)
                }
                temp <- .SD
                MAPPS <- WIMAPPS(WiMapps,temp)
                
                # sort(WiMapps$Masses)
                # MAPPS$
                # temp$masses
                # Intensities <- apply(WiMapps,1,function(x){
                #   sel <- temp$masses>=x[3]&temp$masses<=x[4]
                #   if(any(sel)){
                #     SELE <- max(temp$Int[sel])
                #     return(SELE)
                #   }else{return(as.double(-1))}
                # })
                MAPPS_compiler(MAPPS)
              },.(RT,scan,FAIMS_cv)])
              table(PEPTIDES_init$MatchCount)
              enableJIT(3)
              
              
              
              
              # charge <- which(colnames(PEPTIDES)=="charge")
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
              ExpandFragmentlist(PEPTIDES,precursorInfo,Tselecttemp[Decoy==FALSE],Reverse = F,rname=basename(r))
              ExpandFragmentlist(DECOYPEPTIDES,precursorInfo,Tselecttemp[Decoy==FALSE],Reverse=T,rname=basename(r))
              
              
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
    
    if(ExtractMS1){
      EVENTS <- ms1scans$scanNumber#[which(mz>=lib_precursor[1]&mz<=lib_precursor[2])]
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
        
        XICs <- readXICs(r,masses = mzs,tol = ppm1)# unclear how this function works
        
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


## Extracts Intensities from a rawfile, input: rawfile, InclusionlistPath
Extract_Intensities_Raw <- function(RawFilePath,InclusionlistPath,outpath=RawFilePath,dirout = "PRM_Analyzer_Matches",
                                    parallelExecution=T,ppm1=10,ppm2=10,ChargeDetection="^",useDIA=T,session = NULL,
                                    threads = detectCores(),ScanStringFilterInvert=T,ScanStringFilter=NULL,ExpandingSpecificity=T){
  require(data.table)
  require(rawDiag)
  setwd(RawFilePath)
  #finding RawFiles

  raws <- list.files(pattern = ".raw")
  
  
  # Decoys:
  ST <- fread(paste(InclusionlistPath,"SpectraTable.txt",sep = "/"))
  TT <- fread(paste(InclusionlistPath,"TransitionTable.txt",sep = "/"))
  # ST<- ST[grep("IIGLDQVAGMSETALPGAFK",ST$Sequence)]
  TT <- TT[!is.na(match(TT$SpecID,ST$SpecID))]
  # raws <- paste("../",raws,sep ="/")
  
  TT <- TT[,.SD[!duplicated(Matches),],SpecID]
  
  STdec <- copy(ST)
  TTdec <- copy(TT)
  # STdec[,Sequence := paste(rev(unlist(strsplit(Sequence,""))),collapse = ""),Sequence]
  
  # Reversing Sequences
  DECOYS <- lapply(strsplit(gsub("#.*.","",ST$Sequence),""),rev)
  
  
  DECOYS <- lapply(DECOYS,function(x){
    x <- x
    Carbamido <- unlist(gregexpr("C",x))
    massshift <- 0
    if(all(Carbamido!=-1)){
      massshift <- 57.021464*length(Carbamido)
    }
    mwm <- mw(x,monoisotopic = T)
    mwm <- mwm+massshift
    mwm
  })
  
  
  names(DECOYS) <- ST$SpecID
  
  # Calculating Masses
  # ChargeDetection = "+"
  # ChargeDetection = "^"
  print("Generating Decoy Fragments.")
  
  if(length(session) != 0){
    Max <- TTdec[,1,.(SpecID,Matches)]
    progress <- Progress$new(session,min=0, max=dim(Max)[1])
    on.exit(progress$close())
    progress$set(message = 'Compiling Decoys',
                 detail = "This may take a while",value = 0)
  }
  
  TTdec[,Masses:={
    if(length(session)!=0){
      progress$set(value = .GRP)
    }
    
    cat("\r",.GRP)
    
    temp <- .SD
    GR <- .BY
    # stop()
    pep <- DECOYS[names(DECOYS)==GR$SpecID][[1]]
    pep <- (pep)
    charval <- length(pep)
    ion <- GR$Matches
    ion <- gsub("*","",ion,fixed = T)
    
    charge <- unlist(gregexpr(ChargeDetection,ion,fixed = T))#
    if(ChargeDetection == "+"){
      charge <- substr(ion,charge-2,charge-1)
      
    }
    if(ChargeDetection=="^"){
      charge <- substr(ion,charge+1,charge+1)
      
    }
    charge <- as.numeric(gsub("[\\(acbxyz]","",charge))
    if(is.na(charge)){
      charge <- 1
    }
    if(is.na(charge)){
      "Warning, no charge detected."
      charge=1
    }
    # charge[charge<1] <- 1
    
    position <- (strsplit(ion,"^[abcxyc]")[[1]][2])
    mods <- unlist(strsplit(position,"[-\\(]"))
    mods <- strsplitslot(mods,1,"^",fixed = T)
    position<- as.numeric(mods[1])
    
    additionalMass <- 0
    if(any(grepl("-h20",ion,fixed = T))){
      additionalMass <- additionalMass-18.010565
    }
    if(any(grepl("-NH3",ion,fixed = T))){
      additionalMass <- additionalMass-17.027
    }
    if(any(grepl("+h20",ion,fixed = T))){
      additionalMass <- additionalMass+18.010565
    }
    if(any(grepl("+NH3",ion,fixed = T))){
      additionalMass <- additionalMass+17.027
    }
    
    if(is.na(position)){
      warning("Warning, no position of fragment detected. Returning charval.")
      vec <- charval # double check, not sure what will happen here...
    }else{
      if(any(grepl("[abc]",ion))){
        vec <- 1:position
      }
      if(any(grepl("[xyc]",ion))){
        vec <- ((charval)-position+1):(charval)
      }
    }
    

    REVMW <- sum(pep[vec],na.rm = T)
    REVMW <- (REVMW+charge+additionalMass)/charge
    
    REVMW
    
  },.(SpecID,Matches)]
  if(length(session)!=0){
    progress$close()
    
  }
  
  # plot(TTdec$Masses==TT$Masses)
  # stop()
  #####
  # RawDiag Masses:
  ####
  raws <- paste(RawFilePath,basename(raws),sep = "/")
  if(length(threads)>length(raws)){
    threads <- length(raws)
  }
  if(length(session) != 0){
    Max <- TTdec[,1,.(SpecID,Matches)]
    progress <- Progress$new(session,min=0, max=length(raws))
    on.exit(progress$close())
    progress$set(message = 'Scanning Raw Files...',
                 detail = "This may take a while",value = 0)
  }
  
  if(0){
    # 
    # if(parallelExecution ==T){
    #   raws <<- raws
    #   dirout <<- dirout
    #   cl <- makeCluster(threads)
    #   # PH <- parSapply(cl,uniSeq,function(x){PeptideHydrophobicity(x)})
    #   s1 <-system.time(Cor <- parSapply(cl,raws,RawDiag_Vali_Extracter_DIA,
    #                                     dirout=dirout,outpath=outpath,ST=ST,TT=TT,TTdec=TTdec,ppm1=ppm1,ppm2=ppm2,
    #                                     SystemPath=SystemPath,session=NULL,DIA=useDIA))
    #   stopCluster(cl)
    # }else{
    #   save(raws,RawDiag_Vali_Extracter_DIA,dirout,outpath,ST,TT,TTdec,ppm1,ppm2,SystemPath,useDIA,file="Temp_ExtractorData.rda")
    #   
    #   sapply(raws,RawDiag_Vali_Extracter_DIA,dirout=dirout,outpath=outpath,ST=ST,TT=TT,TTdec=TTdec,ppm1=ppm1,ppm2=ppm2,
    #          SystemPath=SystemPath,DIA=useDIA)
    # }
    # 
  }else{
    CommandFile <- "#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
load(\"Temp_ExtractorData.rda\")
print(raws[as.numeric(args[1])])
try({Vali_Extracter_DIA(raws[as.numeric(args[1])],dirout=dirout,outpath=outpath,ST=ST,TT=TT,TTdec=TTdec,ppm1=ppm1,ppm2=ppm2,
                           SystemPath=SystemPath,DIA=useDIA,descriptionfilter =ScanStringFilter,InvertSearch=ScanStringFilterInvert)})"
    
    raws <- paste(RawFilePath,basename(raws),sep = "/")
    if(!exists("ScanStringFilterInvert")){
      ScanStringFilterInvert <- T
    }
    if(!exists("ScanStringFilter")){
      ScanStringFilter <- NULL
    }else{
      if(ScanStringFilter==""){
        ScanStringFilter <- NULL
      }
    }
    if(length(ScanStringFilter)==0){
      print("STOP NO ScanStringFilter")
    }
    # switch between methods
    Vali_Extracter_DIA <- Rawrr_Vali_Extracter_DIA
    
    save(raws,Vali_Extracter_DIA,dirout,outpath,ST,TT,TTdec,ppm1,ppm2,SystemPath,useDIA,ScanStringFilterInvert,ScanStringFilter,file="Temp_ExtractorData.rda")
    if(length(threads)>length(raws)){
      threads <- length(raws)
    }
    
    initvec <- 1:length(raws)
    write(CommandFile,file= "R_RawFileExtractor_init.R")
    progresslist <- list()
    if(length(session)>0){
      for(i in 1:threads){
        if(length(session) != 0){
          progresslist[[i]] <- Progress$new(session,min=0, max=100)
          progresslist[[i]]$set(detail = paste("Process",i),message="Initiating Scan")
          
          on.exit(progresslist[[i]]$close())
        }
      }
    }
    
    continue <- T
    while(continue){
      
      # Detect active and running session:
      StartedSessions   <- list.files(pattern="SessionStatusMS2",recursive = T)
      StartedSessions_MS1   <- list.files(pattern="SessionStatusMS1",recursive = T)
      
      FinishedSessions  <- list.files(pattern="^RawFile_Finished$",recursive = T)
      Le <- length(StartedSessions)-length(FinishedSessions)
      
      if(Le<threads&length(initvec)>0){
        system(paste("rscript R_RawFileExtractor_init.R",initvec[1]),wait = F,ignore.stdout = T)
        basename(raws[initvec[1]])
        StartedSessions2   <- basename(dirname(list.files(pattern="SessionStatusMS2",recursive = T)))
        while(all(gsub(".raw$","",basename(raws[initvec[1]]))!=StartedSessions2)){
          print("waiting for initiation")
          Sys.sleep(2)
          StartedSessions2   <- basename(dirname(list.files(pattern="SessionStatusMS2",recursive = T)))
        }
        initvec <- initvec[-1]
      }
      # check which rawfiles are active and running 
      DIFF <- is.na(match(basename(dirname(StartedSessions)),basename(dirname(FinishedSessions))))
      if(length(DIFF)>0){
        if(length(session)>0){
          sessioninfo <- strsplit(basename(StartedSessions[DIFF]),"_")
          if(length(StartedSessions_MS1)>0){
            sessioninfo_MS1 <- strsplit(basename(StartedSessions_MS1[DIFF]),"_")
          }
          
          if(length(sessioninfo)>0){
            for(i in 1:length(sessioninfo)){
              try({
                Pro <- progresslist[[i]]
                if(length(StartedSessions_MS1)>0){
                  
                  I <- sessioninfo_MS1[[i]]
                  Add <- "Processing MS1 Scans"
                }else{
                  I <- sessioninfo[[i]]
                  Add <- "Processing MS2 Scans"
                  
                  
                }
                na <- basename(dirname(StartedSessions))[i]
                nac <- nchar(na)
                if(nac>30){
                  na1 <- substr(na,1,30)
                  na2 <- substr(na,31,nac)
                  na <- paste(na1,na2,sep = " ")
                  
                }
                
                Pro$set(detail = na,message=Add,value = as.numeric(I[2])/as.numeric(I[3])*100)
              })
              
            }
            
          }
          
          progress$set(message = 'Scanning Raw Files...',
                       detail = "This may take a while",value = length(FinishedSessions))
          
        }
        
        
        
        
      }
      
      
      if(length(FinishedSessions)>=length(raws)){
        continue = F
      }else{
        Sys.sleep(1)
      }
      
    }
  }
  for(i in 1:length(progresslist)){
    try({progresslist[[i]]$close()})
    try({progress$close()})
  }
  
  if(length(session) != 0){
    progress$set(value = 1)
  }
  print("Finished Extraction")
}
## CompilingDatabase Main Function
CompilingDatabase <- function(mainPath,maxquant = NULL,session = NULL,test = F,SystemPath=SystemPath,threads=1){
  
  # mainPath <<- mainPath
  
  if(length(session) != 0){
    progress <- Progress$new(session,min=0, max=8)
    progress$set(message = 'compiling database',
                 detail = "",value = 0)
  }
  
  
  OutName <- paste(mainPath, "PRM_Analyzer_Matches",sep = "/")
  OutFol <- list.dirs(OutName)
  OutFol <- grep("/PRM_Analyzer_Matches$",OutFol,value = T,invert = T)
  
  For <- unique(unlist(sapply(OutFol,list.files,pattern = "^mz.*.[^-+].txt$",full.names = T)))
  For <- grep("random.txt$",For,invert = T,value = T)
  
  # ReverseHits
  Rev <- unique(unlist(sapply(OutFol,list.files,pattern = "^mz.*.Rev[+].txt$",full.names = T)))
  
  #Rev <- unique(unlist(sapply(OutFol,list.files,pattern = "^mz.*.RevMa.txt$",full.names = T)))
  decoyType = "reverse"
  if(decoyType == "random"){
    Rev <- grep("randomRev[+].txt$",Rev,invert = F,value = T)
  }else{
    Rev <- grep("randomRev[+].txt$",Rev,invert = T,value = T)
    
  }
  setwd(OutName)
  try(dbDisconnect(db),silent = T)
  if(!test){
    unlink("PickyAnalyzer.sqlite",force = T)
  }
  db <- dbConnect(SQLite(),dbname = "PickyAnalyzer.sqlite")
  tempdir <- OutFol 
  TempList <- unique(unlist(sapply(OutFol,list.files,full.name = T,pattern = "^mz")))
  LoopID <<- 0
  if(!test){
    
    try(dbRemoveTable(db,"msmsScans"),silent = T)
    dbSendQuery(conn=db,
                "CREATE TABLE msmsScans
                (
                Raw_File TEXT,
                Scan_Number INT,
                Retention_Time DOUBLE)
                ")
    try(dbRemoveTable(db,"ReverseTable"),silent = T)
    dbSendQuery(conn=db,
                "CREATE TABLE ReverseTable
                (
                RT_min DOUBLE,
                MatchCount INT,
                TransCount INT,
                SCAcut DOUBLE,
                SCAall DOUBLE,
                rawfile TEXT,
                precursor TEXT,
                mz DOUBLE)
                ")
    
    if(length(session) >0){
      progress$set(message = 'compiling database',
                   detail = "Forward Matches",value = 1)
    }
    #
    print("Building Database, Forward Sequences")
    lapply(list(For),function(TempList){
      LoopID <<- LoopID+1
      TempList <<- TempList
      bn <- basename(TempList)
      unib = unique(bn)
      if(length(session)>0){
        progress = Progress$new(session,min=0, max=length(unib))
        on.exit(progress$close())
      }
      
      sapply(1:length(unib),function(b,session,progress){
        
        # b <<- b
        if(length(session)>0){
          progress$set(message = b,value = b)
          
        }
        cat("\r",b,length(unib))
        templ = TempList[unib[b] == bn]
        Na = unique(basename(templ))
        # Combine All TAbles of a precursors
        Outc <- data.frame()
        Reference <- data.frame()
        Na = unique(sapply(templ,function(x){try(readLines(gsub("/mz","/",gsub("txt$","info",x))))}))
        Na <- Na[class(Na)!= "try-error"]
        
        try(Na <- paste("mz",paste(unlist(strsplit(Na,"#_#")),collapse = "_",sep = ""),sep = ""))
        
        for(i in templ){
          
          cl <- class(try(te <- read.csv(i,sep = "\t",stringsAsFactors = F,check.names = F,row.names = NULL),silent = F))
          # if(length(unique(te$rawfile))>22){te <<- te ;print(i);stop()}
          if(cl!= "try-error"){
            
            if(dim(te)[1] > 0){
              te$FDR <- 1
              te$PEP <- 1
              
              te$SCA <- 0
              te$SCAall <- gsub("[","",te$SCAall,fixed = T)
              te$SCAall <- as.numeric(gsub("]","",te$SCAall,fixed = T))
              te$SCAcut <- gsub("[","",te$SCAcut,fixed = T)
              te$SCAcut <- as.numeric(gsub("]","",te$SCAcut,fixed = T))
              Outc <- rbind(Outc,te)
              Reference <- rbind(Reference,te)
            }
            
            
          }
        }
        if(dim(Outc)[1] >0&0){
          # had to be switched off, no compatible with new RawDiag based R Script. COuld be reimplemented, if RawDiag Script includes the Reference information (Spectrum)
          
          Reference <- unique(Reference)
          
          RefOrder <- order(Reference[1,colnames(Reference)[1:match("charge",colnames(Reference))-1]],decreasing = T)
          RefOrder <- c(RefOrder,match("charge",colnames(Reference)):dim(Reference)[2])
          Reference <- Reference[,RefOrder]
          Outc <- Outc[,RefOrder]
        }
        names <- make.unique(make.names(names(Outc)))
        names(names) <- names
        names(Reference) <- names
        #write table into the database
        if(all(dbListTables(db) != Na,na.rm = T)&dim(Outc)[1] >0){
          
          # Initializing Table
          
          SqlTableInit <- GetColumnType(Outc)
          CH <- which(colnames(Outc) == "charge")
          Transitions = colnames(Outc)[1:(CH-1)]
          Rest = colnames(Outc)[(CH):dim(Outc)[2]] 
          
          ChaType = c("INT","DOUBLE","TEXT","INT","DOUBLE","INT","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE")
          ChaType <- ChaType[1:length(Rest)]
          SqlTableInit1 = paste(paste("'",Transitions,"'",sep = ""),"DOUBLE",collapse = ",")
          SqlTableInit2 = paste(paste("'",Rest,"'",sep = ""),ChaType,collapse = ",")
          SqlTableInit = paste(SqlTableInit1,SqlTableInit2,sep = ",")
          dbSendQuery(conn=db,
                      paste(
                        "CREATE TABLE",paste("'",Na,"'",sep = ""),"(",SqlTableInit,")
                        "
                      ))
        }
        if(any(dbListTables(db) == Na,na.rm = T)&dim(Outc)[1] >0){
          
          Outc$rawfile <- basename(Outc$rawfile)
          dbWriteTable(db, name = Na,value = Outc, row.names = FALSE,overwrite = F,append = T)
          # dbWriteTable(db, name = gsub("^mz","Reference",Na),value = Reference, row.names = FALSE,overwrite = F,append = T)
          
        }
        
        
        
        
        
        # if(!all(lengths(tereturn) == 1)){
        #   tereturn <- strsplit(unlist(tereturn),"\t")
        #   
        #   tereturn  <- tereturn[-grep("reference",tereturn,fixed = T)]
        #   
        #   tereturn <- as.data.frame(tereturn)
        #   tereturn <- t(tereturn)
        #   rownames(tereturn) <- NULL
        #   tereturn = tereturn[,match(c("rawfile","index","RT_min"),tereturn[1,])]
        #   colnames(tereturn) <- tereturn[1,]
        #   tereturn <- tereturn[-grep("rawfile",tereturn[,1]),]
        #   if(!is.vector(tereturn)){
        #     write.table(tereturn,msname,sep = "\t",quote = F,row.names = F,append = T)
        #   }else{
        #     write.table(t(as.matrix(tereturn)),msname,sep = "\t",quote = F,row.names = F,append = T,col.names = F)  
        #   }
        # }
        
        
      },session = session,progress = progress)
      
      
      
      
    })
    if(length(session) >0){
      
      progress$set(message = 'compiling database',
                   detail = "Reverse Matches",value = 2)
    }
    print("Building Database, Decoy Sequences")
    
    lapply(list(Rev),function(TempList){
      LoopID <<- LoopID+1
      TempList <<- TempList
      bn <- basename(TempList)
      unib = unique(bn)
      if(length(session)>0){
        progress = Progress$new(session,min=0, max=length(unib))
        on.exit(progress$close())
      }
      
      
      
      # dbListTables(db)
      # dim(te<- dbReadTable(db,"mz509.27182_2_SLHTLFGDK_P02768;P02768-3_ALB_668663.txt"))
      
      sapply(1:length(unib),function(b,session,progress){
        
        if(length(session)>0){
          # print("hump")
          progress$set(message = b,value = b)
          
        }
        cat("\r",b,length(unib))
        templ = TempList[unib[b] == bn]
        Na = unique(basename(templ))
        # Combine All TAbles of a precursors
        Outc <- data.frame()
        Reference <- data.frame()
        Na = unique(sapply(templ,function(x){try(readLines(gsub("/mz","/",gsub("Rev..txt$",".info",x))))}))
        Na <- Na[class(Na)!= "try-error"]
        
        try(Na <- paste("revmz",paste(unlist(strsplit(Na,"#_#")),collapse = "_",sep = ""),sep = ""))
        
        for(i in templ){
          # i <<- i
          # cl <- class(try(te <- fread(i,sep = "\t",stringsAsFactors = F,check.names = F,data.table = F),silent = T))
          # if(cl== "try-error"){
          cl <- class(try(te <- read.csv(i,sep = "\t",stringsAsFactors = F,check.names = F,row.names = NULL),silent = T))
          
          # }
          if(cl!= "try-error"){
            
            if(dim(te)[1] > 0){
              te$FDR <- 1
              te$PEP <- 1
              
              te$SCA <- 0
              te$SCAall <- gsub("[","",te$SCAall,fixed = T)
              te$SCAall <- as.numeric(gsub("]","",te$SCAall,fixed = T))
              te$SCAcut <- gsub("[","",te$SCAcut,fixed = T)
              te$SCAcut <- as.numeric(gsub("]","",te$SCAcut,fixed = T))
              Outc <- rbind(Outc,te[-1,])
              Reference <- rbind(Reference,te[1,])
            }
            
            
          }
        }
        if(dim(Outc)[1] >0&0){
          Reference <- unique(Reference)
          
          RefOrder <- order(Reference[1,colnames(Reference)[1:match("charge",colnames(Reference))-1]],decreasing = T)
          RefOrder <- c(RefOrder,match("charge",colnames(Reference)):dim(Reference)[2])
          Reference <- Reference[,RefOrder]
          Outc <- Outc[,RefOrder]
        }
        #write table into the database
        if(all(dbListTables(db) != Na,na.rm = T)&dim(Outc)[1] >0){
          
          # Initializing Table
          
          SqlTableInit <- GetColumnType(Outc)
          CH <- which(colnames(Outc) == "charge")
          Transitions = colnames(Outc)[1:(CH-1)]
          Rest = colnames(Outc)[(CH):dim(Outc)[2]] 
          
          ChaType = c("INT","DOUBLE","TEXT","INT","DOUBLE","INT","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE")
          ChaType <- ChaType[1:length(Rest)]
          
          SqlTableInit1 = paste(paste("'",Transitions,"'",sep = ""),"DOUBLE",collapse = ",")
          SqlTableInit2 = paste(paste("'",Rest,"'",sep = ""),ChaType,collapse = ",")
          SqlTableInit = paste(SqlTableInit1,SqlTableInit2,sep = ",")
          dbSendQuery(conn=db,
                      paste(
                        "CREATE TABLE",paste("'",Na,"'",sep = ""),"(",SqlTableInit,")
                        "
                      ))
        }
        if(any(dbListTables(db) == Na,na.rm = T)&dim(Outc)[1] >0){
          
          Outc$rawfile <- basename(Outc$rawfile)
          dbWriteTable(db, name = Na,value = Outc, row.names = FALSE,overwrite = F,append = T)
          #dbWriteTable(db, name = gsub("^revmz","revReference",Na),value = Reference, row.names = FALSE,overwrite = F,append = T)
          #  dbReadTable(db,Na)
        }
        
        
        
        
        
        # if(!all(lengths(tereturn) == 1)){
        #   tereturn <- strsplit(unlist(tereturn),"\t")
        #   
        #   tereturn  <- tereturn[-grep("reference",tereturn,fixed = T)]
        #   
        #   tereturn <- as.data.frame(tereturn)
        #   tereturn <- t(tereturn)
        #   rownames(tereturn) <- NULL
        #   tereturn = tereturn[,match(c("rawfile","index","RT_min"),tereturn[1,])]
        #   colnames(tereturn) <- tereturn[1,]
        #   tereturn <- tereturn[-grep("rawfile",tereturn[,1]),]
        #   if(!is.vector(tereturn)){
        #     write.table(tereturn,msname,sep = "\t",quote = F,row.names = F,append = T)
        #   }else{
        #     write.table(t(as.matrix(tereturn)),msname,sep = "\t",quote = F,row.names = F,append = T,col.names = F)  
        #   }
        # }
        
        
      },session = session,progress = progress)
      
      
      
      
    })
    
    
  }
  # grep("ISVLPAYLSYDAPWPVR",dbListTables(db),value = T)
  # temp <- dbReadTable(db,"mz974.0221_2_ISVLPAYLSYDAPWPVR#light_sp|Q10570_CPSF1_40087")
  if(length(session) >0){
    
    progress$set(message = 'FDR',
                 detail = "Scoring",value = 2)
  }
  # test <- dbReadTable(db,"mz530.770195_2_GVSINQFCK#heavy_sp|Q9CQF0_MRPL11_607482")
  # scoring function
  ScoringWrapper_parallelized(db,SystemPath = SystemPath,Parallelh20 = F,session = session,ReScore = F,threads = threads,progress = progress)#requires SystemPath
  # FDR function, returns the model and a Score cutoff at p = 0.01
  if(length(session) >0){
    
    progress$set(message = 'FDR',
                 detail = "Model",value = 2)
  }
  
  print("Estimating FDR")
  ScoreFDR_Loess_md <- ScoreFDR_Loess(db,set_amplifier = 1)
  
  DBtab <- dbListTables(db)
  DBtabFor <- grep("Rev+$",DBtab,value = T,invert = T,fixed = T)
  DBtabFor <- grep("^mz",DBtabFor,value = T)
  #sapply(DBtab,function(x){dim(dbReadTable(db,x))})
  
  compileMSMSSCANS <- function(x,tableName){
    x <- x
    TA <- dbReadTable(db,x)
    if(length(grep("rawfile",TA$charge,fixed = T)) > 0){
      TA  <- TA[-grep("rawfile",TA$charge,fixed = T),]
    }
    
    # TA  <- TA[-grep("rawfile",TA$charge,fixed = T),]
    
    TA = TA[,match(c("rawfile","index","RT_min"),colnames(TA))]
    TA$rawfile = basename(TA$rawfile)
    colnames(TA ) <- c("Raw_File","Scan_Number","Retention_Time")
    dbWriteTable(db, name = tableName,value = TA, row.names = FALSE,overwrite = F,append = T)
  }
  
  if(length(session) >0){
    
    progress$set(message = 'PRM FDR',
                 detail = "Calculating PEP",value = 5)
  }
  
  
  
  
  # Hydrophobicity Based Correction
  if(length(session) >0){
    
    try(progress$set(message = 'Retention Time Alignment',detail = "",value = 6))
  }
  if(length(maxquant) > 0){
    if(file.exists(maxquant)){
      
      try(Fits <- RT_Correction_HP(maxquant,doCompareRT = F))
      
      
    }
  }
  print("Finished RT alignment")
  # Adding Peps to table
  #PeaksFoundCheck
  if(length(session) >0){
    
    try(progress$set(message = 'Updating database',detail = "",value = 6))
  }
  
  if(length(session) > 0){
    progress2 <- Progress$new(session,min=0, max=length(DBtabFor))
    on.exit(progress$close())
  }else{
    progress2 <- NULL
  }
  # Checking TAbles...
  if(length(session) >0){
    
    progress$set(message = 'FDR',
                 detail = "Estimating FDR",value = 2)
  }
  for(it in 1:length(DBtabFor)){
    x <- DBtabFor[it]
    #mz526.74343_2_GDNIYEWR_P51965;P51965-2;P51965-3;Q969T4;Q96LR5_UBE2E1;UBE2E2;UBE2E3_698316 95;P10636-9_MAPT_657300 114
    if(length(session) > 0){
      try(progress2$set(value = it,message = x))
    }
    cat("\r",x,length(DBtabFor[match(x,DBtabFor):length(DBtabFor)]))
    #print(1)
    temp <- dbReadTable(db,x)
    temp$SCA_mz_fdr <- temp$FDR
    temp$SCAfdr <- temp$FDR
    
    
    
    
    #print(4)
    if(exists("Fits")){
      rf_Fits <- Fits$rf_Fits
      template_Fit <- Fits$Template
      
      print("CorrectingFile")
      tempx <- temp
      tempx$RT2 <- NA
      for(r in unique(tempx$rawfile)){
        cat("\r Correcting File",x,"RawFile",r )
        
        sel <- tempx$rawfile == r
        #names(rf_Fits)[names(rf_Fits) == "Beaker_20180919_MN_HS_IMR5_con_1_20180919150059"] <- "Beaker_20180919_MN_HS_IMR5_con_1"
        
        FitsSel <- rf_Fits[names(rf_Fits) == gsub(".raw","",r)]
        if(length(FitsSel) > 0){
          RTHP <- predict(FitsSel[[1]]$md,tempx$RT_min[sel])
          RT <- predict(template_Fit,RTHP)
          # RT <- sapply(tempx$RT_min[sel],function(x){
          #
          #   diff <- abs(x-as.numeric(FitsSel[,2]))
          #   CorrectedRT <- FitsSel[diff == min(diff,na.rm = T),1]
          #   return(CorrectedRT[1])
          # })
          # try(plot(tempx$RT_min[sel],RT))
          tempx$RT2[sel] <- RT
        }
      }
      temp <- tempx
    }
    
    
    dbWriteTable(db, name = x,value = temp, row.names = FALSE,overwrite = T,append = F)
  }
  if(length(session) > 0){
    try( progress2$close())
  }
  if(length(session) >0){
    
    try(progress$set(message = 'Calculating FDR',detail = "DBtabFor",value = 7))
  }
  print("Preparing FDRMatrix")
  
  # db<- dbConnect(SQLite(),"/PickyAnalyzer.sqlite")
  # 
  # dbReadTable(db,"PrecursorSelectionTable")
  Precursors <- grep("^mz",dbListTables(db),invert = F,value = T,fixed = F)
  Ratings <- grep("^RATING",dbListTables(db),value = T)
  Prec <- lapply(Precursors,function(x){
    xsvecdfcor <- NA
    try({
      # x <<- x
      # stop()
      cat("\r",x)
      xs <- dbReadTable(db,x)
      xs <- data.table(xs)
      if(length(xs$DL_Scores)==0){
        xs$DL_Scores <- -99
      }
      if(length(xs$RF_Scores)==0){
        xs$RF_Scores <- 0
      }
      if(length(xs$FDR_DL_all)==0){
        xs$FDR_DL_all <- 1
      }
      
      xsvec <-xs[,.(FDR=min(FDR_DL_all,na.rm = T),
                    DL_Scores=max(DL_Scores,na.rm = T),
                    RF_Scores=max(RF_Scores,na.rm = T),
                    Count=length(rawfile),
                    SCA=max(SCAall)),
                 rawfile]
      RatingsNames <- paste("RATING_PEAKS",gsub("^mz","",x),sep= "")
      sel <- RatingsNames==Ratings
      if(any(sel)){
        RATING <- dbReadTable(db,RatingsNames)
        RATING <- RATING$rating
        
      }else{RATING <- ""}
      xsvec$Rating <- RATING
      NAm<- names(xsvec)
      xsvecdf <- as.data.frame((xsvec))
      M<- match(NAm,c("FDR","DL_Scores","RF_Scores","Count","SCA","Rating","rawfile"))
      # Mna <- M
      # Mna[is.na(Mna)]<- 1
      xsvecdfcor <- xsvecdf[,M]
      xsvecdfcor$Dim <- dim(xs)[1]
    })
    
    # if(any(is.na(M))){
    #   xsvec[,is.na(M)] <- 0
    # }
    # TA <- table(xs$rawfile)
    return((xsvecdfcor))
  })
  names(Prec) <- Precursors
  errorcontrol <- sapply(Prec,function(x){is.na(unlist(x)[1])})
  Prec <- Prec[!errorcontrol]
  # Precdf <- as.data.frame(Prec)
  # Precdf <- t(Precdf)
  # Precdf$names <- rownames(Precdf)
  LEctrl <- table(sapply(Prec,length))
  LEctrl <- LEctrl[LEctrl==max(LEctrl)]
  
  # db <- dbConnect(SQLite(),"PickyAnalyzer.sqlite")
  # dbRemoveTable(db,"PrecursorSelectionTable")
  writtentable <- F
  for(i in 1:length(Prec)){
    tempi <- Prec[[i]]
    tempi <- data.frame(as.list(tempi))
    # names(tempi) <-c("rf","FDR","DL_Scores","RF_Scores","Count","SCA","Rating","Dim") 
    names(tempi) <-c("Rating","rf","FDR","DL_Scores","RF_Scores","Count","SCA","Dim")
    
    tempi$Species <- names(Prec)[i]
    # print(dim(tempi))
    if(dim(tempi)[2] == 9){
      cat("\r",i,"Writing PrecursorSelectionTable")
      
      if(!writtentable){
        # humpi <- dbReadTable(db,"PrecursorSelectionTable")
        writtentable <- dbWriteTable(db,"PrecursorSelectionTable",value = tempi,append = F,overwrite = T)
        
      }else{
        dbWriteTable(db,"PrecursorSelectionTable",value = tempi,append = T)
        
      }
    }else{
      if(length(session)>0){
        showNotification(session,"Problems, when writing PrecursorSelectionTable.")
        print("Problems")
      }else{warning("Problems, when writing PrecursorSelectionTable.")}
    }
    
  }
  # test <- dbReadTable(db,"PrecursorSelectionTable")
  #grep("^mz",dbListTables(db),value =T)det
  
  #pst <- dbReadTable(db,"PrecursorSelectionTable")
  
  try(write("",paste(OutName,"TableCompiled",sep = "/")))
  
  humpi <- try(dbDisconnect(db))
  print(humpi)
  rm(db)
  #file.rename(paste(OutName,"PickyAnalyzer_compiling.sqlite",sep = "/"),paste(OutName,"PickyAnalyzer.sqlite",sep = "/"))
  #print(class(humpis))
  return(NULL)
}
### Column Type/class function, used in compilingDatabase
GetColumnType <- function(Outc){
  CH <- which(colnames(Outc) == "charge")
  Transitions = colnames(Outc)[1:(CH-1)]
  Rest = colnames(Outc)[(CH):dim(Outc)[2]] 
  
  ChaType = c("INT","DOUBLE","TEXT","INT","DOUBLE","INT","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE","DOUBLE")
  SqlTableInit1 = paste(paste("'",Transitions,"'",sep = ""),"DOUBLE",collapse = ",")
  SqlTableInit2 = paste(paste("'",Rest,"'",sep = ""),ChaType,collapse = ",")
  SqlTableInit = paste(SqlTableInit1,SqlTableInit2,sep = ",")
  return(SqlTableInit)
}
### Called in FDR_Calculation
PoissonGLM <- function(x,breaks=30){
  hi <- hist(x,breaks =breaks,plot = F)
  quans <- (hi$breaks[-1]-diff(hi$breaks)/2)
  da <- data.frame(counts = hi$counts,quans = quans)
  gl <- glm(counts ~ quans,data = da,family = quasipoisson(link = "log"))
  # lines(quans,gl$fitted.values,col = 2)
  return(list(glm=gl,x=quans,y=gl$fitted.values,breaks = breaks))
}
### FDR Assessment
FDR_Calculation <- function(forScore,revScore,plotdata = T,byvec=0.1,amplifier = 100,AddMaximalValue = T,usePoissonModel = F,sparvec = 0.1){
  # forScore<<- forScore
  # revScore <<- revScore
  
  if(usePoissonModel){
    
    
    fwd.glm <- PoissonGLM(forScore)
    rev.glm <- PoissonGLM(revScore)
    
    PREDIfwd <- predict(fwd.glm$glm,data.frame(quans = seq(0,1,by = 0.01)),type = "response")*amplifier
    PREDIrev <- predict(rev.glm$glm,data.frame(quans = seq(0,1,by = 0.01)),type = "response")*amplifier
    
    init <- seq(0,1,by = 0.01)
    fwd <- unlist(lapply(1:length(init),function(x){jitter(rep(init[x],round(PREDIfwd[x])))}))
    rev <- unlist(lapply(1:length(init),function(x){jitter(rep(init[x],round(PREDIrev[x])))}))
  }else{
    fwd <-forScore
    rev <- revScore
  }
  fwd <- fwd[!fwd <= 0]
  rev <- rev[!rev<=0]
  names(fwd) <-rep("fwd",length(fwd))
  names(rev) <-rep("rev",length(rev))
  VEC <- c(fwd,rev)
  VEC <- sort(VEC,decreasing = T)
  FDRw <- which(names(VEC)=="rev")
  SCORE <- VEC[FDRw]
  FDR <- c(1:length(FDRw))/FDRw
  
  # pdf("PseudoRocCurve.pdf")
  # plot(SCORE,FDR,type = "l")
  # def.off()
  # abline(Ch=0.01)
  
  
  
  if(plotdata&usePoissonModel){
    try({
      
      hist(forScore,breaks = fwd.glm$breaks,col = 1,xlab = "Score",main = c("Poisson GLM"),xlim =c(0,1))
      hist(revScore,breaks = rev.glm$breaks,add = T,border=2,xlim =c(0,1))
      
      lines(fwd.glm,col = 2)
      lines(rev.glm,col = 3)
    })
  }
  if(max(SCORE,na.rm =T )<3&AddMaximalValue){
    
    seqvec <- seq(max(SCORE),3,by = 0.1)
    seqvec <- rep(seqvec,5)
    SCORE <- c(seqvec,SCORE)
    FDR <- c(rep(0.001,length(seqvec)),FDR)
  }
  
  # smoothScatter(SCORE ,FDR,xlim = c(0,1),ylim= c(0,0.2),xlab = "Score",ylab = "FDR",frame = F,pch = 20)
  # points(AverageScore$Score,AverageScore$FDR,col = 2,type="l",pch = 20)
  # points(FDRMAX$x,FDRMAX$y,type="l",col="blue")
  
  
  
  # correctedFDR <- fdr_MAX(SCORE,FDR)
  
  # Subset Samples, keeping extremes
  # plot(correctedFDR$x,correctedFDR$y,type = "l")
  # or <-order(correctedFDR$x)
  # if(length(or)>10000){
  #   cl<- class(correctedFDR)
  #   
  #   correctedFDR <- lapply(correctedFDR,function(x){x[or]})
  #   
  #  
  #   
  #   endings<-c(1:1000,c((length(or)-1000):length(or)))
  #   sampling <- sample(1:length(or),10000)
  #   sampling <- unique(c(endings,sampling))
  #   correctedFDR <- lapply(correctedFDR,function(x){x[sort(sampling)]})
  #   class(correctedFDR) <- cl
  # }
  
  # sspl <- list(x=SCORE,y = FDR)
  # try({sspl<- smooth.spline(SCORE,FDR)})
  # FDR[20] <- 0.2
  # lines(sspl)
  
  # fdr_MAX(SCORE,FDR)
  # plot(sspl)
  # loess(sspl$)
  
  # plot(sspl,type = "l",xlab = "SCORE",ylab = "FDR")
  
  ## nls smooth 1
  
  # a_start<-8 #param a is the y value when x=0
  # b_start<-2*log(2)/a_start #b is the decay rate
  # nls_result <- nls(FDR~a*exp(-b*SCORE-3.3),start = list(a=a_start,b=b_start))
  # # curve(9*exp(-1.9*x-3.3))
  # pr <- predict(nls_result,data.frame(SCORE=seq(0,1,by = 0.1)),type = "response")
  # if(coefficients(nls_result)[1]<0){
  #   FDRQualityCheck <- F
  # }else{
  #   FDRQualityCheck <- T
  #   
  # }
  # byvec <- 0.1
  AverageScore<- sapply(seq(0,3,by=byvec),function(x){
    # x <<- x
    cat("\r AVERAGE FDR",x)
    c(x+byvec/2,median(jitter(FDR)[SCORE>x&SCORE<(x+byvec)],na.rm = T))
  })
  AverageScore <- data.frame(Score=AverageScore[1,],FDR=AverageScore[2,])
  AverageScore <- AverageScore[!is.na(AverageScore$FDR),]
  correctedFDR <- fdr_MAX(AverageScore$Score,AverageScore$FDR)
  
  # smoothScatter(SCORE ,FDR,xlim = c(0,1),ylim= c(0,0.2),xlab = "Score",ylab = "FDR",frame = F,pch = 20)
  
  
  nls_result <- NULL
  if(plotdata){
    try({
      smoothScatter(SCORE ,FDR,xlim = c(0,1),ylim= c(0,0.2),xlab = "Score",ylab = "FDR",frame = F,pch = 20)
      grid()
      
      abline(h=0.01,lty = "dotted")
      abline(h=0.05,lty = "dotted")
      # lines(sspl,col = 2)
      lines(correctedFDR,col = 3)
      lines(correctedFDR$x,correctedFDR$yold,col = 2)
      # lines(seq(0,1,by = 0.1),pr,col = 3)
      legend("topright",legend=c("average FDR","conservative average FDR"),col = c(2,3),lwd = 1,bty = "n")
      
      # if(!FDRQualityCheck){
      #   legend("top",legend = "Warning, curve is not in line with general assumptions.\n
      #          FDR is corrupted and will not be applied.",text.col = 2,bg ="#FFFFFF99",box.col = NA)
      # }
    })
    
  }
  
  
  
  
  return(list(nls = nls_result,maxFDR = correctedFDR,Score=SCORE,FDR=FDR))
  
  
  
  
  
}
ScoreFDR_Loess <- function(db,Scores = "RF_Scores",set_amplifier=10){
  graphics.off()
  FDR_MD <- NULL
  FDR_MD_All <- NULL
  try({
    pdf("ValiScoreDistributions.pdf",width = 9)
    # random Forest ALL 
    sqliteextractFUN <- function(x){try(dbGetQuery(db,statement = paste(paste("SELECT DL_Scores FROM '",sep = ""),x,"' WHERE MatchCount >= 0",sep = "")),silent = T)}
    forScoreAll <- lapply(grep("^mz",dbListTables(db),value = T),sqliteextractFUN)
    forScoreAll <- forScoreAll[sapply(forScoreAll,class) == "data.frame"]
    revScoreAll <- lapply(grep("^revmz",dbListTables(db),value = T),sqliteextractFUN)
    revScoreAll <- revScoreAll[sapply(revScoreAll,class) == "data.frame"]
    
    FS<- unlist(forScoreAll)
    RS <- unlist(revScoreAll)
    
    if(length(RS)<100000){set_amplifier = set_amplifier}
    if(length(RS)>100000){set_amplifier = 1}
    
    
    
    par(mfrow = c(2,1))
    FDR_RF_All <- FDR_Calculation(forScore = FS,revScore = RS,amplifier = set_amplifier,plotdata = T,byvec=0.01,AddMaximalValue = T)
    # DEEP LEarning all
    sqliteextractFUN <- function(x){try(dbGetQuery(db,statement = paste(paste("SELECT ","DL_scores"," FROM '",sep = ""),x,"' WHERE MatchCount > 0",sep = "")),silent = T)}
    forScoreAll <- lapply(grep("^mz",dbListTables(db),value = T),sqliteextractFUN)
    forScoreAll <- forScoreAll[sapply(forScoreAll,class) == "data.frame"]
    revScoreAll <- lapply(grep("^revmz",dbListTables(db),value = T),sqliteextractFUN)
    revScoreAll <- revScoreAll[sapply(revScoreAll,class) == "data.frame"]
    
    
    
    FS<- unlist(forScoreAll)
    RS <- unlist(revScoreAll)
    
    FDR_DL_All <- FDR_Calculation(forScore = FS,revScore = RS,amplifier = set_amplifier,plotdata = T,byvec = 0.05)
    # # GLM  all
    # sqliteextractFUN <- function(x){try(dbGetQuery(db,statement = paste(paste("SELECT ","GLM_scores"," FROM '",sep = ""),x,"' WHERE MatchCount > 0",sep = "")),silent = T)}
    # forScoreAll <- lapply(grep("^mz",dbListTables(db),value = T),sqliteextractFUN)
    # forScoreAll <- forScoreAll[sapply(forScoreAll,class) == "data.frame"]
    # revScoreAll <- lapply(grep("^revmz",dbListTables(db),value = T),sqliteextractFUN)
    # revScoreAll <- revScoreAll[sapply(revScoreAll,class) == "data.frame"]
    # 
    # FS<- unlist(forScoreAll)
    # RS <- unlist(revScoreAll)
    # 
    # par(mfrow = c(2,1))
    # FDR_GLM_All <- FDR_Calculation(forScore = FS,revScore = RS,amplifier = set_amplifier,plotdata = T)
    # # GBM  all
    # 
    # sqliteextractFUN <- function(x){try(dbGetQuery(db,statement = paste(paste("SELECT ","GBM_scores"," FROM '",sep = ""),x,"' WHERE MatchCount > 0",sep = "")))}
    # forScoreAll <- lapply(grep("^mz",dbListTables(db),value = T),sqliteextractFUN)
    # forScoreAll <- forScoreAll[sapply(forScoreAll,class) == "data.frame"]
    # revScoreAll <- lapply(grep("^revmz",dbListTables(db),value = T),sqliteextractFUN)
    # revScoreAll <- revScoreAll[sapply(revScoreAll,class) == "data.frame"]
    # 
    # FS<- unlist(forScoreAll)
    # RS <- unlist(revScoreAll)
    # 
    # par(mfrow = c(2,1))
    # FDR_GBM_All <- FDR_Calculation(forScore = FS,revScore = RS,amplifier = set_amplifier,plotdata = T)
    # 
    
    
    val <- grep("^mz",dbListTables(db),value = T)
    
    FDR_MD <-   lapply(grep("^mz",dbListTables(db),value = T),function(i){
      cat("\r",i)
      # i <<- i
      par(mfrow= c(1,2))
      mzt <- dbReadTable(db,i)
      # mztrev <- NULL
      # try(mztrev <- dbReadTable(db,gsub("^mz","revmz",i)))
      
      # hist(mzt$RF_Scores,breaks = 100)
      # hist(mztrev$RF_Scores,add = T,border = 2)
      # FDRFU <- FDR_Calculation(forScore = mzt$RF_Scores,revScore = mztrev$RF_Scores,sparvec = 0.1,AddMaximalValue = T)
      FDR <- 1
      FDR_RF_spec<- 1
      FDR_RF_all <- 1
      FDR_DL_all <- 1
      FDR_GLM_all <- 1
      FDR_GBM_all <- 1
      # try({        FDR.spl <- predict(FDRFU$sspl,data.frame(SCORE=mzt$RF_Scores))$y$SCORE})
      # try({        FDR_RF_all <- predict(FDR_RF_All$sspl,data.frame(SCORE=mzt$RF_Scores))$y$SCORE})
      # try({        FDR_DL_all <- predict(FDR_DL_All$sspl,data.frame(SCORE=mzt$DL_Scores))$y$SCORE})
      # try({        FDR_GLM_all <- predict(FDR_GLM_All$sspl,data.frame(SCORE=mzt$GLM_Scores))$y$SCORE})
      # try({        FDR_GBM_all <- predict(FDR_GBM_All$sspl,data.frame(SCORE=mzt$GBM_Scores))$y$SCORE})
      # try({        FDR_RF_spec <- predict(FDRFU$maxFDR,mzt$RF_Scores)})
      
      try({        FDR_RF_all <- predict(FDR_RF_All$maxFDR,mzt$RF_Scores)},silent = T)
      # try({        FDR_DL_all <- predict(FDR_DL_All$maxFDR,mzt$DL_Scores)})
      try({        FDR_DL_all <- predict(FDR_DL_All$maxFDR,mzt$DL_Scores)},silent = T)
      # try({        FDR_GBM_all <- predict(FDR_GBM_All$maxFDR,mzt$GBM_Scores)})
      
      # mzt$FDR_RFspl_peptide <- FDR.spl
      # mzt$FDR_RF_peptide <- FDR_RF_spec
      
      mzt$FDR_RF_all <- FDR_RF_all
      mzt$FDR_DL_all <- FDR_DL_all
      # mzt$FDR_GLM_all <- FDR_GLM_all
      # mzt$FDR_GBM_all <- FDR_GBM_all
      
      dbWriteTable(db,i,mzt,overwrite = T)
      NULL
      #FDRFU
    })
    
  })
  dev.off()
  save(FDR_DL_All,FDR_RF_All,file="FDR_Models.rda")
  
  return(list(FDR_MD,FDR_MD_All))
}

### Predictor of FDR based on Scores
predict.FDRCORRECTED <- function(object,predictionVector,digits = 5,negativeToZero = T){
  # predictionVector <<- predictionVector
  predictionVector <- round(predictionVector,digits)
  unipred <- unique(predictionVector)
  # object <<- object
  unipredResult <- sapply(unipred,function(x){
    # get boundaries
    # x <<- x
    upper <- min(object$x[object$x>=x],na.rm = T)
    
    lower <- max(object$x[object$x<=x],na.rm = T)
    Yval  <- object$y[match(c(upper,lower),object$x)]
    Yval <- Yval[!is.na(Yval)]
    if(length(Yval)>1&upper!=lower){
      xs <- c(upper,lower)
      
      # lm.md <- lm(Yval~xs)
      # pr <- predict(lm.md,data.frame(xs=x))
      
      m <- (Yval[1]-Yval[2])/(upper-lower)
      n <- Yval[1]-m*upper
      pr = x*m+n
    }else{pr =max(Yval) }
    if(negativeToZero){
      pr[pr<0] <- 0
    }
    if(is.na(pr)){
      stop("is.na(pr)")
    }
    pr
  })
  unipredResult <- unipredResult[match(predictionVector,unipred)]
  unipredResult
}
### FDR modifier, keeps always higher value: conservative solution
fdr_MAX <- function(SCORE,FDR){
  ord <- order(SCORE,decreasing = T)
  SCORE <- SCORE[ord]
  FDR <- FDR[ord]
  init_maxFDR <<- 0
  CorrectedFDR <- sapply(FDR,function(x){
    if(x>init_maxFDR){
      init_maxFDR <<- x
    }
    init_maxFDR
  })
  hum <- list(y= CorrectedFDR,yold=FDR,x = SCORE)
  class(hum) <- list("list","FDRCORRECTED")
  return(hum)
}
### Scoringwrapper, main function, used in (compilingdatabase
ScoringWrapper_parallelized <- function(db,mdModel = NULL,SystemPath = NULL,Parallelh20=F,session = NULL,progress=NULL,ReScore=F,threads=3){
  FUN <-try({
    # mdpath <- list.files(paste(SystemPath,"Model",sep = "/"),pattern = "ScoringModel_Full.rda$",full.names =T)
    mdpath <- list.files(paste(SystemPath,"Model",sep = "/"),pattern = "model_",full.names =T)
    mdpath <- list.files(paste(SystemPath,"Model/MojoBackup",sep = "/"),pattern = "model_",full.names =T)
    if(length(mdpath)==0){
      stop("NO MODELS")
    }
    try({if(length(session)>0&length(progress)>0){
      print("Progress Test")
      progress$set(message = 'FDR',
                   detail = "Loading FDR model",value = 2)
    }})
    
    if(require(h2o)){
      h2o.init(nthreads = threads, max_mem_size = '2g', ip = "127.0.0.1", port = 4321)
      # 
      # # localH2O <- h2o.init(nthreads = threads)
      
      # Models.h2o <- lapply(mdpath,function(x){
      #   try(h2o.loadModel(x))
      # })
      Models.h2o <- lapply(mdpath,function(x){
        try(h2o.import_mojo(x))
      })
      names(Models.h2o) <- strsplitslot(basename(mdpath),1,"_")
      Models.h2o <- Models.h2o[sapply(Models.h2o,class)!="try-error"]
      
    }else{
      warning("No package h20 available. Install and run again. Aborting analysis.")
      stop()
    }
    
    
    fi <- grep("^mz|^revmz",dbListTables(db),value = T)
    
    
    graphics.off()
    
    # if(length(session)>0){
    #   progressScoring <- Progress$new(session,min=0, max=length(fi))
    #   on.exit(progressScoring$close())
    # }
    
    ###### new insert ########
    
    print("Extracting Data for Scoreprediction")
    ### Parallel write out of tables ############
    
    dir.create("TempTables")
    if(basename(getwd())=="TempTables"){
      unlink(list.files())
    }
    # cl <- makeCluster(threads,outfile="test.txt",setup_strategy = "sequential")
    cl <- makeCluster(threads,outfile="test.txt",setup_strategy = "sequential")
    
    print("Exporting Models.h2o")
    
    print("Done Exporting Models")
    unlink("./Temp_CheckFile.txt")
    unlink("./REV_Temp_CheckFile.txt")
    Files <- NULL
    RevFiles <- NULL
    cat("\rStarting TempCheckfile")
    try({if(length(session)>0&length(progress)>0){
      progress$set(message = 'FDR',
                   detail = "Generating Featuretable (this will take a while)",value = 2)
    }})
    
    try({
      dbNAME <- db@dbname
      wd <- getwd()
      if(file.exists("Models.h2o")){
        clusterExport(cl, "Models.h2o",envir ="ScoringWrapper_parallelized" )
        
      }
      
      hu <- parLapply(cl,1:length(fi),function(itx,fi,dbNAME=NULL,session=NULL,ReScore=T,SystemPath=NULL,Models.h2o=NULL,db=NULL,wd = NULL){
        library(RSQLite)
        library(data.table)
        library(h2o)
        h2o.init()
        source(paste(SystemPath,"R/EvaluationScript_PRM_sqlite.R",sep = "/"))
        setwd(wd)
        cat("Working on",itx)
        itx <- itx
        mzname <- fi[[itx]]
        # cat("\rExtracting Tables:",round(itx/length(fi)*100,1),"%     ")
        # if(length(session)>0){
        #   progressScoring$set(value = itx)
        # }
        Reverse <- grepl("^revmz",mzname,perl=T)
        
        db2 <- dbConnect(SQLite(),dbNAME)
        mztable <- dbReadTable(db2,mzname)
        dbDisconnect(db2)
        
        mztable <- data.table(mztable)
        
        SEQUENCE <-strsplitslot(mzname,3,"_")
        SEQUENCE <- strsplitslot(SEQUENCE,1,"#",fixed = T)
        SEQUENCE <- strsplitslot(SEQUENCE,1,"#",fixed = T)
        SEQUENCE <- gsub("\\(.*?\\)\\)","",SEQUENCE,perl=T)
        SEQUENCE <- gsub("^_","",SEQUENCE)
        SEQUENCE <- gsub("_.*$","",SEQUENCE)
        mztable$Peptide_Length <- nchar(SEQUENCE)
        
        if(length(unique(mztable$Score))<2|ReScore){
          cat("\rSTARTING RESCORING")
          
          # mztable <- mztable[rawfile!="",]
         
          
          tr <- try({
            mztableSelect <- MakeFeatureTable(mztable,defaultPeptideLength =nchar(SEQUENCE),mzname = mzname,mdpar=Models.h2o$DeepLearning@parameters$x)
          })
          if(class(tr)=="try-error"){
            print("error in MakeFeatureTable. Check variables mztable and mzname in function MakeFeatureTable().")
            # mztable <<- mztable
            # mzname <<- mzname
            stop(tr)
          }
          if(Reverse){
            fwrite(mztableSelect,file = paste("./TempTables/REV_Temp_CheckFile",itx,".txt"),sep = "\t",append = T)
          }else{
            fwrite(mztableSelect,file = paste("./TempTables/Temp_CheckFile",itx,".txt"),sep = "\t",append = T)
          }
        }else{write(x = "Skipping",file = paste("Humpe",as.character(itx),".txt"))}
        
        NULL
      },
      fi=fi,SystemPath=SystemPath,Models.h2o=Models.h2o,dbNAME=dbNAME,wd = wd)
      
      
      Files <- lapply(list.files("./TempTables",pattern="^Temp.CheckFile",full.names = T),fread,sep = "\t",stringsAsFactor=F,integer64 = "double")
      RevFiles <- lapply(list.files("./TempTables",pattern="^REV_Temp.CheckFile",full.names = T),fread,sep = "\t",stringsAsFactor=F,integer64 = "double")
      Files <- rbindlist(Files)
      RevFiles <- rbindlist(RevFiles)
    })
    print("Finished TempCheckfile")
    # setwd("../")
    stopCluster(cl)
    fwrite(RevFiles,file = "REV_Temp_CheckFile.txt",sep = "\t",append = T)
    fwrite(Files,file = "Temp_CheckFile.txt",sep = "\t",append = T)
    
    # Predictin Socres ###########
    
    try({if(length(session)>0&length(progress)>0){
      progress$set(message = 'FDR',
                   detail = "Predicting Scores",value = 2)
    }})
    
    print("Predicting Scores and writing to Database")
    scorelist <- list()
    for(txt in c("Temp_CheckFile.txt","REV_Temp_CheckFile.txt")){
      print(paste("Predicting Scores",txt))
      mztableSelect <- fread(txt,sep = "\t",stringsAsFactors = F)
      if(is.vector(mztableSelect)){
        na <- names(mztableSelect)
        mztableSelect <- t(as.matrix(mztableSelect))
      }
      mztable.h2o <- as.h2o(mztableSelect)
      
      #   
      s2 <-system.time(
        PredictionScores <- sapply(Models.h2o,function(x){
          PREDI <- as.data.frame(h2o.predict(x,mztable.h2o))
          PREDI$predict
        })
      )
      # smoothScatter(mztableSelect$SCAall,PredictionScores$DeepLearning,cex = 1,pch = 20)
      PredictionScores <- data.frame(PredictionScores)
      mztableSelect$RF_Scores <- PredictionScores$DRF
      mztableSelect$DL_Scores <- PredictionScores$DeepLearning
      # smoothScatter(mztableSelect$SCAall,mztableSelect$DL_Scores,cex = 1,pch = 20)
      # smoothScatter(mztableSelect$SCAall,mztableSelect$RF_Scores,cex = 1,pch = 20)
      # smoothScatter(mztableSelect$RF_Scores,mztableSelect$RF_Scores,cex = 1,pch = 20)
      # 
      mztableSelect[,{
        cat("\r Writing",mzname)
        gr<- .BY
        da <- .SD
        temp <- dbReadTable(db,gr$mzname)
        temp <- as.data.table(temp)
        coln <- colnames(temp)
        temp[,c("RF_Scores","DL_Scores"):= {
          tempda <- da[da$rawfile==.BY$rawfile]# rawfile filtering of da
          if(any(duplicated(tempda$scan))){
            warning("Duplicated Scans Found in tempda")
          }
          MatchFun <- match(scan,tempda$scan)
          list(
            tempda$RF_Scores[MatchFun],
            tempda$DL_Scores[MatchFun]
          )
          
        },rawfile]
        
        temp$RF_Scores[is.na(temp$RF_Scores)]<- 0
        temp$DL_Scores[is.na(temp$DL_Scores)]<- 0
        setcolorder(temp,coln)
        # temp <- temp
        dbWriteTable(db,mzname,temp,overwrite = T)
        NULL
      },mzname]
      scorelist <- append(scorelist,list(PredictionScores))
      
      # unlink(txt)
    }
  })
  
  pdf("Scores_density.pdf",width=10,height=5)
  par(mfrow=c(1,2))
  try({
    plot(density(scorelist[[1]]$DeepLearning),main="DL Scores Density",col = "red")
    points(density(scorelist[[2]]$DeepLearning),border = 2,type="l",col = "blue")
    legend("topright",legend=c("fw","rev"),col = c(2,"blue"),lwd = 2)
    
    plot(density(scorelist[[1]]$DRF),main="RF Scores Density",col = "red")
    points(density(scorelist[[2]]$DRF),border = 2,type="l",col = "blue")
    legend("topright",legend=c("fw","rev"),col = c(2,"blue"),lwd = 2)
  })
  dev.off()
  return(FUN)
}
# ScoringWrapper_parallelized(db,SystemPath = SystemPath,Parallelh20 = F,session,ReScore = F,threads = threads)#requires SystemPath

ScoringWrapper <- function(db,mdModel = NULL,SystemPath = NULL,Parallelh20=F,session = NULL,ReScore=F,threads=3){
  FUN <-try({
    # mdpath <- list.files(paste(SystemPath,"Model",sep = "/"),pattern = "ScoringModel_Full.rda$",full.names =T)
    mdpath <- list.files(paste(SystemPath,"Model",sep = "/"),pattern = "model_",full.names =T)
    
    if(require(h2o)){
      h2o.init(nthreads = -1, max_mem_size = '2g', ip = "127.0.0.1", port = 4321)
      
      # localH2O <- h2o.init(nthreads = threads)
      
      # h2o.loadModel(mdpath[[1]])
      Models.h2o <- lapply(mdpath,h2o.loadModel)
      names(Models.h2o) <- strsplitslot(basename(mdpath),1,"_")
      
    }else{
      warning("No package h20 available. Install and run again. Aborting analysis.")
      stop()
    }
    
    # if(class(rf) == "randomForest"){
    #   Features <- rownames(rf$importance)
    #   
    # }else{
    #   rf<- sv
    #   Features <- colnames(sv$SV)
    # }
    
    fi <- grep("^mz|^revmz",dbListTables(db),value = T)
    
    
    graphics.off()
    
    if(length(session)>0){
      progressScoring <- Progress$new(session,min=0, max=length(fi))
      on.exit(progressScoring$close())
    }
    
    pdf("Scoring_plots.pdf")
    sapply(1:length(fi),function(itx){
      # itx <<- itx
      mzname <- fi[[itx]]
      print(mzname)
      if(length(session)>0){
        progressScoring$set(value = itx)
      }
      #mz569.28904_2_HMQNSEIIR#heavy_Q8N3U4;Q8N3U4-2_STAG2_5385025
      plot(mztable$SCAall,mztable$DL_Scores)
      
      mztable <- dbReadTable(db,mzname)
      mztable <- data.table(mztable)
      SEQUENCE <- strsplitslot(mzname,3,"_")
      SEQUENCE <- strsplitslot(SEQUENCE,1,"#",fixed = T)
      SEQUENCE <- gsub("\\(.*?\\)\\)","",SEQUENCE,perl=T)
      SEQUENCE <- gsub("^_","",SEQUENCE)
      SEQUENCE <- gsub("_.*$","",SEQUENCE)
      mztable$Peptide_Length <- nchar(SEQUENCE)
      if(length(unique(mztable$Score))<2|ReScore){
        
        
        mztable <- mztable[rawfile!="",]
        
        mztable <- mztable[,ValiScoringFunction(.SD,.BY),rawfile]
        
        
        
        mztable <- mztable[order(RT_min),]
        mztable <- data.table(mztable)
        mztable <- mztable[order(RT_min),]
        # A2 <- NA
        # A1 <- NA
        if(0){
          PeakFinderError <- try({
            
            try(RTPeaksOri <- PeaksFoundCheck(mztable))# FIND PEAKS WITH originial RT
            try(A1 <- PeaksAligner(mztable,RTPeaksOri))
            # temp2 <- temp
            try(RTPeaksCor <- PeaksFoundCheck(mztable))# FIND PEAKS WITH new RT b-ased on corrected RTs
            if(exists("Fits")){
              try(A2 <- PeaksAligner(mztable,RTPeaksOri))
              
            }else{
              A2 <- A1
            }
            
            
          },silent = T)
          if(class(PeakFinderError) == "try-error"){
            print("WARNING, Problem with PeaksAligner Segment")
            print(PeakFinderError)
          }
          
        }
        
        
        
        
        # mztable$Peaks1 <- A1
        # mztable$Peaks2 <- A2
        
        #### h2o prediction
        TrainSet <- mztable
        # TrainSet$ScaAfter <- NULL
        # TrainSet$ScaAfter2 <- NULL
        # TrainSet$ScaBefore <- NULL
        # TrainSet$ScaBefore2 <- NULL
        # 
        
        TrainSet[,ScaBefore:= {c(SCAall[-1],0)},.(rawfile)]
        TrainSet[,ScaBefore2:= {
          if(length(SCAall)>3){
            ret <-c(SCAall[-c(1:2)],0,0)
          }else{
            ret <- 0
          }
          ret
          
        },.(rawfile)] 
        
        TrainSet[,ScaAfter:= c(0,SCAall[-length(SCAall)]),.(rawfile)]
        TrainSet[,ScaAfter2:= 
                   {
                     if(length(SCAall)>3){
                       ret <- c(0,0,SCAall[-c(length(SCAall):(length(SCAall)-1))])
                     }else{
                       ret <- 0
                     }
                     ret
                     
                   },.(rawfile)]
        TrainSet$MLOGP <- -log10(TrainSet$R_p)
        mdpar <- Models.h2o$DeepLearning@parameters$x
        mdpar_match <- match(mdpar,colnames(TrainSet))
        missing <- mdpar[is.na(mdpar_match)]
        mztableSelect <- TrainSet[,.SD,.SDcols = mdpar_match]
        mztableSelect <- apply(mztableSelect,2,function(x){x[is.na(x)] <- 0;as.numeric(x)})
        
        if(is.vector(mztableSelect)){
          na <- names(mztableSelect)
          mztableSelect <- t(as.matrix(mztableSelect))
        }
        mztable.h2o <- as.h2o(mztableSelect)
        
        
        
        if(Parallelh20){
          cl <- makeCluster(cores <-detectCores())
          cat("Using",cores,"cores")
          # PH <- parSapply(cl,uniSeq,function(x){PeptideHydrophobicity(x)})
          s1 <-system.time(PredictionScores <- parSapply(cl,Models.h2o,function(x,mztable.h2o){
            require(rawDiag)
            require(h2o)
            h2o.init()
            PREDI <- as.data.frame(h2o.predict(x,mztable.h2o))
            PREDI$predict
          },mztable.h2o=mztable.h2o))
          stopCluster(cl)
          
          
          
        }else{
          s2 <-system.time(
            PredictionScores <- sapply(Models.h2o,function(x){
              PREDI <- as.data.frame(h2o.predict(x,mztable.h2o))
              PREDI$predict
            })
          )
        }
        
        
        
        
        PredictionScores <- data.frame(PredictionScores)
        mztable$RF_Scores <- PredictionScores$DRF
        mztable$DL_Scores <- PredictionScores$DeepLearning
        # mztable$GBM_Scores <- PredictionScores$GBM_model_R_1572348889847_41
        # mztable$GLM_Scores <- PredictionScores$GLM_model_R_1572348889847_40
        dbWriteTable(db,mzname,mztable,overwrite = T)
      }
      par(mfrow = c(2,2),mai = c(0.1,0.8,0.4,0.1))
      mztable[,{
        try({
          
          
          plot(RT_min,SCAall,type = "p",pch = 20,cex = 0.2,main = "",ylab = "",axes = F,main = "SCA over time")
          if(!hideMain){
            mtext(rawfile,cex = 0.8,adj = 0,line = 0.5)
          }
          par(new = T)
          plot(RT_min,RF_Scores,type = "p",col = 2,pch = 20,cex = 1,"RF Score over time")
          
          plot(RT_min,SCAall,type = "p",pch = 20,cex = 0.2,main = "",ylab = "",axes = F)
          par(new = T)
          plot(RT_min,DL_Scores,type = "p",col = 3,pch = 20,cex = 1,"DL Score over time")
          # 
          # plot(RT_min,SCAall,type = "p",pch = 20,cex = 0.2,main = "",ylab = "",axes = F)
          # par(new = T)
          # plot(RT_min,GBM_Scores,type = "p",col = 4,pch = 20,cex = 1)
          
          # plot(RT_min,SCAall,type = "p",pch = 20,cex = 0.2,main = "",ylab = "",axes = F)
          # par(new = T)
          # plot(RT_min,GLM_Scores,type = "p",col = 5,pch = 20,cex = 1)
        },silent = T)
        NULL
        
      },rawfile]
      
      
      
    })
    
    
    dev.off()
    #dbDisconnect(db)
    
  })
  
  # db <- dbConnect(SQLite(),"PickyAnalyzer.sqlite")
  # 
  # 
  # forScore <- lapply(grep("^mz",dbListTables(db),value = T),function(x){dbGetQuery(db,statement = paste("SELECT Score FROM '",x,"'",sep = ""))})
  # revScore <- lapply(grep("^revmz",dbListTables(db),value = T),function(x){dbGetQuery(db,statement = paste("SELECT Score FROM '",x,"'",sep = ""))})
  # hist(log(unlist(forScore)),col = 2,breaks = 1000,border = NA,xlab = "log Score")
  # hist(log(unlist(revScore)),col = 3,add = T,breaks = 1000,border = NA)
  # 
  # 
  # dbDisconnect(db)
  return(FUN)
}
## Called in Scoring Wrapper
ValiScoringFunction <- function(tab,group,it=0,parallelization = F,topFragments = 5){
  # tab <<- tab
  
  # group <<- group
  
  cat("\r",it)
  COLUMNS <- grep("charge",names(tab))
  COLUMNS <- COLUMNS-1
  
  mztable.s1 <- tab
  NOTUSED <- apply(mztable.s1[,1:COLUMNS],2,function(x){all(x == -1)})
  NOTUSEDid <- which(!NOTUSED)
  if(length(NOTUSEDid) > topFragments){
    NOTUSEDid <- NOTUSEDid[1:topFragments]
  }
  if(length(NOTUSEDid)>1){
    cross <- combn(NOTUSEDid,2)
    crossSel <- apply(cross,2,function(x){any(is.na(x))})
    cross <- cross[,!crossSel]
    if(is.vector(cross)){
      cross <- as.matrix(cross)
      
    }
  }else{
    cross <- NULL
  }
  
  # print(cross)
  if(length(cross)>0){
    
    cr <- cross
    # print(cross)
    par(mfrow = c(4,1))
    # Cross Correlation between fragments, requires parallelization...
    if(parallelization){#
      cl <- makeCluster(detectCores())
      s1 <-system.time(Cor <- parSapply(cl,as.list(data.frame(cross)),TimeSeriesLocalCor,mztable.s1))
      stopCluster(cl)
      
    }else{
      s2 <- system.time(Cor <- apply(cross,2,TimeSeriesLocalCor,mztable.s1))
      
      if(length(dim(Cor))==0){
        Cor <- t(as.matrix(Cor))
      }
      
    }
    # Cor[is.na(Cor)]<-0
    Corproc <- apply(Cor,1,function(x){RE<- 0;try({RE <- sum(x,na.rm = T)/sum(!is.na(x),na.rm = T)});RE})
    Corproc[is.na(Corproc)] <- 0
    # Corproc <- Corproc+1
    # Corproc <- apply(Cor,1,function(x){median(x,na.rm = T)})
    
    
    di <- c(Corproc,Corproc[length(Corproc)]) -c(Corproc[1],Corproc)
    out <<- 0
    OUT <- sapply(1:length(Corproc),function(x){
      # print(out)
      out <<- di[x]+out
      unlist(out)
    })
    
    Cor[is.na(Cor)] <- 0
    if(is.vector(Cor)){
      Cor <- matrix(Cor,byrow = T,ncol = 3)
    }
    
    if(dim(mztable.s1)[1]>2){
      OUTFU <- apply(Cor,1,function(x){c(median(x,na.rm = T,probs = 0.8),mean(x,na.rm = T),sd(x,na.rm = T),sum(x,na.rm = T))})
      OUTFU <- t(OUTFU)
      OUTFU <- OUTFU+1
      mztable.s1$CorrelationScoreMedian <- c(NA,OUTFU[,1])+1
      mztable.s1$CorrelationScoreMean <- c(NA,OUTFU[,2])+1
      mztable.s1$CorrelationScoreSD<- c(NA,OUTFU[,3])
      mztable.s1$CorrelationScoreSum<- c(NA,OUTFU[,4])+1
      mztable.s1$CorrelationScoreSumSpecial <- c(NA,Corproc)+1
      mztable.s1$CorrelationScoreSumSpecial2 <- c(NA,Corproc*OUT)+1
      
      
      
      # plotTransitions(as.numeric(x$RT_Used),transitions,yl = yl,xlim = xl,TransitionColors = TransitionColors,col = colAll,add = add,frame = F,xlab = xla,ylab = yla,Ppos = Ppos,axes = !blankplot)
      # par(new = T)
      # points(mztable.s1$RT_Used,mztable.s1$CorrelationScoreMean)
    }else{
      mztable.s1$CorrelationScoreMedian <- 1
      
      mztable.s1$CorrelationScoreMean <- 1
      mztable.s1$CorrelationScoreSD<- 0
      mztable.s1$CorrelationScoreSum<-1
      mztable.s1$CorrelationScoreSumSpecial <- 1
      mztable.s1$CorrelationScoreSumSpecial2 <- 1
      
    }
    
    
    mztable.s1[,Score:={
      SCAall*CorrelationScoreSumSpecial2},]
    
  }else{
    # print("HUMPE")
    mztable.s1$CorrelationScoreMedian <- 1
    mztable.s1$CorrelationScoreMean <-1
    mztable.s1$CorrelationScoreSD<- 0
    mztable.s1$CorrelationScoreSum <- 1
    mztable.s1$CorrelationScoreSumSpecial <- 1
    mztable.s1$CorrelationScoreSumSpecial2 <- 1
    
    mztable.s1[,Score:=SCAall,]
    
  }
  
  # mztable.s1.out <<- mztable.s1
  
  mztable.s1 <- as.list(mztable.s1)
  mztable.s1 <- lapply(mztable.s1,as.double)
  
}
### Scoring, called in ValiScoringFunction
TimeSeriesLocalCor <- function(x,mztable.s1,width = 5){
  library(data.table)
  library(gtools)
  FUN <- mztable.s1[,.SD,.SDcols = x]
  FUN[FUN == -1] <- NA
  
  names(FUN) <- c("a","b")
  var_1 <- FUN$a#+0.1
  var_2 <- FUN$b#+0.2
  chg_1 <- diff(var_1)/var_1[-length(var_1)] 
  chg_2 <- diff(var_2)/var_2[-length(var_2)] 
  
  hum <- running(chg_1, chg_2, fun=cor, width=width, by=1, allow.fewer=TRUE, align=c("left"), simplify=TRUE)
  return(abs(hum))
}

## MS1 Compiler input is the working directory with the 
CompileMS <- function(wd,dbname = "VALIdb",folder_PRM_Matches ="PRM_Analyzer_Matches_." ,session = NULL){
  
  db <- dbConnect(SQLite(),dbname)
  
  
  Folders <- list.dirs(dirname(dbname))[-1]
  
  if(length(session) != 0){
    progress <- Progress$new(session,min=0, max=8)
    on.exit(progress$close())
    progress$set(message = 'compiling MS1 Scans',
                 detail = path,value = 0)
  }
  
  uni <- sapply(Folders,function(path){
    # path <<- path
    fi <- list.files(path,pattern = "ms1scans",full.names = T)
  })
  uni <- unique(unlist(uni))
  it <<- 0
  
  hum <- sapply(unique(basename(uni)),function(path){
    path<<- path
    if(length(session) != 0){
      # on.exit(progress$close())
      try({
        progress$set(message = 'compiling MS1 Scans',
                     detail = basename(path),value = it)
        
      })
      
    }
    
    it <<- it+1
    # path <<- path
    
    cat("\r",basename(path),it,length(unique(basename(uni))))
    
    Fi1 <- uni[basename(uni)==path]
    outfi <- c()
    for(i in Fi1){
      tab <- fread(i,sep = "\t",stringsAsFactors = F,showProgress = F)
      outfi <- rbind(outfi,tab)
    }
    if(dim(outfi)[1] == 0){
      
      return(NULL)
    }else{
      print("DATA MATCH")
    }
    file.i <- copy(outfi)
    rm(outfi)
    tab <- which(file.i$`mz1+`!=0)
    tabm <- tab-1
    tabp <- tab+1
    zeros <- setdiff(c(tabm,tabp),tab)
    filter <- sort(unique(c(zeros,tab)))
    
    file.i <- file.i[filter,]
    
    
    # file.i<- fread(pathi,sep = "\t",check.names = F)
    file.i$rawfile <- basename(file.i$rawfile)
    ms1scan <- file.i
    if(length(unique(ms1scan$RT)) == 1){
      ms1scan$RT <- ms1scan$ScanNumber
    }
    ms1scan <- ms1scan[order(RT)]
    # Preparing Ms1Scan
    # ms1scan <- data.table(ms1scan)
    
    Da <- ms1scan[,{.SD},.SDcols = grep("^mz.",colnames(ms1scan),value = T)]
    
    
    cross <- combn(colnames(Da),2)
    DAzero <- apply(Da,1,function(x){all(x == 0)})
    # CorrelationScore
    # s2 <- system.time(Cor <- apply(cross,2,TimeSeriesLocalCor,data.table(Da)))
    # Cor[is.na(Cor)] <- 0
    # if(is.vector(Cor)){
    #   Cor <- matrix(Cor,byrow = T,ncol = 3)
    # }
    
    # OUTFU <- apply(Cor,1,function(x){c(median(x,na.rm = T),mean(x,na.rm = T),sd(x,na.rm = T))})
    # OUTFU <- t(OUTFU)
    # OUTFU <- OUTFU+1
    
    
    ms1scan$MS1_CorrelationScoreMedian <- 1#c(NA,OUTFU[,1])
    ms1scan$MS1_CorrelationScoreMean <- 1#c(NA,OUTFU[,2])
    ms1scan$MS1_CorrelationScoreSD<- 1#c(NA,OUTFU[,3])
    
    # SCA
    # cat("\r",basename(path),"MS1 SCA")
    sequence <- strsplit(gsub("__","_",basename(path)),"_",fixed = T)
    sequence <- unlist(sequence)[3]
    
    AA <-  AAtab()
    
    sequence <- gsub("\\([^()]*\\)", "", gsub("\\([^()]*\\)", "", sequence))
    
    AAtab_OutSequence <<- c()
    sapply(unlist(strsplit(sequence,"")),function(x){
      # x <<- x
      sapply(which(x==AA$AA_L),function(x){
        AAtab_OutSequence <<- rbind(AAtab_OutSequence,AA[x,])
      })
      return(NULL)
    })
    COMPO <- apply(AAtab_OutSequence[,grep("count$",colnames(AAtab_OutSequence))],2,function(x){sum(as.numeric(x))})
    COMPO <- COMPO[COMPO!=0]
    Atoms <- sapply(strsplit(names(COMPO),"_"),function(x){x[1]})
    COMPOSITION <- paste(Atoms,COMPO,sep = "",collapse = "")
    predictPattern <- T
    
    if(require(enviPat)){
      data(isotopes)
      ISOP <- isopattern(isotopes,COMPOSITION)
      
      apply(ISOP[[1]][,-c(1:2)],1,sum)
      ISOP <- ISOP[[1]]
      ISOP <- data.frame(ISOP,stringsAsFactors = F)
      ISOFUN <- aggregate(ISOP$abundance,list(round(ISOP$m.z)),max)
    }else{
      predictPattern <- F
    }
    
    ms1scan <- unique(ms1scan)
    ms1scan[,MS1_SCA:={
      temp <- .SD
      
      
      AD <- temp[,{.SD},.SDcols=grep("^mz.",colnames(temp),value = T)]
      AD <- unlist(AD)
      AD[is.na(AD)] <- 0
      if(length(which(AD!=0))>=3){
        if(predictPattern){
          pattern = ISOFUN$x[1:length(AD)]
          
          
        }else{
          pattern = c(seq(0,1,length.out = length(AD)))
        }
        NormSpecAngle(AD,pattern)
      }else{
        0
      }
      # Predict Pattern:
      
      
    },.(ScanNumber,RT)]
    
    
    
    # ms1scan[,MS1_Score:={
    # MS1_SCA*MS1_CorrelationScoreMedian},]
    # ms1scan$MS1_Score <- ms1scan$MS1_SCA
    ms1scan[,MS1_Count:=apply(.SD,1,function(x){sum(!is.na(x))}),.SDcols=grep("^mz.",colnames(ms1scan)),RT]
    # dbRemoveTable(db,basename(path))
    dbWriteTable(db,ms1scan,name = gsub(".txt$","",gsub(" ","",gsub("__","_",basename(path)),fixed = T)),append = T)
    
    return(NULL)  
  })
  
  if(length(session) != 0){
    progress$close()
  }
  dbDisconnect(db)
  
}



##  NOT USED Calculates Coaccuring peaks across transitions over time xhui= transition table, wi = with in min, dl = rt.accuracy (see cut())
PeaksFoundCheck <- function(xhui,wi = 1,dl = 4){
  xhui <- xhui[order(unfactor(xhui$RT_min)),]
  TransitionsMatrixAll <- xhui[,1:(which(colnames(xhui) == "charge")-1)]
  ch <- sapply(unique(xhui$rawfile),function(rawfile){
    # rawfile <<- rawfile
    rfsel <- xhui$rawfile == rawfile
    TransitionsMatrix <- TransitionsMatrixAll[rfsel,]
    PeaksList <- apply(TransitionsMatrix,2,function(x){unfactor(xhui$RT_min[rfsel])[find_peaks(unfactor(x))]})
    range(unlist(PeaksList))
    x <- PeaksList[[2]]
    Cut <- round(seq(min(xhui$RT_min,na.rm = T),max(xhui$RT_min,na.rm = T),by = wi),1)
    Cut1 <- round(seq(min(xhui$RT_min,na.rm = T)+wi/3,max(xhui$RT_min,na.rm = T),by = wi),1)
    Cut2 <- round(seq(min(xhui$RT_min,na.rm = T)+wi/3*2,max(xhui$RT_min,na.rm = T),by = wi),1)
    
    humpi <- sapply(PeaksList,function(x){
      x <- x
      Hum <-  list(cut(x,Cut,dig.lab = dl),
                   cut(x,Cut1,dig.lab = dl),
                   cut(x,Cut2,dig.lab = dl))
      ta <- table(unlist(Hum))
      ta <- ta[order(names(ta))]
      Range <- ta[ta >0]
      return(names(Range))
      
    })
    ta <- table(unlist(humpi))
    ta <- ta[order(names(ta))]
    OutFu <<- c()
    taMax = ta[ta == max(ta)]
    SAPI<- sapply(1:length(ta),function(x){
      nam <- names(ta)[x]
      nam <- strsplit(gsub("\\(|]","",nam),",")
      nam <- unlist(nam)
      Te <- cbind(nam[1] : nam[2],ta[x])
      OutFu <<- rbind(OutFu,Te)
      return(Te)
    })
    
    
    
    OutFu <- cbind(OutFu,rawfile)
    return(OutFu)
  })
  names(ch) <- unique(xhui$rawfile)
  return(ch)
}
### Used in PeaksFoundCheck
find_peaks <- function (x, m = 3,session = NULL){
  # https://github.com/stas-g/findPeaks
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1],na.rm = T)) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

PeaksAligner <- function(x,PeaksList){
  Peaks = c()
  # par(mfrow = c(3,2))
  for(r in unique(x$rawfile)){
    tp <- PeaksList[names(PeaksList)==r][[1]]
    sel <- x$rawfile == r
    if(is.vector(tp)){
      Peaks <- tp # Not clear what should happen here, should be better evaluateds
      
    }else{
      MaFu <- tp[match(round(x$RT_min[sel],1),as.numeric(tp[,1])),2]
      Peaks[sel ] <- MaFu
      
    }
    
    # plot(x$RT_min[sel],x$MatchCount[sel],type = "l")
    # par(new = T)
    # plot(x$RT_min[sel],MaFu,col = 2,pch = 5)
    # 
    # plot(x$IntensitySum[sel],MaFu)
  }
  return(Peaks)
}


################ UI functions ##############

## Saves Precursor Rating, might be obsolete
SavePrecursorRating <- function(rating,dbtaName,db){
  for(n in dbtaName){
    db <- dbConnect(SQLite(),db)
    dbWriteTable(db,data.frame(rating = rating),name = paste("RATING_",n,sep = ""),overwrite = T)
    dbDisconnect(db)
    # print(paste("RATING_",n,sep = ""))
  }
}
## Browse File Button
BrowseFun <- function(id,session){
  Continue <- T
  library(tcltk)
  if(any(grep("windows",Sys.info(),ignore.case = T),na.rm = T)){
    Ch <- choose.dir()
    if(is.na(Ch)){
      Continue <- F
    }
  }else{
    Ch <- tk_choose.dir()
    if(is.na(Ch)){
      Continue <- F
    }
    
  }
  
  if(Continue&file.exists(Ch)){
    updateTextInput(session,inputId = id,value = Ch)
  }
}
## PEAKS Finding MODUL
RTcandidatesFun <- function(ana,FDRcutoff,RTranges = NULL,supersmooth=F,supersmooth_bw=0.2,Requantify_Priority="DL_Scores"){
  save(ana,FDRcutoff,RTranges,supersmooth,supersmooth_bw,file="RTcandidatesFun_temp.rda")
  if(length(RTranges)!=0){
    if(length(RTrange)!=2){
      RTranges <- range(sapply(ana,function(x){range(x$RT_Used)}),na.rm = T)
    }
    if(any(is.na(RTranges))){
      RTranges <- range(sapply(ana,function(x){range(x$RT_Used)}),na.rm = T)
    }
    if(any(is.na(RTranges))){
      RTranges <- c(0,100000)
    }
    
    RTcandidates <- lapply(ana,function(x){
      x <- data.table(x)
      
      RTsel <- x$RT_Used >= RTrange[1] & x$RT_Used <= RTrange[2]
      
      x[x$FDR<=FDRcutoff&RTsel,.(rawfile,RT_Used,DL_Scores),]
      
    })
  }else{
    RTcandidates <- lapply(ana,function(x){
      x <- data.table(x)
      if(supersmooth){
        x$DL_Scores <- supsmu(x$RT_Used,x$DL_Scores,span = supersmooth_bw)$y
        # plot(x$RT_Used,x$DL_Scores,type="l")
        # points(supsmu(x$RT_Used,x$DL_Scores,bass = 1),col = 2,type="l")
      }
      x[x$FDR<=FDRcutoff,.(rawfile,RT_Used,DL_Scores,IntensitySum),]
      
    })
    
  }
  
  
  
  # Determining best Peaks after FDR CutOff
  # print("Determining best Peak")
  
  # Extracting available RawFiles
  RF <- unique(unlist(sapply(ana,function(x){x$rawfile})))
  # print("COnverting RTcandidates")
  # Converting RTcandidateslist to table
  tempfunCOMBIALL <<- c()
  hum <- lapply(RF,function(rf){
    tempfunCOMBI <<- c()
    tempfun <<-  (lapply(RTcandidates,function(x){
      tempfunCOMBI <<- rbind(tempfunCOMBI,x[x$rawfile == rf,])
    }))
    tempfunCOMBIALL <<- rbind(tempfunCOMBIALL,tempfunCOMBI)
    NULL
  }) 
  
  tempfunCOMBIALL <- data.table(tempfunCOMBIALL)
  # par(mfrow = c(2,4))
  # tempfunCOMBIALL[,plot(RT_Used,GLM_Scores,main =.BY$rawfile) ,rawfile]
  
  # Extracting RT with best Score
  if(length(tempfunCOMBIALL)>0){
    if(F){
      #option: get 90% confidence from rounded RTs
      tempfunCOMBIALL2 <- tempfunCOMBIALL[,.(DL_Scores=quantile(DL_Scores,0.9),IntensitySum=quantile(IntensitySum,0.9)),.(RT_Used=round(RT_Used,1),rawfile)]
    }
    CANDIDATE_RT <- tempfunCOMBIALL[,.SD[DL_Scores==max(DL_Scores,na.rm = T)],rawfile]
    
    if(Requantify_Priority=="DL_Scores"){
      CANDIDATE_RT <- tempfunCOMBIALL[,.SD[DL_Scores==max(DL_Scores,na.rm = T)],rawfile]
    }
    if(Requantify_Priority=="Intensity"){
      CANDIDATE_RT <- tempfunCOMBIALL[,.SD[IntensitySum==max(IntensitySum,na.rm = T)],rawfile]
    }
    
  }else{
    CANDIDATE_RT <-data.table(rawfile = as.vector(RF),RT_Used=NA,DL_Scores = NA)
  }
  CANDIDATE_RT
}
## PLOTTING
plotTransitions <- function(Index,x,yl = NULL,add = F,TransitionColors = T,col = 1,SelectedRetentionTime = NA,Ppos = NULL,addGrid = F,lwdpoints = 1,maximalnumber = 500,plotTransitions=NULL,...){
  
  # x <<- x
  # Index <<- Index
  # yl <<- yl
  if(any(is.infinite(yl))){
    yl <- NULL
  }
  # add <<- add
  TransitionColors <-T
  # lwdpoints <<- lwdpoints
  ylim <- range(as.numeric(unlist(x)),na.rm = T)
  Index <- as.numeric(Index)
  if(length(yl) == 0){
    yl <- ylim
  }
  if(!add){
    # Ppos <<- Ppos
    plot(Index,Index,type = "n",ylim = yl,...)
    if(length(Ppos) > 0){
      rect(xleft = unfactor(Ppos[1]),xright = unfactor(Ppos[2]),ybottom = par()$usr[3],ytop = par()$usr[4],col = "#20202010",border = NA)
      
    }
    
  }
  if(addGrid){
    grid()
  }
  xtemp <- x
  for(i in 1:dim(x)[2]){
    if(length(selected_transitions)==0){
      plot.transition <- T
    }else{
      plot.transition = any(selected_transitions==colnames(x)[i])
      # print("HUI")
    }
    
    if(plot.transition){
      b <- x[,i]
      b[is.na(b)] <- 0
      
      if(TransitionColors){
        # col = rainbow(dim(x)[2])[i]
        coli = col[i]
        lty = 1
      }else{
        lty = i
        coli <- col[1]
      }
      a <- as.numeric(Index)
      b <- as.numeric(b)
      
      if(length(a)>maximalnumber){
        aset  <- cut(a,maximalnumber )
        
        if(length(coli)==length(aset)){
          lwdpointsset <- aggregate(lwdpoints,list(aset),max,na.rm = T)$x
          
        }else{
          lwdpointsset <- lwdpoints
          
        }
        bset <- aggregate(b,list(aset),max,na.rm = T)$x
        aset <- aggregate(a,list(aset),mean,na.rm = T)$x
      }else{
        aset <- a
        bset <- b
        lwdpointsset <- lwdpoints
      }
      
      points(aset,bset,type = "l",col = coli,pch = 20,cex = 0.5,lty = lty,lwd = lwdpointsset)
      
    }
  }
  
  if(!is.na(SelectedRetentionTime)){
    abline(v = SelectedRetentionTime,col = 2)
  }
  
}
plot.analyzed.transition <- function(x,type = "Volcano",yl = NULL,xl = NULL,Ppos = NULL,add = F,FDRCutOff = 0.01,TransitionColors = T,colAll = 2,secPlotType = "Intensity Correlation",p.value = 0.01,RetTime = NULL,DDA = 1:3,PRMonly = T,AddSecondPlot = T,onlySignificant = F,centerPeak = F,colmap=NULL,blankplot = F,FDRType = "FDR",
                                     maximalnumber = 500,selected_transitions=NULL,hideMain=F
                                     ,...){
  xhui <- x
  
  x <- xhui
  x <- x[order(x$RT_Used),]
  xhui <- x
  if(FDRType == "SCAfdr"){
    x$FDR_Choosen <- x$FDR
    xhui$FDR_Choosen <- xhui$FDR
    
  }
  if(FDRType == "FDR"){
    x$FDR_Choosen <- x$FDR
    xhui$FDR_Choosen <- xhui$FDR
    
  }
  # if(FDRType == "SCA_Limit"){
  #   x$FDR_Choosen <- x$SCA_Limit
  #   xhui$FDR_Choosen <- xhui$SCA_Limit
  #   
  # }
  # 
  # if(FDRType == "SCA_mz_fdr"){
  #   x$FDR_Choosen <- x$SCA_mz_fdr
  #   xhui$FDR_Choosen <- xhui$SCA_mz_fdr
  #   
  # }
  # 
  # if(FDRType == "PEP"){
  #   x$FDR_Choosen <- x$PEP
  #   FDR <- FDRCalc(x)
  #   x$FDR_Choosen <- FDR[,2]
  #   x$Score <- FDR[,1]
  # }
  
  # NOT RELEVANT ANYMORE, FDR NOW BASED ON PEP
  
  # x <- x
  # type <<- type
  # yl <<- yl
  # xl <<- xl
  # add <<- add
  # TransitionColors <- TransitionColors
  # colAll <<- colAll
  # secPlotType <<- secPlotType
  print("Preparing Table")
  if(length(x) == 0){
    return(NULL)
  }
  col <- rep(1,length(x$R_p))
  col[x$FDR_Choosen <= FDRCutOff] <- 2
  x$col <- col
  x$IntensitySum <- as.numeric(x$IntensitySum)
  ce <- (x$IntensitySum/max(x$IntensitySum,na.rm = T))
  x$ce <- ce
  # ce[ce < 0.3] <- 0.3
  if(type == "Volcano"){
    pep <- xhui$FDR_Choosen
    pep[pep == 0] <- min(pep[pep!=0],na.rm = T)*0.01
    plot(xhui$SCA,-log10(pep),col = as.numeric(xhui$FDR_Choosen<= FDRCutOff)+1,pch = 20 ,ylab = "-log10 PEP",xlab = "Normalized Spectrum Contrast Angle")
    # x$col[x$col == 2] <- colAll
    # if(!add){
    # plot(x$R,-log10(as.numeric(x$R_p)),ylab = c("-log10 p"),xlab = "R",cex = x$ce,col = x$col,pch = 20,...)
    # }else{
    # points(x$R,-log10(as.numeric(x$R_p)),ylab = c("-log10 p"),xlab = "R",cex = x$ce,col = x$col,pch = 20,...)
    # }
  }else{
    
    
    
    
    if(add){
      par(new = T)
      plot(1,type = "n",xlab = "",ylab = "",axes = F,main = "",yl = yl,xlim = xl)
      
      
      
    }
    if(onlySignificant){
      x$RT_Used[col != 2] <- NA
      xhui$RT_Used[col != 2] <- NA
      
    }
    
    transitions <- x[,1:(grep("charge",colnames(x))-1)]
    transitions$rawfile <- NULL
    if(length(yl) !=2){
      
      yl <- c(0,max(as.numeric(unlist(transitions)),na.rm = T))
    }
    if(length(xl) == 2){
      traset <- transitions[as.numeric(xhui$RT_Used)>=xl[1]&as.numeric(xhui$RT_Used)<= xl[2],]
      traset <- as.numeric(unlist(traset))
      traset <- traset[!is.infinite(traset)]
      mt <- max(traset,na.rm = T)
      
      if(is.infinite(mt)){
        mt <- max(as.numeric(unlist(transitions)),na.rm = T)
      }
      
      if( mt < max(yl,na.rm = T)){
        if(!is.infinite(mt)){
          yl[2] <- mt
        }
      }
    }
    # ggplot(mtcars, aes(wt, mpg)) +
    #   geom_point() +
    #   coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    # ggplot(xhui,aes(RT_Used,b13.H2O)+geom_li)
    # colmap <<- colmap
    if(TransitionColors&length(colmap) > 0){
      colAll <- colmap[match(colnames(transitions),colmap[,1]),2]
    }
    if(blankplot){
      xla = ""
      yla = ""
    }else{
      xla = "time [min]"
      yla = "Intensity"
    }
    # xl = c(2,40)
    # blankplot = F
    # colAll <- 2
    # yl <- c(0,100000000)
    # colAll <- c(1:20)
    # xl<- c(30,40)
    xbackup <- x
    cat("\rStart plotTransitions")
    plotTransitions(Index = as.numeric(x$RT_Used),x = transitions,yl = yl,xlim = xl,
                    TransitionColors = TransitionColors,col = colAll,add = add,
                    frame = F,xlab = xla,ylab = yla,Ppos = Ppos,axes = !blankplot)
    cat("\rFinished plotTransitions")
    
    x <- xbackup
    #abline(v=x$RT_Used[col == 2],lty = "dotted",col = "#99999940")
    
    if(length(RetTime) > 0){
      # points()
      # stop("fjwoeifjoi")
      di <- diff(par("usr")[3:4])
      points(RetTime,rep(par("usr")[3],2)+di*0.01,col = colAll,lwd = 3,type = "l")
    }
    if(AddSecondPlot){
      pax <- par()$usr[1:2]
      # pax <- par()$xaxp[1:2]
      par(new = T)
      if(!TransitionColors){
        x$col[x$col == 2] <- colAll
      }else{
        x$col <- col
      }
      # "FDR","RF_Score","SCAall","SCAcut","MatchCount","Matches","DiffSum","Peaks1","Peaks2"
      if(secPlotType == "R"){
        yv <-  as.numeric(x$R)
        yv[yv < 0] <- 0
        ylab = "Spearman's Correlation"
        
      }
      if(secPlotType == "R_p"){
        yv <-  -log10(as.numeric(x$R))
        yv[yv < 0] <- 0
        ylab = "p-value from Spearman's Correlation"
        
      }
      if(secPlotType == "MatchCount"){
        yv <-  x$MatchCount
        ylab = "# of Fragment Matches"
        
      }
      if(secPlotType == "RF_Scores"){
        yv <-  x$RF_Scores
        ylab = "RF Score"
        
      }
      if(secPlotType == "GLM_Scores"){
        yv <-  x$GLM_Scores
        ylab = "Score"
        
      }
      if(secPlotType == "GBM_Scores"){
        yv <-  x$GBM_Scores
        ylab = "Score"
        
      }
      if(secPlotType == "DL_Scores"){
        yv <-  x$DL_Scores
        ylab = "DL Score"
        
      }
      if(secPlotType == "mz"){
        yv <-  x$mz
        ylab = "m/z"
        
      }
      if(secPlotType == "SCAall"){
        yv <-  x$SCAall
        ylab = "SCA all fragments"
        
      }
      if(secPlotType == "SCAcut"){
        yv <-  x$SCAcut
        ylab = "SCA top 5 fragments"
        
      }
      # secPlotType <- "FDR"
      if(secPlotType == "FDR"){
        yv <-  -x$FDR
        ylab = "FDR"
        
      }
      if(secPlotType == "PEP"){
        yv <-  x$PEP
        yv[yv >1] <- 1
        yv[yv == 0] <- min(yv[yv != 0],na.rm = T)*0.0001
        yv <- -log10(yv)
        ylab = "PEP"
        
      }
      if(secPlotType == "Peaks1"){
        yv <-  x$Peaks
        yv[is.na(yv)] <- 0
        ylab = "Cross Fragement Peaks"
        
        
      }
      if(secPlotType == "Peaks2"){
        yv <-  x$Peaks
        yv[is.na(yv)] <- 0
        ylab = "Cross Fragement Peaks\nafter alignment"
      }
      
      if(secPlotType == "DiffSum"){
        yv <-  -x$DiffSum
        yv[is.na(yv)] <- 100
        ylab = "Delta Mass [Da]"
        
        
      }
      if(secPlotType == "SCA_mz_fdr"){
        yv <-  -log10(x$SCA_mz_fdr)
        
        yv[is.na(yv)] <- 1
        
        ylab = "-log10 FDR"
        
      }
      if(secPlotType == "SCAfdr"){
        # yv <-  -x$SCAfdr
        # yv[yv<=-0.2] <- -0.23
        # yv[is.na(yv)] <- -0.25
        # ylab = "FDR"
        yv <-  -log10(x$SCAfdr)
        
        yv[is.na(yv)] <- 1
        
        ylab = "-log10 FDR"
        
      }
      if(secPlotType == "FDR"){
        # yv <-  -x$SCAfdr
        # yv[yv<=-0.2] <- -0.23
        # yv[is.na(yv)] <- -0.25
        # ylab = "FDR"
        yv <-  -log10(x$FDR)
        
        yv[is.na(yv)] <- 1
        
        ylab = "-log10 FDR"
        
      }
      if(secPlotType == "CorrelationScoreSumSpecial2"){
        yv <-  x$CorrelationScoreSumSpecial2-1
        yv[is.na(yv)] <- 0
        ylab = "CorrelationScore"
        
      }
      if(secPlotType == "CorrelationScoreSumMedian"){
        yv <-  x$CorrelationScoreMedian-1
        yv[is.na(yv)] <- 0
        ylab = "Median CorrelationScore"
        
      }
      if(secPlotType == "CorrelationScoreSumSpecial"){
        yv <-  x$CorrelationScoreSumSpecial-1
        yv[is.na(yv)] <- 0
        ylab = "Median CorrelationScore"
        
      }
      if(secPlotType == "CorrelationScoreMean"){
        yv <-  x$CorrelationScoreMean-1
        yv[is.na(yv)] <- 0
        ylab = "Mean CorrelationScore"
        
      }
      if(secPlotType == "CorrelationScoreSD"){
        yv <-  x$CorrelationScoreSD
        yv[is.na(yv)] <- 0
        ylab = "SD CorelationScore"
        
      }
      
      
      
      if(!exists("yv")){
        matchfun <-  match(secPlotType,colnames(x))
        if(is.na(matchfun)){
          yv <- rep(0,dim(x)[1])
          ylab ="not set"
        }else{
          yv <- x[,matchfun]
          ylab =secPlotType
        }
        
      }
      ra <- rank(yv)
      ra <-ra/max(ra,na.rm = T)
      ra <- round(ra*100)-1
      
      ra[ra <0] <- 0
      ra[nchar(as.character(ra))==1] <- paste(0,ra[nchar(as.character(ra))==1],sep = "")
      col <- x$col
      col[col == 1] <- paste("#999999",ra[col == 1],sep = "")
      col[col == 2] <- paste("#FF0000",ra[col == 2],sep = "")
      if(length(col)==0){
        col <- 1
      }
      # xc<<- x
      
      a <- x$RT_Used
      b <- yv
      if(length(a)>maximalnumber){
        aset  <- cut(a,maximalnumber )
        
        if(length(col)==length(a)){
          # col <<- col
          colset <- aggregate(col,list(aset),function(x){
            xgrep <-grep("#FF00",x,value = T)
            if(length(xgrep)>0){
              # x <<- x
              return(as.character(xgrep[1]))
            }else{
              numi <- as.numeric(substr(x,7,9))
              return(as.character(x[numi==max(numi)][1]))
            }
            
          })$x
          colset <- as.character(colset)
          
        }else{
          colset <- col
          
        }
        # ra <<- ra
        raset <- aggregate(as.numeric(ra),list(aset),max,na.rm = T)$x
        
        bset <- aggregate(b,list(aset),max,na.rm = T)$x
        aset <- aggregate(a,list(aset),mean,na.rm = T)$x
      }else{
        aset <- a
        bset <- b
        raset <- ra
        colset <- col
      }

      plot(aset,bset,axes = F,frame = F,xlab = "",ylab = "",type = "p",col = colset,pch = 2,cex = 0.5*as.numeric(raset)/100,xlim = xl)
      if(!blankplot){
        axis(4,labels = NA,tck = 0.01)
        axis(4,tck = NA,col = NA,line = -2.5)
        
        mtext(ylab,4,line = -3)
        par(new =T)
        # plot(1)
        # axis(4,padj = 0.1)
        plot(as.numeric(xhui$RT_Used),apply(transitions,1,max,na.rm = T),ylim = yl,xlim = xl,type = "n",xlab = "",ylab = "",main = "",frame = F,axes = F)
        par(new = F)
      }
      
    }
  }
}
### MS1
#### Wrapper for plot.ms1scan
plot.ms1scanlist <- function(ms1scan_list,maximalnumber=500){
  plot_ms1scan(ms1scan_list$ms1scan_table,ylab = "MS1 Intensity",xlab = "RT [min]",xl = ms1scan_list$xl,yl = ms1scan_list$yl,add_Second = "MS1_Score")
}
##### MS1 Extraction Functions:
plot_ms1scan <- function(ms1scan_f,xl =NULL,yl = NULL,add_Second = "MS1_Score",...){
  
  # define ranges 
  if(length(xl) ==0|length(xl)!=2){
    xl <- range(ms1scan_f$RT)
  }
  
  Da <- ms1scan_f[,{.SD},.SDcols = grep("^mz.",colnames(ms1scan_f),value = T)]
  ms1scan_f$MS1_Score <- ms1scan_f$MS1_SCA
  try({
    if(add_Second != "none"){
      ms1scan_f[,{
        temp <- unlist(.SD)
        cexvec <- temp/max(temp,na.rm = T)*0.99
        colvec <-as.character(round(cexvec*100))
        colvec[is.na(colvec)] <- 0
        colvec[nchar(colvec) == 1] <- paste("0",colvec[nchar(colvec) == 1],sep = "")
        colvec[colvec=="99"] <- "FF"
        colvec[colvec=="98"] <- "FF"
        colvec[colvec=="97"] <- "FF"
        colvec[colvec=="96"] <- "FF"
        colvec[colvec=="95"] <- "FF"
        
        colvec <- paste("#999999",colvec,sep = "")
        
        plot(RT,unlist(.SD),pch = 3,cex = cexvec,lwd = cexvec,col = as.character(colvec),xlab = "",ylab = "",frame = F,axes = F,xlim = xl)
        axis(4)
        mtext(add_Second,4,line = 2)
        mtext(unique(ms1scan_f$mz)[1],3,cex = 0.8,col = "grey",adj = 0.1,line = -2)
      },.SDcols=add_Second]
      par(new = T)
      
    }
    
  })
  
  
  
  if(length(yl) ==0|length(yl)!=2){
    sel <- ms1scan_f$RT>= min(xl,na.rm = T)&ms1scan_f$RT <= max(xl,na.rm = T)
    yl <- range(Da[which(sel),])
  }
  # plot dummy:
  plot(1,type = "n",xlim = xl,ylim = yl,frame = F,...)
  it <<- 0
  cols <- (terrain.colors(dim(Da)[2]))
  apply(Da,2,function(x){
    it <<- it+1
    x[is.na(x)] <- 0
    
    points(ms1scan_f$RT,x,type = "l",col = cols[it] )
  })
  legend("topright",legend = c(paste("iso",1:5,sep = "")),lwd = 1,col = cols,bty = "n")
  
}

#### Called in plot.ms1scanlist
plot_ms1scan <- function(ms1scan_f,xl =NULL,yl = NULL,add_Second = "MS1_Score",...){
  
  # define ranges 
  if(length(xl) ==0|length(xl)!=2){
    xl <- range(ms1scan_f$RT)
  }
  
  Da <- ms1scan_f[,{.SD},.SDcols = grep("^mz.",colnames(ms1scan_f),value = T)]
  ms1scan_f$MS1_Score <- ms1scan_f$MS1_SCA
  try({
    if(add_Second != "none"){
      ms1scan_f[,{
        temp <- unlist(.SD)
        cexvec <- temp/max(temp,na.rm = T)*0.99
        colvec <-as.character(round(cexvec*100))
        colvec[is.na(colvec)] <- 0
        colvec[nchar(colvec) == 1] <- paste("0",colvec[nchar(colvec) == 1],sep = "")
        colvec[colvec=="99"] <- "FF"
        colvec[colvec=="98"] <- "FF"
        colvec[colvec=="97"] <- "FF"
        colvec[colvec=="96"] <- "FF"
        colvec[colvec=="95"] <- "FF"
        
        colvec <- paste("#999999",colvec,sep = "")
        
        plot(RT,unlist(.SD),pch = 3,cex = cexvec,lwd = cexvec,col = as.character(colvec),xlab = "",ylab = "",frame = F,axes = F,xlim = xl)
        axis(4)
        mtext(add_Second,4,line = 2)
        mtext(unique(ms1scan_f$mz)[1],3,cex = 0.8,col = "grey",adj = 0.1,line = -2)
      },.SDcols=add_Second]
      par(new = T)
      
    }
    
  })
  
  
  
  if(length(yl) ==0|length(yl)!=2){
    sel <- ms1scan_f$RT>= min(xl,na.rm = T)&ms1scan_f$RT <= max(xl,na.rm = T)
    yl <- range(Da[which(sel),])
  }
  # plot dummy:
  plot(1,type = "n",xlim = xl,ylim = yl,frame = F,...)
  it <<- 0
  cols <- (terrain.colors(dim(Da)[2]))
  apply(Da,2,function(x){
    it <<- it+1
    x[is.na(x)] <- 0
    
    points(ms1scan_f$RT,x,type = "l",col = cols[it] )
  })
  legend("topright",legend = c(paste("iso",1:5,sep = "")),lwd = 1,col = cols,bty = "n")
  
}

plot.TransplotList <- function(x,SimultaneousMassShiftvec= T,maximalnumber=500,selected_transitions=NULL,specificFragments=NULL){
  
  
  if(length(specificFragments)>0){
    com <- x$colmap
    com[,2] <- "#33333310"
    com[match(specificFragments,com[,1]),2] <- "#FF0000"
    x$colmap <- com
  }
  plot.analyzed.transition(x=x$x,
                           type = x$type,
                           frame = x$frame,
                           FDRCutOff = x$FDRCutOff,
                           Ppos = (x$Ppos),
                           colmap = x$colmap,
                           xl =x$xl,
                           yl = x$yl,
                           TransitionColors = x$TransitionColors,
                           col = x$col,
                           add = x$add,
                           secPlotType = x$secPlotType,
                           p.value = x$p.value,
                           RetTime = x$RetTime,
                           PRMonly = x$PRMonly,
                           AddSecondPlot = all(x$AddSecondPlot),
                           onlySignificant = x$onlySignificant,
                           xlab = x$xlab,
                           ylab = x$ylab,
                           blankplot = x$blankplot,
                           SelectedMass=x$SelectedMass,
                           maximalnumber = maximalnumber,selected_transitions=selected_transitions)
  
  cat("\rEnded TracePlot")
  # x <<- TransplotList
  if(SimultaneousMassShiftvec){
    try({
      if(x$namesAna == x$SelectedPrecursor){
        col = 2
      }else{
        col = 1
      }
      mtext(strsplit(x$namesAna,"_")[[1]][1],side = 3,cex = 0.8,line = -0.5,col = col,adj = 0.1,xpd = NA)
    })
  }
  if(setallrf){
    rfall.i.initialized <- rfall.i
    rfall.i <- rfall
    requiredPlots <- length(rfall.i) * length(ana)
    
    
    
    # par(mai= c(0.1,0.1,0.1,0),mfcol = c(ceiling(requiredPlots^0.5),ceiling(requiredPlots^0.5)))
    CenterPeakAlternative <- T
    
  }else{
    CenterPeakAlternative <-F
  }
  
  if(dim(x$SelectedMass)[1] != 0 &!x$frame){#
    try({
      
      print("Decision1")
      startpos <- grepl("^Start",colnames(x$SelectedMass))
      endpos <- grepl("^End",colnames(x$SelectedMass))
      
      apply(x$SelectedMass,1,function(x){
        abline(v=c(x[c(which(startpos),which(endpos))]),lty = "dotted",col = "blue")
        return(NULL)
      })
      rt <- RT_Extractor(x$SelectedMass)
      abline(v = rt ,lty = "dashed",col = "#00009980")
      text(rt,par()$usr[4],"predicted RT",cex = 0.6,col = "blue",srt = 90,pos = 2,xpd = NA)
      
      pr <- par()$usr
      yd<-pr[3]-diff(pr[3:4])*0.025
      # Peaks <<- Peaks
      points(x$Peaks,rep(yd,length(x$Peaks)),pch = 2,xpd = NA,col ="grey",cex = 1)
      
      try(points(x$PeaksSel,rep(yd,length(x$PeaksSel)),pch = 17,xpd = NA,col =2,cex = 1))
      
      
      legend(pr[2],pr[4],x$colmap[,1],col = x$colmap[,2],xpd = NA,cex = 0.5,pch = 20,lty = "solid",bty = "n")
      
      # Writing Rectangle of predicted Peak:
    })
    
  }else{
    print("Decision2")
    
    char <- dist(par()$usr[1:2])/dist(par("cxy"))
    rfi <- x$rf
    if(nchar(x$rf) > char){
      rfi <- c(substr(x$rf,1,char),substr(x$rf,char+1,nchar(x$rf)))
      
    }
    
    mtext(paste(rfi,collapse = "\n"),3,cex = 0.8,adj = 1,line = -length(rfi),xpd =T )
    # mtext(rf,3,cex = 0.8,adj = 1,line = -2)
    mtext(round(x$xl[1],1),2,col = "grey",line = -1)
    mtext(round(x$xl[2]),4,col = "grey",line = -1)
    
    
  }
}
#### Called by plot.TransplotList
RT_Extractor <- function(ilt){
  RTpredicted <- ilt$Start..min.+((ilt$End..min.-ilt$Start..min.)/2)
  RTpredicted[RTpredicted<0] <- 0
  return(RTpredicted)
}

### Protein Stuff
ProteinQuan <- function(PeaksTable,name="",ValuesAcrossSampleThreshold = 1,log10 = T){
  # PeaksTable <<- PeaksTable
  # name <<- name
  # ValuesAcrossSampleThreshold <<- ValuesAcrossSampleThreshold
  RT <- PeaksTable$Q1+(PeaksTable$Q1+PeaksTable$Q2)/2
  
  PeaksConvMaC <- PeaksTable[,1:which("Q1" == colnames(PeaksTable))-1]
  
  PeaksConvMaC[is.na(PeaksConvMaC)] <-0
  rows <- PeaksTable$rawfile
  
  # TransSel 
  SEL <- apply(PeaksConvMaC,2,function(x){length(x[x> 0])>=ValuesAcrossSampleThreshold})
  PeaksConvMaC <- PeaksConvMaC[,SEL]
  PeaksConvMaC <- apply(PeaksConvMaC,2,function(x){
    x <- as.numeric(x)
    x[x <= 0] <- NA
    if(log10){
      x <- log10(x)
    }
    if(length(x) == 2){
      x <- x/median(x,na.rm = T)
      
    }else{
      x <- (scale(x))
      
    }
    return(x)
    
  })
  
  PeaksConvMaCM <- apply(PeaksConvMaC,1,function(x){median(as.numeric(x),na.rm = T)})
  PeaksConvMaCsd <- apply(PeaksConvMaC,1,function(x){sd(as.numeric(x),na.rm = T)})
  
  names(PeaksConvMaCM) <- rows
  LEN <- apply(PeaksConvMaC,2,function(x){length(x[!is.na(x)])})
  
  PeaksRawF <- (list(M = PeaksConvMaCM,SD = PeaksConvMaCsd,coln = as.character(rows),all = PeaksConvMaC[,LEN >= ValuesAcrossSampleThreshold],name =  name,fdr=PeaksTable$SCAfdr,RT=RT))
  
  
  
  
  
  return(PeaksRawF)
}
plotPeptideList <- function(PeptideList,colmap = NULL,alpha = 0.01){
  # PeptideList <<-PeptideList
  # colmap <<- colmap
  for(tempi in 1:length(PeptideList)){
    tempd <- PeptideList[tempi]
    for(hmpi in  1:length(tempd)){
      # tempd <<- tempd
      M <- tempd[[hmpi]]$M
      S <- tempd[[hmpi]]$SD
      A <- tempd[[hmpi]]$all
      fdr <- tempd[[hmpi]]$fdr
      if(length(fdr) == 0){
        colvec = "turquoise4"
      }else{
        colvec = c("grey","turquoise4")[as.numeric(fdr<=alpha)+1]
      }
      hu <- barplot2(M,plot.ci = T,ci.l = M-S,ci.u = M+S,las = 1,main = "",las = 2,ylab = "relative Intensity",border = NA,col = colvec)
      # NameTags <- unlist(strsplit(names(PeptideList),"__"))
      # H1 <- paste("Seq:",NameTags[2])
      # H2 <- paste(gsub("mz","mz:",NameTags[1]),"Time:",names(tempd)[hmpi])
      # H3 <- paste("Info:",NameTags[3])
      # H4 <- paste("Sample:",NameTags[4])
      # mtext(paste(H1,H2,H3,H4,sep = "\n"),3,line = 2.5,adj = 0,cex = 0.8,outer = F)
      # mtext(paste(unlist(strsplit(names(PeptideList)[tempi],"_")),names(tempd)[hmpi],sep = "\n"),3,line = 2,adj = 0)
      
      if(dim(A)[2] > 1){
        Cols <- colmap[match(colnames(A),colmap[,1]),2]
        sapply(1:dim(A)[2],function(x){
          
          coll = gsub("..$","80",rainbow(dim(A)[2])[x])
          points(jitter(hu[,1]),A[,x],col = gsub("FF$","80",Cols[x]),type = "o",pch = 20,xpd = NA)
          
        })
      }
      
      
      
    }
  }
  return(hu)
}


# SQLITE helper functions ##########
## SQLite Read Table
dbread <- function(x,db = "./PickyAnalyzer.sqlite"){
  if(!file.exists(paste(dirname(db),"TableCompiled",sep = "/"))){
    return(NULL)
  }
  db <- dbConnect(SQLite(),db)
  temp <- dbReadTable(db,x)
  dbDisconnect(db)
  return(temp)
}
dbwrite <- function(x,name,db = "./PickyAnalyzer.sqlite",...){
  if(!file.exists(paste(dirname(db),"TableCompiled",sep = "/"))){
    return(NULL)
  }
  
  db <- dbConnect(SQLite(),db)
  temp <- dbWriteTable(db,name = name,value = x,...)
  dbDisconnect(db)
  return(temp)
}
## SQLite List Tables
dblistT<- function(db = "./PickyAnalyzer.sqlite"){
  fp <- paste(dirname(db),"TableCompiled",sep = "/")
  if(!file.exists(fp)){
    warning("Database not compiled")
    return(NULL)
  }
  
  db <- dbConnect(SQLite(),db)
  temp <- dbListTables(db)
  dbDisconnect(db)
  return(temp)
}
## SQLITE Quesry helper
dbqueryT<- function(db = "./PickyAnalyzer.sqlite",query = NULL){
  if(!file.exists(paste(dirname(db),"TableCompiled",sep = "/"))){
    print("Database not compiled")
    return(NULL)
  }
  
  db <- dbConnect(SQLite(),db)
  temp <- dbSendQuery(db,query)
  dbDisconnect(db)
  return(temp)
}




# ???? #############


## ?
FDRtoManyRaw <- function(subdat){
  FDR <- matrix(NA,dim(subdat)[1],2)
  for(rfi in unique(subdat$rawfile)){
    FDRTemp <- FDRCalc(subdat[selectedrf <-subdat$rawfile == rfi,])
    FDR[selectedrf,] <- FDRTemp
  }
  return(FDR)
}


# Retention Time Correction ########
## RT average HP corrected, input: maxquant run (evidence.txt)
RT_Correction_HP <- function(maxquant,span = 1,doCompareRT = F){
  print("Reading TABLE")
  evi <- fread(paste(maxquant,"evidence.txt",sep = "/"),sep = "\t",stringsAsFactors = F,data.table = F)
  cat("\r Calculating HP")
  S <- unique(evi$Sequence[evi$Modifications == "Unmodified"])
  H <- sapply(S,PeptideHydrophobicity)
  evi$HP <- H[match(evi$Sequence,S)]
  evi$HP2 <- NA
  RA <-  aggregate(evi$HP2,list(evi$`Raw file`),range,na.rm = T)
  diff <- apply(RA$x,1,diff)
  selected <- RA[diff == max(diff),1]
  
  evicount <- table(evi$`Raw file`[evi$Modifications == "Unmodified"])
  selected <- names(evicount)[max(evicount)==evicount]
  tempx <- evi[evi$`Raw file`==selected,]
  a <- tempx$`Retention time`
  b <- tempx$HP
  # plot(a,b)
  # TemplateHP <- loess(unfactor(b)~unfactor(a),span = span)
  ## points(TemplateHP$x,TemplateHP$fitted,col = 2,cex = 0.5)
  #Template <- loess(TemplateHP$x~TemplateHP$fitted)
  Template <- loess(unfactor(a)~unfactor(b),span = 1,control = loess.control(surface="direct"))
  cat("\r Correcting HP")
  print("Analyzing Raw-Files")
  
  pdf("RetentionTimeCorrection.pdf")
  # evi <<- evi
  evi$RT2 <- NA
  RT2 <<- c()
  Fits <- lapply(unique(evi$`Raw file`),function(x){
    # x <<- x
    cat("\r",x)
    tempx <- evi[evi$`Raw file`==x,]
    # o <- order(a)
    # a <- a[o]
    # b <- b[o]
    a <- tempx$`Retention time`-tempx$`Retention time calibration`
    b <- tempx$HP
    if(length(a) < 1000){
      span <- 1
    }
    hui <- loess(unfactor(b)~unfactor(a),span = span,loess.control(surface="direct"))
    HP2 <- predict(hui,a)
    RTn <- predict(Template,HP2)
    plot(a,b,pch = 20,cex = 0.5)
    points(hui$x,hui$fitted,type = "p",col = 2)
    PredictedRT <- predict(Template,hui$fitted)
    points(PredictedRT,hui$fitted,col = 4)
    # plot(RT_Correction_HP,a)
    # evi$HP2[evi$`Raw file`==x] <- hui$fitted
    RT2[evi$`Raw file`==x] <<- RTn
    # df <- data.frame(PredictedRT = PredictedRT,RT= hui$x,rawfile = x,seq = "Modified_Sequence")
    return(list(df = data.frame(PredictedRT = PredictedRT,RT= hui$x),md = hui))
    
  })
  dev.off()
  print("FINISHED LOOP")
  # plot(evi$`Retention time`,RT2)
  names(Fits) <- unique(evi$`Raw file`)
  evi$RT2 <- RT2
  if(doCompareRT){
    try(hm <- CompareRT(evi))
  }
  
  #system("open RetentionTimeCorrection.pdf")
  return(list(rf_Fits = Fits,Template = Template))
}
## RT average HP correction, called in RT_Correction_HP
CompareRT <- function(evi){
  TempEvi <- paste(evi$`Modified sequence`,evi$`Raw file`)
  TempEviU <- unique(TempEvi)
  TempEviU <- strsplitslot(TempEviU,1," ",fixed = T)
  hm <- table(TempEviU)
  sel <- evi$Modifications == "Unmodified"
  
  evix <- evi[sel,]
  evix$OldRT  <-  evix$`Retention time`-evix$`Retention time calibration`
  
  AGGFUN <- aggregate(subset(evix,select = c("Retention time","RT2","OldRT")),list(evix$`Modified sequence`),sd,na.rm = T)
  AGGFUNMean <- aggregate(subset(evix,select = c("Retention time")),list(evix$`Modified sequence`),mean,na.rm = T)
  AGGFUNMean <- AGGFUNMean[order(AGGFUNMean$`Retention time`),]
  AGGFUN <- AGGFUN[match(AGGFUNMean$Group.1,AGGFUN$Group.1),]
  # PreSelect <- sapply(unique(evix$`Modified sequence`),function(ms){
  #   ms <- ms
  #   tempx <- evix[evix$`Modified sequence` == ms,]
  #   tempx$OldRT <- tempx$`Retention time`-tempx$`Retention time calibration`
  #   # plot(tempx$`Retention time`,tempx$OldRT)
  #   Cr <- sd(tempx$`Retention time`,na.rm = T)
  #   Cr2 <- sd(tempx$RT2)
  #   Cr3 <- sd(tempx$OldRT)
  #   return(c(median(tempx$`Retention time`),RT = Cr,RT2 = Cr2,RT3 = Cr3))
  # })
  # PreSelect2 <- t(PreSelect)
  pdf("RT_Alignment.pdf",width = 10)
  hist(AGGFUN$`Retention time`,col = "#00990060",border = NA,breaks = 1000,xlab = "RT sd of Seq across Raw Files",main = "RT Alignment")
  hist(AGGFUN$RT2,add = T,col = "#99000060",border = NA,breaks = 1000)
  
  legend("topright",c("MQ RT sd","Hydrophobicity Corrected RT sd"),title = "sd",fill = c("#00990060","#99000060"))
  plot(AGGFUNMean$`Retention time`,AGGFUN$`Retention time`,type = "l",ylab = "sd RT",xlab = "time in min",col = "#00990060")
  points(AGGFUNMean$`Retention time`,AGGFUN$RT2,type = "l",col = "#99000060")
  RelEl <- round(AGGFUNMean$`Retention time`/max(AGGFUNMean$`Retention time`)*100)
  col <- colorRampPalette(c("blue","red"))(100)
  col <- paste(col,99,sep = "")
  plot(log2(AGGFUN$`Retention time`),log2(AGGFUN$RT2),pch = 20,cex = 0.5,col =col[RelEl] ,xlab = "MQ RT sd " ,ylab = "HP RT sd")
  grid()
  abline(0,1,h = 0,v = 0)
  legend("bottomright",legend = c("Begin","End"),title = "Gradient",fill = c("blue","red"))
  dev.off()
  system("open RT_Alignment.pdf")
  
  
}


# NEX'''#########







############## Peak Detection ################

## Called By Detect Peak
ThresholdingAlgo <- function(y,lag,threshold,influence) {
  #lag lag: the lag parameter determines how much your data will be smoothed and how adaptive the algorithm is to changes in the long-term average of the data. The more stationary your data is, the more lags you should include (this should improve the robustness of the algorithm). If your data contains time-varying trends, you should consider how quickly you want the algorithm to adapt to these trends. I.e., if you put lag at 10, it takes 10 'periods' before the algorithm's treshold is adjusted to any systematic changes in the long-term average. So choose the lag parameter based on the trending behavior of your data and how adaptive you want the algorithm to be.
  # #influence: this parameter determines the influence of signals on the algorithm's detection threshold. If put at 0, signals have no influence on the threshold, such that future signals are detected based on a threshold that is calculated with a mean and standard deviation that is not influenced by past signals. Another way to think about this is that if you put the influence at 0, you implicitly assume stationarity (i.e. no matter how many signals there are, the time series always returns to the same average over the long term). If this is not the case, you should put the influence parameter somewhere between 0 and 1, depending on the extent to which signals can systematically influence the time-varying trend of the data. E.g., if signals lead to a structural break of the long-term average of the time series, the influence parameter should be put high (close to 1) so the threshold can adjust to these changes quickly.
  # # threshold: the threshold parameter is the number of standard deviations from the moving mean above which the algorithm will classify a new datapoint as being a signal. For example, if a new datapoint is 4.0 standard deviations above the moving mean and the threshold parameter is set as 3.5, the algorithm will identify the datapoint as a signal. This parameter should be set based on how many signals you expect. For example, if your data is normally distributed, a threshold (or: z-score) of 3.5 corresponds to a signaling probability of 0.00047 (from this table), which implies that you expect a signal once every 2128 datapoints (1/0.00047). The threshold therefore directly influences how sensitive the algorithm is and thereby also how often the algorithm signals. Examine your own data and determine a sensible threshold that makes the algorithm signal when you want it to (some trial-and-error might be needed here to get to a good threshold for your purpose).
  ## Taken from https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data
  
  y <- y
  if(all(is.na(y))){
    print("Warning: ThresholdingAlgo, y object all(is.na(y)) == T.")
    y[is.na(y)] <- 0
  }
  y[y==0] <- abs(jitter(y[y==0],amount =100))
  #lagi <<- lag
  #threshold <<- threshold
  #influence <<- influence
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag],na.rm = T)
  stdFilter[lag] <- sd(y[0:lag],na.rm = T)
  stdFilter[is.na(stdFilter)] <- 0
  avgFilter[is.na(avgFilter)] <- 0
  if(lag>length(y)){
    lag <- 2
  }
  # if(length(y)<3){
  #   print("To less data points for ThresholdingAlgo")
  # }
  
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i],na.rm = T)
    stdFilter[i] <- sd(filteredY[(i-lag):i],na.rm = T)
  }
  LI <- list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter,lag = lag, threshold =threshold,influence = influence,y = y)
  class(LI) <- c("ThresholdingAlgo","list")
  return(LI)
}
## Called by DetectPeakWrapper
infoquantileParser <- function(infoquantile){
  # baef<<- infoquantile
  # infoquantile <- baef
  if(dim(infoquantile)[1] > 0){
    if(length(infoquantile$FDR_Choosen) == 0){
      infoquantile$FDR_Choosen <- infoquantile$SCAfdr
    }
    sel1 <- infoquantile$FDR_Choosen == min(infoquantile$FDR_Choosen,na.rm = T)
    sel1[is.na(sel1)] <- F
    infoquantile <- infoquantile[sel1,]
    infoquantile <- infoquantile[infoquantile$RF_Scores == max(infoquantile$RF_Scores,na.rm = T),]
    infoquantile <- infoquantile[infoquantile$SCAall == max(infoquantile$SCAall,na.rm = T),]
    infoquantile <- subset(infoquantile,select = c("RF_Scores","DL_Scores","FDR","FDR_RF_all","FDR_DL_all"))
    
    # infoquantile <- data.frame(RF_Scores=infoquantile$PEP[1],SCA = infoquantile$SCA[1],MatchCount = infoquantile$MatchCount[1],SCAfdr = infoquantile$SCAfdr,SCA_mz_fdr = infoquantile$SCA_mz_fdr)
  }else{
    infoquantile <- data.frame(RF_Scores=NA,DL_Scores = NA,FDR = NA,FDR_RF_all = NA,FDR_DL_all=NA)
    
  }
  return(infoquantile)
}
## Requantify Related, generic in principle, split Transitions and Info from Table
SplitTransitionInfo <- function(x,ColumnSplitName="charge"){
  Charge <- which(colnames(x)==ColumnSplitName)
  if(length(Charge)==1){
    x <-as.data.frame(x)
    Transitions <- x[,1:(Charge-1)]
    Info <- x[,(Charge):dim(x)[2]]
    if(length(Transitions$rawfile)>0){
      Info$rawfile <- Transitions$rawfile
      Transitions$rawfile <- NULL
    }
    liout <- list(Transitions=Transitions,Info=Info)
    
  }else{
    liout <-  list(Transitions=NULL,Info=NULL)
  }
  liout
}
## Detects Peak, Main Function
# DetectPeak <- function(rt_pep,Peakwidth,transitions,RT,lag = 5,threshold = 5,influence = 0.01,presetQuantiles = NULL,session = NULL){
#   rt_pep     <<- unique(rt_pep)
#   Peakwidth  <<- Peakwidth
#   transitions <<- transitions
#   RT.DetectPeak <<- RT
#   transitions.DetectPeak <- data.frame(transitions)
#   
#   ord_vec <- order(RT)
#   transitions <- transitions[ord_vec,]
#   transitions$rawfile <- NULL
#   RT <- RT[ord_vec]
#   transitions[transitions == -1] <- 0
#   rt_pep <- median(unique(rt_pep))
#   
#   
#   RT_range <- c(rt_pep-Peakwidth/2,rt_pep+Peakwidth/2)
#   xl<-RT_range
#   sel <- which(xl[1] <=RT&xl[2] >=RT)
#   Seli <- which(abs(RT-rt_pep) == min(abs(RT-rt_pep),na.rm = T))
#   
#   if(length(presetQuantiles) != 2&length(Seli) > 0){
# 
#     win <- sapply(1:dim(transitions)[2],function(xi){
#       x <- transitions[,xi]
#      
#       check <- try(xthresh <- ThresholdingAlgo(x,lag,threshold,influence)$signals)
#       if(class(check)=="try-error"){
#         x <- rep(0,length(x))
#       }else{
#         x <- xthresh
#         x[x<0] <- 0
#       }
#       x <- x
#       x[is.na(x)] <- 0
#       if(any(x[Seli]== 1,na.rm = T)){
#         tempx <- x[sel]
#         
#         po <- which(!is.na(match(sel,Seli)))
#         if(length(po) ==0){
#           return(rep(NA,2))
#         }
#         tempxRight <- tempx[min(po):length(tempx)]
#         tempxRightlimit <- min(which(tempxRight <=0),na.rm = T)+min(po)
#         if(is.infinite(tempxRightlimit)){
#           tempxRightlimit <- 1
#         }
#         
#         tempxLeft <- tempx[1:max(po)]
#         tempxLeftlimit <- max(which(tempxLeft <= 0))
#         if(is.infinite(tempxLeftlimit)){
#           tempxLeftlimit <- length(tempxLeft)
#         }
#         
#         ELUWIN <-RT[sel][c(tempxLeftlimit,tempxRightlimit)]
#         
#         # abline(v = ELUWIN,col = "#00990060")
#         return(ELUWIN)
#       }else{return(rep(NA,2))}
#     })
#     qua <- quantile(unlist(win),probs = c(0.1,0.9),na.rm = T)
#   }else{
#     qua <- presetQuantiles
#     win <- NA
#   }
#   wise <- RT>= min(qua,na.rm = T)&RT <=max(qua,na.rm = T)
#   if(!all(wise)){
#     if(length(session) !=0){
#       showNotification(session,"No spectra at this RT window (DetectPeak).")
#       
#     }
#   }
#   traPeak <- data.frame(transitions[wise,])
#   quaPeak <- RT[wise]
#   
#   XIC <- sapply(1:dim(traPeak)[2],function(x){trapz(quaPeak,traPeak[,x])})
#   names(XIC) <- colnames(traPeak)
#   INT <- apply(traPeak,2,max,na.rm = T)
#   return(list(quantile= qua,all  = win,XIC = XIC,intensity = INT))
# }
Neighbour_Mean_imputation <- function(rawda){
  rawdanew <- apply(rawda,1,function(xu){
    xu <- xu
    NAs <- which(!is.na(xu))
    if(any(is.na(xu))){
      repi <- sapply(1:length(xu),function(y){
        v1 <- xu[y]
        y <- y
        if(is.na(v1)){
          mi <- max(NAs[NAs<y])
          ma <- min(NAs[NAs>y])
          v1 <- mean(xu[c(mi,ma)],na.rm = T)
        }
        v1
        
      })
      
    }else{
      repi <- xu
    }
  })
  rawdanew <- t(rawdanew)
  rawda <- data.table(rawdanew)
  rawda
}
PeakDetection <- function(transitions,RTs,movavg_window=5,movavg_type="t",QualityPosition=NA,summary_bw=0.1,ZeroSmooth=T,Diff_Smooth=T,supersmooth_I=F,supersmooth_bw=0.2){
  transitionstemp <<- transitions
  # movavg_type <<- movavg_type
  # ZeroSmooth <<- ZeroSmooth
  # Diff_Smooth <<- Diff_Smooth
  # RTs <<- RTs
  # RTs <- RT.DetectPeak
  RTs <<- RTs
  plot(RTs,transitionstemp$a2,type="l")
  plot(RTs,transitionstemp$y8,type="l")
  
  RT.DetectPeak <- RTs 
  if(ZeroSmooth){
    transitionstemp <- apply(transitionstemp,2,function(x){
      x <- x
      xi <- x
      le <- sapply(1:length(x),function(i){
        p <- F
        try(p<- (x[i]==0&x[(i-1)]!=0)&(x[i]==0&x[(i+1)]!=0))
        if(length(p)==0){p <- F}
        p
      })
      
      x[le] <- NA
      x
    })
    transitionstemp <- t(transitionstemp)
    transitionstemp <- Neighbour_Mean_imputation(transitionstemp)
  }
  
  if(Diff_Smooth){
    transitionstemp <- apply(transitionstemp,1,function(x){
      x[abs(diff(x))<100] <- NA
      x
    })
    transitionstemp <- t(transitionstemp)
    transitionstemp <- Neighbour_Mean_imputation(transitionstemp)
  }
  temp <- lapply(1:dim(transitionstemp)[2],function(counter){
    x <- as.data.frame(transitionstemp)[,counter]
    
    # plot(RTs,x,type="l")
    # try({
    #   if(sd(diff(diff(x[!is.na(x)],na.rm = T)))>10000){
    #     movavg_window <- 10
    #     print("movavg_window",movavg_window)
    #     # x[x==0] <- NA
    #     
    #   }
    # })
    # 
    # xm<- movavg(x,movavg_window,movavg_type,n = 1)
    # xm <- smooth(x,c("3RS3R", "3RSS", "3RSR", "3R", "3", "S")[2],twiceit = T)
    # rt_pep
    if(supersmooth_I){
      xm <- supsmu(RTs,x,span=supersmooth_bw)$y
      
    }else{
      xm <- x
    }
    # points(xm,typeC="l",col=2,type="l")
    te <- data.table(findpeaks(xm,nups = 1,ndowns = 1,minpeakheight = 1000))
  })
  
  # if findpeaks with nups and ndowns is too conservative
  if(all(lengths(temp)==0)){
    temp <- lapply(1:dim(transitionstemp)[2],function(counter){
      x <- as.data.frame(transitionstemp)[,counter]
      
      
      if(supersmooth_I){
        xm <- supsmu(RTs,x,span=supersmooth_bw)$y
        
      }else{
        xm <- x
      }
      # points(xm,typeC="l",col=2,type="l")
      te <- data.table(findpeaks(xm,nups = 0,ndowns = 0,minpeakheight = 1000))
    })
    
  }
  unique(temp)
  temp <- lapply(1:length(temp),function(i){
    temptemp <- temp[[i]]
    if(dim(temptemp)[1]>0){
      temptemp$Counter <- i
    }else{return(NULL)}
    temptemp
  })
  temp <- rbindlist(temp)
  if(length(temp)==0){
    tfu <- data.table(V2=NA,s=NA,start=NA,end=NA)
  }else{
    tfu <- temp[,{
      data.table(s=sum(V1/(Counter),na.rm = T),
                 start=(quantile(RTs[V3],probs=0.1,na.rm = T)),
                 end=(quantile(RTs[V4],probs=0.9,na.rm = T)))
      
    },.(V2=round(RTs[V2],1))]
  }
  
  
  
  
  PeakSummary <- NA
  d <- NA
  list(Peaks=tfu,PeakSummary=PeakSummary,PeakDensity=d)
}

DetectPeak <- function(rt_pep,Peakwidth,transitions,RT,
                       scores=NULL,      
                       lag = 5,threshold = 5,influence = 0.01,#please check if obsolete
                       presetQuantiles = NULL,# defines manual set borders
                       session = NULL,# for shiny
                       MaxPeakWidth = 1/60*30, # for manual cut of detected Peaks
                       MinPeakWidth=1/60*10,# if Peaks ar too small, they will be symmetrically widened to this window
                       supersmooth_I_set =F, # apply smoother before detection
                       supersmooth_bw_set=0.1, 
                       ApplyMaximumWidth=T,
                       FDR=NULL,
                       alpha=0.01
){
  rt_pep     <<- unique(rt_pep)
  Peakwidth  <<- Peakwidth
  transitions <<- transitions
  scores <<- scores
  RT <<- RT
  RT.DetectPeak <- RT
  # putting in RT order:
  ord_vec <- order(RT)
  transitions <- transitions[ord_vec,]
  RT <- RT[ord_vec]
  RT.DetectPeak <- RT.DetectPeak[ord_vec]
  if(length(scores)==length(ord_vec)){
    scores <- scores[ord_vec]
  }
  if(length(FDR)==length(ord_vec)){
    FDR <- FDR[ord_vec]
  }
  
  
  
  transitions[transitions == -1] <- 0
  transitions$rawfile <- NULL
  # print(paste(Sys.time(),"Running Detect Peak"))
  
  
  
  if(0){
    transitions.DetectPeak <- data.frame(transitions)
    transitions.DetectPeak_temp <- transitions.DetectPeak
    transitions.DetectPeak_temp$RT <- RT.DetectPeak
    transitions.DetectPeak_temp_long <- melt(transitions.DetectPeak_temp,id.vars="RT")
    g <- ggplot(transitions.DetectPeak_temp_long,aes(x=RT,y=value,group=variable,color=variable))+geom_line()
    # Smoothing:
    
  }
  
  
  rt_pep <- median(unique(rt_pep))
  
  
  RT_range <- c(rt_pep-Peakwidth/2,rt_pep+Peakwidth/2)
  xl<-RT_range
  sel <- which(xl[1] <=RT&xl[2] >=RT)
  Seli <- which(abs(RT-rt_pep) == min(abs(RT-rt_pep),na.rm = T))
  
  if(length(presetQuantiles) != 2&length(Seli) > 0){
    rm("Peaks")
    # print("Peak Detection")
    Peaks <<- PeakDetection(transitions,RT.DetectPeak,movavg_window=5,movavg_type="t",ZeroSmooth = F,Diff_Smooth = F,supersmooth_I = supersmooth_I_set,supersmooth_bw = supersmooth_bw_set)
    # Peaks <<- PeakDetection(transitions,RT.DetectPeak,movavg_window=5,movavg_type="t",ZeroSmooth = F,Diff_Smooth = F,supersmooth_I = F)
    # scores <<- scores
    # plot(RT.DetectPeak,transitionstemp$)
    # RT.DetectPeak <<- RT.DetectPeak
    # FDR <<- FDR
    # Adding ScoreInformation:
    ######
    # plot(RT.DetectPeak,scores)
    PS <- Peaks$Peaks
    if(length(PS$s)>0){
      PS$MaxIntensity <- PS$s
      PS$peak <- PS$V2
    }
    # adding FDR  and score information:
    PS[,score:=.({
      gr <- .BY
      sel <- RT.DetectPeak>=gr$start&RT.DetectPeak<=gr$end
      if(sum(sel,na.rm = T)==0){
        o <- (0)
      }else{
        sc <- sort(scores[sel],decreasing=T)
        if(length(sc)>3){
          sc <- sc[1:3]
        }
        # o <- sum(sc,na.rm = T)/length(sc)
        o <- max(sc,na.rm = T)
      }
      o
    }),.(start,end,MaxIntensity,peak)]
    if(length(FDR)>0){
      PS[,c("FDRaverage","FDRmin"):={
        gr <- .BY
        sel <- RT.DetectPeak>=gr$start&RT.DetectPeak<=gr$end
        if(sum(sel,na.rm = T)==0){
          o <- (1)
          o_min <- 1
        }else{
          sc <- sort(FDR[sel],decreasing=T)
          # if(length(sc)>3){
          #   sc <- sc[1:3]
          # }
          o <- sum(sc,na.rm = T)/length(sc)
          o_min <- min(FDR,na.rm = T)
        }
        list(o,o_min)
      },.(start,end,MaxIntensity,peak)]
      
    }
    
    
    # PS$RTdiff <- abs(PS$peak-rt_pep)
    q1 <- quantile(PS$score,0.8,na.rm = T)
    q2 <- quantile(PS$MaxIntensity,0.8,na.rm = T)
    q3 <- quantile(PS$RTdiff,0.2,na.rm = T)
    
    # step 1 filter on best included Scores:
    PS <- PS[score>=q1]
    ObligatoricFDR <- T
    if(length(FDR)>0){
      if(ObligatoricFDR&length(presetQuantiles)==0){
        PS <- PS[FDRmin<alpha]
        
      }
    }
    
    qset <- 0.8
    while(diff(range(PS$score))>0.1){
      qset <- qset+0.1
      q1 <- quantile(PS$score,0.8,na.rm = T)
      PS <- PS[score>=q1]
    }
    
    # step on 
    # step 2 filter on best Intensities:
    if(dim(PS)[1]>1&any(PS$MaxIntensity>=q2)){
      PS <- PS[MaxIntensity>=q2]
    }
    # step 3 filter on most close Match :
    # qua  <- PS[RTdiff<=q3]
    qua <- PS
    qua <- qua[MaxIntensity==max(MaxIntensity,na.rm = T)]
    qua <- qua[,.(quantile(start,0.1,na.rm = T),quantile(end,0.9,na.rm = T),median(peak,na.rm = T))]
    qua <- unlist(qua[1,])[1:3]
    
    if(ApplyMaximumWidth){
      PeakBoundaries_Algo <- "CutOff"
      if(PeakBoundaries_Algo=="CutOff"){
        # qua <<- qua
        if(all(!is.na(qua))){
          if(diff(qua[1:2])>MaxPeakWidth){
            # checking to the left
            LeftSize <- abs(diff(qua[c(1,3)]))
            if(LeftSize>MaxPeakWidth/2){
              qua[1] <- qua[3]-MaxPeakWidth/2
            }
            # checking to the right
            RightSize <- abs(diff(qua[c(2,3)]))
            if(RightSize>MaxPeakWidth/2){
              qua[2] <- qua[3]+MaxPeakWidth/2
            }
            
          }
          
        }
        #   qua <- qua[,.(quantile(start,0.1),quantile(end,0.9),median(peak,na.rm = T))]
        
      }
      if(all(!is.na(qua))){
        if(diff(range(qua[1:2]))<MinPeakWidth){
          if(abs(diff(range(qua[c(1,3)])))<MinPeakWidth/2){
            qua[1] <- qua[3]-MinPeakWidth/2
          }
          if(abs(diff(range(qua[c(3,2)])))<MinPeakWidth/2){
            qua[2] <- qua[3]+MinPeakWidth/2
          }
        }
        
      }
    }
    
    
    
    # if(PeakBoundaries_Algo=="lm"){
    #   lm_border_leftright<- function(hum,qua,transitions){
    #     # qua is a vector with three values: leftborder, rightborder and peak
    #     left_lm_data <- hum[hum$RT>=qua[1]&hum$RT<=qua[3],]
    #     right_lm_data <- hum[hum$RT>=qua[3]&hum$RT<=qua[2],]
    #     
    #     co <- coefficients(lm(left_lm_data$value~left_lm_data$RT))
    #     if(length(co)==2){
    #       borderLeft <- (0-co[1])/co[2]
    #     }else{
    #       borderLeft <- qua[1]
    #     }
    #     co <- coefficients(lm(right_lm_data$value~right_lm_data$RT))
    #     if(length(co)==2){
    #       borderRight <- (0-co[1])/co[2]
    #     }else{
    #       borderRight <- qua[1]
    #     }
    #     
    #     lm_border <- function(right_lm_data,type="right",transitions){
    #       lmTemp <- lapply(unique(right_lm_data$RT),function(rti){
    #         rti <- rti
    #         if(type=="right"){
    #           right_lm_data_temp <- right_lm_data[right_lm_data$RT<=rti,]
    #           
    #         }else{
    #           right_lm_data_temp <- right_lm_data[right_lm_data$RT>=rti,]
    #           
    #         }
    #         co <- coefficients(lmfun <- lm(right_lm_data_temp$value~right_lm_data_temp$RT))
    #         if(length(co)==2){
    #           borderRight <- (0-co[1])/co[2]
    #         }else{
    #           borderRight <- NA
    #         }
    #         list(coefficients=co,borderRight=borderRight,residuals = summary(lmfun)$residuals)
    #       })
    #       resicheck <- (sapply(lmTemp,function(x){c(sum(abs(x$residuals))/max(transitions),length(x$residuals))}))
    #       # plot(unique(right_lm_data$RT),resicheck[1,])
    #       tempqua <- sapply(lmTemp[which(resicheck[1,]<quantile(resicheck[1,],0.2))],function(x){x$borderRight})
    #       borderRight <- median(tempqua,na.rm = T)
    #       borderRight
    #     }
    #     borderRight_lm <- lm_border(right_lm_data,"right",hum$value)
    #     borderLeft_lm <- lm_border(left_lm_data,"left",hum$value)
    #     if(is.na(borderRight_lm)){
    #       borderRight_lm <- borderRight
    #     }
    #     if(is.na(borderLeft_lm)){
    #       borderLeft_lm <- borderLeft
    #     }
    #     list(qua=c(borderLeft_lm,borderRight_lm))
    #   }
    #   
    #   transitionsTemp <- transitions
    #   transitionsTemp$RT <- RT.DetectPeak
    #   hum <- melt(transitionsTemp,id.vars="RT")
    #   qua <- lm_border_leftright(hum,qua)$qua
    # }
    win <- NA
    
    
    
    
    
  }else{
    m <- list("Manual")
    qua <- presetQuantiles
    win <- NA
    Peaks <- NA
  }
  wise <- RT>= min(qua,na.rm = T)&RT <=max(qua,na.rm = T)
  if(!all(wise)){
    if(length(session) !=0){
      showNotification(session,"No spectra at this RT window (DetectPeak).")
      
    }
  }
  traPeak <- data.frame(transitions[wise,])
  quaPeak <- RT[wise]
  
  XIC <- sapply(1:dim(traPeak)[2],function(x){trapz(quaPeak,traPeak[,x])})
  names(XIC) <- colnames(traPeak)
  INT <- apply(traPeak,2,max,na.rm = T)
  return(list(quantile= qua,all  = win,XIC = XIC,intensity = INT,Peaks=Peaks))
}

DetectPeak_v2_old <- function(rt_pep,Peakwidth,transitions,RT,estimatedPeakDistance=0.5,BootstrapRounds=20,minpeakheight=0.5,...){
  save(rt_pep,Peakwidth,transitions,RT,estimatedPeakDistance,BootstrapRounds,minpeakheight,file="DetectPeak_V2_TempNameSpace.rda")
  RT.DetectPeak <- RT
  
  #   # estimatedPeakDistance <- 0.5
  # BootstrapRounds <- 20
  # minpeakheight <- 0.5
  
  transitions.DetectPeak <- data.frame(transitions)
  transitions.DetectPeak <- apply(transitions.DetectPeak,2,function(x){scale(x)+abs(min(scale(x)))})
  plot(transitions.DetectPeak[,1])
  
  # option: correlation cutoff
  # apply(cor(transitions.DetectPeak_smooth),2,median) > 0.5
  # apply(cor(transitions.DetectPeak),2,median) > 0.5
  transitions.DetectPeak_smooth <- apply(transitions.DetectPeak,2,savgol,fl=151,forder=4,dorder=0)
  transitions.DetectPeak_smooth_strong <- apply(transitions.DetectPeak,2,savgol,fl=501,forder=4,dorder=0)
  # transitions.DetectPeak_smooth_strong[transitions.DetectPeak_smooth_strong<0] <- 0
  # transitions.DetectPeak_smooth[transitions.DetectPeak_smooth<0] <- 0
  ZeroRemove <- colSums(transitions.DetectPeak_smooth_strong,na.rm = T)
  transitions.DetectPeak_smooth_strong <- transitions.DetectPeak_smooth_strong[,ZeroRemove>0]
  
  iti<-0
  trendPeak_average <- NULL
  trend <- apply(transitions.DetectPeak_smooth_strong[,apply(cor(transitions.DetectPeak_smooth_strong),2,median) > 0.5],1,median,na.rm = T)
  
  trend2 <- c(rep(min(trend,na.rm = T)*0.9,9),min(trend,na.rm = T)*0.95,trend,min(trend,na.rm = T)*0.95,rep(min(trend,na.rm = T)*0.9,9))
  thresholdVal <- 0.1
  while(length(trendPeak_average)==0){
    iti <- iti+1
    
    # trend <- trend+abs(min(trend,na.rm = T))+0.001
    
    # plot(trend2)
    # trendPeak_average <- findpeaks(,nups = 1,ndowns = 1,minpeakheight = 0.5,minpeakdistance = 3,threshold = 0.1)
    
    trendPeak_average <- findpeaks(trend2,nups = 1,ndowns = 1,minpeakheight = minpeakheight,minpeakdistance = 2,threshold = thresholdVal)
    thresholdVal <- thresholdVal*0.75
    minpeakheight <- minpeakheight*0.9
    if(length(trendPeak_average)>0){
      trendPeak_average[,2:4]  <- trendPeak_average[,2:4]  - 10
      trendPeak_average[trendPeak_average<=0] <- 1
      trendPeak_average[trendPeak_average>length(RT.DetectPeak)] <- length(RT.DetectPeak)
      
    }
    if(iti>50){
      trendPeak <- matrix(NA,1,4)
      break()
    }
  }
  # stop()
  if(length(trendPeak_average)==0){
    BootStrapPeak <- lapply(1:BootstrapRounds,function(x){
      cat("\r",x)
      # transitions.DetectPeak_smooth <- apply(transitions.DetectPeak,2,savgol,fl=151,forder=4,dorder=0)
      # transitions.DetectPeak_smooth[transitions.DetectPeak_smooth<0] <- 0
      Sampled <- transitions.DetectPeak_smooth[,sample(1:dim(transitions.DetectPeak_smooth)[2],replace = T)]
      trend <- apply(Sampled[,apply(cor(Sampled),2,median) > 0.5],1,median,na.rm = T)
      thresholdVal <- 0.1
      trendPeak <- findpeaks(trend,nups = 1,ndowns = 1,minpeakheight = minpeakheight,minpeakdistance = 2,threshold = thresholdVal)
      iti<-0
      while(length(trendPeak)==0){
        iti <- iti+1
        trend2 <- c(rep(0,10),trend,rep(0,10))
        thresholdVal <- thresholdVal*0.75
        trendPeak <- findpeaks(trend2,nups = 0,ndowns = 1,minpeakheight = minpeakheight,minpeakdistance = 2,threshold = thresholdVal)
        if(length(trendPeak)>0){
          trendPeak[,2:4] <- trendPeak[,2:4]-10
          trendPeak[trendPeak<=0] <- 1
        }
        if(iti>50){
          trendPeak <- matrix(NA,1,4)
          break()
        }
      }
      
      
      colnames(trendPeak) <- c("Height","Peak","Start","End")
      data.table(trendPeak)
    })
    BootStrapPeak2 <- lapply(1:BootstrapRounds,function(x){
      cat("\r",x)
      # transitions.DetectPeak_smooth <- apply(transitions.DetectPeak,2,savgol,fl=151,forder=4,dorder=0)
      # transitions.DetectPeak_smooth[transitions.DetectPeak_smooth<0] <- 0
      Sampled <- transitions.DetectPeak_smooth_strong[,sample(1:dim(transitions.DetectPeak_smooth_strong)[2],replace = T)]
      trend <- apply(Sampled[,apply(cor(Sampled),2,median) > 0.5],1,median,na.rm = T)
      thresholdVal <- 0.1
      trendPeak <- findpeaks(trend,nups = 1,ndowns = 1,minpeakheight = minpeakheight,minpeakdistance = 2,threshold = thresholdVal)
      iti<-0
      while(length(trendPeak)==0){
        iti <- iti+1
        trend2 <- c(rep(0,10),trend,rep(0,10))
        thresholdVal <- thresholdVal*0.75
        trendPeak <- findpeaks(trend2,nups = 0,ndowns = 1,minpeakheight = minpeakheight,minpeakdistance = 2,threshold = thresholdVal)
        if(length(trendPeak)>0){
          trendPeak[,2:4] <- trendPeak[,2:4]-10
          trendPeak[trendPeak<=0] <- 1
        }
        if(iti>50){
          trendPeak <- matrix(NA,1,4)
          break()
        }
      }
      
      
      colnames(trendPeak) <- c("Height","Peak","Start","End")
      data.table(trendPeak)
    })
    
    # apply(BootStrapPeak,2,range)
    print("hu")
    BootStrapPeak <- rbindlist(BootStrapPeak)
    BootStrapPeak2 <- rbindlist(BootStrapPeak2)
    BootStrapPeak <- rbind(BootStrapPeak,BootStrapPeak2)
    BootStrapPeak <- BootStrapPeak[!is.na(BootStrapPeak$Height)]
  }else{
    BootStrapPeak <- data.table(trendPeak_average)
    colnames(BootStrapPeak) <-     c("Height","Peak","Start","End")
    
  }
  # BootStrapPeak$Peak <- RT.DetectPeak[BootStrapPeak$Peak]
  # BootStrapPeak$Start <- RT.DetectPeak[BootStrapPeak$Start]
  # BootStrapPeak$End <- RT.DetectPeak[BoCotStrapPeak$End]
  # dim_return <- NULL
  BootStrapPeak$PeakDiff <-abs(RT.DetectPeak[BootStrapPeak$Peak]-rt_pep)
  print("hu")
  qcut <- 0
  while(min(BootStrapPeak$PeakDiff)>estimatedPeakDistance){
    qcut <- qcut+0.1
    cat("\r",qcut)
    estimatedPeakDistance <- quantile(BootStrapPeak$PeakDiff,qcut)
  }
  BootStrapPeak <- BootStrapPeak[BootStrapPeak$PeakDiff <= estimatedPeakDistance,]
  
  print("hu")
  
  # hist(BootStrapPeak)
  
  # PEAKSFUN <- BootStrapPeak[PeakDiff<estimatedPeakDistance,{
  #   te <<- .SD
  #   te$Height <- NULL
  #   te_quantiles <- lapply(1:dim(te)[2],function(x){
  #     list(c(quantile(as.data.frame(te)[,x],probs=c(0.1,0.5,0.9),na.rm = T),sd(as.data.frame(te)[,x],na.rm = T)))
  #   })
  #   te_quantiles_c <- rbindlist(te_quantiles)
  #   te_quantiles_c <-t(te_quantiles_c)
  #   colnames(te_quantiles_c) <- sapply(colnames(te),function(x){paste(x,c(0.1,0.5,0.9,"SD"),sep = "_")})
  #   data.table(te_quantiles_c)
  # },round(Peak,1)]
  BootStrapPeakResults <- apply(BootStrapPeak,2,quantile,probs=c(0.1,0.5,0.9))
  BootStrapPeakResults_vec <- as.vector(BootStrapPeakResults)
  names(BootStrapPeakResults_vec) <- sapply(colnames(BootStrapPeakResults),function(x){paste(x,rownames(BootStrapPeakResults),sep = "_")})
  BootStrapPeakResults_vec <- as.data.frame(t(as.data.frame(BootStrapPeakResults_vec)))
  
  # quantify
  wise <- c(NA,NA)
  try(wise <- round(BootStrapPeakResults_vec$`Start_10%`:BootStrapPeakResults_vec$`End_90%`))
  
  
  pl <- F
  if(pl){
    gs <- lapply(list((transitions.DetectPeak_smooth),transitions.DetectPeak),function(input){
      transitions.DetectPeak_temp <<- data.frame(input)
      transitions.DetectPeak_temp$RT <- RT.DetectPeak
      transitions.DetectPeak_temp_long <- melt(as.data.table(transitions.DetectPeak_temp),id.vars="RT")
      library(ggplot2)
      library(pracma)
      g <- ggplot(transitions.DetectPeak_temp_long,aes(x=RT,y=value,group=variable,color=variable))+geom_line()
      
      g
    })
  }
  
  # medianfun <- apply(BootStrapPeak[BootStrapPeak$PeakDiff<estimatedPeakDistance,],2,median)
  # TempMeanSmooth <- savgol(TempMean,11,forder= 4,dorder=0)
  # Smoothing:
  # library(gridExtra)
  # grid.arrange(gs[[1]]+xlim(c(5,8.5))+geom_vline(xintercept = c((BootStrapPeak$Start),BootStrapPeak$Peak,(BootStrapPeak$End))))
  # BootStrapPeak <<- BootStrapPeak
  
  
  traPeak <- data.frame(transitions[wise,])
  quaPeak <- RT.DetectPeak[wise]
  
  XIC <- sapply(1:dim(traPeak)[2],function(x){trapz(quaPeak,traPeak[,x])})
  names(XIC) <- colnames(traPeak)
  INT <- apply(traPeak,2,max,na.rm = T)
  
  colnames(BootStrapPeak) <- c("V2","s","start","end","PeakDiff")
  Peaks <- list(Peaks,PeaksSummary=NA,PeakDensity=NA)
  
  
  
  gs[[1]]+geom_vline(xintercept = RT.DetectPeak[range(wise)])+geom_vline(xintercept = RT.DetectPeak[trendPeak_average[1,3:4]],color = "blue")
  
  # cuti <- wise
  
  return(list(quantile= RT.DetectPeak[range(wise)],all  = NA,XIC = XIC,intensity = INT,Peaks=Peaks))
  
  
}

DetectPeak_v2_v1 <- function(rt_pep,Peakwidth,transitions,RT,estimatedPeakDistance=1,RelaxationRounds=20,BootstrapRounds=20,minpeakheight=0.1,pl=F,smoothingFactor=51,presetQuantiles=NULL,ApplyMaximumWidth=T,Score = NULL,Identifier=NULL,...){
  RT.DetectPeak <- RT
  # load("/Users/henno/Documents/Skripte/R-Functions/Selbach-Functions/DataAnalysis/20220209_MVB_Rawfiles_library_bio_benchmarking/DetectPeak_V2_TempNameSpace.rda")
  save(rt_pep,Peakwidth,transitions,smoothingFactor,
       RT.DetectPeak,estimatedPeakDistance,BootstrapRounds,
       minpeakheight,presetQuantiles,ApplyMaximumWidth,Score,
       file=paste("DetectPeak_V2_TempNameSpace",Identifier,".rda",sep = ""))
  # scaling data
  
    trendPeak_average <- NULL
    
    # scaling data
    # stop()
    thresholdVal <- 0.1
    if(dim(transitions)[1]>5){
      if(length(presetQuantiles)!=2){
        
        
        # transitions <- Neighbour_Mean_imputation(transitions)
        class(transitions) <- "data.frame"
        transitions.DetectPeak <- as.data.frame(transitions)
        ZeroRemove <- colSums(transitions,na.rm = T)
        
        transitions.DetectPeak <- apply(transitions.DetectPeak,2,function(x){scale(x)})
        # attempt to filter based on summed absolutec differences across each trace, idea is, that fluctating data is causing higher values
        diffs <- apply(transitions.DetectPeak,2,function(x){sum(abs(diff(x[!is.na(x)])))})
        diffs <- diffs/dim(transitions.DetectPeak)[1]
        CorInit <- F
        corThresh <- 0.75
        presetwarning <- options()$warn
        options(warn=-1)
        while(sum(CorInit)<3){
          try({
            CorInit1 <- sapply(1:dim(transitions)[2],function(x){cor(transitions[,1],transitions[,x])})>corThresh
            CorInit2 <- sapply(1:dim(transitions)[2],function(x){cor(transitions[,2],transitions[,x])})>corThresh
            CorInit3 <- sapply(1:dim(transitions)[2],function(x){cor(transitions[,3],transitions[,x])})>corThresh
            
            CorInit1[is.na(CorInit1)]  <- FALSE
            CorInit2[is.na(CorInit2)]  <- FALSE
            CorInit3[is.na(CorInit3)]  <- FALSE
            CorFinal <- CorInit1|CorInit2|CorInit3
          },silent = T)
          
          # CorInit <- apply(cor(transitions[,ZeroRemove>0]),2,median,na.rm = T)>corThresh
          corThresh <- corThresh-0.01
          if(corThresh<0.2){
            break()
          }
        }
        options(warn=presetwarning)
        
        # CorFinal <- rep(F,rep(length(ZeroRemove)))
        # CorFinal[ZeroRemove>0] <- CorInit
        # if(any(transitions.DetectPeak)<0){
        #   transitions.DetectPeak <- transitions.DetectPeak+abs(min(transitions.DetectPeak))
        # }
        # smoothing strong or weak
        if(1){
          smoothingFactor <- 11
          minpeakheight <- 0.1
          # estimatedPeakDistance <- 1
        }
        repeatit <- T
        it <- 0
        thresholdVal <- 0.1
        sumvecCutoff <- 0.3
        diffThreshold <- 1# switched off, set to 0.1 to take action. It can cause problems
        
        RTdiffCollect <- c()
        trendraw <- apply(transitions.DetectPeak,1,function(x){M <- NA;try({M <- median(x[x!=0],na.rm = T)});M})
        trendPeak_trendPeak_average_raw <- findpeaks(trendraw,nups = 1,ndowns = 1,minpeakheight = 0,minpeakdistance = 1,threshold = 0)
        RTdiff <- abs(rt_pep-RT.DetectPeak[trendPeak_trendPeak_average_raw[,2]])
        trendPeak_trendPeak_average_raw <- as.data.frame(trendPeak_trendPeak_average_raw)
        trendPeak_trendPeak_average_raw <- trendPeak_trendPeak_average_raw[RTdiff<= estimatedPeakDistance,]
        # plot(RT.DetectPeak,transitions.DetectPeak[,1],type = "l")
        # points(RT.DetectPeak,Score,col = 2)
        # abline(v=rt_pep,col = 3)
        # points(RT.DetectPeak,transitions.DetectPeak[,1],col = 4,type = "l")
        # trendPeak_trendPeak_average_raw <- trendPeak_trendPeak_average_raw[trendPeak_trendPeak_average_raw[,1]>quantile(range(trendraw),0.1),]
        if(dim(trendPeak_trendPeak_average_raw)[1]!=1){
          while(repeatit){
            it <- it+1
            cat("\rRound",it)
            transitions.DetectPeak_smooth <- apply(transitions.DetectPeak,2,function(trend){
              trend <- trend
              trend2 <- c(rep(NA,10),trend,rep(NA,10))
              trend2raw <- trend2
              prwarnsettings <- options()$warn
              options(warn=-1)
              try({
                
                sdfun <- sd(trend[trend<=quantile(trend,0.3,na.rm = T)])
                if(is.na(sdfun)){
                  sdfun <- 0.1
                }
                if(sdfun==0){
                  sdfun <- 0.1
                }
                
                adv <-c( jitter(rep(min(trend,na.rm = T)-diff(range(trend,na.rm = T))*0.02,9),amount = ),min(trend,na.rm = T)-diff(range(trend,na.rm = T))*0.01)
                
                trend2 <- c(adv,trend,rev(adv))
                if(!all(is.na(trend2))){
                  trend2[is.na(trend2)] <- min(trend2,na.rm = T)
                  
                }
              },silent =T)
              
              try({
                
                trend2 <- savgol(trend2,fl=smoothingFactor,forder=4,dorder=0)
                if(any(is.na(trend2))){
                  trend2[is.na(trend2)] <- jitter(min(trend))
                }
                trend2[trend2<min(trend)] <- jitter(min(trend))
                
              },silent = T)
              options(warn=prwarnsettings)
              # plot(trend2raw/trend2)
              # par(new=T)
              # plot(trend2raw,type = "l")
              # points(trend2,col = 2,type = "l")
              trend2
            })
            
            transitions.DetectPeak_smooth <- apply(transitions.DetectPeak_smooth,2,scale)
            transitions.DetectPeak_smooth[is.infinite(transitions.DetectPeak_smooth)] <- NA
            # plot(transitions.DetectPeak_smooth[,2],type = "l")
            # points(transitions.DetectPeak[,2],col = 2,type = "l")
            
            # transitions.DetectPeak_smooth[transitions.DetectPeak_smooth<min(transitions.DetectPeak,na.rm = T)] <- min(transitions.DetectPeak,na.rm = T)
            
            transitions.DetectPeak_smooth <- as.data.frame(transitions.DetectPeak_smooth)
            transitions.DetectPeak_smooth <- transitions.DetectPeak_smooth[,ZeroRemove>0&CorFinal]
            if(is.vector(transitions.DetectPeak_smooth)){
              transitions.DetectPeak_smooth <- (data.frame(nclue=transitions.DetectPeak_smooth))
            }
            
            iti<-0
            #Correlation as a potential quality metric
            # corres <- NA
            # try(    corres <- apply(cor(transitions.DetectPeak_smooth),2,median))
            
            RTdiff <- 0
            trendPeak_average <- NULL
            
            try({
              # Smoothing data
              # if(!is.matrix(transitions.DetectPeak_smooth)){
              #   save(rt_pep,Peakwidth,transitions,RT.DetectPeak,estimatedPeakDistance,BootstrapRounds,minpeakheight,presetQuantiles,ApplyMaximumWidth,file="DetectPeak_V2_TempNameSpace_Error.rda")
              #   
              #   stop()
              # }
              trend <- apply(transitions.DetectPeak_smooth,1,function(x){M <- NA;try({M <- median(x[x!=0],na.rm = T)});M})
              sdfuncollect <- c(sdfuncollect,sumvec <- sum(abs(trendraw-trend[-c(c(1:10),(length(trend)-9):length(trend))]),na.rm = T)/length(trend))
              trend[trend<min(trendraw)] <- trendraw[trend<min(trendraw)]
              # plot(trendraw,type= "l")
              #add Trendbefore and After the data
              trend2 <- trend #c
              # points(trend[-(1:10)],type = "l",col = 2)
              #finding peaks
              trendPeak_average <- findpeaks(trend2,nups = 1,ndowns = 1,minpeakheight = minpeakheight,minpeakdistance = 2,threshold = thresholdVal)
              
              
              if(dim(trendPeak_average)[1]>0){
                trendPeak_average[,2:4] <- trendPeak_average[,2:4]-10
                trendPeak_average[trendPeak_average<=0] <- 1
                trendPeak_average[trendPeak_average>length(RT.DetectPeak)] <- length(RT.DetectPeak)
                if(length(trendPeak_average)>0){
                  correctedPosition <- trendPeak_average[,2]-10
                  correctedPosition[correctedPosition<=0] <- 1
                  correctedPosition[correctedPosition>length(RT.DetectPeak)] <- length(RT.DetectPeak)
                  RTdiff <- abs(rt_pep-RT.DetectPeak[correctedPosition])
                  if(length(RTdiff)==0){
                    RTdiff <- rep(0,dim(trendPeak_average)[1])
                  }
                }
              }else{
                trendPeak_average <- NULL
              }
              
              # RTdiffCollect <- c(RTdiffCollect,min(RTdiff))
              # Determining distance from best scoring RT position
              
            })
            
            if((length(trendPeak_average)==0|all(RTdiff>estimatedPeakDistance))&sumvec<sumvecCutoff){
              # less restrictions if nothing passed the threshholds, simultanteous approach could become problematic, but seems to be fine at the moment
              minpeakheight <- minpeakheight-0.01
              if(minpeakheight<0){
                minpeakheight <- 0
              }
              thresholdVal <- thresholdVal*0.9
              smoothingFactor <- smoothingFactor+4
              diffThreshold <- diffThreshold+0.02
            }else{
              #stopping loop
              repeatit <- F
              # Preparing Data
              if(length(trendPeak_average)>0){
                trendPeak_average <- data.frame(trendPeak_average)
                colnames(trendPeak_average) <- c("Height","Peak","Start","End")
                trendPeak_average$Diff <- abs(RTdiff)
                # Filtering for minimal distance
                
                trendPeak_averagetemp <- trendPeak_average[trendPeak_average$Diff <=estimatedPeakDistance*2,]
                if(length(Score)>0){
                  sc <-  Score[trendPeak_averagetemp$Peak]
                  trendPeak_averagetemp <- trendPeak_averagetemp[sc==max(sc),]
                }
                if(dim(trendPeak_averagetemp)[1]==0){
                  trendPeak_averagetemp <- trendPeak_average[trendPeak_average$Diff ==min(trendPeak_average$Diff),]
                  # trendPeak_averagetemp <- data.frame(Height=NA,Peak=NA,"Start"=NA,"End"=NA,Diff=NA)
                }
                # Replacing artificialy zscore height with real summed intensities
                trendPeak_averagetemp$Height <- rowSums(transitions[trendPeak_averagetemp$Peak,])
                # Final Filtering on summed intensities
                trendPeak_average <- trendPeak_averagetemp[trendPeak_averagetemp$Height==max(trendPeak_averagetemp$Height),]#[1,]#apply(trendPeak_averagetemp,2,median)
                
              }
              
            }
            if(it==RelaxationRounds){
              break()
            }
            
            # stop()
          }
          
        }else{
          trendPeak_average <- trendPeak_trendPeak_average_raw
        }
        # stop()
        
        # plot(RT.DetectPeak,trend)
        # abline(v=RT.DetectPeak[trendPeak_average[,2]])
        if(length(trendPeak_average)==0){
          if(length(rt_pep)>0){
            trendPeak_average <- data.frame(Height=NA,Peak=rt_pep,"Start"=rt_pep-estimatedPeakDistance/2,"End"=rt_pep+estimatedPeakDistance/2,Diff=0)
          }else{
            print("Warning:")
            print(paste("rt_pep is empty, returning NULL"))
            trendPeak_average <- NULL
          }
        }
        wise <- unlist(trendPeak_average[1,3:4])
      }else{
        
        wise <- which(RT.DetectPeak>=presetQuantiles[1]&RT.DetectPeak<=presetQuantiles[2])
        
      }
    }else{
      wise <- c(NA,NA)
      trendPeak_average <- NULL
      
    }
    
    
    
    
    traPeak <- data.frame(transitions[wise,])
    quaPeak <- RT.DetectPeak[wise]
    
    XIC <- sapply(1:dim(traPeak)[2],function(x){trapz(quaPeak,traPeak[,x])})
    names(XIC) <- colnames(traPeak)
    INT <- apply(traPeak,2,max,na.rm = T)
    
    Peaks <- list(trendPeak_average,PeaksSummary=NA,PeakDensity=NA)
    
    pl <- F
    if(pl){
      gs <- lapply(list(transitions),function(input){
        transitions.DetectPeak_temp <<- data.frame(input)
        transitions.DetectPeak_temp$RT <- RT.DetectPeak
        transitions.DetectPeak_temp_long <- melt(as.data.table(transitions.DetectPeak_temp),id.vars="RT")
        library(ggplot2)
        library(pracma)
        g <- ggplot(transitions.DetectPeak_temp_long,aes(x=RT,y=value,group=variable,color=variable))+geom_line()+theme(legend.position = "none")
        
        g
      })
      plot(gs[[1]]+geom_vline(xintercept = RT.DetectPeak[range(wise)]))
    }
  # 
  # 
  RTrange <- RT.DetectPeak[range(wise)]
  if(all(!is.na(RTrange))){
    if(diff(RTrange)<estimatedPeakDistance/3&ApplyMaximumWidth){
      try(showNotification("Adding default PeakDistance"),silent = T)
      addit <- estimatedPeakDistance/3-abs(diff(RTrange))/2
      RTrange <- RTrange+c(-addit,addit)
    }  
  }
  
  return(list(quantile= RTrange,all  = NA,XIC = XIC,intensity = INT,Peaks=Peaks, QualityMeasure=list(minpeakheight=minpeakheight,
                                                                                                     thresholdVal=thresholdVal,
                                                                                                     smoothingFactor=smoothingFactor
                                                                                                     # R_M=median(corres,na.rm = T),
                                                                                                     # R_sigma=sd(corres,na.rm = T),
                                                                                                     # R_n = length(corres[!is.na(corres)])
                                                                                                     )))
  
}
DetectPeak_v2 <- function(rt_pep,Peakwidth,transitions,RT,estimatedPeakDistance=1,
                          RelaxationRounds=20,BootstrapRounds=20,minpeakheight=0.1,
                          pl=F,smoothingFactor=51,presetQuantiles=NULL,ApplyMaximumWidth=T,Score = NULL,Identifier=NULL,...){
  RT.DetectPeak <- RT
  Identifier <- NULL
  save(rt_pep,Peakwidth,transitions,smoothingFactor,
       RT.DetectPeak,estimatedPeakDistance,BootstrapRounds,
       minpeakheight,presetQuantiles,ApplyMaximumWidth,Score,
       file=paste("DetectPeak_V2_TempNameSpace",Identifier,".rda",sep = ""))
  # scaling data
  corThresh <- 0.75 # initial required correlation Threshold
  smoothingFactor <- 11 # initial smoothin value
  minpeakheight <- 0.1 # minimal Peak Intensity after scaling
  thresholdVal <- 0.1 # minimal intensity of peaks (after scaling)
  
  # estimatedPeakDistance <- 1
  
  if(dim(transitions)[1]>=3){
    if(length(presetQuantiles)!=2){
      
      
      # transitions <- Neighbour_Mean_imputation(transitions)
      class(transitions) <- "data.frame"
      transitions.DetectPeak <- as.data.frame(transitions)
      ZeroRemove <- colSums(transitions,na.rm = T)
      
      transitions.DetectPeak <- apply(transitions.DetectPeak,2,function(x){scale(x)})
      # attempt to filter based on summed absolutec differences across each trace, idea is, that fluctating data is causing higher values
      diffs <- apply(transitions.DetectPeak,2,function(x){sum(abs(diff(x[!is.na(x)])))})
      diffs <- diffs/dim(transitions.DetectPeak)[1]
      CorInit <- F
      presetwarning <- options()$warn
      options(warn=-1)
      while(sum(CorInit)<3){
        try({
          CorInit1 <- sapply(1:dim(transitions)[2],function(x){cor(transitions[,1],transitions[,x])})>corThresh
          CorInit2 <- sapply(1:dim(transitions)[2],function(x){cor(transitions[,2],transitions[,x])})>corThresh
          CorInit3 <- sapply(1:dim(transitions)[2],function(x){cor(transitions[,3],transitions[,x])})>corThresh
          
          CorInit1[is.na(CorInit1)]  <- FALSE
          CorInit2[is.na(CorInit2)]  <- FALSE
          CorInit3[is.na(CorInit3)]  <- FALSE
          CorFinal <- CorInit1|CorInit2|CorInit3
        },silent = T)
        
        # CorInit <- apply(cor(transitions[,ZeroRemove>0]),2,median,na.rm = T)>corThresh
        corThresh <- corThresh-0.01
        if(corThresh<0.2){
          break()
        }
      }
      options(warn=presetwarning)
      
      # CorFinal <- rep(F,rep(length(ZeroRemove)))
      # CorFinal[ZeroRemove>0] <- CorInit
      # if(any(transitions.DetectPeak)<0){
      #   transitions.DetectPeak <- transitions.DetectPeak+abs(min(transitions.DetectPeak))
      # }
      # smoothing strong or weak
      
      repeatit <- T
      it <- 0
      thresholdVal <- 0.1
      sumvecCutoff <- 0.3
      diffThreshold <- 1# switched off, set to 0.1 to take action. It can cause problems
      
      trendraw <- apply(transitions.DetectPeak,1,function(x){M <- NA;try({M <- median(x[x!=0],na.rm = T)});M})
      trendPeak_trendPeak_average_raw <- findpeaks(trendraw,nups = 1,ndowns = 1,minpeakheight = 0,minpeakdistance = 1,threshold = 0)
      RTdiff <- abs(rt_pep-RT.DetectPeak[trendPeak_trendPeak_average_raw[,2]])
      trendPeak_trendPeak_average_raw <- as.data.frame(trendPeak_trendPeak_average_raw)
      trendPeak_trendPeak_average_raw <- trendPeak_trendPeak_average_raw[RTdiff<= estimatedPeakDistance,]
      # plot(RT.DetectPeak,transitions.DetectPeak[,1],type = "l")
      # points(RT.DetectPeak,Score,col = 2)
      # abline(v=rt_pep,col = 3)
      # points(RT.DetectPeak,transitions.DetectPeak[,1],col = 4,type = "l")
      # trendPeak_trendPeak_average_raw <- trendPeak_trendPeak_average_raw[trendPeak_trendPeak_average_raw[,1]>quantile(range(trendraw),0.1),]
      if(dim(trendPeak_trendPeak_average_raw)[1]!=1){
        while(repeatit){
          it <- it+1
          cat("\rRound",it)
          transitions.DetectPeak_smooth <- apply(transitions.DetectPeak,2,function(trend){
            trend <- trend
            trend2 <- c(rep(NA,10),trend,rep(NA,10))
            trend2raw <- trend2
            prwarnsettings <- options()$warn
            options(warn=-1)
            try({
              
              sdfun <- sd(trend[trend<=quantile(trend,0.3,na.rm = T)])
              if(is.na(sdfun)){
                sdfun <- 0.1
              }
              if(sdfun==0){
                sdfun <- 0.1
              }
              
              adv <-c( jitter(rep(min(trend,na.rm = T)-diff(range(trend,na.rm = T))*0.02,9),amount = ),min(trend,na.rm = T)-diff(range(trend,na.rm = T))*0.01)
              
              trend2 <- c(adv,trend,rev(adv))
              if(!all(is.na(trend2))){
                trend2[is.na(trend2)] <- min(trend2,na.rm = T)
                
              }
            },silent =T)
            
            try({
              
              trend2 <- savgol(trend2,fl=smoothingFactor,forder=4,dorder=0)
              if(any(is.na(trend2))){
                trend2[is.na(trend2)] <- jitter(min(trend))
              }
              trend2[trend2<min(trend)] <- jitter(min(trend))
              
            },silent = T)
            options(warn=prwarnsettings)
            # plot(trend2raw/trend2)
            # par(new=T)
            # plot(trend2raw,type = "l")
            # points(trend2,col = 2,type = "l")
            trend2
          })
          
          transitions.DetectPeak_smooth <- apply(transitions.DetectPeak_smooth,2,scale)
          transitions.DetectPeak_smooth[is.infinite(transitions.DetectPeak_smooth)] <- NA
          
          transitions.DetectPeak_smooth <- as.data.frame(transitions.DetectPeak_smooth)
          transitions.DetectPeak_smooth <- transitions.DetectPeak_smooth[,ZeroRemove>0&CorFinal]
          if(is.vector(transitions.DetectPeak_smooth)){
            transitions.DetectPeak_smooth <- (data.frame(nclue=transitions.DetectPeak_smooth))
          }
          
          iti<-0
          #Correlation as a potential quality metric
          # corres <- NA
          # try(    corres <- apply(cor(transitions.DetectPeak_smooth),2,median))
          
          RTdiff <- 0
          trendPeak_average <- NULL
          
          ({
            # Smoothing data
            
            trend <- apply(transitions.DetectPeak_smooth,1,function(x){M <- NA;try({M <- median(x[x!=0],na.rm = T)});M})
            trend[trend<min(trendraw,na.rm = T)] <- min(trendraw,na.rm = T)
            # plot(trendraw,type= "l")
            #add Trendbefore and After the data
            trend2 <- trend #c
            trendPeak_average <- findpeaks(trend2,nups = 1,ndowns = 1,minpeakheight = minpeakheight,minpeakdistance = 2,threshold = thresholdVal)
            if(length(trendPeak_average)>0){
              check2 <- dim(trendPeak_average)[1]>0 
              check1 <- T
            }else{
              check2 <- F
              check1 <- F
            }
            if(check1&check2){
              trendPeak_average[,2:4] <- trendPeak_average[,2:4]-10
              trendPeak_average[trendPeak_average<=0] <- 1
              trendPeak_average[trendPeak_average>length(RT.DetectPeak)] <- length(RT.DetectPeak)
              if(length(trendPeak_average)>0){
                correctedPosition <- trendPeak_average[,2]-10
                correctedPosition[correctedPosition<=0] <- 1
                correctedPosition[correctedPosition>length(RT.DetectPeak)] <- length(RT.DetectPeak)
                RTdiff <- abs(rt_pep-RT.DetectPeak[correctedPosition])
                if(length(RTdiff)==0){
                  RTdiff <- rep(0,dim(trendPeak_average)[1])
                }
              }
            }else{
              trendPeak_average <- NULL
            }
            
            # RTdiffCollect <- c(RTdiffCollect,min(RTdiff))
            # Determining distance from best scoring RT position
            
          })
          sumvec <- sum(abs(trendraw-trend[-c(c(1:10),(length(trend)-9):length(trend))]),na.rm = T)/length(trend)
          
          if((length(trendPeak_average)==0|all(RTdiff>estimatedPeakDistance))&sumvec<sumvecCutoff){
            # less restrictions if nothing passed the threshholds, simultanteous approach could become problematic, but seems to be fine at the moment
            minpeakheight <- minpeakheight-0.01
            if(minpeakheight<0){
              minpeakheight <- 0
            }
            thresholdVal <- thresholdVal*0.9
            smoothingFactor <- smoothingFactor+4
            diffThreshold <- diffThreshold+0.02
          }else{
            #stopping loop
            repeatit <- F
            # Preparing Data
            if(length(trendPeak_average)>0){
              trendPeak_average <- data.frame(trendPeak_average)
              colnames(trendPeak_average) <- c("Height","Peak","Start","End")
              trendPeak_average$Diff <- abs(RTdiff)
              # Filtering for minimal distance
              
              trendPeak_averagetemp <- trendPeak_average[trendPeak_average$Diff <=estimatedPeakDistance*2,]
              if(length(Score)>0){
                sc <-  Score[trendPeak_averagetemp$Peak]
                trendPeak_averagetemp <- trendPeak_averagetemp[sc==max(sc),]
              }
              if(dim(trendPeak_averagetemp)[1]==0){
                trendPeak_averagetemp <- trendPeak_average[trendPeak_average$Diff ==min(trendPeak_average$Diff),]
                # trendPeak_averagetemp <- data.frame(Height=NA,Peak=NA,"Start"=NA,"End"=NA,Diff=NA)
              }
              # Replacing artificialy zscore height with real summed intensities
              trendPeak_averagetemp$Height <- rowSums(transitions[trendPeak_averagetemp$Peak,])
              # Final Filtering on summed intensities
              trendPeak_average <- trendPeak_averagetemp[trendPeak_averagetemp$Height==max(trendPeak_averagetemp$Height),]#[1,]#apply(trendPeak_averagetemp,2,median)
              
            }
            
          }
          if(it==RelaxationRounds){
            break()
          }
          
          # stop()
        }
        
      }else{
        trendPeak_average <- trendPeak_trendPeak_average_raw
      }
      # stop()
      
      
      if(length(trendPeak_average)==0){
        if(length(rt_pep)>0){
          trendPeak_average <- data.frame(Height=NA,Peak=rt_pep,"Start"=rt_pep-estimatedPeakDistance/2,"End"=rt_pep+estimatedPeakDistance/2,Diff=0)
        }else{
          print("Warning:")
          print(paste("rt_pep is empty, returning NULL"))
          trendPeak_average <- NULL
        }
      }
      wise <- unlist(trendPeak_average[1,3:4])
    }else{
      
      wise <- which(RT.DetectPeak>=presetQuantiles[1]&RT.DetectPeak<=presetQuantiles[2])
      
    }
    trendPeak_average <- NULL
    
  }else{
    wise <- c(NA,NA)
    trendPeak_average <- NULL
    
  }
  
  
  
  
  traPeak <- data.frame(transitions[wise,])
  quaPeak <- RT.DetectPeak[wise]
  
  XIC <- sapply(1:dim(traPeak)[2],function(x){trapz(quaPeak,traPeak[,x])})
  names(XIC) <- colnames(traPeak)
  prewarnsettings <- options()$warn
  options(warn=-1)
  INT <- apply(traPeak,2,max,na.rm = T)
  options(warn=prewarnsettings)
  Peaks <- list(trendPeak_average,PeaksSummary=NA,PeakDensity=NA)
  
  pl <- F
  if(pl){
    gs <- lapply(list(transitions),function(input){
      transitions.DetectPeak_temp <<- data.frame(input)
      transitions.DetectPeak_temp$RT <- RT.DetectPeak
      transitions.DetectPeak_temp_long <- melt(as.data.table(transitions.DetectPeak_temp),id.vars="RT")
      library(ggplot2)
      library(pracma)
      g <- ggplot(transitions.DetectPeak_temp_long,aes(x=RT,y=value,group=variable,color=variable))+geom_line()+theme(legend.position = "none")
      
      g
    })
    plot(gs[[1]]+geom_vline(xintercept = RT.DetectPeak[range(wise)]))
  }
  # 
  # 
  RTrange <- RT.DetectPeak[range(wise)]
  if(all(!is.na(RTrange))){
    if(diff(RTrange)<estimatedPeakDistance/3&ApplyMaximumWidth){
      try(showNotification("Adding default PeakDistance"),silent = T)
      addit <- estimatedPeakDistance/3-abs(diff(RTrange))/2
      RTrange <- RTrange+c(-addit,addit)
    }  
  }
  
  return(list(quantile= RTrange,all  = NA,XIC = XIC,intensity = INT,Peaks=Peaks, QualityMeasure=list(minpeakheight=minpeakheight,
                                                                                                     thresholdVal=thresholdVal,
                                                                                                     smoothingFactor=smoothingFactor,
                                                                                                     thresholdVal=thresholdVal
                                                                                                     # R_M=median(corres,na.rm = T),
                                                                                                     # R_sigma=sd(corres,na.rm = T),
                                                                                                     # R_n = length(corres[!is.na(corres)])
  )))
}


############## Export ###################
dbtaNameExport <- function(ana,dbpath,pw){
  
  dbtaName <- paste(gsub("^mz","PEAKS",ana),sep = "_")
  return(dbtaName)
}
CondenseNames <- function(uniRaw){
  TabSplit <- strsplit(uniRaw,"")
  LE <- min(lengths(TabSplit))
  StartT <- sapply(TabSplit,function(x){x[1:LE]})
  EndT <- sapply(TabSplit,function(x){x[(length(x)-(LE-1)):length(x)]})
  StartT <- apply(StartT,1,function(x){length(unique(x)) == 1})
  FrontCut <- min(which(!StartT))
  EndT <- apply(EndT,1,function(x){length(unique(x)) == 1})
  EndCut <- max(which(!EndT))
  uniRawCut <- sapply(TabSplit,function(x){paste(x[FrontCut : EndCut],collapse = "")})
  return(uniRawCut)
}
# Wrapper function for automated Peak Detection During export
# Wrapper function for automated Peak Detection During export
# DetectPeakWrapper <- function(
#   ana,#TransitionLists
#   CANDIDATE_RT,# RT candidates (based on best FDR after CutOff)
#   dbp,#Database path
#   RetentionTimeWindow=5,
#   QType= "Intensities",
#   RT_BASED_onbestScore=T,
#   Reanalysis = F,
#   ApplyMaximumWidth=T,
#   MinPeakWidth=0.1,
#   MaxPeakWidth=0.5,
#   supersmooth_I_set=F,
#   supersmooth_bw_set=0.1,
#   Requantify_Priority="DL_Scores",
#   session = NULL
# ){
#   print("PeakWrapper")
#   save(dbp,ana,CANDIDATE_RT,RetentionTimeWindow,QType,RT_BASED_onbestScore,Reanalysis,ApplyMaximumWidth,MinPeakWidth,MaxPeakWidth,
#        supersmooth_I_set,supersmooth_bw_set,Requantify_Priority,file="temp.DetectPeakWrapper.rda")
#   # First DPLIST and Writing it to Database for Int and XIC
#   if(length(session)>0){
#     le <-length(ana)*length(unique(ana[[1]]$rawfile))
#     
#     progress <- Progress$new(session,min=0, max=le)
#     on.exit(progress$close())
#     progress$set(message = 'Compiling scans',
#                    detail = "almost done",value = 0)
#     progressit <<- 0
#     
#   }
#   DPlist <- lapply(1:length(ana),function(it){ it<-it
#   # it
#   x <- ana[[it]]
#   dbtaNameVec <- dbtaName(ana,dbp)[it]
#   # print(dbtaNameVec)
#   if(all(unlist(dblistT(db =dbp)) != dbtaNameVec)|Reanalysis){
#     x <- data.table(x)
#     temp <- x[rawfile==rawfile[1]]
#     # Going through rawfile by rawfile
#     DP <- x[,{
#       if(length(session)>0){
#         
#         progressit <<- progressit+1
#         progress$set(value=progressit)
#       }
#       # cat("\r",.GRP,.BY$rawfile)
#       temp <- .SD
#       grp <- .BY
#       temp2 <- temp
#       Candidate <- CANDIDATE_RT[CANDIDATE_RT$rawfile == grp$rawfile]
#       RTset <- median(Candidate$RT_Used,na.rm=T)
#       RTsetWindow <- RetentionTimeWindow/2
# 
#       temp <- temp[RT_Used>(RTset-RTsetWindow)&RT_Used<(RTset+RTsetWindow),.SD,]
#       if(dim(temp)[1]==0){
#         temp <- temp2
#         cleandata <- 1
#         #stop("PROBLEM WITH PeakFUN")
#       }else{
#         cleandata <- 0
#       }
#       
#       transidf <- temp[,.SD,.SDcols=1:(which(names(temp) =="charge")-1)]
#       transidf <- as.data.frame(transidf)
#       transidf[transidf==-1]<- 0
#       # print("DetectPeak Start")
#       rm("PeakDetected")
#       try({PeakDetected <- DetectPeak_v2(rt_pep = Candidate$RT_Used,
#                                       Peakwidth = RetentionTimeWindow,
#                                       transitions = transidf,
#                                       RT = temp$RT_Used,
#                                       scores=temp$DL_Scores,
#                                       FDR=temp$FDR,
#                                       RelaxationRounds=20,
#                                       MinPeakWidth = MinPeakWidth,
#                                       BootstrapRounds = 20,
#                                       MaxPeakWidth = MaxPeakWidth,supersmooth_I_set = supersmooth_I_set,supersmooth_bw_set = supersmooth_bw_set,ApplyMaximumWidth = T
#       )})
#       if(!exists("PeakDetected")){
#         return(NULL)
#       }
#       #$quantile
#       # print(PeakDetected)
#       # print("DetectPeak Finished")
#       
#       # Compiling Stuff:
#       infoquantile <- temp[temp$RT_min>= min(PeakDetected$quantile,na.rm = T)&temp$RT_min<=max(PeakDetected$quantile,na.rm= T),]
#       infoquantile <- infoquantileParser(infoquantile)
#       infoquantile <-rbind(infoquantile,infoquantile)
#       infoquantile <- apply(infoquantile,2,as.double)
#       
#       PeakDetected$all <- NULL     
#       QuantitationType <- c("XIC","Intensities")
#       Q <- PeakDetected$quantile
#       if(length(Q) == 0){
#         Q = as.double(c(NA,NA))
#       }
#       
#       IntensityStuff <- rbind(PeakDetected$XIC,PeakDetected$intensity)
#       IntensityStuff <- data.frame(IntensityStuff)
#       IntensityStuff$Q1 <- Q[1]
#       IntensityStuff$Q2 <- Q[2]
#       IntensityStuff$Q1align <- Q[1]
#       IntensityStuff$Q2align <- Q[2]
#       IntensityStuff$Peaks <- as.double(NA)
#       IntensityStuff$QuantitationType <- QuantitationType
#       OutTable <- cbind(IntensityStuff,infoquantile)
#       
#       if(dim(temp)[1]==0){
#         temp <- temp2
#         cleandata <- 1
#         #stop("PROBLEM WITH PeakFUN")
#       }else{
#         cleandata <- 0
#       }
#       
#       OutTable
#       
#     },rawfile]
#     DP[is.na(Q1)&Q2==1,c("Q1","Q2","Q1align","Q2align"):=NA]
#     setcolorder(DP,c(names(DP)[-1],"rawfile"))
#     
#     # DP_PeaksFun<- DP
#     dbsuccess <- dbwrite(x = DP,name = dbtaNameVec,db =dbp,overwrite = T)
#     # print("end Writing Peakfun Peaksfun")
#     
#   }else{
#     # print("Found Peak File")
#     DP <- dbread(x = dbtaNameVec,db =dbp)
#   }
#   # DP <- DP[DP$QuantitationType == QType,]# double check, if this is necessary
#   DP
#   })
#   names(DPlist) <- names(ana)
#   
#   # fu <- (sapply(DPlist,function(x){x <<- x;(cbind(x$Q1,x$Q2,x$rawfile))}))
#   DPlistBackup <- DPlist
#   
#   
#   
#   if(RT_BASED_onbestScore&length(DPlist)>1&Requantify_Priority!="none"){
#     print("Reassign")
#     # Extracting Windows
#     DPlistDT <<- lapply(1:length(DPlist),function(x){
#       tempu <- DPlist[[x]];
#       tempu$Precursor<-names(DPlist)[x];
#       tempu <- data.table(tempu,stringsAsFactors=F)
#       tempu[is.na(Q1)&Q2==1,c("Q1","Q2","Q1align","Q2align"):=NA]
#       tempu
#     })
#     DPDT <<- rbindlist(DPlistDT)
#     # Deciding on precursors per rawfile, which need to be requantified
#     Reassign <- DPDT[,{
#       temp <- .SD
#       temp$FDR[is.na(temp$FDR)] <- 10
#       
#       if(Requantify_Priority=="DL_Scores"){
#         sel <- temp$DL_Scores==max(temp$DL_Scores,na.rm = T)
#       }
#       if(Requantify_Priority=="FDR"){
#         sel <- FDR==max(FDR,na.rm = T)
#       }
#       if(Requantify_Priority=="Intensity"){
#         MaxFun <- as.numeric(apply(tempSelect <- temp[QuantitationType=="Intensities",],1,function(x){max(as.numeric(x)[1:(which(names(temp)=="Q1")-1)],na.rm = T)}))
#         sel <- Precursor==tempSelect$Precursor[MaxFun==max(MaxFun)]
#         
#       }
#       if(Requantify_Priority=="Light"){
#         
#         mz <- as.numeric(gsub("mz","",sapply(strsplit(temp$Precursor,"_"),"[[",1)))
#         sel <- mz==min(mz)
#       }
#       if(Requantify_Priority=="Heavy"){
#         mz <- as.numeric(gsub("mz","",sapply(strsplit(temp$Precursor,"_"),"[[",1)))
#         sel <- mz==max(mz)
#       }
#       Quantiles <- temp[sel,.(Q1,Q2)]
#       
#       Quantiles <- unique(Quantiles)
#       
#       if(all(is.na(Quantiles))){
#         # print("NULL")
#         print(paste("NULL for", .BY$rawfile))
#         Quantiles <- NULL
#       }else{
#         # print("FUN")
#         NeedsAnUpdate <- temp[sel,.(Precursor)]
#         Quantiles <- Quantiles
#         Precursor <- NeedsAnUpdate$Precursor
#         Quantiles <- cbind(Quantiles,Precursor)
#         Quantiles <- unique(Quantiles)
#       }
#       # print(length(Quantiles))
#       # print(Quantiles)
#       as.list(Quantiles)
#       
#     },.(rawfile)]
#     
#     Reassign <- unique(Reassign)
#     Reassign <- Reassign[!is.na(Reassign$Precursor),]
#     if(dim(Reassign)[1]>0){
#      
#       sapply(1:dim(Reassign)[1],function(its){
#         tempReassing <- Reassign[its,]
#         gr <- tempReassing$Precursor
#         
#         itlist <- which(names(ana)!=gr)
#         if(length(itlist)>0){
#           for(it in itlist){
#             tempana <- ana[[it]]
#             dbtaNameVec <- dbtaName(ana,dbp)[it]
#             DP <- dbread(dbtaNameVec,dbp)
#             DP <- data.table(DP,stringsAsFactors=F)
#             tempana <- data.table(tempana,stringsAsFactors=F)
#             Corrected <- tempana[!is.na(match(rawfile,tempReassing$rawfile)),{
#               # cat("\r",.GRP,.BY$rawfile)
#               tempana_rf <- .SD
#               grp <- .BY
#               if(dim(tempana_rf)[1]==0){
#                 OutTable <- NULL
#               }else{
#                 # stop()
#                 # Extract transitions
#                 transidf <- tempana_rf[,.SD,.SDcols=1:(which(names(tempana) =="charge")-1)]
#                 transidf <- as.data.frame(transidf)
#                 trans <<- SplitTransitionInfo(tempana_rf)
#                 trans$Transitions[trans$Transitions==-1]<-0
#                 # Specify RT range
#                 rangeinit <<- tempReassing[rawfile==grp$rawfile,.(Q1,Q2)]
#                 if(dim(rangeinit)[1]>0){
#                   
#                   RANGE <- range(rangeinit)
#                   # print("DetectPeak Start")
#                   PeakDetected <- list(quantile=c(NA,NA),XIC=NA,intensity=NA)
#                   try({PeakDetected <- DetectPeak_v2(rt_pep = mean(RANGE),
#                                                   Peakwidth = diff(RANGE),
#                                                   transitions = trans$Transitions,
#                                                   RT = trans$Info$RT_Used,
#                                                   presetQuantiles = RANGE,
#                                                   FDR=NULL,
#                                                   scores=trans$Info$DL_Scores,
#                                                   MinPeakWidth = MinPeakWidth,
#                                                   MaxPeakWidth = MaxPeakWidth,supersmooth_I_set  = supersmooth_I_set,supersmooth_bw_set = supersmooth_bw_set
#                   )})#$quantile
#                   temp <- trans$Info
#                   infoquantile <- temp[temp$RT_min>= min(PeakDetected$quantile,na.rm = T)&temp$RT_min<=max(PeakDetected$quantile,na.rm= T),]
#                   infoquantile <- infoquantileParser(infoquantile)
#                   infoquantile <-rbind(infoquantile,infoquantile)
#                   infoquantile <- apply(infoquantile,2,as.double)
#                   
#                   PeakDetected$all <- NULL     
#                   QuantitationType <- c("XIC","Intensities")
#                   Q <- PeakDetected$quantile
#                   if(length(Q) == 0){
#                     Q = as.double(c(NA,NA))
#                   }
#                   
#                   IntensityStuff <- rbind(PeakDetected$XIC,PeakDetected$intensity)
#                   IntensityStuff <- data.frame(IntensityStuff)
#                   IntensityStuff$Q1 <- Q[1]
#                   IntensityStuff$Q2 <- Q[2]
#                   IntensityStuff$Q1align <- Q[1]
#                   IntensityStuff$Q2align <- Q[2]
#                   IntensityStuff$Peaks <- as.double(NA)
#                   IntensityStuff$QuantitationType <- QuantitationType
#                   OutTable <- cbind(IntensityStuff,infoquantile)
#                 }else{
#                   OutTable <- NULL
#                 }
#                 
#                 # if(cleandata){
#                 #   OutTable$Q1 <- as.double(NA)
#                 #   OutTable$Q2 <- as.double(NA)
#                 #   OutTable$Q1align <- as.double(NA)
#                 #   OutTable$Q2align <- as.double(NA)
#                 #   
#                 #   OutTable[,1:which(colnames(OutTable)=="Q1")]<- as.double(NA)
#                 #   
#                 # }else{
#                 #   OutTable <-   as.list(OutTable)
#                 #   
#                 # }
#                 # OutTable <<- OutTable
#               }
#               
#               OutTable
#             },rawfile]
#             DPnochange <- DP[is.na(match(rawfile,Corrected$rawfile)),]
#             if(dim(Corrected)[1]>0){
#               DPnew <- rbind(DPnochange,Corrected)
#               
#               cat("\rrewriting",dbtaNameVec)
#               
#               if(dim(DPnew)[1]==dim(DP)[1]){
#                 dbwrite(DPnew,dbtaNameVec,dbp,overwrite = T)
#               }
#             }
#            
#           }
#           # DPnew
#         }else{
#           NULL
#         }
#         NULL
#       })
#       
#       dbtaNameVec <- dbtaName(ana,dbp)
#       DPlist <- lapply(dbtaNameVec,function(x){DP <- dbread(x,dbp);   #   DP <- DP[DP$QuantitationType == QType,];
#       DP
#       
#       
#       })
#     }
#   }
#   
#   DPlist
# }
DetectPeakWrapper <- function(
  ana,#TransitionLists
  CANDIDATE_RT,# RT candidates (based on best FDR after CutOff)
  dbp,#Database path
  RetentionTimeWindow=5,
  QType= "Intensities",
  RT_BASED_onbestScore=T,
  Reanalysis = F,
  ApplyMaximumWidth=T,
  MinPeakWidth=0.1,
  MaxPeakWidth=0.5,
  supersmooth_I_set=F,
  supersmooth_bw_set=0.1,
  Requantify_Priority="DL_Scores",
  session = NULL,
  parallel_db=NULL
){
  if(length(parallel_db)>0){
    dbp_write <- as.character(parallel_db)
    if(!file.exists(dbp_write)){
      library(RSQLite)
      db <- dbConnect(SQLite(),dbp_write)
      dbDisconnect(db)
    }
  }else{
    dbp_write <- dbp
  }
  print("PeakWrapper")
  try({save(dbp,ana,CANDIDATE_RT,dbp_write,RetentionTimeWindow,QType,RT_BASED_onbestScore,Reanalysis,ApplyMaximumWidth,MinPeakWidth,MaxPeakWidth,
       supersmooth_I_set,supersmooth_bw_set,Requantify_Priority,file="temp.DetectPeakWrapper.rda")})
  # First DPLIST and Writing it to Database for Int and XIC
  if(length(session)>0){
    le <-length(ana)*length(unique(ana[[1]]$rawfile))
    
    progress <- Progress$new(session,min=0, max=le)
    on.exit(progress$close())
    progress$set(message = 'Compiling scans',
                 detail = "almost done",value = 0)
    progressit <<- 0
    
  }
  DPlist <- lapply(1:length(ana),function(it){ it<-it
  # it
  it <<- it
  x <- ana[[it]]
  dbtaNameVec <- dbtaName(ana,dbp)[it]
  # print(dbtaNameVec)
  if(all(unlist(dblistT(db =dbp)) != dbtaNameVec)|Reanalysis){
    x <- data.table(x)
    temp <- x[rawfile==rawfile[1]]
    grp <- list(rawfile=x$rawfile[1])
    # Going through rawfile by rawfile
    DP <- x[,{
      if(length(session)>0){
        
        progressit <<- progressit+1
        progress$set(value=progressit)
      }
      # cat("\r",.GRP,.BY$rawfile)
      temp <- .SD
      grp <<- .BY
      temp2 <<- temp
      Candidate <- CANDIDATE_RT[CANDIDATE_RT$rawfile == grp$rawfile]
      RTset <- median(Candidate$RT_Used,na.rm=T)
      RTsetWindow <- RetentionTimeWindow/2
      
      temp <- temp[RT_Used>(RTset-RTsetWindow)&RT_Used<(RTset+RTsetWindow),.SD,]
      if(dim(temp)[1]==0){
        temp <- temp2
        cleandata <- 1
        #stop("PROBLEM WITH PeakFUN")
      }else{
        cleandata <- 0
      }
      
      transidf <- temp[,.SD,.SDcols=1:(which(names(temp) =="charge")-1)]
      transidf <- as.data.frame(transidf)
      transidf[transidf==-1]<- 0
      # print("DetectPeak Start")
      # try({rm(PeakDetected)},silent = T)
      PeakDetected <- list(quantile=c(NA,NA),XIC=NA,intensity=NA)
      
      ({PeakDetected <- DetectPeak_v2(rt_pep = Candidate$RT_Used,
                                         Peakwidth = RetentionTimeWindow,
                                         transitions = transidf,
                                         RT = temp$RT_Used,
                                         scores=temp$DL_Scores,
                                         FDR=temp$FDR,
                                         RelaxationRounds=20,
                                         MinPeakWidth = MinPeakWidth,
                                         BootstrapRounds = 20,
                                         MaxPeakWidth = MaxPeakWidth,
                                         supersmooth_I_set = supersmooth_I_set,
                                         supersmooth_bw_set = supersmooth_bw_set,
                                         ApplyMaximumWidth = T,
                                         Score=temp$DL_Scores,
                                         Identifier=grp$rawfile,
                                          info=temp
      )})
      if(!exists("PeakDetected")){
        OutTable <- NULL
      }else{
        infoquantile <- temp[temp$RT_min>= min(PeakDetected$quantile,na.rm = T)&temp$RT_min<=max(PeakDetected$quantile,na.rm= T),]
        infoquantile <- infoquantileParser(infoquantile)
        infoquantile <-rbind(infoquantile,infoquantile)
        infoquantile <- apply(infoquantile,2,as.double)
        
        PeakDetected$all <- NULL     
        QuantitationType <- c("XIC","Intensities")
        Q <- PeakDetected$quantile
        if(length(Q) == 0){
          Q = as.double(c(NA,NA))
        }
        
        IntensityStuff <- rbind(PeakDetected$XIC,PeakDetected$intensity)
        IntensityStuff <- data.frame(IntensityStuff)
        IntensityStuff$Q1 <- Q[1]
        IntensityStuff$Q2 <- Q[2]
        IntensityStuff$Q1align <- Q[1]
        IntensityStuff$Q2align <- Q[2]
        IntensityStuff$Peaks <- as.double(NA)
        IntensityStuff$QuantitationType <- QuantitationType
        OutTable <- cbind(IntensityStuff,infoquantile)
        
        if(dim(temp)[1]==0){
          temp <- temp2
          cleandata <- 1
          #stop("PROBLEM WITH PeakFUN")
        }else{
          cleandata <- 0
        }
        
        
      }
      #$quantile
      # print(PeakDetected)
      # print("DetectPeak Finished")
      
      # Compiling Stuff:
      OutTable
      
    },rawfile]
    DP[is.na(Q1)&Q2==1,c("Q1","Q2","Q1align","Q2align"):=NA]
    setcolorder(DP,c(names(DP)[-1],"rawfile"))
    
    # DP_PeaksFun<- DP
    dbsuccess <- dbwrite(x = DP,name = dbtaNameVec,db =dbp_write,overwrite = T)
    # print("end Writing Peakfun Peaksfun")
    
  }else{
    print("Found Peak File")
    DP <- dbread(x = dbtaNameVec,db =dbp)
  }
  # DP <- DP[DP$QuantitationType == QType,]# double check, if this is necessary
  DP
  })
names(DPlist) <- names(ana)
  
  # fu <- (sapply(DPlist,function(x){x <<- x;(cbind(x$Q1,x$Q2,x$rawfile))}))
  DPlistBackup <- DPlist
  
  
  
  if(RT_BASED_onbestScore&length(DPlist)>1&Requantify_Priority!="none"){
    print("Reassign")
    # Extracting Windows
    DPlistDT <<- lapply(1:length(DPlist),function(x){
      tempu <- DPlist[[x]];
      tempu$Precursor<-names(DPlist)[x];
      tempu <- data.table(tempu,stringsAsFactors=F)
      tempu[is.na(Q1)&Q2==1,c("Q1","Q2","Q1align","Q2align"):=NA]
      tempu
    })
    DPDT <<- rbindlist(DPlistDT)
    # Deciding on precursors per rawfile, which need to be requantified
    Reassign <- DPDT[,{
      temp <- .SD
      temp$FDR[is.na(temp$FDR)] <- 10
      
      if(Requantify_Priority=="DL_Scores"){
        sel <- temp$DL_Scores==max(temp$DL_Scores,na.rm = T)
      }
      if(Requantify_Priority=="FDR"){
        sel <- FDR==max(FDR,na.rm = T)
      }
      if(Requantify_Priority=="Intensity"){
        MaxFun <- as.numeric(apply(tempSelect <- temp[QuantitationType=="Intensities",],1,function(x){max(as.numeric(x)[1:(which(names(temp)=="Q1")-1)],na.rm = T)}))
        sel <- Precursor==tempSelect$Precursor[MaxFun==max(MaxFun)]
        
      }
      if(Requantify_Priority=="Light"){
        
        mz <- as.numeric(gsub("mz","",sapply(strsplit(temp$Precursor,"_"),"[[",1)))
        sel <- mz==min(mz)
      }
      if(Requantify_Priority=="Heavy"){
        mz <- as.numeric(gsub("mz","",sapply(strsplit(temp$Precursor,"_"),"[[",1)))
        sel <- mz==max(mz)
      }
      Quantiles <- temp[sel,.(Q1,Q2)]
      
      Quantiles <- unique(Quantiles)
      
      if(all(is.na(Quantiles))){
        # print("NULL")
        print(paste("NULL for", .BY$rawfile))
        Quantiles <- NULL
      }else{
        # print("FUN")
        NeedsAnUpdate <- temp[sel,.(Precursor)]
        Quantiles <- Quantiles
        Precursor <- NeedsAnUpdate$Precursor
        Quantiles <- cbind(Quantiles,Precursor)
        Quantiles <- unique(Quantiles)
      }
      # print(length(Quantiles))
      # print(Quantiles)
      as.list(Quantiles)
      
    },.(rawfile)]
    
    Reassign <- unique(Reassign)
    Reassign <- Reassign[!is.na(Reassign$Precursor),]
    if(dim(Reassign)[1]>0){
      
      sapply(1:dim(Reassign)[1],function(its){
        tempReassing <- Reassign[its,]
        gr <- tempReassing$Precursor
        
        itlist <- which(names(ana)!=gr)
        if(length(itlist)>0){
          for(it in itlist){
            tempana <- ana[[it]]
            dbtaNameVec <- dbtaName(ana,dbp)[it]
            DP <- dbread(dbtaNameVec,dbp)
            DP <- data.table(DP,stringsAsFactors=F)
            tempana <- data.table(tempana,stringsAsFactors=F)
            Corrected <- tempana[!is.na(match(rawfile,tempReassing$rawfile)),{
              # cat("\r",.GRP,.BY$rawfile)
              tempana_rf <- .SD
              grp <- .BY
              if(dim(tempana_rf)[1]==0){
                OutTable <- NULL
              }else{
                # stop()
                # Extract transitions
                transidf <- tempana_rf[,.SD,.SDcols=1:(which(names(tempana) =="charge")-1)]
                transidf <- as.data.frame(transidf)
                trans <<- SplitTransitionInfo(tempana_rf)
                trans$Transitions[trans$Transitions==-1]<-0
                # Specify RT range
                rangeinit <<- tempReassing[rawfile==grp$rawfile,.(Q1,Q2)]
                if(dim(rangeinit)[1]>0){
                  
                  RANGE <- range(rangeinit)
                  # print("DetectPeak Start")
                  PeakDetected <- list(quantile=c(NA,NA),XIC=NA,intensity=NA)
                  try({PeakDetected <- DetectPeak_v2(rt_pep = mean(RANGE),
                                                     Peakwidth = diff(RANGE),
                                                     transitions = trans$Transitions,
                                                     RT = trans$Info$RT_Used,
                                                     presetQuantiles = RANGE,
                                                     FDR=NULL,
                                                     scores=trans$Info$DL_Scores,
                                                     MinPeakWidth = MinPeakWidth,
                                                     MaxPeakWidth = MaxPeakWidth,supersmooth_I_set  = supersmooth_I_set,supersmooth_bw_set = supersmooth_bw_set,
                                                     Score=tempana_rf$DL_Scores,
                                                     Identifier = grp$rawfile,
                                                     info=trans$Info
                                                     
                  )})#$quantile
                  temp <- trans$Info
                  infoquantile <- temp[temp$RT_min>= min(PeakDetected$quantile,na.rm = T)&temp$RT_min<=max(PeakDetected$quantile,na.rm= T),]
                  infoquantile <- infoquantileParser(infoquantile)
                  infoquantile <-rbind(infoquantile,infoquantile)
                  infoquantile <- apply(infoquantile,2,as.double)
                  
                  PeakDetected$all <- NULL     
                  QuantitationType <- c("XIC","Intensities")
                  Q <- PeakDetected$quantile
                  if(length(Q) == 0){
                    Q = as.double(c(NA,NA))
                  }
                  
                  IntensityStuff <- rbind(PeakDetected$XIC,PeakDetected$intensity)
                  IntensityStuff <- data.frame(IntensityStuff)
                  IntensityStuff$Q1 <- Q[1]
                  IntensityStuff$Q2 <- Q[2]
                  IntensityStuff$Q1align <- Q[1]
                  IntensityStuff$Q2align <- Q[2]
                  IntensityStuff$Peaks <- as.double(NA)
                  IntensityStuff$QuantitationType <- QuantitationType
                  OutTable <- cbind(IntensityStuff,infoquantile)
                }else{
                  OutTable <- NULL
                }
                
                # if(cleandata){
                #   OutTable$Q1 <- as.double(NA)
                #   OutTable$Q2 <- as.double(NA)
                #   OutTable$Q1align <- as.double(NA)
                #   OutTable$Q2align <- as.double(NA)
                #   
                #   OutTable[,1:which(colnames(OutTable)=="Q1")]<- as.double(NA)
                #   
                # }else{
                #   OutTable <-   as.list(OutTable)
                #   
                # }
                # OutTable <<- OutTable
              }
              
              OutTable
            },rawfile]
            DPnochange <- DP[is.na(match(rawfile,Corrected$rawfile)),]
            if(dim(Corrected)[1]>0){
              DPnew <- rbind(DPnochange,Corrected)
              
              cat("\rrewriting",dbtaNameVec)
              
              if(dim(DPnew)[1]==dim(DP)[1]){
                dbwrite(DPnew,dbtaNameVec,dbp_write,overwrite = T)
              }
            }
            
          }
          # DPnew
        }else{
          NULL
        }
        NULL
      })
      
      dbtaNameVec <- dbtaName(ana,dbp)
      DPlist <- lapply(dbtaNameVec,function(x){DP <- dbread(x,dbp);   #   DP <- DP[DP$QuantitationType == QType,];
      DP
      
      
      })
    }
  }
  
  DPlist
}

## Export Ratio Calculate Function
erl.smooth <- function(x,y,smooth = F){
  ep.data <- data.frame(x=x,y=y)
  if(smooth){
    ep.data <- ep.data[order(ep.data$x),]
    ep.data <- unique(ep.data)
    ep.rle <- rle(ep.data$y)
    stair.midpoints <- cumsum(ep.rle$lengths) - floor(ep.rle$lengths/1)
    ep.data.sm <- ep.data[stair.midpoints,]
  }else{
    ep.data.sm <- ep.data
  }
  
  return(list(smooth=ep.data.sm,ori=ep.data))
}
FWHM.finder <- function(ep.data, mu.index){
  peak.height <- ep.data$y[mu.index]
  fxn.for.roots <- ep.data$y - peak.height/2
  indices <- 1:nrow(ep.data)
  root.indices <- which(diff(sign(fxn.for.roots))!=0)
  tmp <- c(root.indices,mu.index) %>% sort
  tmp2 <- which(tmp == mu.index)
  first.root <- root.indices[tmp2 -1]
  second.root <- root.indices[tmp2]
  HWHM1 <- ep.data$x[mu.index] - ep.data$x[first.root]
  HWHM2 <- ep.data$x[second.root] - ep.data$x[mu.index]
  FWHM <- HWHM2 + HWHM1
  FWHM2 = 2*min(c(HWHM1,HWHM2))
  return(list(HWHM1 = HWHM1,HWHM2 = HWHM2,FWHM = FWHM,FWHM2 = FWHM2))
}
# last change: 20200213
PeaksProcessingLM <- function(peaks,Remove_y1 =T,Plot=T,progressObj=NULL,outpath="./export"){
  CollectFun <<- c() # backupcollection of peaks
  if(all(is.na(peaks$mz))){
    peaks$mz <- as.numeric(gsub("mz","",peaks$Precursor))
  }
  
  graphics.off()
  RA <<- c()
  

  pdf(paste(outpath,"/",Sys.getpid(),"PeakView.pdf",sep = ""),width = 7)
  GRPmax  <-dim(peaks[,1,.(Sequence=Sequence,rawfile=rawfile,charge = charge,Gene=Gene)])[1]
  
  RatioCalc <-peaks[,{
    # Plot <- T
    if(length(progressObj)>0){
      progressObj$set(message = "Working on ratio estimation.",
                             detail = paste(length(ID),"Processes running."),value = .GRP/GRPmax)
    }else{
      li <- list.files(outpath,pattern=paste("RatioProcessing",Sys.getpid(),sep = "_"),full.names = T)
      if(length(li)>0){
        sapply(li,unlink)
      }
      write("",paste(outpath,paste("RatioProcessing",Sys.getpid(),.GRP,GRPmax,sep = "_"),sep = "/"))
    }
    
    cat("\r PeakView of",.GRP)
    
    
    hm <- .SD
    hmTemp <- hm
    # stop()
    if(Remove_y1){
      hm <- hm[variable!= "y1"]
      
    }
    gr <- .BY
    # Obtain Label combinations:
    hm$Label[is.na(hm$Label)] <- "light"
    
    
    if(length(unique(hm$Label))<2){
      
      UNI <- c(unique(hm$Label),"no_pair")
    }else{
      UNI <- unique(hm$Label)
      
    }
    UNI <- sort(unique(hm$mz),decreasing = T)

        if(length(UNI) >1){
      CombnTab <- combn(UNI,2)
      
      #align RT:
      # optional not yet done
      par(mfrow = c(4,4),mai = c(0.4,0.4,0.5,0.1),mgp = c(1.3,0.4,0),pointsize = 5)
      
      HUM <- hm[,{
        hmsel <- .SD
        cat("\r",.GRP)
        
        if(.GRP==1){
          started <<- 1
        }
        
        # hmsel <- .SD
        
        Ratios <- apply(CombnTab,2,function(sel){
          # sel <<- sel
          
          sel <- sort(sel,decreasing = T)
          s1 <- hmsel[mz==sel[2],] # label one # Light
          s2 <- hmsel[mz==sel[1],] # label two # Heavy
          RatioInfo <- paste(unique(s2$PrecursorInfo)," #/# ",unique(s1$PrecursorInfo),sep = "",collapse = "")
          if(dim(s1)[1] == 0){
            return(NULL)
          }
          if(dim(s2)[1] == 0){
            return(NULL)
          }
          
          AddZeros <- function(x,basevalue = 0,time = F,timeref=x){
            x <- x
            if(time){
              SF <- sd(diff(timeref[!is.na(timeref)]),na.rm = T)
              if(is.na(SF)){
                SF <- 0.1
              }
              c(min(x,na.rm = T)-SF,x,max(x,na.rm = T)+SF)
              
            }else{c(basevalue,x,basevalue)}
            
            
          }
          s2x <- AddZeros(s2$RT_min,time = T,timeref = c(s2$RT_min,s1$RT_min))
          s2y <- AddZeros(s2$value)
          s1x <- AddZeros(s1$RT_min,time = T)
          s1y <- AddZeros(s1$value)
          if(Plot){
            plot(s1x,s1y,xlim = range(c(s1x,s2x),na.rm = T),ylim = range(c(s1y,s2y)),na.rm = T,
                 xlab = "Time",ylab = "Fragment Intensity",main = "",frame = F,sub ="")
            if(started){
              mtext(hui<<- paste(paste(unlist(gr),collapse = " ")),3,line = 2,xpd = NA,adj = 0,cex = 0.7)
              # graphics.off()
              # system("open PeakView.pdf")
              started <<- 0
            }
            points(s2x,s2y,col = 2)
            legend("topright",legend = c(s1$Label[1],s2$Label[1]),col = c(1,2),pch = 20,lty = "solid",cex = 0.4,bty = "n")
            mtext(unlist(.BY)[1],3,line = 0.1,cex = 0.4)
            
          }
          # smoothing
          OutList <<- list(s1x,s1y,s2x,s2y)
          # commonrange <- range(c(OutList[[1]],OutList[[3]]),na.rm = T)
          
          commonrange <- range(c(s1x,s2x),na.rm = T)
          
          r1 <- range(s1x)
          r2 <- range(s2x)
          # intersect(seq(r1[1],r1[2],length.out = 1000),seq(r2[1],r2[2],length.out = 1000))
          
          # fit1 <- Fit.FunG(OutList[[1]],OutList[[2]],rangeval = commonrange)
          # fit2 <- Fit.FunG(OutList[[3]],OutList[[4]],rangeval = commonrange)
          # s1x <- OutList[[1]]
          # s1y <- OutList[[2]]
          # s2x <- OutList[[3]]
          # s2y <- OutList[[4]]
          # 
          
          
          SimpleFit <- Fun.Fit.lm(s1x,s1y,s2x,s2y,labelnames = c(sel[2],sel[1]))
          mtext(unlist(.BY)[1],3,line = 0,cex = 0.4)
          
          # fit1 <- Fit.FunG(s1x,s1y,rangeval = commonrange,iterations=100)
          # fit2 <- Fit.FunG(s2x,s2y,rangeval = commonrange,iterations=100)
          # 
          # if(Plot){
          #   lines(fit1$dffit)
          #   lines(fit2$dffit,col = 2)
          # }
          
          x <- SimpleFit$MappedIntensities$x
          y <- SimpleFit$MappedIntensities$y
          
          # x <- fit1$dffit$y
          # y <- fit2$dffit$y
          if(Plot&length(x)>0){
            # plot(c(0,x),c(0,y),xlab = paste("Intensity",s1$Label[1]),ylab = paste("Intensity",s2$Label[1]),frame = F)
            
          }
          # intercept <- 0
          # ra <- y/x
          
          LI <- list(x= as.double(x),
                     y= as.double(y),
                     Ratio= SimpleFit$RatioFit$coefficients,
                     IntensityRatio=max(s2y,na.rm = T)/max(s1y,na.rm = T),
                     "Relative Residuals"= sd(abs(SimpleFit$RatioFit$residuals)),
                     R.Pearson=as.double(0),
                     R.pvalue=as.double(1),
                     Label1=as.character(s1$Label[1]),
                     Label2=as.character(s2$Label[1]),
                     
                     Info1=as.character(s1$PrecursorInfo[1]),
                     Info2=as.character(s2$PrecursorInfo[1]),
                     Counts=as.double(1),
                     RatioType="NearestRatio",
                     RT1=as.double(SimpleFit$PeakRT),
                     RT2=as.double(SimpleFit$PeakRT),
                     RT1_Sigma = as.double(NA),
                     RT2_Sigma =as.double(NA)
          )
          
          if(0){
            if(abs(fit1$MaxDifferencePerc)>10|abs(fit2$MaxDifferencePerc)>10){
              RatioType <- "Intensity"
              
              LI <- list(x= as.double(x),
                         y= as.double(y),
                         Ratio= fit2$max/fit1$max,
                         "Relative Residuals"= as.double(-100),
                         R.Pearson=as.double(0),
                         R.pvalue=as.double(1),
                         Label1=as.character(s1$PrecursorInfo[1]),
                         Label2=as.character(s2$PrecursorInfo[1]),
                         Counts=as.double(1),
                         RatioType=RatioType,
                         RT1=as.double(NA),
                         RT2=as.double(NA),
                         RT1_Sigma = as.double(NA),
                         RT2_Sigma =as.double(NA)
              )
              
            }else{
              # y[y<0] <- 0
              # x[x<0] <- 0
              fit <- list(coefficients =  NA,residuals = NA)
              
              try(fit <- lm(I(y - intercept) ~ 0 + x),silent = F)
              try(sdresid <- as.double(sd(abs(rstandard(fit)),na.rm = T)))
              
              
              co <- list(estimate = NA,p.value = NA)
              try(co <- cor.test(x,y),silent = T)
              if(Plot){
                try(abline(0,fit$coefficients))
                
              }
              RatioType <- "Slope"
              LI<- list(x=as.double(x),
                        y=as.double(y),
                        Ratio=as.double(fit$coefficients),
                        "Relative Residuals"= as.double(sdresid),
                        R.Pearson=as.double(co$estimate),
                        R.pvalue=as.double(co$p.value),
                        Label1=as.character(s1$PrecursorInfo[1]),
                        Label2=as.character(s2$PrecursorInfo[1]),
                        Counts=as.double(length(ra[!is.na(ra)])),
                        RatioType = RatioType,
                        RT1=fit1$PeakRT,
                        RT2=fit2$PeakRT,
                        RT1_Sigma = fit1$PeakSigma,
                        RT2_Sigma =fit2$PeakSigma
                        
              )
              
              
              
            }
            
          }
          
          # LIFun <<- LI
          LI
          
        })
        as.list(rbindlist(Ratios))
      },variable]
      
      # HUM <<- HUM
      # HUMbackup <<- HUM
      
      # HUM <- HUMbackup
      HUM <- data.table(HUM)
      # stop()
      HUM$rf <- gr$rawfile
      
      cat("Removing Bad Fits")
      if(dim(HUM)[1]>0){
        # try(    HUM <- HUM[R.pvalue < 0.05&(R.Pearson)>0.5,,],silent = T)
        
        # outlier removal:
        cat("Removing Outliers")
        lr <- log2(HUM$Ratio)
        lr[is.infinite(lr)] <- NA
        
        slr <- sd(lr,na.rm = T)*2
        mlr <- mean(lr,na.rm = T)
        lrfun <- lr>=(mlr-slr)&lr<=(mlr+slr)
        table(lrfun)
        HUM <- HUM[lrfun,]
      }
      
      
      if(dim(HUM)[1]==0#|all(is.na(HUM$Ratio))
      ){
        RatioOut2 <- NULL
        # stop()
      }else{
        RatioOut2 <-HUM[,{
          # Outlier Selection within sd:
          Mean <- mean(Ratio,na.rm = T)
          CVs <- sd(Ratio,na.rm = T)
          SEL <- (Ratio>=(Mean-CVs*2))&(Ratio<=(Mean+CVs*2))
          SEL[is.na(SEL)] <- F
          # SEL2[is.na(SEL2)] <- F
          SEL2 <- SEL
          SEL2 <- y>0&x>0
          if(all(!SEL)){
            SEL <- rep(T,length(Ratio))
          }
          fit <- list(coefficients =  NA,residuals = NA)
          
          try(fit <- lm(I(y[SEL&SEL2] - 0) ~ 0 + x[SEL&SEL2]),silent = T)
          # try(sdresid <- as.double(sd(abs(rstandard(fit)),na.rm = T)))
          # Calculating Confidence Intervals of the slope:
          CONF <- as.double(rep(NA,2))
          try({
            CONF   <- confint(fit, level=0.9)
          })
          
          # CO <- fit$coefficients[1]
          # CO <-
          
          slopeCV <- NA
          try(slopeCV <- as.vector(vcov(fit))[1],silent = T)
          colvec <- as.factor(variable)
          colvec <- rainbow(max(as.numeric(colvec)))[as.numeric(colvec)]
          colvec[!SEL] <- "grey"
          par(mfrow = c(1,1),mai = c(0.8,0.8,0.5,0.1),mgp = c(1.5,0.5,0))
          plot(x,y,col = colvec,pch = 20,xlab = "Intensity Label 1",ylab = "Intensity Label 2",frame = F)
          grid()
          mtext(paste(unlist(.BY)[1],paste(unlist(gr),collapse = "\n")),3,line = 0.1,xpd = NA,adj = 0,cex = 0.5)
          try(abline(0,fit$coefficients),silent = T)
          
          cols <- unique(cbind(variable,colvec))
          legend("bottomright",cols[,1],col = cols[,2],cex = 0.3,pch = 20,bty = "n")
          co <- list(estimate = NA,p.value = NA)
          try(co <- cor.test(x,y),silent = T)
          list(Ratio = as.double(fit$coefficients),
               "Relative.Residuals" = as.double(sd(fit$residuals/median(y,na.rm = T),na.rm = T)),
               R.Pearson = as.double(co$estimate),
               R.pvalue = as.double(co$p.value),
               IntensityRatio=NA,
               Counts=length(x),
               Info = unlist(.BY)[1],
               RatioCI_lo=CONF[1],
               RatioCI_up=CONF[2],
               Log10_Fraction_CI=diff(log10(c(CONF[1],CONF[2])))/log10(as.double(fit$coefficients)),
               RatioCI_diff= diff(as.vector(CONF)),
               RatioVariance = as.double(slopeCV),
               IntensityRatioCounts = sum(RatioType == "Intensity"),SlopeRatioCounts = sum(RatioType == "Intensity"))
          
        },.(Label1,Label2)]
        
      }
      
      # TempCollect<<- rbind(TempCollect,RatioOut2)
      
      
      RatioOut3 <- as.list(RatioOut2)
      
    }else{
      RatioOut3 <- NULL
    }
    OrderFun <- hm[,unique(variable),variable]$variable
    hm$variable <- factor(hm$variable,OrderFun,OrderFun)
    if(length(RatioOut3)>0){
      ra <- log10(HUM$IntensityRatio)
      sel <- is.na(ra)|is.infinite(ra)
      # 
      # plot(HUM$`Relative Residuals`,log10(HUM$IntensityRatio))
      
      RatioOut3$IntensityRatio <- 10^mean(ra[!sel],na.rm = T)
      
      RatioOut3$Label1 <- as.character(RatioOut3$Label1)
      RatioOut3$Label2 <- as.character(RatioOut3$Label2)
      RatioOut3$Info <- as.character(RatioOut3$Info)
      
    }else{
      # g <- ggplot(hm,aes(RT_min,value,group=variable,color=variable))+geom_line()+ylab("Intensity")+xlab("Intensity [min]")
    }
    
    
    as.list(RatioOut3)
  },.(Sequence=Sequence,rawfile=rawfile,charge = charge,Gene=Gene)]
  
  graphics.off()
  # system("open PeakView.pdf")
  
  
  
  
  
  head(RatioCalc)
  RatioCalc$Info <- RA
  RatioCalc
}
# Ratios <- PeaksProcessingLM(peaks)

Fit.Fun <- function(x,y,rangeval = NULL,add.edges=T,iterations = 100){
  xtemp <- x
  ytemp <- y
  if(add.edges){
    if(length(rangeval)!=2){
      xr <- range(x,na.rm = T)
      
    }else{
      xr <- rangeval
    }
    fr <- seq(min(xr)-diff(xr)/2,min(xr),length.out = 5)
    la <- seq(max(xr),max(xr)+diff(xr)/2,length.out = 5)
    
    x <- c(fr,x,la)
    y <- c(rep(0,length(fr)),y,rep(0,length(la)))
  }
  # sm <- erl.smooth(x,y)
  ep.data <- data.frame(x=x,y=y)
  gamma.max <- which(max(ep.data$y)==ep.data$y)
  FWHM <- lapply(gamma.max, FWHM.finder, ep.data = ep.data)
  #sd of peaks
  gamma.sigma <- unlist(sapply(FWHM, '[', 'FWHM2'))/2.355
  if(is.na(gamma.sigma[1])){
    gamma.sigma <- sd(xr)
    warning("No Sigma, too less data points. Using sd of RT ranges.")
  }
  if(gamma.sigma[1]==0){
    gamma.sigma <- sd(xr)
    warning("No Sigma, too less data points. Using sd of RT ranges.")
  }
  peak.heights <- ep.data$y[gamma.max]
  duplis<- duplicated(x)&duplicated(y)
  # x <- x[!duplis]
  # y <- y[!duplis]
  continueLoop <- T
  gamma.sigmavec <- gamma.sigma[1]
  gamma.sigmavecFrac <- gamma.sigmavec*0.1
  itloop <- 0
  while(continueLoop){
    itloop<- itloop+1
    trcl <- class(try(fit <- nls(y ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2))),
                                 data = unique(ep.data),
                                 start = list(mean1 = ep.data$x[gamma.max][1],
                                              sigma1 = gamma.sigmavec,
                                              C1 = peak.heights[1]),
                                 algorithm = "port",
                                 control = list(maxiter = iterations, warnOnly=T))))
    if(trcl=="try-error"){
      gamma.sigmavec <- gamma.sigmavec+gamma.sigmavecFrac*itloop
      remove <- y!=max(y)|y!=0
      # x <- x[remove]
      # y <- y[remove]
      FitQuality <- "failed"
    }else{
      FitQuality <- "worked"
      continueLoop<- F
    }
  }
  
  predictRange <- seq(min(ep.data$x),max(ep.data$x),length.out = 50)
  
  dffit <- data.frame(x=seq(min(ep.data$x),max(ep.data$x),length.out = 50))
  dffit$y <- predict(fit, newdata=dffit)
  maxFit <- predict(fit,newdata = data.frame(x=ep.data$x[gamma.max][1]))
  maxDifference <- (maxFit-peak.heights[1])/peak.heights[1]*100
  # plot(x,y,ylim = c(0,8000))
  # lines(dffit)
  
  return(list(fit = fit,dffit = dffit,max = max(y,na.rm = T),Fit = FitQuality,MaxDifferencePerc = maxDifference,PeakSigma = gamma.sigmavec,PeakRT=ep.data$x[gamma.max][1]))
}
# Fit.FunG <- cmpfun(Fit.Fun)
# Creating Pairs
Fun.Fit.lm <- function(s1x,s1y,s2x,s2y,labelnames = NULL,...){
  # plot(s2x,s2y,type="b",ylim = c(0,60000))
  # points(s1x,s1y,type="b",col= 2)
  
  if(length(s1x)>length(s2x)){
    Refx <- s2x
    Refy <- s2y
    Checkx <- s1x
    Checky <- s1y
    format = 1
  }else{
    Refx <- s1x
    Refy <- s1y
    Checkx <- s2x
    Checky <- s2y
    format = 0
  }
  
  
  VAL <- lapply(1:length(Refx),function(x){
    cat("\r",x)
    xsel <- Refx[x]
    xdif <- Checkx-xsel
    xdifp <- xdif>=0
    xdifm <- xdif<=0
    if(any(xdifp)&any(xdifm)){
      sel <- c(max(which(xdifm)),min(which(xdifp)))
      Y_M <-Checky[sel]
      X_M <- Checkx[sel]
      fit <- lm(y~x,data.frame(x=X_M,y=Y_M))
      FitVal <- predict(fit,newdata = data.frame(x=xsel))
      if(format){
        LI <- list(x=FitVal,y=Refy[x])
        
      }else{    
        LI <- list(x=Refy[x],y=FitVal)
      }
      
      return(LI)
      
    }else{return(NULL)}
    
    
  })
  VAL <- VAL[lengths(VAL)>0]
  
  VAL <- rbindlist(VAL)
  VAL <- unique(VAL)
  if(length(VAL$x)>0){
    
    plot(VAL$x,VAL$y,xlim = c(0,max(c(VAL$x,VAL$y))),ylim =c(0,max(c(VAL$x,VAL$y))),xlab = paste(labelnames[1],"Intensity"),ylab = paste(labelnames[2],"Intensity"),...)
    fit <- lm(y~x+0,VAL)
    abline(0,1,lty = "dotted")
    points(seq(0,8000,length.out = 10),seq(0,8000,length.out = 10)*fit$coefficients,type="l")
    
  }else{fit <- list(x=as.double(NA),y=as.double(NA),coefficients=as.double(NA),residuals=as.double(NA))}
  
  
  list(RatioFit=fit,MappedIntensities =VAL,max = NA,Fit = "worked",MaxDifferencePerc=0,PeakSigma=0,PeakRT=median(Refx))
}
# lapply(DPlist,function(x){x$Q1})
# one <- DPlist[[1]]
# one <- one[order(one$rawfile),]
# two <- DPlist[[2]]
# two <- two[order(one$rawfile),]
# plot(one$Q1,two$Q1)
# any(one$Q1==two$Q1,na.rm = T)
# unique(one$Q1,two$Q1)
############## Not implmented here #############
ConvertMSMStoPicky <- function(mspath,IL = NULL,sep = ";",ppm1 = 100,  UseBestScoringPeptide = T){
  # IL <<- IL
  # mspath <<- mspath
  if(length(IL) == 0){
    path <- dirname(mspath)
  }else{
    path <- dirname(IL)
  }
  
  ms <- read.csv(mspath,sep = "\t",stringsAsFactors = F,check.names = F)
  ms <- ms[ms$Intensities!="",]
  # ms[ms$Sequence == "MTKQPAKK",]$`m/z`
  
  if(length(IL)> 0){
    il <- read.csv(IL,sep = sep)
    SelMZ <- cbind(il$Mass..m.z.,il$Species)
    
    
  }else{
    SelMZ <- cbind(ms$`m/z`,ms$Charge)
  }
  
  
  if(UseBestScoringPeptide){
    Sel <- apply(unique(SelMZ),1,function(x){
      # x <<-x
      ppmwindow <- x[1]/(1/ppm1*1000000)
      ppmwindow <- c(x[1]-ppmwindow/2,x[1]+ppmwindow/2)
      if(!is.na(x[2])){
        which(ms$`m/z` >=min(ppmwindow) & ms$`m/z` <= max(ppmwindow)& ms$Charge == x[2])
        
      }else{
        which(ms$`m/z` >=min(ppmwindow) & ms$`m/z` <= max(ppmwindow))
        
      }
    })
    
  }else{
    Sel <- 1:dim(ms)[1]
    
    
  }
  
  names(Sel) <- unique(apply(unique(SelMZ),1,paste,collapse = "#"))
  lapply(Sel,function(x){unique(ms$`Modified sequence`[x])})
  
  # il[lengths(Sel) ==0,]
  # 466.2733
  
  Sel <- Sel[lengths(Sel) != 0]
  id <<- 0
  if(length(Sel) > 0){
    
    colnames <- c("Sequence",	"Modified.sequence",	"Charge",	"Fragmentation",	"Mass.analyzer",	"mz",	"Retention.time",	"Score",	"ETD.identification.type",	"CollisionEnergy",	"SpecID",	"Proteins",	"GeneSymbol",	"Hydrophobicity",	"PredictionRange",	"Start",	"End",	"ProteinColor",	"ID",	"ProteinsTemp","Matches","Intensities","Masses")
    colnamesCheck <- gsub("."," ",colnames,fixed = T) 
    colnamesCheck <- gsub("mz","m/z",colnamesCheck,fixed = T)
    if(all(colnamesCheck != "Modified sequence")){
      ms$`Modified sequence` <- ms$Sequence
    }
    for(i in colnamesCheck){
      if(all(i!= colnames(ms))){
        ms <- cbind(ms,NA)
        colnames(ms)[dim(ms)[2]] <- i
      }
    }
    RET <- sapply(Sel,function(x){
      id <<- id+1
      # x <<- x
      temp <- ms[x,]
      if(dim(temp)[1] > 1){
        xout <-c() 
        if(all(is.na(temp$`Modified sequence`))){
          
          temp$`Modified sequence` <- temp$Sequence  
        }
        for(i in unique(temp$`Modified sequence`)){
          tempx <- temp[temp$`Modified sequence` == i,]
          tempx <- tempx[tempx$PEP == min(tempx$PEP,na.rm = T),]
          LE <- lengths(strsplit(tempx$Matches,";"))
          tempx <- tempx[LE == max(LE,na.rm = T),]
          tempx <- tempx[tempx$Score == max(tempx$Score,na.rm = T),]
          if(length(tempx$`Intensity coverage`) > 0){
            tempx <- tempx[tempx$`Intensity coverage` == max(tempx$`Intensity coverage`,na.rm = T),]
            
          }
          xout <- rbind(xout,tempx)
        }
        temp <- xout[1,]
      }
      checknames <- c("Sequence","Modified sequence","Charge","Type","m/z","Retention time","Score","ETD identification type","Proteins","Gene Names")
      for(c in checknames){
        if(all(colnames(temp) != c)){
          temp <- cbind(temp,NA)
          colnames(temp)[dim(temp)[2]] <- c
        }
      }
      xout <- cbind(temp$Sequence,temp$`Modified sequence`,temp$Charge,temp$Type,"",temp$`m/z`,temp$`Retention time`,temp$Score,temp$`ETD identification type`,"",id,temp$Proteins,temp$`Gene Names`,"","",temp$`Retention time`,temp$`Retention time`,"","","",temp$Matches,temp$Intensities,temp$Masses)
      
      
      return(t(xout))
      
    })
    RET <- t(RET)
    colnames(RET)<- colnames
    
    RET <- data.frame(RET,stringsAsFactors = F)
    RET$Modified.sequence[is.na(RET$Modified.sequence)] <- as.character(RET$Sequence[is.na(RET$Modified.sequence)])
    
    Matches <- unlist(MatchesList <- strsplit(as.character(RET$Matches),";"  ))
    Intensities <- unlist(strsplit(as.character(RET$Intensities),";"  ))
    Masses  <- unlist(strsplit(as.character(RET$Masses),";"  ))
    SpecID <- unlist(sapply(1:length(MatchesList),function(x){rep(RET$SpecID[x],lengths(MatchesList)[x])}))
    SpecType <-2
    Intensities_RAW <- Intensities
    tt <- cbind(SpecID,Masses,Intensities,Matches,Intensities_RAW,SpecType)
    wd <- getwd()
    setwd(path)
    
    if(length(IL)== 0){
      Inclusionlist <- cbind(RET$mz,NA,NA,RET$Charge,"Positive",NA,NA,26,NA,paste(RET$Modified.sequence,RET$Proteins,RET$GeneSymbol,sep = "#"))
      colnames(Inclusionlist) <- c("Mass [m/z]","Formula [M]","Species","CS [z]","Polarity","Start [min]","End [min]","(N)CE","MSX ID","Comment")
      write.table(Inclusionlist,"InclusionList.txt",sep = "\t",quote = F,row.names =F)
    }
    write.table(tt,"TransitionTable.txt",sep = "\t",quote = F,row.names = F)
    write.table(RET,"SpectraTable.txt",sep = "\t",quote = F,row.names = F)
    setwd(wd)
    return(path)
  }else{print("Empty Selection, skipping conversion")}
}
SilacConverter <-function(ilpath,ShiftTableIni){
  setwd(ilpath)
  if(!file.exists("original/SpectraTable.txt")){
    dir.create("original")
    file.rename("SpectraTable.txt","original/SpectraTable.txt")
    file.rename("TransitionTable.txt","original/TransitionTable.txt")
    
  }
  st <- read.csv("original/SpectraTable.txt",sep = "\t",stringsAsFactors = F)
  tt <- read.csv("original/TransitionTable.txt",sep = "\t",stringsAsFactors = F)
  st$Matches <- tolower(st$Matches)
  tt$Matches <- tolower(tt$Matches)
  
  stshifted <- grep("#",st$Sequence,fixed = T)
  if(length(stshifted)>0){
    print("Warning, found shifts")
    stshiftedSPECID <- st$SpecID[stshifted]
    
    st <- st[-stshifted,]
    ttshifted <- unique(unlist(sapply(stshiftedSPECID,function(x){
      which(x == tt$SpecID)
    })))
    tt <- tt[-ttshifted,]
  }
  stOri <- st
  
  
  # 1. Convert Spectra TABLE
  #ShiftTable <- data.table(mass = c(8.01419),AA = c("K"))
  #ShiftTable <- data.table(mass = c(4.0251069836),AA = c("K"))
  
  ShiftGrepIni <- lapply(ShiftTableIni$AA,gregexpr,st$Sequence)
  names(ShiftGrepIni) <- ShiftTableIni$TYPE
  
  tt_ShiftOut <- c()
  st_ShiftOut <- c()
  for(shifttype in unique(ShiftTableIni$TYPE)){
    ShiftTable <- ShiftTableIni[ShiftTableIni$TYPE == shifttype,]
    ShiftGrep <- ShiftGrepIni[ShiftTableIni$TYPE == shifttype]
    
    ttshift <- tt
    st <- stOri
    for(x in 1:dim(st)[1]){
      tempx <- st[x,]
      cat("\r",x)
      grepResults   <- lapply(ShiftGrep,function(y){y[[x]][y[[x]] != -1]})
      ShiftMass     <- sum(lengths(grepResults)*ShiftTable$mass,na.rm = T)
      ShiftMZ       <- ShiftMass/tempx$Charge
      st$mz[x]      <- tempx$mz+ShiftMZ
      SHIFTS <- names(grepResults)[lengths(grepResults)>0]
      if(length(SHIFTS) > 0){
        st$Sequence[x] <- paste(st$Sequence[x],paste(SHIFTS,collapse = "",sep = ""),sep = "#")
        st$Modified_sequence[x] <- gsub("_$",paste("#",SHIFTS,"_",collapse = "",sep = ""),st$Modified_sequence[x])
        
      }
      # Correct TransitionTABLE
      se <- tempx$Sequence
      forw <- seq(1:nchar(se))
      revw <- rev(forw)
      Kpos <- ShiftGrep$K[[x]]
      Kpos <- Kpos[Kpos != -1]
      Rpos <- ShiftGrep$R[[x]]
      Rpos <- Rpos[Rpos != -1]
      
      Pos <- lapply(ShiftGrep,function(y){y[[x]]})
      posvec <- revw
      FragmentShift <- function(posvec,ShiftTable,reverse = F){
        sapply(1:length(posvec),function(i){
          
          
          LE <- posvec[1:i]
          LEmapped <- lapply(Pos,function(y){
            # y <<- y
            # if(reverse){
            # y[y!=-1] <- nchar(se)-y[y!=-1]+1
            # }
            m <- match(LE,y)
            length(m[!is.na(m)])
          })
          SU <- lapply(1:dim(ShiftTable)[1],function(x){
            sum(LEmapped[[x]]*ShiftTable$mass[x],na.rm = T)
          })
          return(sum(unlist(SU),na.rm = T))
          
        })
        
        
      }
      
      abcshift <- FragmentShift(forw,ShiftTable)
      
      xyzshift <- FragmentShift(revw,ShiftTable,reverse = T)
      
      abc <- lapply(c("a","b","c"),function(x){
        c(paste(x,1:nchar(se),sep = ""))
      })
      xyz <- lapply(c("x","y","z"),function(x){
        c(paste(x,1:nchar(se),sep = ""))
      })
      
      #xyz ions
      ttshiftx <- ttshift[ttshift$SpecID == tempx$SpecID,]
      FragmentShift <- sapply(ttshiftx$Matches,function(y){
        # y <<- y
        y <- strsplitslot(y,1,"[-(]")
        typeselect <- sapply(c("^a|^b|^c","^x|^y|^z"),function(x){
          grep(x,y)
        })
        shift1 = 0
        if(lengths(typeselect)[1]==1){
          yl <- lapply(abc,function(z){which(z==y)})
          if(all(lengths(yl) == 0)){
            shift1 = 0
          }else{
            yl <- unlist(yl[lengths(yl) != 0])
            shift1 <- abcshift[yl]
          }
        }
        if(lengths(typeselect)[2]==1){
          yl <- lapply(xyz,function(z){which(z==y)})
          if(all(lengths(yl) == 0)){
            shift1 = 0
          }else{
            yl <- unlist(yl[lengths(yl) != 0])
            shift1 <- xyzshift[yl]
          }
        }
        return(shift1)
      })
      
      gr <- gregexpr(".[+]",ttshiftx$Matches)
      charge <- sapply(1:length(gr),function(g){
        y <- gr[[g]]
        # y <<- y
        y[is.na(y)] <- -1
        if(y!= "-1"){
          return(as.numeric(substr(ttshiftx$Matches[g],(y),(y))))
        }else{
          return(1)
        }
      })
      charge[is.na(charge)] <- 1
      
      ttshift$Masses[ttshift$SpecID == tempx$SpecID] <- ttshiftx$Masses+FragmentShift/charge
    }
    barplot(table(round(ttshift$Masses-c(tt$Masses),2)),las = 2)
    
    st$SpecID <- st$SpecID +max(st$SpecID)
    ttshift$SpecID <- ttshift$SpecID+max(ttshift$SpecID)
    
    tt_ShiftOut <- rbind(tt_ShiftOut,ttshift)
    st_ShiftOut <- rbind(st_ShiftOut,st) 
    # stop()
    
  }
  
  # stop()
  # OutputTable --------
  tt$Shift <- 0
  tt_ShiftOut$Shift <- tt_ShiftOut$Masses-tt$Masses
  ttall <- rbind(tt,tt_ShiftOut)
  stall <- rbind(stOri,st_ShiftOut)
  
  # stop()F
  # dir.create("original")
  # setwd("Label_Added")
  write.table(stall,"SpectraTable.txt",sep = "\t",quote = F,row.names = F)
  write.table(ttall,"TransitionTable.txt",sep = "\t",quote = F,row.names = F)
  setwd("../")
  
}
#ShiftTableIni <- data.table(mass = c(8.01419,10.0082,-8.01419,-10.0082),AA = c("K","R"),TYPE = c("heavy","heavy","light","light"))

# Transfers Identification between PickyAnalyzer Databases
transferPeaks_db <- function(db1,db2){
  # db1 <<- db1
  # db2 <<- db2
  db1 <- dbConnect(SQLite(),db1)
  db2 <- dbConnect(SQLite(),db2)
  
  PEAKS <- grep("^PEAKS",dbListTables(db1),value = T)
  if(length(PEAKS)>0){
    sapply(PEAKS,function(x){
      temp <- dbReadTable(db1,x)
      dbWriteTable(db2,x,temp,overwrite = T)
    })
  }
  dbDisconnect(db1)
  dbDisconnect(db2)
  
}

# plots transtable with GGplot (is significantly slower)
TransitionGGplot <- function(subDat,
                             ILtemp,#RetentionTeim
                             Ppos,#PEAK
                             secPlotType,xl,yl,PEAKS=NULL){
  
  transitionNames <- 1:(which(names(subDat) =="charge")-1)
  transitionNames <- names(subDat)[transitionNames]
  transitionNames <- transitionNames[transitionNames!="rawfile"]
  
  # RetTime = c(ILtemp$Start,ILtemp$End)
  # subDat <- subDat[,,]
  
  p1 <- ggplot(subDat,aes_string(x="RT_Used",y=transitionNames[1]))+xlab("RT")+ylab("Intensity")+xlim(xl)+ylim(yl)
  
  p1 <- p1+geom_rect(aes(xmin = Ppos[1], xmax = Ppos[2], ymin = -Inf, ymax = Inf),
                     fill = "lightgrey", alpha = 0.03)
  #### adding Lines of different Transitions
  for(i in transitionNames){
    col <- TransCole[TransCole[,1]==i,2]
    if(length(col)==0){
      col <- "grey"
    }
    p1 <- p1+geom_line(aes_string(y=i),color  = col)
  }
  # overlay otherplot
  Selected <- subDat[secPlotType][,1]
  minSelected <- min(Selected,na.rm =T)
  yl2 <- range(Selected)
  # normalizing to axis
  Selected <- Selected-minSelected
  maxSelected <- max(Selected)
  Selected <- Selected/maxSelected*max(yl)
  # cex scale
  cexOfPoints <-  Selected/max(Selected)
  cexOfPoints <- cexOfPoints
  cexOfPoints[cexOfPoints<0.01]<- 0.001
  # adding to plot
  # colorrange<- colorRampPalette(c("lightgrey","black"))(100)
  col<- rep("black",length(cexOfPoints))
  col[subDat$FDR< FDRcutoff]<-"red"
  
  p1 <- p1+geom_point(aes(y=Selected,size = cexOfPoints,alpha = cexOfPoints*0.2),shape = 2,color = col)+scale_size(range = c(0, 1))+theme() 
  
  p1<- p1+scale_y_continuous(sec.axis=sec_axis(~.*maxSelected/max(yl2)/maxSelected-minSelected))
  
  # p1 <- p1+geom_vline(xintercept = as.numeric(RetTime))
  p1 <- p1+theme_classic()+theme(legend.position = "none")
  if(length(PEAKS)!=0){
    PEAKS <- PEAKS[!is.na(PEAKS)]
    for(i in PEAKS){
      p1 <- p1+geom_point(aes(x=i,y=-diff(yl)*0.02),shape=2)
    }
    # Peaks <<- Peaks
    # points(Peaks,rep(yd,length(Peaks)),pch = 2,xpd = NA,col ="grey",cex = 1)
  }
  p1+ annotate("text",  x=Inf, y = Inf, label = paste("mz:",as.character(round(unique(subDat$mz),1))), vjust=1, hjust=1)
}

MakeFeatureTable <- function(mztable,defaultPeptideLength=NULL,mzname="unknown",mdpar=c("UNKNOWN")){
  library(h2o)
  h2o.init()
  mztable[,id:=1:dim(mztable)[1]]
  mztable <- mztable[,.SD[order(RT_min)],rawfile]
  
  mztable <- mztable[,ValiScoringFunction(.SD,.BY),rawfile]
  # mztable <- mztable[,ValiScoringFunction(.SD,.BY),]
  
  
  mztable <- data.table(mztable)
  
  
  #### h2o prediction
  TrainSet <- mztable
  # TrainSet$ScaAfter <- NULL
  # TrainSet$ScaAfter2 <- NULL
  # TrainSet$ScaBefore <- NULL
  # TrainSet$ScaBefore2 <- NULL
  # 
  
  TrainSet[,ScaBefore:= {c(SCAall[-1],0)},.(rawfile)]
  TrainSet[,ScaBefore2:= {
    if(length(SCAall)>3){
      ret <-c(SCAall[-c(1:2)],0,0)
    }else{
      ret <- 0
    }
    ret
    
  },.(rawfile)] 
  
  TrainSet[,ScaAfter:= c(0,SCAall[-length(SCAall)]),.(rawfile)]
  TrainSet[,ScaAfter2:= 
             {
               if(length(SCAall)>3){
                 ret <- c(0,0,SCAall[-c(length(SCAall):(length(SCAall)-1))])
               }else{
                 ret <- 0
               }
               ret
               
             },.(rawfile)]
  TrainSet$MLOGP <- -log10(TrainSet$R_p)
  if(length(defaultPeptideLength)>0){
    print("Adding PeptideLength")
    TrainSet$Peptide_Length <- defaultPeptideLength
  }
  # evaluating:
  
  
  mdpar_match <- match(mdpar,colnames(TrainSet))
  missing <- mdpar[is.na(mdpar_match)]
  print(missing)
  
  if(length(missing)>0){
    print(paste("ReturningNUll, missing",missing))
    return(NULL)
  }
  mztableSelect <- TrainSet[,.SD,.SDcols = mdpar_match]
  mztableSelect$scan <- TrainSet$scan
  mztableSelect$rawfile <- TrainSet$rawfile
  
  mztableSelect$mzname <- mzname
  mztableSelect <- mztableSelect[order(TrainSet$id)]
  mztableSelect
}
# Expanding Spectrum Specificity: #########
ExpandSpectrumSpecificity <- function(InclusionlistPath){
  print("Finding specific fragments")
  if(file.exists(paste(InclusionlistPath,"SpectraTableOriginal.txt",sep = "/"))){
    ST <- fread(paste(InclusionlistPath,"SpectraTableOriginal.txt",sep = "/"),sep = "\t")
    
  }else{
    ST <- fread(paste(InclusionlistPath,"SpectraTable.txt",sep = "/"),sep = "\t")
    
  }
  if(file.exists(paste(InclusionlistPath,"TransitionTableOriginal.txt",sep = "/"))){
    TT <- fread(paste(InclusionlistPath,"TransitionTableOriginal.txt",sep = "/"),sep ="\t")
    
  }else{
    TT <- fread(paste(InclusionlistPath,"TransitionTable.txt",sep = "/"),sep = "\t")
    
  }
  tempfun <- environment()
  tempfun$TToptimizedTemp  <- TT
  tempfun$TToptimizedTemp$Specific <- F
  
  
  # Retentionszeit # optional# potential parameter. Only Sequence implemented at the moment
  PairComparisson <- c("Mass","ExactMass","Sequence")[3]
  if(PairComparisson=="Mass"){
    CombiCheck <- combn(ST$SpecID,2)
    Pair<- apply(CombiCheck,2,function(x,Type="ppm",ppmset=10,WindowTolerance=1.3){
      x <<- x
      
      if(Type=="Da"){
        WithinTolerance <- abs(diff(ST$mz[x]))<(WindowTolerance/2)
      }
      if(Type=="ppm"){
        xm <- sort(ST$mz[x])
        wi <- ppmWindow(xm[2],ppm = ppmset)
        WithinTolerance <- xm[1]>=wi[1]&xm[1]<=wi[2]
      }
      WithinTolerance
    })
    Isobaric <- CombiCheck[,Pair]
    
  }
  
  if(PairComparisson=="Sequence"){
    STFUN <- ST[,.(modSeq=unique(`Modified sequence`),SpecID=unique(SpecID)),Sequence]
    STFUN$StrippedSequence <- sapply(strsplit(STFUN$Sequence,"#"),"[[",1)
    Pairs <- STFUN[,{
      if(length(unique(SpecID))>1){
        re <- combn(unique(SpecID),2)
        re2 <<- re
        re <- as.list(data.frame(t(re)))
        refu <- apply(re2,2,function(x){
          x <- x
          # print(x)
          x <- x
          xf <- STFUN[match(x,SpecID)]$modSeq
          fufi <- fun_pairwisecomparison(xf[1],xf[2])
          fufi2 <- fun_pairwisecomparison(xf[2],xf[1])
          if(length(fufi)>0|length(fufi2)>0){
            da <- as.data.frame(rbind(cbind(fufi,xf[1]),cbind(fufi2,xf[2])))
            names(da) <- c("fragments","seq")
          }else{
            da <- NULL
          }
          
          da
        })
        refu2 <<- refu
        if(length(refu)>0){
          # print(refu)
          # print(lengths(refu))
          # print("rbindlistStart")
          goodfragments <- rbindlist(refu)
          # print("rbindlistStopp")
          
          goodfragments[,{
            # print("HU")
            temp <- .SD
            grp <- .BY
            for(x in ST[`Modified sequence`==grp$seq,]$SpecID){
              sel1 <<- !is.na(match(TT$SpecID,x))
              sel2 <<- !is.na(match(TT[sel1,]$Matches,temp$fragments))
              
              tempfun$TToptimizedTemp$Specific[sel1][sel2] <- T
            }
            # print("HU3")
            NULL
          },seq]
          # print("OtherStuff")
        }
        
        
      }else{
        re <- NULL
      }
      # re <<- re
      re
    },Sequence]
 
    TToptimized <- tempfun$TToptimizedTemp
    TToptimized <- TToptimized[Specific==TRUE,]
    
    SToptimized <- ST[!is.na(match(ST$SpecID,TToptimized$SpecID))]

  }
  
  #rewriting Library:
  # fwrite(ST,paste(InclusionlistPath,"SpectraTableOriginal.txt",sep = "/"),sep = "\t")
  # fwrite(TT,paste(InclusionlistPath,"TransitionTableOriginal.txt",sep = "/"),sep = "\t")
  fwrite(SToptimized,paste(InclusionlistPath,"SpectraTable_specific.txt",sep = "/"),sep = "\t")
  fwrite(TToptimized,paste(InclusionlistPath,"TransitionTable_specific.txt",sep = "/"),sep = "\t")
  dbp <- paste(InclusionlistPath,"PRM_Analyzer_Matches/PickyAnalyzer.sqlite",sep = "/")
  if(file.exists(dbp)){
    try({
      db <- dbConnect(SQLite(),dbp)
      dbWriteTable(db,"ST_Specific",SToptimized,overwrite=T)
      dbWriteTable(db,"TT_Specific",TToptimized,overwrite=T)
      dbl <- dbListTables(db)
      dbla <- grep("^mz",dbl,value = T)
      dblatemp <- rbindlist(lapply(strsplit(dbla,"_"),function(x){(data.frame(mz=gsub("mz","",x[1]),z=x[2],se=x[3],id=x[length(x)]))}))
      dblatemp$FastaName <- dbla
      # venn(list(dblatemp$id[grep("Specific",dblatemp$se,invert=T)],unique(TToptimized$SpecID)))
      
      dblatemp[,{
        grp <- .BY
        cat("\r",grp$id)
        TTfound <- TToptimized[TToptimized$SpecID==grp$id]
        if(dim(TTfound)[1]>0){
          dbWriteTable(db,gsub("^mz","SpecificFragments#",grp$FastaName),TTfound,overwrite=T)
        }
        NULL
      },.(id,FastaName)]
     
    })

    
  }
  print("Finished finding specific fragments.")
  
}

# pairwise comparison function, spits out unique fragments. 
fun_pairwisecomparison <- function(Spectratable_seq1, Spectratable_seq2,SearchMod){
  # Author: Mirjam Van Bentum 2021, modified by Henrik Zauber
  # syntax Spectratable_seq1 _IADPEHDHTGFLTEY(Phospho (STY))VATR_
  #read in characteristics sequences
  
  if(grepl("Acetyl|Oxidation",Spectratable_seq1)|grepl("Acetyl|Oxidation",Spectratable_seq2)){
    return(NULL)
  }
  
  Spectratable_seq1 <- gsub("\\(Phospho \\(STY\\)\\)", 
                            "p", Spectratable_seq1)
  seq1 <- gsub("_", "", Spectratable_seq1)
  
  loc_1 <- regexpr("[STY]p", seq1)[[1]][1]
  aa1 <- substr(seq1, loc_1, 
                loc_1 +1)
  
  Spectratable_seq2 <- gsub("\\(Phospho \\(STY\\)\\)", 
                            "p", Spectratable_seq2)
  seq2 <- gsub("_", "", Spectratable_seq2)
  
  loc_2 <- regexpr("[STY]p", seq2)[[1]][1]
  aa2 <- substr(seq1, loc_2, 
                loc_2 +1)
  
  loc_high <- max(loc_1, loc_2)
  loc_low <- min(loc_1, loc_2)
  
  seq_len <- nchar(gsub("p", "", seq2))
  
  #find unique y ions: 
  ymin <- seq_len - loc_high + 1 
  ymax <- seq_len - loc_low
  y_unique <- paste0("y", ymin:ymax)
  
  #find unique b ions
  bmin <- loc_low
  bmax <- loc_high -1
  b_unique <- paste0("b", bmin:bmax)
  
  fragmentvector <- c(y_unique, b_unique)
  
  # add other fragment variants
  # fragment naming from MQ
  fragmentvector <- unique(c(fragmentvector,
                             paste0(fragmentvector, "*"),
                             paste0(fragmentvector, "-NH3"),
                             paste0(fragmentvector, "-H2O"),
                             paste0(fragmentvector, "(2+)"),
                             gsub("b", "a", fragmentvector)))
  
  # pY ion unique identifier? 
  if((aa1 == "Yp"|aa2 == "Yp") & aa1 != aa2){
    fragmentvector <- c(fragmentvector, "pY")
    return(fragmentvector)
  }
  
  return(fragmentvector) }
