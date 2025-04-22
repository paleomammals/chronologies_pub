checkAndInstall <- function(required.packages) {
  uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(uninstalled.packages) > 0) install.packages(uninstalled.packages)
  sapply(required.packages, function(x) eval(parse(text = paste("require(",x,")"))))
  if (any(required.packages == "oxcAAR")) {
    quickSetupOxcal(path = paste0(c(.libPaths()[1],"/oxcAAR/"),collapse = ""))}
}

checkAndInstall(required.packages = c("oxcAAR","neotoma2","here","httr","jsonlite","dplyr","rlist"))

i_am("workflow/scripts/functions.R")

#save
save.image(here(".RData"))

assembleDateString <- function(x) {
  paste0("R_Date('",as.numeric(x["geochronid"]),"[",x["labnumber"],"]',",x["age"],",",
         mean(as.numeric(c(x["errorolder"],x["erroryounger"])),na.rm = T),");",collapse = "")
}

assembleCDateString <- function(x) {
  paste0("C_Date('",x["geochronid"],"[",x["labnumber"],", ",x["geochrontypeid"],"]',",x["age"],",",
         mean(as.numeric(c(x["errorolder"],x["erroryounger"])),na.rm = T),");",collapse = "")
}

unitDates <- function(data,auid="all"){
  if (auid == "all") auid <- unique(data$analysisunitid)
  dateString <- vector(mode = "list",length = length(auid))
  for (i in 1:length(auid)) {
    rdates <- subset(data,data$analysisunitid %in% auid[i] & 
                       (data$agetypeid == "14C yr BP" | data$agetypeid == "Radiocarbon years BP" | data$agetypeid == "4"))
    rdates <- rdates[rev(order(rdates$age)),]
    if (nrow(rdates) > 1) {
      dateString[i] <- paste0(apply(rdates,1,assembleDateString),collapse = "\n      ")
    } else if (nrow(rdates) == 1) {
      dateString[i] <- assembleDateString(rdates)  
    }
  }
  codeString <- vector(mode = "list",length = length(auid))
  for (i in 1:length(auid)) {
      codeString[i] <- paste0("\n     Boundary(\"beginAnalysisUnit\");\n     Phase(\"Analysis Unit ", 
                              auid[i], 
                              "\")\n     {\n      ", 
                              dateString[i], 
                              "\n      Date(\"Event\");\n     };\n     Boundary(\"endAnalysisUnit\");",
                              collapse = "")
    }
  reso <- max(c(10,10^(floor(log10(max(c(data$errorolder,data$erroryounger),na.rm=T))))/10))
  code <- paste0("Options()\n   {\n    BCAD=FALSE;\n    Curve(\"IntCal20\",\"IntCal20.14c\");",
                 paste0("\n    Resolution=",reso,";",collapse=""),
                 "\n   };\n   Plot()\n   {\n    Sequence(\"Collection Unit ",
                 data$collectionunitid[1],
                 "\")\n    {",
                 paste0(codeString,collapse = "\n"),
                 "\n    };\n   };\n  ",
                 collapse = "")
  execute <- executeOxcalScript(code)
  result <- parseOxcalOutput(readOxcalOutput(execute),only.R_Date = F)
  return(result)
}

unitDates.code <- function(data, ordered){
  auid <- unique(data$analysisunitid)
  #returns only code, not modeling results
  #does NOT contain collection beginning and ending boundaries
  #for use in complex models
  dateString <- vector(mode = "list",length = length(auid))
  for (i in 1:length(auid)) {
    rdates <- subset(data,data$analysisunitid %in% auid[i] & 
                       (data$agetypeid == "14C yr BP" | data$agetypeid == "Radiocarbon years BP" | data$agetypeid == "4"))
    rdates <- rdates[rev(order(rdates$age)),]
    if (nrow(rdates) > 1) {
      dateString[i] <- paste0(apply(rdates,1,assembleDateString),collapse = "\n      ")
    } else if (nrow(rdates) == 1) {
      dateString[i] <- assembleDateString(rdates)  
    }
  }
  codeString <- vector(mode = "list",length = length(auid))
  for (i in 1:length(auid)) {
    codeString[i] <- paste0("\n    ",
                            if (ordered == F) {
                              paste0("Sequence(\"AU ",auid[i],"\")\n   {\n    ",collapse = "")
                            },
                            "Boundary(beginAnalysisUnit);\n    Phase(\"Analysis Unit ",auid[i],"\")\n     {\n      ", 
                            dateString[i], 
                            "\n      Date(\"Event\");\n     };",
                            "\n    Boundary(endAnalysisUnit);",
                            if (ordered == F) {
                              "\n    };"
                            },
                            collapse = "")
  }
  return(codeString)
}

trapezoid <- function(p){
  #returns area under a trapezoid
  #p must have arguments in order: x1, y1, x2, y2
  return(0.5 * abs(p[2] + p[4]) * abs(p[1] - p[3]))
}

bc2bp <- function(x) {
  #adjusts for oxcAAR BC/BP date reporting bug
  #subtracts 1950 and returns date as a positive value
  return(round(1950 - x))
}

medianProbDist <- function(xy, datename = "dates", probname = "probabilities", quantile=0.5) {
  #new version using cumsum, much more efficient
  #converts result to BP
  #xy <- x$posterior_probabilities
  xy$cum <- unlist(cumsum(xy[probname]))
  median <- min(which(xy$cum >= quantile*max(xy$cum)))
  result <- mean(xy[c(median,median - 1),datename])
  result <- bc2bp(result)
  return(result)
}

extractChronology <- function(full){
  #wrapper for medianProbDist that takes unmodified result of unitDates()
  if (is.null(full[[1]])) {
    result <- data.frame(begin.bound = NA,begin.event = NA,middle = NA,end.event = NA,end.bound = NA)
  } else {
    begin.bound <- medianProbDist(full$begin$posterior_probabilities)
    begin.event <- medianProbDist(full$Event$posterior_probabilities,quantile = 0.3173)
    middle <- medianProbDist(full$Event$posterior_probabilities)
    end.event <- medianProbDist(full$Event$posterior_probabilities,quantile = 0.6827)
    end.bound <- medianProbDist(full$end$posterior_probabilities)
    result <- data.frame(begin.bound, begin.event, middle, end.event, end.bound)
  }
  return(result)
}

extractChronology.multiple <- function(full) {
#wrapper for medianProbDist that takes unmodified result of oxcalResults()
#  chron <- data.frame(collunitid = substr(names(full)[length(full)],17,nchar(names(full)[length(full)])),
#                      chronologies.ageboundolder = medianProbDist(full$begin$posterior_probabilities),
#                      chronologies.ageboundyounger = medianProbDist(full$end$posterior_probabilities))
  AUages <- data.frame(analysisunitid = rep(as.numeric(sapply(grep("^Analysis",names(full),value = T),
                      function(x) substr(x,15,nchar(x)))),2))
  AUages$chronology.agemodel <- c(rep("event",length(AUages$analysisunitid)/2),
                                      rep("bounds",length(AUages$analysisunitid)/2))
  AUages$sampleages.age <- c(sapply(full[grepl("Event$",names(full))],
                                    function(x) medianProbDist(x$posterior_probabilities,quantile = 0.5)),
                             rep(NA,length(which(AUages$chronology.agemodel == "bounds"))))
  AUages$sampleages.ageolder <- c(sapply(full[grepl("Event",names(full))],
                                function(x) medianProbDist(x$posterior_probabilities,quantile = 0.3173)),
                              sapply(full[grepl("begin",names(full))],
                               function(x) medianProbDist(x$posterior_probabilities,quantile = 0.5)))
  AUages$sampleages.ageyounger <- c(sapply(full[grepl("Event",names(full))],
                                function(x) medianProbDist(x$posterior_probabilities,quantile = 0.6827)),
                              sapply(full[grepl("end",names(full))],
                                function(x) medianProbDist(x$posterior_probabilities,quantile = 0.5)))
  controls <- data.frame(geochronid = as.numeric(sapply(grep("^[0-9]",names(full),value = T),
                                                        function(x) strsplit(x,"\\[")[[1]][1])))
  controls$geochronology.labnumber <- sapply(grep("^[0-9]",names(full),value = T), 
                                             function(x) strsplit(strsplit(x,split = "\\[")[[1]],split = "\\]")[[2]])
  controls$chroncontrols.age <- sapply(full[grepl("^[0-9]",names(full))],
                       function(x) medianProbDist(x$posterior_probabilities,quantile = 0.5))
  controls$chroncontrols.agelimitolder <- sapply(full[grepl("^[0-9]",names(full))],
                        function(x) medianProbDist(x$posterior_probabilities,quantile = 0.0455))
  controls$chroncontrols.agelimityounger <- sapply(full[grepl("^[0-9]",names(full))],
                       function(x) medianProbDist(x$posterior_probabilities,quantile = 0.9545))
  return(list(sampleages = AUages,chroncontrols = controls))
}

cleandates <- function(data){
  drop <- which(is.na(data$age) | is.na(data$errorolder) | is.na(data$erroryounger) |
                data$age < 10 | data$errorolder == 0 |
                data$infinite == TRUE | data$human == TRUE | data$rejected == TRUE)
  if (length(drop) > 0) {result <- data[-drop,]} else {result <- data}
  return(result)
}

dateranges <- function(collid,data) {
  coll_dates <- subset(data,data$collectionunitid == collid)
  coll_dates <- cleandates(coll_dates)
  print(paste0(c("Site ",coll_dates$siteid[1],", collection ",collid,": ", nrow(coll_dates)," dates"),collapse = ""))
  if (nrow(coll_dates) == 0) return()
  oxcalResult <- unitDates(coll_dates)
  chronology <- extractChronology(oxcalResult)
  AUage <- data.frame(
    analysisunitid = rep(coll_dates$analysisunitid[1],2),
    chronology.agemodel = c("event","bounds"),
    sampleages.age = c(chronology$middle,NA),
    sampleages.ageolder = c(chronology$begin.event,chronology$begin.bound),
    sampleages.ageyounger = c(chronology$end.event,chronology$end.bound)
  )
  controls <- data.frame(unname(t(sapply(grep("^[0-9]",names(oxcalResult),value = T),
                                         function(x) strsplit(x,"\\[|\\]")[[1]][1:2]))))
  colnames(controls) <- c("geochronid","geochronology.labnumber")
  controls$chroncontrols.age <- sapply(oxcalResult[grepl("^[0-9]",names(oxcalResult))],
                                       function(x) medianProbDist(x$posterior_probabilities,quantile = 0.5))
  controls$chroncontrols.agelimitolder <- sapply(oxcalResult[grepl("^[0-9]",names(oxcalResult))],
                                                 function(x) medianProbDist(x$posterior_probabilities,quantile = 0.0455))
  controls$chroncontrols.agelimityounger <- sapply(oxcalResult[grepl("^[0-9]",names(oxcalResult))],
                                                   function(x) medianProbDist(x$posterior_probabilities,quantile = 0.9545))
  return(list(sampleages = AUage,chroncontrols = controls))
}

dateranges.multi <- function(data, ordered = "unknown") {
  #if units are ordered and non-NA, sort them
  if (ordered == "unknown") {
    if (is.null(data[[1]]$auorder) | any(is.na(data[[1]]$auorder))) {
      ordered <- F
      data[[1]]$analysisunitid <- factor(data[[1]]$analysisunitid, ordered = F)
    } else {
      ordered <- T
      order <- distinct(data[[1]][,c("analysisunitid","auorder")])
      order <- order[order(order$auorder),]
      data[[1]]$analysisunitid <- factor(data[[1]]$analysisunitid, 
                                         ordered = T, levels = order$analysisunitid)
      data[[1]] <- data[[1]][order(data[[1]]$analysisunitid),]
    }
  }
  data[[1]] <- cleandates(data[[1]])
  print(paste0(c("Site ",data[[1]]$siteid[1],", collection ",data[[1]]$collectionunitid[1],": ", 
                 nrow(data[[1]])," dates, ",length(unique(data[[1]]$analysisunitid)),
                 " analysis units, ordered = ",ordered),collapse = ""))
  if (nrow(data[[1]]) == 0) return()
  #assemble oxcal code
  reso <- max(c(10,10^(floor(log10(max(c(data[[1]]$errorolder,data[[1]]$erroryounger),na.rm=T))))/10))
  if (ordered == T) {
    code <- paste0("Options()\n {\n  BCAD=FALSE;\n  Curve(\"IntCal20\",\"IntCal20.14c\");",
                   paste0("\n  Resolution=",reso,";\n };",collapse=""),
                   "\n Plot()\n {\n  Sequence(\"Collection Unit ", names(data),"\")",
                   "\n  {",
                   paste0(unitDates.code(data[[1]], ordered = ordered),collapse = ""),
                   "\n  };\n };",collapse = "")
  } else if (ordered == F) {
    code <- paste0("Options()\n {\n  BCAD=FALSE;\n  Curve(\"IntCal20\",\"IntCal20.14c\");",
                   paste0("\n  Resolution=",reso,";\n };",collapse=""),
                   "\n Plot()\n {\n  Phase(\"Collection Unit ", names(data),"\")",
                   "\n  {",
                   paste0(unitDates.code(data[[1]], ordered = ordered),collapse = ""),
                   "\n  };\n };",collapse = "")
  }
  #run oxcal code and parse results
  execute <- executeOxcalScript(code)
  oxcalResult <- parseOxcalOutput(readOxcalOutput(execute),only.R_Date = F)
  result <- extractChronology.multiple(oxcalResult)
  return(result)
}

IsSmallerOrEqual <- function(a,b) {   
  if (   class(all.equal(a, b)) == "logical" && (a < b | all.equal(a, b))) { return(TRUE)
  } else if (a < b) { return(TRUE)
  } else { return(FALSE) }
}
IsBiggerOrEqual <- function(a,b) {
  if (   class(all.equal(a, b)) == "logical" && (a > b | all.equal(a, b))) { return(TRUE)
  } else if (a > b) { return(TRUE)
  } else { return(FALSE) }
}
IsEqual <- function(a,b) {
  if (   class(all.equal(a, b)) == "logical" ) { return(TRUE)
  } else { return(FALSE) }
}
checkLocation <- function(x){
  if (IsBiggerOrEqual(as.numeric(x["longitude"]),as.numeric(x["longitudeeast"])) 
     & IsSmallerOrEqual(as.numeric(x["longitude"]),as.numeric(x["longitudewest"])) 
     & IsSmallerOrEqual(as.numeric(x["latitude"]),as.numeric(x["latitudenorth"]))
     & IsBiggerOrEqual(as.numeric(x["latitude"]),as.numeric(x["latitudesouth"]))) {
    return(TRUE)} else{return(FALSE)}
}

getDistance <- function(x){
  #  a<-round(abs(mean(c(as.numeric(x["latitudesouth"]),as.numeric(x["latitudenorth"])),na.rm=T)-as.numeric(x["latitude"]))*110.574,5)
  #  b<-round(abs(cos(mean(c(as.numeric(x["longitudeeast"]),as.numeric(x["longitudewest"])),na.rm=T))-cos(as.numeric(x["longitude"])))*111.320,5)
  a <- abs(min(as.numeric(x["latitudenorth"]) - as.numeric(x["latitude"]),as.numeric(x["latitudenorth"]) - as.numeric(x["latitude"])))
  b <- abs(min(as.numeric(x["longitudeeast"]) - as.numeric(x["longitude"]),as.numeric(x["longitudewest"]) - as.numeric(x["longitude"])))
  return(c(a,b))
}

uncalMatch <- function(caldate, curve = intcal20, errordirection){
  bounds <- c(max(curve[which(curve$calbp <= caldate),"uncalbp"]),
            min(curve[which(curve$calbp >= caldate),"uncalbp"]))
  sigma <- mean(curve[which.max(curve$calbp <= caldate),"sigma"],
              curve[which.min(curve$calbp >= caldate) - 1,"sigma"])
  if (errordirection == "+") {
    return(round(mean(bounds + sigma),-1))
  } else if (errordirection == "-") {
    return(round(mean(bounds - sigma),-1)) 
  } else if (errordirection == "0") {
    return(round(mean(bounds),-1))
  }
}

getagelimits.bycoll <- function(obj,i,j){
  cu <- obj[[i]]$collunits[[j]]$collectionunitid[[1]]
  datelims <- as.data.frame(chronologies(obj[[i]]$collunits[[j]]))[,c("ageboundolder","ageboundyounger")]
  return(data.frame(collectionunitid = cu,
                    agemax = max(datelims$ageboundolder),
                    agemin = min(datelims$ageboundyounger)))
}

findTop <- function(rowno, tab, IDcol = 1, topIDcol = 3){
  #this function cannot be used with apply()
  curr <- tab[rowno,IDcol]; top <- tab[rowno,topIDcol]
  if (curr != top) {
    while (curr != top) {
      curr <- tab[tab[,IDcol] == top,topIDcol]
      top <- tab[curr,topIDcol]
    }}
  return(curr)
}

getvariables.bydataset <- function(id) {
  url <- paste0(c("https://tilia.neotomadb.org/api/?method=ti.getdatasetvariables&_datasetid=",
                id),collapse = "")
  obj <- fromJSON(getURL(url))
  vals <- obj$data$taxoncode
  vals <- na.omit(smalltaxa[match(vals,smalltaxa$taxoncode),"taxonid"])
  return(data.frame(datasetid = rep(id,length(vals)),taxonid = vals))
}

checkIntervals <- function(x,int = intervals.14Cybp){ 
  result <- as.numeric(x["agemax"]) >= int$end & 
    as.numeric(x["agemin"]) <= int$start
  names(result) <- int$name
  return(result)
}

uniqueByColumn <- function(x) {apply(x,2,function(x) length(unique(x)))}

renumber <- function(x) {rownames(x) <- 1:nrow(x); return(x)}

get_from_tilia <- function(values, params, meth) {
  require(jsonlite);require(httr)
  paramstring <- paste0("&_",params,"=",values,collapse = "")
  url <- paste0(c("https://tilia.neotomadb.org/api/?method=ti.", meth, paramstring), collapse = "")
  obj <- try(GET(url), silent = TRUE)
  return(fromJSON(content(obj,as = "text")))
}

get_analysisunits <- function(siteid) {
  require(neotoma2)
  colls <- unique(getids(get_sites(siteid))$collunitid)
  aus <- unname(unlist(sapply(colls,function(x) get_from_tilia(x,"collectionunitid","getanalysisunit")$data$analysisunitid)))
  au.names <- data.frame(aus,as.data.frame(t(sapply(aus, function(x) get_from_tilia(x,"analunitid","getanalysisunitbyid")$data))))
  return(au.names)
}

makeOxcalCode <- function(data, auid = "all"){
  data <- cleandates(data)
  if (auid == "all") auid <- unique(data$analysisunitid)
  dateString <- vector(mode = "list",length = length(auid))
  for (i in 1:length(auid)) {
    rdates <- subset(data,data$analysisunitid == auid[i] & data$agetypeid == "Radiocarbon years BP")
    cdates <- subset(data,data$analysisunitid == auid[i] & data$agetypeid != "Radiocarbon years BP")
    if (nrow(rdates) > 1) {
      dateString[i] <- paste0(apply(rdates,1,assembleDateString),collapse = "\n      ")
    } else if (nrow(rdates) == 1) {
      dateString[i] <- assembleDateString(rdates)  
    }
    #C_Date does not work in oxcAAR
    if (nrow(cdates) > 0) {
      if (nrow(cdates) > 1) {
        cdateString <- paste0(apply(cdates,1,assembleCDateString),collapse = "\n      ")
      } else { 
        cdateString <- assembleCDateString(cdates)
      }
      if (nrow(rdates) > 0) {
        dateString[i] <- paste0(c(dateString[i],cdateString),collapse = "\n      ")
      } else {
        dateString[i] <- cdateString
      }
    }
  }
  codeString <- vector(mode = "list",length = length(auid))
  for (i in 1:length(auid)) {
    codeString[i] <- paste0("\n    Sequence(\"AU ",auid[i],"\")\n   {\n    ",
                            "Boundary(beginAnalysisUnit);\n    Phase(\"Analysis Unit ",auid[i],"\")\n     {\n      ", 
                            dateString[i], 
                            "\n      Date(\"Event\");\n     };",
                            "\n    Boundary(endAnalysisUnit);\n    };",
                            collapse = "")
  }
  reso <- max(c(10,10^(floor(log10(max(c(data$errorolder,data$erroryounger),na.rm=T))))/10))
  code <- paste0("Options()\n {\n  BCAD=FALSE;\n  Curve(\"IntCal20\",\"IntCal20.14c\");",
                 paste0("\n  Resolution=",reso,";",collapse=""),
                 "\n };\n Plot()\n {\n  Phase(\"Collection Unit ", unique(data$collectionunitid),"\")",
                 "\n  {",
                 paste0(codeString,collapse = ""),
                 "\n  };\n };",collapse = "")
  return(code)
}

extractChronology.Table <- function(data, idtype = "geochronid") {
  L <- split(data,data$index)
  names(L) <- apply(distinct(data[,1:3])[,2:3],1,paste0,collapse = " ")
  dateNames <- data.frame(original = grep("[0-9]",names(L),value = T))
  controls <- data.frame(id = as.numeric(sapply(dateNames$original,
                 function(x) strsplit(strsplit(x,"_Date ")[[1]],"\\[")[[2]][1])))
  rownames(controls) <- 1:nrow(controls)
  controls$geochronology.labnumber <- sapply(dateNames$original, 
                                             function(x) gsub("_"," ",strsplit(strsplit(x,split = "\\[")[[1]],split = "\\]")[[2]]))
  controls$chroncontrols.age <- sapply(L[grepl("[0-9]",names(L))],
                                       function(x) medianProbDist(subset(x,x$type == "posterior"),quantile = 0.5, 
                                                                  datename = "value", probname = "probability"))
  controls$chroncontrols.agelimitolder <- sapply(L[grepl("[0-9]",names(L))],
                                                 function(x) medianProbDist(subset(x,x$type == "posterior"),quantile = 0.0455, 
                                                                            datename = "value", probname = "probability"))
  controls$chroncontrols.agelimityounger <- sapply(L[grepl("[0-9]",names(L))],
                                                   function(x) medianProbDist(subset(x,x$type == "posterior"),quantile = 0.9545, 
                                                                              datename = "value", probname = "probability"))

  AUages <- data.frame(analysisunitid = rep(sort(unique(subset(alldates,alldates[[idtype]] %in% controls$id)$analysisunitid)),2))
  AUages$chronology.agemodel <- c(rep("event",length(AUages$analysisunitid)/2),
                                  rep("bounds",length(AUages$analysisunitid)/2))
  AUages$sampleages.age <- c(sapply(L[grepl("^Calculate Event",names(L))],
                                    function(x) medianProbDist(subset(x,x$type == "posterior"),quantile = 0.5, 
                                                               datename = "value", probname = "probability")),
                             rep("NA",length(AUages$analysisunitid)/2))
  AUages$sampleages.ageolder <- c(sapply(L[grepl("^Calculate Event",names(L))],
                                         function(x) medianProbDist(subset(x,x$type == "posterior"),quantile = 0.3173, 
                                                                    datename = "value", probname = "probability")),
                                  sapply(L[grepl("begin",names(L))],
                                         function(x) medianProbDist(subset(x,x$type == "posterior"),quantile = 0.5, 
                                                                    datename = "value", probname = "probability")))
  AUages$sampleages.ageyounger <- c(sapply(L[grepl("^Calculate Event",names(L))],
                                           function(x) medianProbDist(subset(x,x$type == "posterior"),quantile = 0.6827, 
                                                                      datename = "value", probname = "probability")),
                                    sapply(L[grepl("end",names(L))],
                                           function(x) medianProbDist(subset(x,x$type == "posterior"),quantile = 0.5, 
                                                                      datename = "value", probname = "probability")))
  colnames(controls)[1] <- idtype
  result <- list(sampleages = AUages,chroncontrols = controls)
  return(result) 
}

duplicates <- function(x){
  return(which(x%in%x[duplicated(x)]))
}

join_sampleages <- function(coll) {
  require(dplyr)
  result <- left_join(get_from_tilia(coll,"collunitid","getanalysisunitsbycollunitid")$data[,1:2],
            subset(newchron.all,newchron.all$collectionunitid==coll)[,2:6])
  result <- cbind(rep(coll,nrow(result)),result)
  return(result)
}


number_aus <- function(data) {
  #analysis units in data must be in correct order
  tt <- data.frame(distinct(data[,c("analysisunitid","analysisunitname")]),
                   number=1:length(unique(data$analysisunitname)))
  result <- data
  result$order <- tt$number[match(result$analysisunitid,tt$analysisunitid)]
  return(result)
}

trim_unitages <- function(data) {
  if(is.null(data$order)) {print("no order provided");break}
  result <- subset(data,!is.na(data$chronology.agemodel))
  result$trimolder <- result$trimyounger <- NA
  for (i in 1:nrow(data)) {
    result[i,]$trimyounger <- max(subset(data,data$order <= i)$sampleages.ageyounger)
    result[i,]$trimolder <- min(subset(data,data$order >= i)$sampleages.ageolder)
  }
  return(result)
}

interpolate_ages <- function(data) {
  data <- number_aus(data)
  data <- split(data,data$chronology.agemodel)
}

checkRecent <- function(x,age="1 month") {
  if (max(as.Date(x$recdatemodified), na.rm = T) > 
      seq(Sys.Date(), length = 2, by = paste0(c("-",age),collapse=""))[2]) {
    return("TRUE")
  } else return("FALSE")
}

loadGeochronTable <- function() {
  i_am("workflow/scripts/functions.R")
  if (file.exists("workflow/data/ndb_geochronology.csv")) {
    geochronology.table <- read.csv("workflow/data/ndb_geochronology.csv")
  }
  if (exists("geochronology.table")) {if (checkRecent(geochronology.table) == F) {rm(geochronology.table)}}
  if (!exists("geochronology.table")) {
    print("Downloading table 'geochronology'...")
    geochronology.table <- distinct(get_table("geochronology", limit = 99999999))
    write.csv(geochronology.table, "workflow/data/ndb_geochronology.csv", row.names = F)
  }
  return(geochronology.table)
}
