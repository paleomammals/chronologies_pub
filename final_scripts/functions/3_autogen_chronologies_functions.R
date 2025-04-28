checkAndInstall <- function(required.packages) {
  uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(uninstalled.packages) > 0) install.packages(uninstalled.packages)
  sapply(required.packages, function(x) eval(parse(text = paste("require(",x,")"))))
  if (any(required.packages == "oxcAAR")) {
    quickSetupOxcal(path = paste0(c(.libPaths()[1],"/oxcAAR/"),collapse = ""))}
}

checkAndInstall(required.packages = c("oxcAAR","dplyr","here","rlist"))

cleandates <- function(data){
  drop <- which(is.na(data$age) | is.na(data$errorolder) | is.na(data$erroryounger) |
                  data$age < 10 | data$errorolder == 0 |
                  data$infinite == TRUE | data$human == TRUE | data$rejected == TRUE)
  if (length(drop) > 0) {result <- data[-drop,]} else {result <- data}
  return(result)
}

assembleDateString <- function(x) {
  paste0("R_Date('",as.numeric(x["geochronid"]),"[",x["labnumber"],"]',",x["age"],",",
         mean(as.numeric(c(x["errorolder"],x["erroryounger"])),na.rm = T),");",collapse = "")
}

unitDates <- function(data,auid="all"){
  if (auid == "all") auid <- unique(data$analysisunitid)
  dateString <- vector(mode = "list",length = length(auid))
  for (i in 1:length(auid)) {
    rdates <- subset(data,data$analysisunitid %in% auid[i] &
                       rdates <- rdates[rev(order(rdates$age)),])
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

medianProbDist <- function(xy, datename = "dates", probname = "probabilities", quantile = 0.5) {
  #new version using cumsum, much more efficient
  #converts result from BC to BP
  #xy <- x$posterior_probabilities
  xy$cum <- unlist(cumsum(xy[probname]))
  median <- min(which(xy$cum >= quantile*max(xy$cum)))
  result <- mean(xy[c(median,median - 1),datename])
  result <- round(1950 - result)
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

unitDates.code <- function(data, ordered){
  auid <- unique(data$analysisunitid)
  #returns only code, not modeling results
  #does NOT contain collection beginning and ending boundaries
  #for use in complex models
  dateString <- vector(mode = "list",length = length(auid))
  for (i in 1:length(auid)) {
    rdates <- subset(data,data$analysisunitid == auid[i] & 
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

extractChronology.multiple <- function(full) {
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
