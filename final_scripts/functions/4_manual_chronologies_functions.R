checkAndInstall <- function(required.packages) {
  uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(uninstalled.packages) > 0) install.packages(uninstalled.packages)
  sapply(required.packages, function(x) eval(parse(text = paste("require(",x,")"))))
  if (any(required.packages == "oxcAAR")) {
    quickSetupOxcal(path = paste0(c(.libPaths()[1],"/oxcAAR/"),collapse = ""))}
}

checkAndInstall(required.packages = c("oxcAAR","neotoma2","here","httr","jsonlite","dplyr","rlist"))

assembleDateString <- function(x) {
  paste0("R_Date('",as.numeric(x["geochronid"]),"[",x["labnumber"],"]',",x["age"],",",
         mean(as.numeric(c(x["errorolder"],x["erroryounger"])),na.rm = T),");",collapse = "")
}

assembleCDateString <- function(x) {
  paste0("C_Date('",x["geochronid"],"[",x["labnumber"],", ",x["geochrontypeid"],"]',",x["age"],",",
         mean(as.numeric(c(x["errorolder"],x["erroryounger"])),na.rm = T),");",collapse = "")
}

cleandates <- function(data){
  drop <- which(is.na(data$age) | is.na(data$errorolder) | is.na(data$erroryounger) |
                  data$age < 50 | data$errorolder == 0 |
                  data$infinite == TRUE | data$human == TRUE | data$rejected == TRUE)
  if (length(drop) > 0) {result <- data[-drop,]} else {result <- data}
  return(result)
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

extractChronology.Table <- function(data, idtype = "geochronid") {
  require(dplyr)
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
                                                 function(x) medianProbDist(subset(x,x$type == "posterior"),quantile = 0.0455,                                                 datename = "value", probname = "probability"))
  controls$chroncontrols.agelimityounger <- sapply(L[grepl("[0-9]",names(L))],
                                                   function(x) medianProbDist(subset(x,x$type == "posterior"),quantile = 0.9545,                                                 datename = "value", probname = "probability"))
  AUages <- data.frame(analysisunitid = rep(sort(unique(subset(alldates,as.numeric(alldates[[idtype]]) %in% as.numeric(controls$id))$analysisunitid)),2))
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
