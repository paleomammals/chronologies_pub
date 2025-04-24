checkAndInstall <- function(required.packages) {
  uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(uninstalled.packages) > 0) install.packages(uninstalled.packages)
  sapply(required.packages, function(x) eval(parse(text = paste("require(",x,")"))))
  if (any(required.packages == "oxcAAR")) {
    quickSetupOxcal(path = paste0(c(.libPaths()[1],"/oxcAAR/"),collapse = ""))}
}

checkAndInstall(required.packages = c("neotoma2","here","dplyr","rlist"))

getbounds <- function(collmod) {
  min <- min(collmod["sampleages.ageyounger"])
  max <- max(collmod["sampleages.ageolder"])
  return(c(max,min))
}

dfbounds <- function(coll) {
  bounds <- data.frame(t(do.call(cbind,tapply(coll,coll$chronology.agemodel,getbounds))))
  result <- data.frame(collectionunitid = coll$collectionunitid[1],
                       agemodel = rownames(bounds), 
                       bounds, row.names = NULL)
  colnames(result)[3:4] <- c("ageboundolder","ageboundyounger")
  return(result)
}
