checkAndInstall <- function(required.packages) {
  uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(uninstalled.packages) > 0) install.packages(uninstalled.packages)
  sapply(required.packages, function(x) eval(parse(text = paste("require(",x,")"))))
  if (any(required.packages == "oxcAAR")) {
    quickSetupOxcal(path = paste0(c(.libPaths()[1],"/oxcAAR/"),collapse = ""))}
}

checkAndInstall(required.packages = c("neotoma2","here","rlist"))

checkRecent <- function(x) {
  if (max(as.Date(x$recdatemodified), na.rm = T) > seq(Sys.Date(), length = 2, by = "-1 month")[2]) {
    return("TRUE")
  } else return("FALSE")
}

cleandates <- function(data){
  drop <- which(is.na(data$age) | is.na(data$errorolder) | is.na(data$erroryounger) |
                  data$age < 10 | data$errorolder == 0 |
                  data$infinite == TRUE | data$human == TRUE | data$rejected == TRUE)
  if (length(drop) > 0) {result <- data[-drop,]} else {result <- data}
  return(result)
}

get_from_tilia <- function(values, params, meth) {
  require(jsonlite);require(httr)
  paramstring <- paste0("&_",params,"=",values,collapse = "")
  url <- paste0(c("https://tilia.neotomadb.org/api/?method=ti.", meth, paramstring), collapse = "")
  obj <- try(GET(url), silent = TRUE)
  return(fromJSON(content(obj,as = "text")))
}