checkAndInstall <- function(required.packages) {
  uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(uninstalled.packages) > 0) install.packages(uninstalled.packages)
  sapply(required.packages, function(x) eval(parse(text = paste("require(",x,")"))))
  if (any(required.packages == "oxcAAR")) {
    quickSetupOxcal(path = paste0(c(.libPaths()[1],"/oxcAAR/"),collapse = ""))}
}
get_from_tilia <- function(values, params, meth) {
  require(jsonlite);require(httr)
  paramstring <- paste0("&_",params,"=",values,collapse = "")
  url <- paste0(c("https://tilia.neotomadb.org/api/?method=ti.", meth, paramstring), collapse = "")
  obj <- try(GET(url), silent = TRUE)
  return(fromJSON(content(obj,as = "text")))
}
cleandates <- function(data){
  drop <- which(is.na(data$age) | is.na(data$errorolder) | is.na(data$erroryounger) |
                  data$age < 50 | data$errorolder == 0 |
                  data$infinite == TRUE | data$human == TRUE | data$rejected == TRUE)
  if (length(drop) > 0) {result <- data[-drop,]} else {result <- data}
  return(result)
}
cleandates_report <- function(data){
  drop_counts <- data.frame(infinite_or_modern = length(which(data$age < 50  | data$infinite == TRUE)),
                            noage_or_error = length(which(is.na(data$age) | 
                                                            is.na(data$errorolder) | 
                                                            is.na(data$erroryounger) | 
                                                            data$errorolder == 0)),
                            rejected = length(which(data$rejected == TRUE)),
                            human = length(which(data$human == TRUE))
                  )
  return(drop_counts)
}
