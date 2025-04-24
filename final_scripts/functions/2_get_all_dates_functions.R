checkAndInstall <- function(required.packages) {
  uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(uninstalled.packages) > 0) install.packages(uninstalled.packages)
  sapply(required.packages, function(x) eval(parse(text = paste("require(",x,")"))))
  if (any(required.packages == "oxcAAR")) {
    quickSetupOxcal(path = paste0(c(.libPaths()[1],"/oxcAAR/"),collapse = ""))}
}

checkAndInstall(required.packages = c("neotoma2","here","rlist"))

get_from_tilia <- function(values, params, meth) {
  require(jsonlite);require(httr)
  paramstring <- paste0("&_",params,"=",values,collapse = "")
  url <- paste0(c("https://tilia.neotomadb.org/api/?method=ti.", meth, paramstring), collapse = "")
  obj <- try(GET(url), silent = TRUE)
  return(fromJSON(content(obj,as = "text")))
}