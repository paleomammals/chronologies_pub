checkAndInstall <- function(required.packages) {
  uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(uninstalled.packages) > 0) install.packages(uninstalled.packages)
  sapply(required.packages, function(x) eval(parse(text = paste("require(",x,")"))))
  if (any(required.packages == "oxcAAR")) {
    quickSetupOxcal(path = paste0(c(.libPaths()[1],"/oxcAAR/"),collapse = ""))}
}

checkAndInstall(required.packages = c("oxcAAR","neotoma2","here","httr","jsonlite","dplyr","rlist"))

get_sampleids <- function(vert_and_chron) {
  require(httr);require(jsonlite);require(dplyr);require(rlist)
  samples.ti <- vector(mode = "list",length = nrow(vert_and_chron))
  print(paste0(c("Downloading geochron sample data for",nrow(vert_and_chron),"datasets..."),collapse = " "))
  for (i in 1:nrow(vert_and_chron)) {
    temp <- get_from_tilia(values = vert_and_chron$datasetid[i], 
                           params = "datasetid", 
                           meth = "getgeochronanalysisunitsbydatasetid")$data
    if (length(temp) == 0) { temp <- data.frame(sampleid = NA,collectionunitid = NA,
                                                analysisunitid = NA,analysisunitname = NA,
                                                depth = NA,thickness = NA) }
    samples.ti[[i]] <- data.frame(datasetid = vert_and_chron$datasetid[i], temp)
  }
  samples.ids <- distinct(list.stack(samples.ti)[, c("datasetid", "sampleid", "analysisunitid")])
  samples.ids <- subset(samples.ids,!is.na(samples.ids$sampleid))
  print(paste0(c("Downloaded data for",nrow(samples.ids),"samples."),collapse = " "))
  return(samples.ids)
}

get_geochrons <- function(vert_and_chron) {
  require(httr);require(jsonlite);require(dplyr);require(rlist)
  chron <- subset(vert_and_chron,vert_and_chron$datasettype == "geochron")
  geochrons.ti <- vector(mode = "list",length = nrow(chron))
  print(paste0(c("Downloading dates from",nrow(chron),"datasets..."),collapse = " "))
  for (i in 1:nrow(chron)) {
    temp <- get_from_tilia(values = chron$datasetid[i], 
                           params = "datasetid", 
                           meth = "getgeochronbydatasetid")$data
    if (length(temp) == 0) { 
      temp <- data.frame(geochronid = NA, geochrontypeid = NA, geochrontype = NA, 
                         agetype = NA, depth = NA, thickness = NA, analysisunitid = NA, 
                         analysisunitname = NA, age = NA, errorolder = NA, erroryounger = NA,
                         infinite = NA, labnumber = NA, materialdated = NA, notes = NA) }
    geochrons.ti[[i]] <- data.frame(datasetid = chron$datasetid[i], temp)
  }
  geochrons.ti <- distinct(list.stack(geochrons.ti)[, c(1,2,4,5,8:16)])
  geochrons.ti <- subset(geochrons.ti,!is.na(geochrons.ti$geochronid))
  print(paste0(c("Downloaded",nrow(geochrons.ti),"dates."),collapse = " "))
  return(geochrons.ti)
}

get_colldata <- function(vert_and_chron){
  require(httr);require(jsonlite);require(dplyr);require(rlist)
  colls <- unique(vert_and_chron$collectionunitid)
  collections.ti <- vector(mode = "list",length = length(colls))
  print(paste0(c("Downloading data for",length(colls),"collections..."),collapse = " "))
  for (i in 1:length(colls)) {
    temp <- get_from_tilia(values = colls[i], 
                           params = "collectionunitid", 
                           meth = "getcollunitbyid")$data
    if (length(temp) == 0) { 
      temp <- data.frame(collectionunitid = NA,handle = NA,siteid = NA,colltypeid = NA,
                         depenvtid = NA,collunitname = NA,colldate = NA,colldevice = NA,
                         gpslatitude = NA,gpslongitude = NA,gpsaltitude = NA,gpserror = NA,
                         waterdepth = NA,substrateid = NA,slopeaspect = NA,slopeangle = NA,
                         location = NA,notes = NA) }
    collections.ti[[i]] <- data.frame(datasetid = colls[i], temp)
  }
  collections.ti <- distinct(list.stack(collections.ti)[,c("collectionunitid","depenvtid")])
  print(paste0(c("Downloaded data for",nrow(collections.ti),"collections."),collapse = " "))
  return(collections.ti)
}

get_smallvert_datasets <- function(vert_and_chron, site.geochrons,datasets) {
  require(httr);require(jsonlite);require(dplyr);require(rlist)
  if (!file.exists(here("data/smallmammaltaxa.csv"))) {
    source(here("final_scripts/functions/filter_taxa.R"))}
  smalltaxa <- read.csv(here("data/smallmammaltaxa.csv"))
  vert.datasets <- unique(intersect(subset(vert_and_chron,
                                           vert_and_chron$collectionunitid %in% site.geochrons$collectionunitid)[,"datasetid"],
                                    subset(datasets,datasets$datasettype == "vertebrate fauna")$datasetid))
  temp <- vector(mode = "list",length = length(vert.datasets))
  print(paste0(c("Downloading taxon data for",length(vert.datasets),"datasets..."),collapse = " "))
  for (i in 1:length(vert.datasets)) {
    verts.temp <- try(get_from_tilia(values = vert.datasets[i],
                                     params = "datasetid",
                                     meth = "getdatasetvariables")$data) 
    vals <- smalltaxa[na.omit(match(verts.temp[,"taxoncode"],smalltaxa[,"taxoncode"])), "taxonid"]
    temp[[i]] <- data.frame(datasetid = rep(vert.datasets[i],length(vals)),
                            taxonid = vals)
  }
  taxonids <- distinct(rlist::list.stack(temp))
  smallmamm.colls <- unique(subset(vert_and_chron, 
                                   vert_and_chron$datasetid %in% taxonids$datasetid)$collectionunitid)
  print(paste0(c("Found",length(smallmamm.colls),"datasets with small mammals."),collapse = " "))
  return(smallmamm.colls)
}

checkRecent <- function(x) {
  if (max(as.Date(x$recdatemodified), na.rm = T) > seq(Sys.Date(), length = 2, by = "-1 month")[2]) {
    return("TRUE")
  } else return("FALSE")
}

loadGeochronTable <- function() {
  if (file.exists(here("data/ndb_geochronology.csv"))) {
    geochronology.table <- read.csv(here("data/ndb_geochronology.csv"))
  }
  if (exists("geochronology.table")) {
    if (checkRecent(geochronology.table) == F) {rm(geochronology.table)}}
  if (!exists("geochronology.table")) {
    print("Downloading table 'geochronology'...")
    geochronology.table <- distinct(get_table("geochronology", limit = 99999999))
    write.csv(geochronology.table, here("data/ndb_geochronology.csv"), row.names = F)
  }
  return(geochronology.table)
}

get_site_geochrons <- function(siteid, geochronology.table = geochronology.table, agelimit = 30000) {
  require(neotoma2);require(here);require(dplyr)
  site_data <- get_sites(siteid, all_data = T)
  if (is.null(site_data)) {
    print("Site not found, returning NA"); return(NA)  }
  siteinfo <- show(site_data)
  datasets <- as.data.frame(datasets(site_data))
  ids <- inner_join(getids(site_data),datasets[,c("datasetid","datasettype")], by = join_by(datasetid))
  vert_and_chron <- distinct(subset(ids, ids$datasettype %in% c("geochronologic", "vertebrate fauna")))
  #rename column and convert to integer
  if (nrow(vert_and_chron) == 0) {
    print("No datasets found, returning NA"); return(NA)  }
  colnames(vert_and_chron)[which(colnames(vert_and_chron) == "collunitid")] <- "collectionunitid"
  for (i in 1:3) vert_and_chron[,i] <- as.numeric(vert_and_chron[,i])
  vert_and_chron$datasettype <- factor(vert_and_chron$datasettype)
  #rename the dataset types
  levels(vert_and_chron$datasettype) <- c("geochron","vertebrate")
  #use Tilia calls to assemble list of samples
  sampleids <- get_sampleids(vert_and_chron)
  geochrons.ti <- get_geochrons(vert_and_chron)
  collections.ti <- get_colldata(vert_and_chron)
  #download geochronology table if necessary; save it externally
  geochronology.table <- loadGeochronTable()
  #merge all together
  a <- distinct(inner_join(subset(vert_and_chron,vert_and_chron$datasettype == "geochron"),
                           distinct(sampleids[,c(1,3)]),by = "datasetid"))
  b <- distinct(right_join(a, geochrons.ti, by = "analysisunitid"))
  c <- distinct(inner_join(b,data.frame(siteid = as.numeric(siteinfo$siteid),siteinfo = siteinfo$sitename),
                           by = "siteid"))
  d <- distinct(inner_join(c, collections.ti, by = "collectionunitid"))
  e <- distinct(left_join(d,distinct(geochronology.table[,c("geochronid","sampleid","delta13c")]),
                          by = "geochronid"))
  site.geochrons <- e[,c("siteid", "collectionunitid", "sampleid", "datasetid.x", 
                         "analysisunitid", "geochronid", "geochrontype", "agetype", 
                         "age", "errorolder", "erroryounger", "infinite", 
                         "delta13c", "labnumber", "materialdated", "notes", 
                         "depenvtid", "siteinfo","analysisunitname")]
  colnames(site.geochrons) <- c("siteid", "collectionunitid", "sampleid", "datasetid", 
                                "analysisunitid", "geochronid", "geochrontypeid", "agetypeid", 
                                "age", "errorolder", "erroryounger", "infinite", 
                                "delta13c", "labnumber", "materialdated", "notes", 
                                "depenvtid", "sitename","analysisunitname")
  print(paste(nrow(site.geochrons), "dates found from", 
              length(unique(site.geochrons$collectionunitid)), "collections", collapse = " "))
  #filter by taxa and age 
  ##filter site.geochrons to only <30ka
  colls.agelimit <- subset(site.geochrons, site.geochrons$age <= agelimit | 
                         (site.geochrons$age - site.geochrons$erroryounger) <= agelimit)$collectionunitid
  result <- subset(site.geochrons, site.geochrons$collectionunitid %in% colls.agelimit) 
  if (nrow(result) == 0) {
    print("Collection minimum age > specified age limit, returning NA"); return(NA)  }
  ##get list of small mammal taxa
  smalltaxa <- try(read.csv(here("data/smallmammaltaxa.csv")))
  if (!exists("smalltaxa")) {
    source(here("scripts/filter_taxa.R"))
  }
  smallmamm.colls <- get_smallvert_datasets(vert_and_chron, site.geochrons,datasets)
  result <- subset(result, result$collectionunitid %in% smallmamm.colls)
  result <- distinct(result)
  if (nrow(result) == 0) {
    print("No small mammals in collection, returning NA"); return(NA)  }
  result$source <- "literature search"
  
  #tag rejected or human remains
  result$rejected <- FALSE
  result[grepl("[Dd]ate rejected",result$notes) | 
           grepl("[Dd]ate too",result$notes) | 
           grepl("[Cc]ontamination probable",result$notes),"rejected"] <- TRUE
  result$human <- FALSE
  result[grepl("Teeth, Homo",result$materialdated) | 
           grepl("Bone, Homo",result$materialdated) | 
           grepl("human bone",result$notes) |
           grepl("Homo sapiens bone",result$materialdated),"human"] <- TRUE
  
  #assign category by collection
  result.clean <- cleandates(result)
  if (nrow(result.clean) == 0) {
    print("No usable dates; returning NA"); return(NA)  }
  collids.cat <- data.frame(collid = unique(result.clean$collectionunitid))
  collids.cat$category <- NA
  for (i in 1:length(collids.cat$collid)) {
    dd <- which(result.clean$collectionunitid == collids.cat$collid[i])
    if (any(result.clean[dd,]$agetypeid != "Radiocarbon years BP")) {
      if (any(result.clean[dd,]$agetypeid == "Radiocarbon years BP")) {
        collids.cat[i,"category"] <- 3 } else { collids.cat[i,"category"] <- 4 }
    } else {
      if (length(unique(result.clean[dd,]$analysisunitid)) == 1) {
        collids.cat[i,"category"] <- 1 } else { collids.cat[i,"category"] <- 2 }
    }
  }
  result$category <- collids.cat$category[match(result$collectionunitid,collids.cat$collid)]
  return(result)
}

get_from_tilia <- function(values, params, meth) {
  require(jsonlite);require(httr)
  paramstring <- paste0("&_",params,"=",values,collapse = "")
  url <- paste0(c("https://tilia.neotomadb.org/api/?method=ti.", meth, paramstring), collapse = "")
  obj <- try(GET(url), silent = TRUE)
  return(fromJSON(content(obj,as = "text")))
}