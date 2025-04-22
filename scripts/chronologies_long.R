#New workflow as of April 2024. Meant to run remotely on pinnacles cluster.
#Takes simple (single analysis unit) collections
#Returns age models for dates and collections
#Requires "dates_for_analysis_smallmamm.csv" and "functions.R" in same directory.
#Load data, functions, and libraries
require(here)
i_am("workflow/scripts/chronologies_long.R")
require(oxcAAR);quickSetupOxcal()
source(here("workflow/scripts/functions.R"))

#import dates list
#dates <- read.csv(here("workflow/data/dates_for_analysis_smallmamm.csv"))
dates <- read.csv(here("workflow/data/alldates.csv"))
#Get uncalibrated dates with one analysis unit
singleAU_uncal <- subset(dates,dates$category %in% c(1,2,7) & 
                           dates$infinite == FALSE &
                           dates$rejected == FALSE &
                           dates$age > 5 &
                           !is.na(dates$errorolder) & 
                           !is.na(dates$erroryounger) & 
                           !is.na(dates$age))

#Pick 5 random collection unit id labels
#colls <- unique(singleAU_uncal$collectionunitid)[runif(5,min = 0,max = length(unique(singleAU_uncal$collectionunitid)))]
#or a single collection
#colls <- 5603
#or the whole list (but check first)
colls <- unique(singleAU_uncal$collectionunitid)
#minus the short ones
colls <- subset(colls,sapply(colls,function(x) nrow(subset(singleAU_uncal,singleAU_uncal$collectionunitid == x))) >= 25)
#minus the ones we already calculated
done <- as.numeric(unlist(sapply(dir(here("workflow/outputs/sampleages")),
                  function(x) return(strsplit(x,".csv")))))
colls <- subset(colls, !colls %in% done)
#Return values for Neotoma  
#test <- sapply(colls, dateranges, data = singleAU_uncal, simplify = F)
#names(test) <- colls
#chronologies <- rlist::list.stack(sapply(test,function(x) x["chronology"]))
#write.csv(chronologies,here("workflow/outputs/new_sampleages_long.csv"),row.names = F)
for (i in 1:length(colls)) {
  result <- dateranges(colls[i],data = singleAU_uncal)
  write.csv(result$sampleages,here("workflow/outputs/sampleages/",paste0(colls[i],".csv",collapse = "")),row.names = F)
  write.csv(result$chroncontrols,here("workflow/outputs/chroncontrols/",paste0(colls[i],".csv",collapse = "")),row.names = F)
}
