#This script recursively filters the taxa table down to small mammals (Rodentia, Lagomorpha, Eulipotyphla), returns the result as table "smalltaxa", and saves a copy as "smallmammaltaxa.csv".

i_am("final_scripts/functions/filter_taxa.R")

#download entire taxa table and filter out non-mammals
if (!exists("taxontable")) {
  taxontable <- get_table("taxa",limit = 9999999)
  taxontable <- taxontable[taxontable$taxagroupid == "MAM",]
}

#target clades: Rodentia, Lagomorpha, Lipotyphla
#recursively find all taxa downstream from these three
toptaxaID <- taxontable[which(taxontable$taxonname == "Rodentia" |
                                taxontable$taxonname == "Lagomorpha" |
                                taxontable$taxonname == "Lipotyphla" |
                                taxontable$taxonname == "Dasypodidae" |
                                taxontable$taxonname == "Didelphidae" ),"taxonid"]
alltaxaID <- toptaxaID
for (i in 1:length(toptaxaID)) {
  a <- toptaxaID[i]
  b <- c(a,taxontable[taxontable$highertaxonid %in% toptaxaID[i],"taxonid"])
  while (length(a) < length(b)) {
    a <- b
    b <- unique(c(b,taxontable[taxontable$highertaxonid %in% a,"taxonid"])) 
  }
  alltaxaID <- unique(c(alltaxaID,b))
}
#store as "smalltaxa"
smalltaxa <- taxontable[taxontable$taxonid %in% alltaxaID,]

#save result and clear unnecessary objects
rm(taxontable,toptaxaID,alltaxaID,a,b)
write.csv(smalltaxa[,c(1:12)],file = here("data/smallmammaltaxa.csv"))
