rcmethods <- data.frame(geochronid = alldates$geochronid,radiocarbonmethod = NA)
I <- 1
for (i in I:nrow(rcmethods)) {
  res <- try(get_from_tilia(values = rcmethods$geochronid[i], 
                            params = "geochronidlist", meth = "getradiocarbonbygeochronid"))
  if (length(res$data) > 0) rcmethods$radiocarbonmethod[i] <- res$data$radiocarbonmethod
}
alldates$radiocarbonmethod <- rcmethods[match(alldates$geochronid,rcmethods$geochronid),]$radiocarbonmethod
table(alldates$radiocarbonmethod)
write.csv(alldates,file = "alldates_method.csv",row.names = F)
