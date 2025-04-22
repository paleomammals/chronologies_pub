These are the publication versions of Val's scripts for this project. They are intended to be run on computers that are not set up exactly like Val's and eventually for archiving at publication. 

Scripts in this folder will run on any computer with the right R libraries installed, no need to install Postgres. They operate by calling the Neotoma 2.0 API using library "neotoma2", and where necessary the Neotoma Tilia API using library "RCurl".

Use "collections_loc_time_depenvt.Rmd" to produce the list of collections with age bins and depositional environment in "collectionpoints.csv". (The script "getsites.R" is currently broken due to the new API update, so don't use it or "sites_for_sim_API.Rmd" which depends on it.)

Of the files that it outputs, the ones you probably want are as follows:
* collectionpoints.csv (table for taphonomy model that includes all small mammal collections with _chronologies_ now)
* dates_for_analysis_smallmamm.csv (all geochronology dates <30k 14Cybp from collections that have small mammal taxa, from sites in North America)
    * dates_for_analysis_allvert.csv (same as previoius but without filtering to small mammals only)
* smallmammaltaxa.csv (version of Neotoma's table "taxa" filtered to only taxa in Rodentia, Lagomorpha, and Eulipotyphla, equivalent to the "ecologicalgroups" RODE, LAGO, and SORI in Neotoma's "samples" table)
