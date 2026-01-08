# FIFE_RProject_public
Scripts and results for the Fire Injury and Fecundity (FIFE) project, led 
by PhD student Nina E. Venuti at the University of California, Davis, with 
mentorship from Dr. Andrew M. Latimer and Dr. Derek J. N. Young.//


This repository is organized into three main folders: **1) data**, 
**2) scripts**, and **3) results**.

**1) The data folder** is empty, and serves as a receptacle for the data
files stored in the FIFE project's data repository, published on Dryad
[permanent link to come]. Individuals interested in running the scripts
for the FIFE project should fork the GitHub repository ("FIFE_RProject_
public"), navigate to the Dryad repository [link to come], download the
.csv and .txt files contained therein, and store those files in the
the data folder within their local (forked) copy of the GitHub repository.
All references to the data directory within FIFE project scripts
are relative; scripts should run easily once data files are downloaded
from Dryad and deposited in the data folder.

**2) The scripts folder** contains a series of scripts, ordered from
01 through 07, that document the steps taken to clean and prepare the
project data for analysis (scripts 01 and 02), model the effects of
fire injury on tree fecundity in each of the project's focal taxa
(white fir, yellow pine, and sugar pine) (script 03, labeled according 
to taxon code), and summarize and present the results of these models
via model summary tables, and taxon- and sample year-specific results
figures and tables (scripts 04, 05, 06, and 07, labeled according to taxon
code, as appropriate). Running these scripts in order will allow interested 
parties to both trace decisions made during FIFE data analysis, and
reproduce final project results. It is recommended that interested parties
review script 01, but begin running code with script 02, which utilizes
the full, cleaned dataset available to the public via the Dryad repository. 
Running script 01 relies on access to the original, raw data files (which 
contain some errors, and are thus not made available to the public).

**3) The results folder** is divided into two sub-folders, containing
final model results figures and tables formatted for publication:

  - The **figures** folder contains a series of .png and .pdf 
    files displaying the effects of fire injury (and, in some
    cases, tree size) on seed cone production in each focal 
    taxon and sample year. Some of these figures are composite
    "panel plots," formatted for publication.
        
  - The **tables** folder contains a series of .csv files
    displaying key model summary metrics, as well as taxon-
    and sample year-specific information on seed cone production
    at particular values of fire injury and tree size (to aid in
    results interpretation).
    


