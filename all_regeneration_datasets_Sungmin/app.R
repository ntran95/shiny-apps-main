library(shiny)
library(cowplot)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(shinyWidgets)
library(dplyr)


# =================================================================== Preamble


# ===================================================================== Server


# ========================================================================= UI


# =============================================================================
shinyApp(ui = ui, server = server)


# ==== Command line tools
# Upload app to shinyaqpps.io (must be directory that contains app.R)
# start R session

# Deploy to shinyapps.io
# if you have repo issues with bioconductor

# options(repos = BiocManager::repositories())
# getOption("repos")


# rsconnect::deployApp('/Volumes/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/all_regeneration_datasets_Sungmin', account = 'piotrowskilab')
# rsconnect::deployApp('/n/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/all_regeneration_datasets_Sungmin', account = 'piotrowskilab')

# Execute app locally
# options(shiny.reactlog=TRUE, shiny.fullstacktrace = TRUE); shiny::runApp('/Volumes/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/all_regeneration_datasets_Sungmin/app.R')

# Logs
# rsconnect::showLogs(account = 'piotrowskilab', appName = 'all_regeneration_datasets_Sungmin')

# available.packages(contriburl = "https://cran.rstudio.com")