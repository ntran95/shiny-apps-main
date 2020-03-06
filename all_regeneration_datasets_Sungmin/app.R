library(shiny)
library(cowplot)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(shinyWidgets)
library(dplyr)
library(devtools)

# Git branch
branch <- "split-app"

# =================================================================== Preamble
source_url(paste0("https://github.com/diazdc/shiny-apps-main/",
  "blob/", branch,"/all_regeneration_datasets_Sungmin/preamble.R"))

# ===================================================================== Server
source_url(paste0("https://github.com/diazdc/shiny-apps-main/",
  "blob/", branch,"/all_regeneration_datasets_Sungmin/server.R"))

# ========================================================================= UI
source_url(paste0("https://github.com/diazdc/shiny-apps-main/",
  "blob/", branch,"/all_regeneration_datasets_Sungmin/ui.R"))

# =============================================================================
shinyApp(ui = ui, server = server)


# ==== Deploy/execute tools
# options(repos = BiocManager::repositories())
# getOption("repos")

if (FALSE) {
  # Deploy local
  rsconnect::deployApp(paste0("/Volumes/projects/ddiaz/Analysis/",
    "Scripts/rsconnect/shinyapps.io/all_regeneration_datasets_Sungmin"),
    account = "piotrowskilab")
  # Server
  rsconnect::deployApp(paste0("/n/projects/ddiaz/Analysis/",
    "Scripts/rsconnect/shinyapps.io/all_regeneration_datasets_Sungmin"),
    account = "piotrowskilab")

  #Execute app locally
  options(shiny.reactlog=TRUE, shiny.fullstacktrace = TRUE)
  shiny::runApp(paste0("/Volumes/projects/ddiaz/Analysis/",
    "Scripts/rsconnect/shinyapps.io/all_regeneration_datasets_Sungmin/app.R"))

  # Logs
  rsconnect::showLogs(account = 'piotrowskilab',
    appName = 'all_regeneration_datasets_Sungmin')
}
# available.packages(contriburl = "https://cran.rstudio.com")