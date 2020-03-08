library(shiny)
library(cowplot)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(shinyWidgets)
library(dplyr)
library(rsconnect)

# Git branch
branch <- "split-app"
app_name <- "test_split_regen"

# =================================================================== Preamble
source(paste0("https://raw.githubusercontent.com/diazdc/shiny-apps-main/", branch, "/", app_name, "/preamble.R"))


# ===================================================================== Server
source(paste0("https://raw.githubusercontent.com/diazdc/shiny-apps-main/",
  branch, "/", app_name, "/server.R"))


# ========================================================================= UI
source(paste0("https://raw.githubusercontent.com/diazdc/shiny-apps-main/",
  branch, "/", app_name, "/ui.R"))


# =============================================================================
shinyApp(ui = ui, server = server)


# ========== Deploy/execute tools
# options(repos = BiocManager::repositories())
# getOption("repos")

if (FALSE) {

  # Deploy local
  rsconnect::deployApp(paste0("/Volumes/projects/ddiaz/Analysis/",
    "Scripts/rsconnect/shinyapps.io/", app_name),
    account = "piotrowskilab")
  
  # Server
  rsconnect::deployApp(paste0("/n/projects/ddiaz/Analysis/",
    "Scripts/rsconnect/shinyapps.io/", app_name),
    account = "piotrowskilab")

  #Execute app locally
  options(shiny.reactlog = TRUE, shiny.fullstacktrace = TRUE)
  shiny::runApp(paste0("/Volumes/projects/ddiaz/Analysis/",
    "Scripts/rsconnect/shinyapps.io/", app_name, "/app.R"))

  # Logs
  rsconnect::showLogs(account = 'piotrowskilab',
    appName = app_name)
}
# available.packages(contriburl = "https://cran.rstudio.com")