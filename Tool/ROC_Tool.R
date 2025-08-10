# This file is a part of the https://github.com/BruhZul/roc-and-iso-roc-tool repository.
#
# The tool is free software and can be redistributed and/or modified under the 
# GNU General Public License v3.0.
#
# This tool is distributed but without any warranty or implied warranty of
# merchantability or fitness for a particular purpose. See the LICENSE file for
# more details, or refer to gnu.org/licenses for more details.


#install packages
packages_list = c("shiny","ggplot2","shinyjs","RColorBrewer","xml2")
not_installed_packages = packages_list[!(packages_list %in% installed.packages()[,"Package"])]
if(length(not_installed_packages)){
  cat("The following packages has been installed:", not_installed_packages,"\n")
  install.packages(not_installed_packages)
}

library(shiny)
library(ggplot2)
library(iRRA)
library(shinyjs)
library(RColorBrewer)
library(xml2)

# Load needed functions.
source("./utils/build_iso_cost.R")
source("./utils/build_iso_curves.R")
source("./utils/compute_AUC_iso_curves.R")
source("./utils/compute_RRA_iso_curves.R")

# Function that orders ROC points (needed to work with ggplot2)
sort_roc_data <- function(roc_data) {
  roc_data <- roc_data[order(roc_data$FPR, roc_data$TPR), ]
  return(roc_data)
}

# generate ROC curve from csv file.
generate_roc_curve <- function(data) {
  if (!all(sapply(data, is.numeric))) {
   stop("")
  }
  
  if (!("FPR" %in% names(data)) || !("TPR" %in% names(data))) {
    stop("CSV file must contain the columns 'FPR' e 'TPR'.")
  }
  
  points <- list(x_values = data$FPR, y_values = data$TPR)
  return(points)
}

# add iso_curves
add_iso_curves <- function(rho, metrics, metric_values, precision = 0.01, curve_type = "AUC") {
  if (is.null(metrics)) return(data.frame())
  if (!is.vector(metrics)) metrics <- c(metrics)
  
  iso_df <- data.frame()
  group_id <- 1
  
  for (metric in metrics) {
    iso_curves <- if (curve_type == "AUC") {
      build_iso_curves(metric, metric_values, rho, precision)
    } else {
      build_iso_curves(metric, metric_values, rho, precision)
    }
    
    for (points in iso_curves) {
      if (!is.null(points) && length(points$x) > 1 && length(points$y) > 1) {
        iso_df <- rbind(iso_df, data.frame(FPR = points$x, TPR = points$y, group = group_id, CurveType = curve_type, Metric = metric))
        group_id <- group_id + 1
      }
    }
  }
  
  if (nrow(iso_df) == 0) {
    warning("No curve has been generated.")
  }
  
  return(iso_df)
}


# iso-cost curves
add_iso_curves_cost <- function(rho, cost_range, cFN, cFP) {
  if (is.null(cost_range)) return(data.frame()) 
  
  iso_df <- data.frame()
  group_id <- 1
  
  iso_curves <- build_iso_norm_cost_curve(rho, cost_range, cFN, cFP) 
  
  for (points in iso_curves) {
    if (!is.null(points) && length(points$x) > 1 && length(points$y) > 1) {
      valid_indices <- which(points$x >= 0 & points$x <= 1 & points$y >= 0 & points$y <= 1) 
      if (length(valid_indices) > 1) {  
        iso_df <- rbind(iso_df,data.frame(FPR = points$x[valid_indices], TPR = points$y[valid_indices], group = group_id, CurveType = "Cost"))
        group_id <- group_id + 1
      }
    }
  }
  
  if (nrow(iso_df) == 0) {
    warning("No iso-cost curve has been generated.")
  }
  
  return(iso_df)
}


# compute AUC
calculate_auc <- function(x_values, y_values) {
  if (x_values[1] == 1) {
    x_values <- rev(x_values)
    y_values <- rev(y_values)
  }
  
  auc <- sum(diff(x_values) * (head(y_values, -1) + tail(y_values, -1)) / 2)
  
  return(auc)
}

# compute RoI
calculate_rra_region <- function(roc_data, rho) {
  rra_compuation = suppressWarnings(rra(roc_data$FPR, roc_data$TPR, 1, (1-rho)/rho, 
                       recall = TRUE, fallout = TRUE, print = FALSE, plot = FALSE))
  
  combined_fpr <- unique(roc_data$FPR)
  combined_fpr <- sort(combined_fpr)
  upper_tpr <- sapply(combined_fpr, function(fpr) {
    tpr_values <- roc_data$TPR[roc_data$FPR <= fpr]
    max(tpr_values, na.rm = TRUE)
  })
  valid_indices <- combined_fpr <= rho & upper_tpr >= rho
  rra_region <- data.frame(FPR = unlist(rra_compuation$under_x),
                           TPR = unlist(rra_compuation$under_y))
  
  if (any(valid_indices)) {
    first_point <- data.frame(
      FPR = combined_fpr[which(valid_indices)[1]],
      TPR = upper_tpr[which(valid_indices)[1]]
    )
    last_point <- data.frame(
      FPR = combined_fpr[which(valid_indices)[length(which(valid_indices))]],
      TPR = upper_tpr[which(valid_indices)[length(which(valid_indices))]]
    )
  } else {
    first_point <- data.frame(FPR = NA, TPR = NA)
    last_point <- data.frame(FPR = NA, TPR = NA)
  }
  
  # save list
  list(
    rra_region = rra_region,
    first_point = first_point,
    last_point = last_point
  )
}


# UI
ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$script(src = "clickHandler.js")),
  titlePanel("Calculate ROC Curve"),
  
  fluidRow(
    # Left column for menu and inputs
    column(
      width = 4,
      
      actionButton("menu_button", "Menu", class = "btn-primary"),
      actionButton("upload_config_button", "Open Project", class = "btn-info"),
      actionButton("reset_button", "Reset", class = "btn-danger"),
      actionButton("help_button", "Help", class = "btn-secondary"),
      
      
      hidden(
        div(
          id = "menu_fields",
          fileInput("file_upload", "Upload CSV files", multiple = TRUE, accept = ".csv"),
          
          actionButton("rename_button", "Rename ROC Curves", class = "btn-secondary"),
          br(),
          hidden(
            div(
              id = "rename_fields",
              uiOutput("rename_inputs"), 
              actionButton("apply_rename_button", "Save Changes", class = "btn-success")
            )
          ),
          br(),
          
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            numericInput("AP", "AP:", value = 1, width = "100px"),
            numericInput("AN", "AN:", value = 1, width = "100px")
          ),
          
          br(),
          
          actionButton("isoCurves_button", "Iso-Curves", class = "btn-secondary", width = "100px"),
          br(),
          hidden(
            div(
              id = "isoCurves_fields",
              uiOutput("isoCurves_inputs"), 
              selectInput("iso_metric", "Iso-Curve Metric:",
                          choices = c("TPR", "FPR", "MCC", "BA", "Gmean", "GM", "D2H", "precision", 
                                      "F-score", "NPV", "TNR", "NM", "Markedness"), 
                          selected = "", multiple = TRUE),
              numericInput("iso_interval", "Interval Iso-Curve:", value = 0.1, min = 0.01, step = 0.01),
              
              selectInput("iso_metric_auc", "Metric for Iso-Curve based on AUC:",
                          choices = c("TPR", "FPR", "MCC", "BA", "Gmean", "GM", "D2H", "precision", 
                                      "F-score", "NPV", "TNR", "NM", "Markedness"), 
                          selected = "", multiple = TRUE),
              selectInput("iso_metric_rra", "Metric for Iso-Curve based on RRA:",
                          choices = c("TPR", "FPR", "MCC", "BA", "Gmean", "GM", "D2H", "precision", 
                                      "F-score", "NPV", "TNR", "NM", "Markedness"), 
                          selected = "", multiple = TRUE)
            )
          ),
          br(),
          
          actionButton("cost_button", "Cost", class = "btn-secondary", width = "100px"),
          br(),
          hidden(
            div(
              id = "cost_fields",
              uiOutput("cost_inputs"), 
              numericInput("cFN", "FN (cFN):", value = 1, min = 0, step = 0.1),
              numericInput("cFP", "FP (cFP):", value = 1, min = 0, step = 0.1),
              numericInput("iso_interval_cost", "Interval Iso-Curves Cost:", value = 0.1, min = 0.01, step = 0.01)
            )
          ),
          br(),
          
          numericInput("plot_size_input", "Plot Size (px):", value = 700, min = 500, step = 50, width = "100px"),
          br(),
          
          actionButton("show_hide_button", "Show/Hide Lines", class = "btn-secondary"),
          br(),
          hidden(
            div(
              id = "show_hide_fields",
              checkboxInput("show_diagonal", "Show Diagonal (FPR = TPR)", FALSE),
              checkboxInput("show_roi", "Show RoI", TRUE),
              #checkboxInput("hide_iso_curves", "Show Iso-Curve", TRUE),
              checkboxInput("hide_iso_curve_cost", "Show Iso-Curve (Cost)", FALSE),
              checkboxInput("show_threshold_points", "Show RoI Threshold Points", TRUE),
              uiOutput("curve_visibility_controls"),
            )
          ),
          
          #uiOutput("curve_visibility_controls"),
          br(),
          
          actionButton("submit_button", "Submit", class = "btn-success"),
          downloadButton("save_config_button", "Save Project", class = "btn-warning")
          
        )
      ),
    ),
    
    # Right column. Plot, legend, output.
    column(
      width = 8,
      
      hidden(
        div(
          id = "upload_fields",
          fileInput("file_upload_config", "Load Project", accept = ".rds"),
        )
      ),
      
      mainPanel(
        plotOutput("roc_plot", width = "auto", height = "auto"),
        verbatimTextOutput("auc_output"),
        #actionButton("save_partial_xhtml_button", "Save Partial XML Plot", class = "btn-secondary"),
        #downloadButton("download_partial_xml", "Download Partial XML", class = "btn-primary"),
        hidden(downloadButton("download_plot_button", "Export Plot and Data", class = "btn-warning"))
      )
    )
  )
)




# Server
server <- function(input, output, session) {
  roc_data_list <- reactiveVal(list())    
  auc_rra_values <- reactiveVal(list()) 
  file_names_list <- reactiveVal(list())
  iso_values = reactiveVal(list())
  last_processed_file = reactiveVal(NULL)
  pressed_reset = reactiveVal(FALSE)
  
  observeEvent(input$menu_button, {
    shinyjs::toggle("menu_fields") 
  })
  
  observeEvent(input$upload_config_button, {
    shinyjs::toggle("upload_fields")
  })
  
  observeEvent(input$show_hide_button, {
    shinyjs::toggle("show_hide_fields")
  })
  
  observeEvent(input$rename_button, {
    shinyjs::toggle("rename_fields")
  })
  
  observeEvent(input$isoCurves_button, {
    shinyjs::toggle("isoCurves_fields")
  })
  
  observeEvent(input$cost_button, {
    shinyjs::toggle("cost_fields")
  })
  
  # For renaming curves
  output$rename_inputs <- renderUI({
    file_names <- file_names_list()
    if (length(file_names) == 0) {
      return(tags$p("No curves to rename."))
    }
    
    lapply(seq_along(file_names), function(i) {
      textInput(paste0("curve_name_", i), label = paste("Rename Curve", i), value = file_names[[i]])
    })
  })
  # For renaming curves
  observeEvent(input$apply_rename_button, {
    file_names <- file_names_list()
    if (length(file_names) == 0) return(NULL)
    
    new_names <- sapply(seq_along(file_names), function(i) {
      input[[paste0("curve_name_", i)]]
    })
  
    file_names_list(new_names)
    
    output$roc_plot <- renderPlot({
      generate_roc_plot()
    }, width = function() { input$plot_size_input }, height = function() { input$plot_size_input })
    
    showNotification("Curve names updated successfully.", type = "message")
  })
  
  
  # save and download configuration
  output$save_config_button <- downloadHandler(
    filename = function() {
      paste("project_", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) { 
      # Use reactive variables
      roc_data_list_val <- roc_data_list()
      auc_rra_values_val <- auc_rra_values()
      file_for_roc = file_names_list()
      
      # Create configuartion file
      config <- list(
        AP = input$AP,
        AN = input$AN,
        show_roi = input$show_roi,
        iso_metric = input$iso_metric,
        iso_interval = input$iso_interval,
        iso_metric_auc = input$iso_metric_auc,
        iso_metric_rra = input$iso_metric_rra,
        cFN = input$cFN,
        cFP = input$cFP,
        iso_interval_cost = input$iso_interval_cost,
        hide_iso_curves = input$hide_iso_curves,
        hide_iso_curve_cost = input$hide_iso_curve_cost,
        show_threshold_points = input$show_threshold_points,
        plot_size_input = input$plot_size_input,
        show_diagonal = input$show_diagonal,
        roc_data_list = roc_data_list_val,
        auc_rra_values = auc_rra_values_val,
        file_names = file_for_roc
      )
      
      # Save file
      saveRDS(config, file)
      showNotification("Configuration saved successfully.", type = "message")
    }
  )
  
  
  # Load configuration
  observeEvent(input$file_upload_config, {
    req(input$file_upload_config)
    
    file <- input$file_upload_config$datapath
    
    if (file != "") {
      config <- readRDS(file) 
      
      updateNumericInput(session, "AP", value = config$AP)
      updateNumericInput(session, "AN", value = config$AN)
      updateCheckboxInput(session, "show_roi", value = config$show_roi)
      updateSelectInput(session, "iso_metric", selected = config$iso_metric)
      updateNumericInput(session, "iso_interval", value = config$iso_interval)
      updateSelectInput(session, "iso_metric_auc", selected = config$iso_metric_auc)
      updateSelectInput(session, "iso_metric_rra", selected = config$iso_metric_rra)
      updateNumericInput(session, "cFN", value = config$cFN)
      updateNumericInput(session, "cFP", value = config$cFP)
      updateNumericInput(session, "iso_interval_cost", value = config$iso_interval_cost)
      updateCheckboxInput(session, "hide_iso_curves", value = config$hide_iso_curves)
      updateCheckboxInput(session, "show_threshold_points", value = config$show_threshold_points)
      updateCheckboxInput(session, "hide_iso_curve_cost", value = config$hide_iso_curve_cost)
      updateNumericInput(session, "plot_size_input", value = config$plot_size_input)
      updateCheckboxInput(session, "show_diagonal", value = config$show_diagonal)
      
      if (!is.null(config$roc_data_list)) {
        roc_data_list_l <- config$roc_data_list
        auc_rra_list <- list()
        roc_files = config$file_names
        
        for (i in seq_along(roc_data_list_l)) {
          roc_data <- roc_data_list_l[[i]]
          
          roc_data <- roc_data[order(roc_data$FPR, roc_data$TPR), ]   
          
          auc_value <- calculate_auc(roc_data$FPR, roc_data$TPR)
          
          rho <- config$AP / (config$AP + config$AN)
          rra_value <- suppressWarnings(rra(roc_data$FPR, roc_data$TPR, config$AP, config$AN, 
                                            fallout = TRUE, recall = TRUE, plot = FALSE, print = FALSE))
          
          auc_rra_list[[i]] <- list(auc = auc_value, rra = rra_value$rra) 
          
          roc_data_list_l[[i]] <- roc_data  
        }
        
        roc_data_list(roc_data_list_l)
        auc_rra_values(auc_rra_list)
        file_names_list(roc_files)
        
        num_curves <- length(roc_data_list_l)
        showNotification(paste("Load", num_curves, "ROC Curve."), type = "message")
        
        shinyjs::click("submit_button")
      }
    }
  })
  
  
  observeEvent(input$submit_button, {
    req(input$file_upload) 
    show("roc_plot")
    show("download_plot_button")
    
    file_paths <- input$file_upload$datapath
    
    if(!identical(file_paths,last_processed_file())){
      last_processed_file(file_paths)
      AP <- input$AP
      AN <- input$AN
      N <- AP + AN
      rho <- AP / N
      
      file_names <- input$file_upload$name
      for (i in seq_along(file_paths)) {
        file <- file_paths[i]
        file_name <-  substr(file_names[i], 1, nchar(file_names[i])-4)
        
        roc_data <- read.csv(file)
        
        # needed, otherwise ggplot does not plot the ROC curve correctly
        roc_data <- roc_data[order(roc_data$FPR, roc_data$TPR), ] 
        
        if (!"Thresholds" %in% colnames(roc_data)) {
          roc_data$Thresholds <- rep("Not Provided", nrow(roc_data))
        }
        
        contains_names = "Name" %in% colnames(roc_data)
        
        if(contains_names){
          curve_names = unique(roc_data$Name)
        }else{
          curve_names = file_name
        }
        
        for(name in curve_names){
          if(contains_names){
            roc_name = roc_data[roc_data$Name == name,]
          }else{
            roc_name = roc_data
          }
          
          roc_points <- data.frame(FPR = roc_name$FPR, TPR = roc_name$TPR, Thresholds = roc_name$Thresholds)
          
          auc_value <- calculate_auc(roc_points$FPR, roc_points$TPR)  
          
          rra_region <- calculate_rra_region(roc_points, rho) 
          
          auc_rra <- suppressWarnings(rra(roc_points$FPR, roc_points$TPR, AP, AN, 
                                          fallout = TRUE, recall = TRUE, plot = FALSE, print = FALSE))
          
          roc_data_list(append(roc_data_list(), list(roc_points)))
          auc_rra_values(append(auc_rra_values(), list(list(auc = auc_value, rra = auc_rra$rra))))
          file_names_list(append(file_names_list(), name))
        }
      }
    }
  })
  
  
  
  observeEvent(input$save_partial_xhtml_button, {
    save_dir <- tempdir() 
    xhtml_filename <- "roc_data_partial.xml"
    xhtml_path <- file.path(save_dir, xhtml_filename)
    
    tryCatch({
      auc_rra_list <- auc_rra_values()
      file_names <- file_names_list()
      new_xhtml_content <- generate_xhtml_content(auc_rra_list, file_names)
      
      if (!is.null(new_xhtml_content)) {
        writeLines(new_xhtml_content, xhtml_path)
        
        file_state$xhtml_content <- new_xhtml_content 
        file_state$xhtml_path <- xhtml_path
        
        showNotification("Partial XML saved successfully. File available for download.", type = "message")
      } else {
        showNotification("Error generating XML content.", type = "error")
      }
    }, error = function(e) {
      showNotification("Error saving partial XML.", type = "error")
    })
  })
  
  output$download_partial_xml <- downloadHandler(
    filename = function() {
      paste0("partial_roc_", Sys.Date(), ".xml")
    },
    content = function(file) {
      if (!is.null(file_state$xhtml_path) && file.exists(file_state$xhtml_path)) {
        file.copy(file_state$xhtml_path, file)
      } else {
        stop("Partial XML file not available.")
      }
    },
    contentType = "application/xml"
  )
  
  
  #salva e scarica plot
  file_state <- reactiveValues(
    xhtml_path = NULL,
    xhtml_content = NULL
  )
  
  output$download_plot_button <- downloadHandler(
    filename = function() {
      paste("roc_plot_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      temp_dir <- tempdir()
      png_filename <- "roc_plot.png"
      xhtml_filename <- "roc_data.xhtml"
      zip_filename <- "roc_plot.zip"
      
      tryCatch({
        plot <- generate_roc_plot()
        if (is.null(plot)) stop("Plot is not available.")
        ggsave(filename = png_filename, plot = plot, device = "png", width = 8, height = 6)
      }, error = function(e) {
        showNotification(paste("Error saving PNG file:", e$message), type = "error")
        return(NULL)
      })
      
      png_base64 <- tryCatch({
        png_path = png_filename
        png_data <- readBin(png_path, what = "raw", n = file.info(png_path)$size)
        base64enc::base64encode(png_data)
      }, error = function(e) {
        showNotification("Error encoding PNG file.", type = "error")
        return(NULL)
      })
      
      auc_rra_list <- auc_rra_values()
      file_names <- file_names_list()
      new_xhtml_content <- generate_xhtml_content(auc_rra_list, file_names)
      
      if (!is.null(new_xhtml_content)) {
        xhtml_path = xhtml_filename
        writeLines(new_xhtml_content, xhtml_path)
        file_state$xhtml_path <- xhtml_path
      }
      
      tryCatch({
        zip(zipfile = zip_filename, files = c(png_filename, file_state$xhtml_path))
        file.copy(zip_filename, file, overwrite = TRUE)
      }, error = function(e) {
        showNotification("Error creating ZIP file.", type = "error")
        return(NULL)
      })
      
      showNotification("Download completed successfully.", type = "message")
    },
    contentType = "application/zip"
  )
  
  generate_xhtml_content <- function(auc_rra_list, file_names) {
    auc_metric <- input$iso_metric_auc
    rra_metric <- input$iso_metric_rra
    rho <- input$AP / (input$AP + input$AN)
    
    roc_curves = roc_data_list()
    
    ctr = 1
    
    tryCatch({
      content <- '<?xml version="1.0" encoding="UTF-8"?>\n'
      content <- paste0(content, '<!DOCTYPE html>\n<html xmlns="http://www.w3.org/1999/xhtml">\n')
      content <- paste0(content, "<head>\n<title>ROC Data</title>\n")
      content <- paste0(content, "<style>\nimg { display: none; max-width: 800px; max-height: 600px; }\n</style>\n")
      content <- paste0(content, "<script>\nfunction toggleImage(id) {\n")
      content <- paste0(content, "  var img = document.getElementById(id);\n")
      content <- paste0(content, "  img.style.display = (img.style.display === 'none') ? 'block' : 'none';\n")
      content <- paste0(content, "}\n</script>\n</head>\n<body>\n")
      
      for (i in seq_along(auc_rra_list)) {
        curve_id <- paste0("curve", i)
        content <- paste0(content, "<div>\n")
        content <- paste0(content, "<p>Curve Name: ", file_names[[i]], "</p>\n")
        content <- paste0(content, "<p>AUC: ", round(auc_rra_list[[i]]$auc, 3), "</p>\n")
        content <- paste0(content, "<p>RRA: ", round(auc_rra_list[[i]]$rra, 3), "</p>\n")
        
        # duplicated code (bandaid fix for now)
        curve <- roc_curves[[ctr]]
        rra_result <- calculate_rra_region(curve, rho)
        first_point <- rra_result$first_point
        last_point <- rra_result$last_point

        if(curve$Thresholds[1] == "Not Provided" || is.na(first_point$FPR)){
          min_thresh <- "Not Provided"
          max_thresh <- "Not Provided"
        }else{
          min_thresh <- curve$Thresholds[which.min(abs(curve$FPR - first_point$FPR))]
          max_thresh <- curve$Thresholds[which.min(abs(curve$FPR - last_point$FPR))]
        }
        content <- paste0(content, "<p>First point in RoI: FPR = ", 
                          round(first_point$FPR, 3), ", TPR = ", round(first_point$TPR, 3),
                          ifelse(min_thresh == "Not Provided", "", 
                                 paste0(", Threshold = ", round(min_thresh, 3))), "</p>\n")
        content <- paste0(content, "<p>Last point in RoI: FPR = ", 
                          round(last_point$FPR, 3), ", TPR = ", round(last_point$TPR, 3),
                          ifelse(max_thresh == "Not Provided", "", 
                                 paste0(", Threshold = ", round(max_thresh, 3))), "</p>\n")
        ctr = ctr + 1
        for(metr in auc_metric){
          iso_auc_value <- compute_iso_value_from_auc(metr, auc_rra_list[[i]]$auc, rho)
          content <- paste0(content,"<p>Iso-",metr," AUC:", round(iso_auc_value, 3), "</p>\n")
        }
        
        for(metr in rra_metric){
          iso_rra_value <- compute_iso_value_from_rra(metr, auc_rra_list[[i]]$rra, rho)
          content <- paste0(content,"<p>Iso-", metr," RRA:", round(iso_rra_value, 3), "</p>\n") 
        }
        content <- paste0(content, "</div>\n")
      }
      content <- paste0(content, "<div>\n")
      content <- paste0(content, '<p><a href="#" onclick="toggleImage(\'', "plot", '\')">Show/Hide Plot</a></p>\n')
      content  = paste0(content, "<img id='plot' src='./roc_plot.png'></img>")
      content <- paste0(content, "</div>\n")
      
      content <- paste0(content, "</body>\n</html>")
      content
    }, error = function(e) {
      print(e)
      NULL
    })
  }
  
  # dynamic checkboxes
  output$curve_visibility_controls <- renderUI({
    req(file_names_list())
    file_names <- file_names_list()
    
    if (length(file_names) == 0) {
      return(tags$p("No curves to display."))
    }
    
    controls <- list()
    #qui
    iso_metric = input$iso_metric
    for(iso in iso_metric){
      controls[[length(controls) + 1]] <- checkboxInput(paste0("isocurve_",iso), 
                                                        label = paste0("Show Iso-",iso," Curves"), 
                                                        value = TRUE) 
    }
    auc_metric <- input$iso_metric_auc
    rra_metric <- input$iso_metric_rra
    
    for (i in seq_along(file_names)) {
      curve_id <- paste0("curve_", i)
      controls[[length(controls) + 1]] <- checkboxInput(curve_id, 
                                                        label = paste("Show", file_names[[i]]), 
                                                        value = TRUE) 
      
      for(metr in auc_metric){
        controls[[length(controls) + 1]] <- checkboxInput(paste0(curve_id,"_",metr, "_auc"), 
                                                          label = paste0("Show Iso-",metr," Curves based on AUC for", file_names[[i]]), 
                                                          value = TRUE)
      }
      
      for(metr in rra_metric){
        controls[[length(controls) + 1]] <- checkboxInput(paste0(curve_id,"_",metr, "_rra"), 
                                                          label = paste0("Show Iso-",metr," Curves based on RRA for", file_names[[i]]), 
                                                          value = TRUE) 
      }
    }
    do.call(tagList, controls)
  })
  outputOptions(output, "curve_visibility_controls", suspendWhenHidden = FALSE)
  
  
  observeEvent(input$reset_button, {
    is_reset_pressed = pressed_reset()
    if(!is_reset_pressed){
      showNotification("To reset the tool, press Reset again", type = "error")
    }else{
      roc_data_list(list())
      auc_rra_values(list())
      file_names_list(list())
      
      updateNumericInput(session, "AP", value = 1)
      updateNumericInput(session, "AN", value = 1)
      updateCheckboxInput(session, "show_roi", value = TRUE)
      updateCheckboxInput(session, "hide_iso_curves", value = FALSE)
      updateNumericInput(session, "plot_size_input", value = 700)
      updateCheckboxInput(session, "show_diagonal", value = FALSE)
      updateSelectInput(session, "iso_metric", selected = NULL)
      
      hide("roc_plot")
      hide("download_plot_button")
      
      showNotification("All data and settings have been reset.", type = "message")
    }
    
    is_reset_pressed = !is_reset_pressed
    pressed_reset(is_reset_pressed)
  })
  

  observeEvent(input$help_button, {
    system(paste0('open "', "./doc/Documentation.pdf", '"'))
    showModal(modalDialog(
      title = "Help",
      pre(includeText("help.txt")),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  #plot
  generate_roc_plot <- function() {
    req(roc_data_list())
    
    if (length(roc_data_list()) == 0) {
      return(ggplot() +
               labs(title = "ROC Curve",
                    x = "FPR",
                    y = "TPR") +
               theme_minimal())
    }
    
    plot_size <- input$plot_size_input
    rho <- input$AP / (input$AP + input$AN)
    
    num_curves <- length(roc_data_list())
    colors <- if (num_curves < 3) {
      c("red", "blue", "green")[1:num_curves]
    } else {
      brewer.pal(min(num_curves, 9), "Set1")
    }
    
    iso_colors <- suppressWarnings(brewer.pal(min(length(input$iso_metric), 9), "Dark2")) 
    
    iso_legend_labels <- lapply(input$iso_metric, function(metric) paste0("Iso-", metric, " Curve"))
    names(iso_legend_labels) <- unlist(iso_legend_labels)
    
    p <- ggplot() +
      labs(title = "ROC Curve",
           x = "FPR",
           y = "TPR") +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      coord_fixed(ratio = 1)
    
    roc_data <- roc_data_list()
    auc_rra <- auc_rra_values()
    file_names <- file_names_list()
    shown_file_names = c()
    shown_colors = c()
    
    if (length(file_names) != length(roc_data)) {
      file_names <- paste("Curve", seq_along(roc_data))
    }
    
    linetype_labels <- list(
      "TPR Line" = "TPR Line",
      "FPR Line" = "FPR Line",
      "Diagonal" = "Diagonal",
      "Iso-Curves" = paste("Iso-Curves -", paste(input$iso_metric, collapse = ", ")),
      "Cost" = "dashed" 
    )
    
    
    iso_curves_data <- list()
    
    for (i in seq_along(roc_data)) {
      roc_points <- roc_data[[i]][order(roc_data[[i]]$FPR, roc_data[[i]]$TPR), ]
      roc_points$curve_id <- as.factor(i)
      
      curve_id <- paste0("curve_", i)
      show_curve <- input[[curve_id]]
      show_curve = ifelse(is.null(show_curve),FALSE,show_curve)
      auc_metric <- input$iso_metric_auc
      rra_metric <- input$iso_metric_rra
      
      curve_color <- colors[i]
      rra_region <- calculate_rra_region(roc_points, rho)
      
      if (input$show_roi) {
        p <- p + geom_hline(aes(yintercept = rho, linetype = "TPR RoI Line"), color = "pink", size = 1)
        p <- p + geom_vline(aes(xintercept = rho, linetype = "FPR RoI Line"), color = "pink", size = 1)
      }
      if (input$show_diagonal) {
        p <- p + geom_abline(aes(intercept=0, slope=1, linetype = "Diagonal"), color = "lightgrey", size = 1, alpha=0.6)
      }
      
      #iso-curve
      if (length(roc_data) > 0) { 
        metrics <- input$iso_metric
        metric_values <- seq(0, 1, by = input$iso_interval)
        for (metric_idx in seq_along(metrics)) {
          iso_metric <- metrics[metric_idx]
          if(ifelse(is.null(input[[paste0("isocurve_",iso_metric)]]), FALSE, input[[paste0("isocurve_",iso_metric)]])){
            iso_curves <- add_iso_curves(rho, iso_metric, metric_values)
            iso_color <- iso_colors[metric_idx]
            iso_label <- paste0("Iso-", iso_metric, " Curve")
            
            if (nrow(iso_curves) > 0) {
              iso_curves$legend_label <- iso_label
              p <- p + geom_path(data = iso_curves, aes(x = FPR, y = TPR, group = group, linetype = legend_label), color = iso_color, size = 0.8)
            }
          }
        }
      }
      
      if (input$hide_iso_curve_cost) {
        iso_cost_data <- add_iso_curves_cost(rho, seq(0, 1, by = input$iso_interval_cost), cFN = input$cFN, cFP = input$cFP)  
        if (nrow(iso_cost_data) > 0) {
          p <- p + 
            geom_path(data = iso_cost_data, aes(x = FPR, y = TPR, group = group, linetype = CurveType),
                      color = "purple", size = 0.5) +
            scale_linetype_manual(values = c("Cost" = "dashed"))
        }
      }
      
      
      if (show_curve) {   
        shown_file_names = c(shown_file_names, file_names[i])
        shown_colors = c(shown_colors, colors[i])
        p <- p + geom_line(data = roc_points, aes(x = FPR, y = TPR, color = curve_id), size = 1)
        
        rra_result <- calculate_rra_region(roc_points, rho)
        rra_region <- rra_result$rra_region
        first_point <- rra_result$first_point
        last_point <- rra_result$last_point
        
        if(input$show_roi){
          if (!is.null(rra_region) && length(rra_region) > 0) {
            p <- p + geom_polygon(data = rra_region, aes(x = FPR, y = TPR, fill = "RoI Region"), alpha = 0.2)
          }
        }
        
        if(input$show_threshold_points){
          if (!is.na(first_point$FPR)) {
            p <- p + geom_point(data = first_point, aes(x = FPR, y = TPR), color = "black", size = 2.7)
          }
          if (!is.na(last_point$FPR)) {
            p <- p + geom_point(data = last_point, aes(x = FPR, y = TPR), color = "black", size = 3.5)
          }
        }
        
        for(metr in auc_metric){
          if(ifelse(is.null(input[[paste0(curve_id,"_",metr, "_auc")]]), FALSE, input[[paste0(curve_id,"_",metr, "_auc")]])){
            auc_label <- paste0("Iso-",metr," Curve AUC ",file_names[i])
            iso_curves_auc <- add_iso_curves(rho, metr, 
                                             compute_iso_value_from_auc(metr, auc_rra[[i]]$auc, rho), 
                                             precision = input$iso_interval, curve_type = "AUC")
            if (nrow(iso_curves_auc) > 0) {
              iso_curves_auc$legend_label <- auc_label
              p <- p + geom_line(data = iso_curves_auc, aes(x = FPR, y = TPR, linetype = legend_label), 
                                 color = curve_color, size = 1)
              linetype_labels[[auc_label]] <- auc_label
            } 
          }
        }
        
        for(metr in rra_metric){
          if(ifelse(is.null(input[[paste0(curve_id,"_",metr, "_rra")]]), FALSE, input[[paste0(curve_id,"_",metr, "_rra")]])){
            rra_label <- paste0("Iso-",metr," Curve RRA ",file_names[i])
            iso_curves_rra <- add_iso_curves(rho, metr, 
                                             compute_iso_value_from_rra(metr, auc_rra[[i]]$rra, rho), 
                                             precision = input$iso_interval, curve_type = "RRA")
            if (nrow(iso_curves_rra) > 0) {
              iso_curves_rra$legend_label <- rra_label
              p <- p + geom_line(data = iso_curves_rra, aes(x = FPR, y = TPR, linetype = legend_label), 
                                 color = curve_color, size = 1)
              linetype_labels[[rra_label]] <- rra_label
            }  
          }
        }
      }
    }
    
    
    p <- p + scale_color_manual(name = "ROC Curves", values = shown_colors, labels = shown_file_names) +
      scale_linetype_manual(name = "Lines", 
                            values = c(
                              "TPR RoI Line" = "dashed", 
                              "FPR RoI Line" = "dashed", 
                              "Diagonal" = "solid", 
                              
                              setNames(rep("dotted", length(iso_legend_labels)), names(iso_legend_labels)),
                              
                              setNames(rep("twodash", length(linetype_labels[-c(1:5)][grepl(" AUC ", linetype_labels[-c(1:5)], fixed = TRUE)])),
                                       names(linetype_labels[-c(1:5)][grepl(" AUC ", linetype_labels[-c(1:5)], fixed = TRUE)])),
                              
                              setNames(rep("dotdash", length(linetype_labels[-c(1:5)][grepl(" RRA ", linetype_labels[-c(1:5)], fixed = TRUE)])),
                                       names(linetype_labels[-c(1:5)][grepl(" RRA ", linetype_labels[-c(1:5)], fixed = TRUE)])),
                              
                              "Cost" = "dashed")) +
      scale_fill_manual(name = "Region", values = c("RoI Region" = "yellow")) +
      theme(legend.position = "right", legend.text = element_text(size = 10))
      #theme(legend.position = "right", legend.text = element_text(size = 10), legend.key.size = unit(1, 'cm'))
    
    return(p)
  }
  
  
  
  
  
  output$roc_plot <- renderPlot({
    if (length(roc_data_list()) == 0) {
      return(NULL)
    }
    
    generate_roc_plot()
  }, width = function() { input$plot_size_input }, height = function() { input$plot_size_input })
  
  
  
  # output AUC, RRA e Thresholds of ROI
  output$auc_output <- renderPrint({
    req(auc_rra_values(), roc_data_list())
    auc_rra_list <- auc_rra_values()
    file_names <- file_names_list()
    rho <- input$AP / (input$AP + input$AN)
    
    auc_metric <- input$iso_metric_auc
    rra_metric <- input$iso_metric_rra
    
    for (i in seq_along(auc_rra_list)) {
      curve_id <- paste0("curve_", i)
      show_curve <- input[[curve_id]]
      if(!(show_curve)){ 
        next
      }
      curve <- roc_data_list()[[i]]
      auc_rra <- auc_rra_values()[[i]]
      
      rra_result <- calculate_rra_region(curve, rho)
      first_point <- rra_result$first_point
      last_point <- rra_result$last_point
      
      if(curve$Thresholds[1] == "Not Provided" || is.na(first_point$FPR)){
        min_thresh <- "Not Provided"
        max_thresh <- "Not Provided"
      }else{
        min_thresh <- curve$Thresholds[which.min(abs(curve$FPR - first_point$FPR))]
        max_thresh <- curve$Thresholds[which.min(abs(curve$FPR - last_point$FPR))]
      }
      
      curve_name <- ifelse(length(file_names) >= i, file_names[[i]], paste("Curve", i))
      cat("ROC Curve -", curve_name, ":\n")
      cat("  AUC:", round(auc_rra_list[[i]]$auc, 3), "\n")
      cat("  RRA:", round(auc_rra_list[[i]]$rra, 3), "\n")
      
      for (metr in auc_metric) {
        if(ifelse(is.null(input[[paste0(curve_id,"_",metr, "_auc")]]), FALSE, input[[paste0(curve_id,"_",metr, "_auc")]])){
          iso_auc_value <- compute_iso_value_from_auc(metr, auc_rra_list[[i]]$auc, rho)
          cat("    Iso-", metr, " AUC:", round(iso_auc_value, 3), "\n", sep = "")
        }
      }
      
      for (metr in rra_metric) {
        if(ifelse(is.null(input[[paste0(curve_id,"_",metr, "_rra")]]), FALSE, input[[paste0(curve_id,"_",metr, "_rra")]])){
          iso_rra_value <- compute_iso_value_from_rra(metr, auc_rra_list[[i]]$rra, rho)
          cat("    Iso-", metr, " RRA:", round(iso_rra_value, 3), "\n", sep = "")
        }
      }
      
      if(input$show_threshold_points){
        cat("  First Point in ROI:\n")
        cat("    FPR:", round(first_point$FPR, 3), "TPR:", round(first_point$TPR, 3), "\n")
        if(min_thresh != "Not Provided"){
          cat("    Threshold:", ifelse(is.na(min_thresh), "NA", round(min_thresh, 3)), "\n")
        }
        
        cat("  Last Point in ROI:\n")
        cat("    FPR:", round(last_point$FPR, 3), "TPR:", round(last_point$TPR, 3), "\n")
        if(max_thresh != "Not Provided"){
          cat("    Threshold:", ifelse(is.na(max_thresh), "NA", round(max_thresh, 3)), "\n")
        }
      }
      cat("\n")
    }
  })
  
  
}


# Run app
shinyApp(ui = ui, server = server)


