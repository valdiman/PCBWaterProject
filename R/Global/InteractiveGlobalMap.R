# Code to create interactive maps of PCB concentrations for all locations.
# Data on the plots are aggregated per week to help better visualizing the data

# Install packages
install.packages("dplyr")
install.packages("ggplot2")
install.packages("leaflet")
install.packages("shiny")
install.packages("tidyverse")

# Load libraries
{
  library(dplyr)
  library(ggplot2)
  library(leaflet)
  library(shiny)
  library(stringr) # str_detect
}

# Shiny app ---------------------------------------------------------------
# Define the Shiny UI
ui <- fluidPage(
  titlePanel("PCB Water Concentration Data Visualization"),
  leafletOutput("map"),
  splitLayout(
    tableOutput("data"),
    plotOutput("plot", height = "300px", width = "700px"),
    cellWidths = c("35%", "65%")
  ),
  verbatimTextOutput("plot_text")
)

# Define the Shiny server
server <- function(input, output, session) {
  # Data in pg/L
  wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")
  
  # Data preparation
  # Remove samples (rows) with total PCBs = 0
  wdc.1 <- wdc[!(rowSums(wdc[, c(14:117)], na.rm = TRUE) == 0), ]
  # Calculate total PCB
  tpcb <- rowSums(wdc.1[, c(14:117)], na.rm = TRUE)
  # Create data.frame
  wdc.2 <- data.frame(
    SiteID = with(wdc.1, SiteID),
    date = with(wdc.1, as.Date(SampleDate, format = "%m/%d/%y")),
    Latitude = with(wdc.1, as.numeric(Latitude)),
    Longitude = with(wdc.1, as.numeric(Longitude)),
    tPCB = as.numeric(tpcb)
  )
  
  # Render the map
  output$map <- renderLeaflet({
    leaflet(wdc.2) %>%
      addTiles() %>%
      addMarkers(
        lng = ~Longitude,
        lat = ~Latitude,
        label = ~as.character(SiteID),
        labelOptions = labelOptions(noHide = TRUE)
      )
  })
  
  # Render the table
  output$data <- renderTable({
    if (!is.null(input$map_marker_click)) {
      siteid <- input$map_marker_click$id
      filtered_data <- subset(wdc.2, SiteID == siteid)[, c("SiteID", "date", "tPCB")]
      filtered_data$date <- format(as.Date(filtered_data$date, format = "%m/%d/%y"), "%m-%d-%Y")
      colnames(filtered_data)[3] <- paste("\u03A3", "PCB ", "(ng/L)", sep = "")
      
      # Sort the data by date
      filtered_data <- filtered_data[order(as.Date(filtered_data$date, format = "%m-%d-%Y")), ]
      
      # Get the number of samples
      num_samples <- nrow(filtered_data)
      
      colnames(filtered_data)[1] <- paste("SiteID (n =", num_samples, ")")
      
      return(filtered_data)
    } else {
      return(data.frame())
    }
  })
  
  output$plot <- renderPlot({
    if (!is.null(input$map_marker_click)) {
      siteid <- input$map_marker_click$id
      filtered_data <- subset(wdc.2, SiteID == siteid)
      
      if (nrow(filtered_data) == 0) {
        # No data available for the selected SiteID
        return(NULL)
      }
      
      # Aggregate data by week and calculate the average of PCB values
      filtered_data$week <- cut(filtered_data$date, breaks = "week")
      data_agg <- aggregate(tPCB ~ week, data = filtered_data, mean)
      
      num_values <- nrow(data_agg)
      width <- ifelse(num_values > 5, 0.8, 0.2 + (num_values * 0.1))
      
      p <- ggplot(data_agg, aes(x = week, y = tPCB)) +
        geom_col(fill = "steelblue", width = width) +
        labs(x = NULL, y = paste("\u03A3", "PCB (ng/L)", sep = "")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 12),  # Increase the font size of y-axis numbers
          axis.title.y = element_text(size = 14),  # Increase the font size of y-axis title
          axis.text = element_text(size = 8)
        )
      
      # Check if y-values are too large for linear scale
      max_value <- max(data_agg$tPCB)
      if (max_value >= 50000) {
        p <- p + scale_y_log10()
      }
      
      return(p)
    } else {
      return(NULL)
    }
  })
  
  output$plot_text <- renderPrint({
    if (!is.null(input$map_marker_click)) {
      cat("Plots are showing the aggregated data per week.\n")
      cat("If the maximum tPCB is too large (>80000 ng/L), the y-axis changes to log10 scale.\n")
      cat("Source:")
    }
  })  
  
  observe({
    leafletProxy("map", data = wdc.2) %>%
      clearMarkers() %>%
      addMarkers(
        lng = ~Longitude,
        lat = ~Latitude,
        layerId = ~SiteID,
        popup = ~paste(
          "SiteID: ", SiteID, "<br>",
          "Latitude: ", Latitude, "<br>",
          "Longitude: ", Longitude, "<br>",
          "Number of Samples: ", as.character(table(wdc.2$SiteID)[as.character(SiteID)])
        )
      )
  })
}

# Run the Shiny app
shinyApp(ui, server)
