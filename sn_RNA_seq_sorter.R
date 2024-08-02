# Load required R packages
library(shiny)
library(shinythemes)
library(readxl)
library(dplyr)
library(openxlsx)

#Enables large files to be processed
options(shiny.maxRequestSize = 30*1024^2)  # 30 MB

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  "snRNAseq Data Sorter",
                  tabPanel("Process Data File",
                           sidebarPanel(
                             tags$h3("Select File:"),
                             fileInput("file1", "Choose Excel File",
                                       accept = c(".xlsx", ".xls", ".xltx")),
                             actionButton("process", "Process Data"),
                             br(),
                             br(),
                             downloadButton("download", "Download Processed Data")
                           ),
                           
                           mainPanel(
                             textOutput("message")
                           ) # mainPanel
                           
                  ), 
                  tabPanel("User Instructions", 
                           HTML("<p>Must remove spaces from category names (e.g. Change 'ALS - Borderline' to 'ALS-BL') prior to processing</p>"),
                           HTML("<p>Input file should resemble the following format:</p>"),
                           tags$img(src = "https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEjC3jMAUjKx2PARssYWEzBDdunFA575JdYv84PxXL68Y-2IYnh75WKh7p7lrjJLVteIjGle1xqw-N-kI4rWjuoLngom08-H8f8Qd2-k6NfxLPiosqxPnaI0rozkx_M04QYypSzSZ3HmlsMkjMOyieazzBvkku-m9JTme2P3w5eZcKTJEVZtdV5hRQ8hyphenhyphenh8/w627-h496/sample_file.png", height = "300px", width = "400px")
                  ) 
                )
) # end of fluidPage


# Define server function  
server <- function(input, output, session) {
  
# Initiate reactive values to store data
  data <- reactiveValues(
    excel_data = NULL,
    processed_data = NULL
)
  
# Reactive value to store processing message
  processing_message <- reactiveVal("")
  
# Read uploaded Excel file
  observeEvent(input$file1, {
    req(input$file1)  # Ensure file is uploaded
    data$excel_data <- read_excel(input$file1$datapath)
    output$message <- renderText({
      paste("File", input$file1$name, "uploaded successfully.")
    })
  })

# Process data when button is clicked
  observeEvent(input$process, {
    req(data$excel_data)  # Ensure data is loaded
    processed_data <- data$excel_data
    
    # Processing code
    
    #Identify unique items in each category
    unique_cell_types <- sort(unique(processed_data$Cell_Type))
    unique_tissue_types <- sort(unique(processed_data$Tissue))
    unique_diagnois_types <- unique(processed_data$Diagnosis)
    
    #Determine all possible combinations of unique values
    all_combinations <- expand.grid(
      unique_diagnois_types = unique_diagnois_types,
      unique_cell_types = unique_cell_types,
      unique_tissue_types = unique_tissue_types,
      stringsAsFactors = FALSE
    )
    
#Initialize empty list for organized data to later be stored in
organized_data <- list()

#Add each category as a data frame in the organized_data list
#Add columns for gene lists and gene counts to the organized_data file
for (i in 1:nrow(all_combinations)) {
  names <- as.character(all_combinations[i, ])
  new_df <- data.frame(
    'Log2FC_0.25' = character(1),
    'Count_0.25' = 0,
    'Upregulated_0.25' = character(1),
    'Count_0.25_Upregulated' = 0,
    'Downregulated_0.25' = character(1),
    'Count_0.25_Downregulated' = 0,
    'Log2FC_1' = character(1),
    'Count_1' = 0,
    'Upregulated_1' = character(1),
    'Count_1_Upregulated' = 0,
    'Downregulated_1' = character(1),
    'Count_1_Downregulated' = 0,
    'Log2FC_1.5' = character(1),
    'Count_1.5' = 0,
    'Upregulated_1.5' = character(1),
    'Count_1.5_Upregulated' = 0,
    'Downregulated_1.5' = character(1),
    'Count_1.5_Downregulated' = 0,
    'Log2FC_2' = character(1),
    'Count_2' = 0,
    'Upregulated_2' = character(1),
    'Count_2_Upregulated' = 0,
    'Downregulated_2' = character(1),
    'Count_2_Downregulated' = 0,
    stringsAsFactors = FALSE
  )
  organized_data[[paste(names, collapse = ' ', sep = '')]] <- new_df
}

# Filter for significant values and store in data_signif variable
data_signif <- processed_data %>% filter(padj < 0.05)

# Processing loop to fill in columns of organized_data for each category
for (i in seq_along(organized_data)) {
  df <- organized_data[[i]]
  df_name <- names(organized_data[i])
  words <- strsplit(df_name, " ")[[1]]
  df_diagnosis <- words[1]
  df_cell <- words[2]
  df_tissue <- words[3]
  
  for (j in 1:nrow(data_signif)) {
    gene_j <- data_signif$Gene[j]
    fold_change_j <- data_signif$log2FoldChange[j]
    diagnosis_j <- data_signif$Diagnosis[j]
    tissue_j <- data_signif$Tissue[j]
    cell_type_j <- data_signif$Cell_Type[j]
    
    # Fill in log2FC > 0.25 genes
    if ( abs(fold_change_j) > 0.25 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) ) 
      
    {
      
      # Index of first empty row ("") in 'data' column
      first_empty_index_0.25 <- which(df$Log2FC_0.25 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_0.25, 'Log2FC_0.25'] <- gene_j
      
      logFC0.25 <- TRUE
      
      df[1, 'Count_0.25'] <- df[1, 'Count_0.25'] + 1
      
    }
    
    else {
      logFC0.25 <- FALSE
    }
    
    if ( fold_change_j > 0.25 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
      
    {
      # Index of first empty row ("") in 'data' column
      first_empty_index_upreg_0.25 <- which(df$Upregulated_0.25 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_upreg_0.25, 'Upregulated_0.25'] <- gene_j
      
      df[1, 'Count_0.25_Upregulated'] <- df[1, 'Count_0.25_Upregulated'] + 1
      
    }
    
    if ( fold_change_j < -0.25 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
      
    {
      # Index of first empty row ("") in 'data' column
      first_empty_index_downreg_0.25 <- which(df$Downregulated_0.25 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_downreg_0.25, 'Downregulated_0.25'] <- gene_j
      
      df[1, 'Count_0.25_Downregulated'] <- df[1, 'Count_0.25_Downregulated'] + 1
      
    }
    
    # Fill in log2FC > 1 genes
    if ( abs(fold_change_j) > 1 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
    {
      
      # Index of first empty row ("") in 'data' column
      first_empty_index_1 <- which(df$Log2FC_1 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_1, 'Log2FC_1'] <- gene_j
      logFC1 <- TRUE
      df[1, 'Count_1'] <- df[1, 'Count_1'] + 1
      
    }
    
    else {
      logFC1 <- FALSE
    }
    
    if ( fold_change_j > 1 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
      
    {
      # Index of first empty row ("") in 'data' column
      first_empty_index_upreg_1 <- which(df$Upregulated_1 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_upreg_1, 'Upregulated_1'] <- gene_j
      
      df[1, 'Count_1_Upregulated'] <- df[1, 'Count_1_Upregulated'] + 1
      
    }
    
    if ( fold_change_j < -1 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
      
    {
      # Index of first empty row ("") in 'data' column
      first_empty_index_downreg_1 <- which(df$Downregulated_1 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_downreg_1, 'Downregulated_1'] <- gene_j
      
      df[1, 'Count_1_Downregulated'] <- df[1, 'Count_1_Downregulated'] + 1
      
    }
    
    # Fill in log2FC > 1.5 genes
    if ( abs(fold_change_j) > 1.5 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
    {
      
      # Index of first empty row ("") in 'data' column
      first_empty_index_1.5 <- which(df$Log2FC_1.5 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_1.5, 'Log2FC_1.5'] <- gene_j
      
      logFC1.5 <- TRUE
      
      df[1, 'Count_1.5'] <- df[1, 'Count_1.5'] + 1
      
    }
    
    else {
      logFC1.5 <- FALSE
    }
    
    if ( fold_change_j > 1.5 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
      
    {
      # Index of first empty row ("") in 'data' column
      first_empty_index_upreg_1.5 <- which(df$Upregulated_1.5 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_upreg_1.5, 'Upregulated_1.5'] <- gene_j
      
      df[1, 'Count_1.5_Upregulated'] <- df[1, 'Count_1.5_Upregulated'] + 1
      
    }
    
    if ( fold_change_j < -1.5 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
      
    {
      # Index of first empty row ("") in 'data' column
      first_empty_index_downreg_1.5 <- which(df$Downregulated_1.5 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_downreg_1.5, 'Downregulated_1.5'] <- gene_j
      
      df[1, 'Count_1.5_Downregulated'] <- df[1, 'Count_1.5_Downregulated'] + 1
      
    }
    
    # Fill in log2FC > 2 genes
    if ( abs(fold_change_j) > 2 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
      
    {
      
      # Index of first empty row ("") in 'data' column
      first_empty_index_2 <- which(df$Log2FC_2 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_2, 'Log2FC_2'] <- gene_j
      
      logFC2 <- TRUE
      
      df[1, 'Count_2'] <- df[1, 'Count_2'] + 1
      
    }
    
    else {
      logFC2 <- FALSE
    }
    
    if ( fold_change_j > 2 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
      
    {
      # Index of first empty row ("") in 'data' column
      first_empty_index_upreg_2 <- which(df$Upregulated_2 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_upreg_2, 'Upregulated_2'] <- gene_j
      
      df[1, 'Count_2_Upregulated'] <- df[1, 'Count_2_Upregulated'] + 1
      
    }
    
    if ( fold_change_j < -2 &
         (df_diagnosis == diagnosis_j) &
         (df_cell == cell_type_j) &
         (df_tissue == tissue_j) )
      
    {
      # Index of first empty row ("") in 'data' column
      first_empty_index_downreg_2 <- which(df$Downregulated_2 == "")[1]
      
      # Fill the empty row in log2FC column with Gene value
      df[first_empty_index_downreg_2, 'Downregulated_2'] <- gene_j
      
      df[1, 'Count_2_Downregulated'] <- df[1, 'Count_2_Downregulated'] + 1
      
    }
    
    #If a gene got added to any category...
    if ( logFC0.25 == TRUE | logFC1 == TRUE | logFC1.5 == TRUE | logFC2 == TRUE ) 
      
    {
      # Add a new row to the data frame df
      df <- rbind(df, new_df)
      
      # Update the data frame in organized_data with the modified df
      organized_data[[i]] <- df
      
    }
  }
}

# Save processed data in reactiveValues
data$processed_data <- organized_data

# Display message once processing is complete
processing_message("Data processed successfully!")

# Return the processing message to be displayed
output$message <- renderText({
  processing_message() })
  
  })

# Download processed data when button is pressed
  output$download <- downloadHandler(
    filename = function() {
      "processed_scRNAseq_data.xlsx"
    },
    content = function(file) {
      # Save organized_data to multi-worksheet Excel file
      wb <- openxlsx::createWorkbook()
      for (df_name in names(data$processed_data)) {
        openxlsx::addWorksheet(wb, sheetName = df_name)
        openxlsx::writeData(wb, sheet = df_name, x = data$processed_data[[df_name]], 
                            startCol = 1, startRow = 1)
      }
      openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
    }
  )
  
} # server


# Create Shiny object
shinyApp(ui, server)