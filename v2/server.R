# Goals for v2
# Utilize the whole expression dataset
# Add pearson and spearman correlation
# Add heatmaps for summarizing results?

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

server <- function(input, output){
  
  
  options(shiny.maxRequestSize=50*1024^2) 
  
  
  output$brushtop5 <- renderTable({
    
    validate(
      need(input$brushtop5, "Draw a rectangle around data points for further information")
    )
    
    brushedPoints(top5_df_brush, input$brushtop5)
  }, striped = T)
  
  
  shinyjs::disable("download_res")
  shinyjs::disable("download_top5")
  
  observeEvent(analyzed_df(), {
    shinyjs::enable("download_res")
    shinyjs::enable("download_top5")
  })
  
  
  # output$downloadData_ui <- renderUI({
  #    req(analyzed_df())
  #    downloadButton("downloadData", "Download results", class="down_but")
  #  })  #THIS ALSO WORKS TO SHOW BUTTON AFTER ANALYSIS
  
  
  ################################################################################################################################
  # Define conditional dynamic file upload prompted when user selects "Custom" as reference data
  output$ui_sel_ref <- renderUI ({
    
    if (input$sel_reference == "Custom"){
      fileInput("ref_file", "Upload custom reference file",
                multiple = F, 
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"))
    }
  })
  
  
  output$ui_sel_ref_annot <- renderUI ({
    
    if (input$sel_reference == "Custom"){
      fileInput("annot_file", "Upload custom annotation file",
                multiple = F, 
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"))
    }
  })
  
  
  output$sample_data_file <- renderImage({
    
    list(src = "data/cluster_expr_IMG.png",
         alt = "Sample SCseq data",
         width=500)
  
  }, deleteFile = F)
  
  output$sample_reference_file <- renderImage({
    
    list(src = "data/ref_data_IMG.png",
         alt = "Sample reference data",
         width=500)
    
  }, deleteFile = F)
  
  output$sample_annotation_file <- renderImage({
    
    list(src = "data/ref_annot_IMG.png",
         alt = "Sample annotation data",
         width=500)
    
  }, deleteFile = F)
  
  ################################################################################################################################
  # Read uploaded differential expression file
  expr_data <- reactive({
    
    inFile <- input$expr_file
    
    
    
    if(is.null(inFile) & input$example_data == T){
      
      
      
      dat <- readRDS("data/til_scseq_exprs_mean_subset.rds")
      
      req(input$run)
      
      dat
      
      
      
    } else {
      
      validate(
        need(input$expr_file != "", "Please upload a data set or use example data")
      )
      
      # Make sure the file type is correct
      validate(
        need(tools::file_ext(inFile$name) %in% c(
          'text/csv',
          'text/comma-separated-values',
          'text/plain',
          'csv'
        ), "Wrong File Format. File needs to be a .csv file."))
      

      
      dat <- read.csv(inFile$datapath, check.names=TRUE, strip.white = TRUE, stringsAsFactors = F)

      
      # Make sure the column names are proper for correct subsetting
      validate(
        need(
          {if({sum(grepl(pattern = "gene", colnames(dat), ignore.case=T)) == 1}) TRUE else FALSE},
          "Formatting error: Make sure your dataset contains a column named 'gene' (capitalization is not important, duplicate gene column is not allowed)"
      )
      )
      
      
      
      req(input$run)
      
      gene_column <<- grep("gene", colnames(dat), ignore.case = T, value = T)
      
      dat[,gene_column] <- tolower(dat[,gene_column])
      
      dat
      
    }
    
  }) # close expr_data reactive object
  
  
  ################################################################################################################################
  # Read reference dataset
  
  ref_data <- reactive({
    
    if(input$sel_reference == "ImmGen"){
      
      reference <- readRDS("data/immgen_combined.rds")
      
      ref_gene_column <<- grep("gene", colnames(reference), ignore.case = T, value = T)
      
    } else if (input$sel_reference == "Custom"){
      
      in_refFile <- input$ref_file
      
      reference <- read.csv(in_refFile$datapath, check.names=FALSE, strip.white = TRUE, stringsAsFactors = F)
      
      ref_gene_column <<- grep("gene", colnames(reference), ignore.case = T, value = T)
      
      reference[,ref_gene_column] <- tolower(reference[,ref_gene_column])
      
      
    } else {NULL}
    
    
    
    reference
    
    
  })
  
  ################################################################################################################################
  # Read immgen annotation file for explanations of cell types
  
  
  reference_annotation <- reactive({
    
    if(input$sel_reference == "ImmGen"){
      
      ref_annotation <- readRDS("data/imm_annot.rds")
      ref_annotation
      
    } else if(input$sel_reference == "Custom"){
      
      annotFile <- input$annot_file
      
      ref_annotation <- read.csv(annotFile$datapath, check.names=FALSE, strip.white = TRUE, stringsAsFactors = F)
      ref_annotation
      
    }
  })
  
  ################################################################################################################################
  # Define a reactive cluster object that will store cluster information
  clusters <- reactive({
    
    
    gtools::mixedsort(
      levels(
        as.factor(
          colnames(expr_data())[!grepl("gene", colnames(expr_data()), ignore.case = T)]
          )
        )
      )
    
  }) # close clusters reactive object
  
  
  
  ################################################################################################################################
   # Compare expr_data against reference file
  analyzed_df <- reactive({
    

    
   
    dat_genes <- expr_data()[gene_column] %>% pull() %>% as.character
    ref_genes <- ref_data()[ref_gene_column] %>% pull() %>% as.character
    
    common_genes <- intersect(dat_genes, ref_genes)
    
    trim_dat <- expr_data() %>%
      filter(!!rlang::sym(gene_column) %in% common_genes) %>%
      arrange_(.dots=gene_column) %>%
      select_(.dots = paste0("-", gene_column))
    
    trim_ref <- ref_data() %>%
      filter(!!rlang::sym(ref_gene_column) %in% common_genes) %>%
      arrange_(.dots=ref_gene_column) %>%
      select_(.dots = paste0("-", ref_gene_column))
    
    
    
    master_df <- data.frame()
    
    withProgress(message = 'Analysis in progress', value = 0, {
      
      if(input$comp_method == "Spearman") comp_method = "spearman" else comp_method = "pearson"
      
      for (i in clusters()) {
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(clusters()), detail = paste("Analyzing cluster", i))
        
        
        cor_df <- cor(trim_dat[i], trim_ref, method = comp_method)
        
        
        df <- data.frame(cor_coefficient = cor_df[1,])
        
        df <- rownames_to_column(df, var="reference_id")
        
        # FIX THIS LATER AFTER UPDATING THE IMMGEN/CUSTOM REFERENCE UPLOADING CODE ABOVE SECTION
        
        if(input$sel_reference == "ImmGen"){
          
          df <- left_join(df, reference_annotation(), by=c("reference_id" = "short_name"))
          
          
          
          
          
        } else if (input$sel_reference == "Custom" & !is.null(input$annot_file)){
          
          df <- left_join(df, reference_annotation(), by=c("reference_id" = "short_name"))
          
          
        } else if(input$sel_reference == "Custom" & is.null(input$annot_file)){
          
          df$reference_cell_type <- rep("Upload annotation file", dim(ref_data())[2]-1)
          df$short_name <- colnames(ref_data())[!colnames(ref_data()) %in% ref_gene_column]
          df$long_name <- rep("Upload annotation file", dim(ref_data())[2]-1)
          df$description <- rep("Upload annotation file", dim(ref_data())[2]-1)
          
        }
        
        
        
        
        df$cluster <- i
        
        # Add confidence-of-prediction calculations here and append to the df
        # Calculate the mean and standard deviation of the aggregate scores per reference cell type
        mean_cor_coeff <- mean(df$cor_coefficient)
        cor_coeff_sd <- sd(df$cor_coefficient)
        
        # Calculate the distance of the identity score from population mean (how many std devs apart?)
        df$z_score <- (df$cor_coefficient - mean_cor_coeff)/cor_coeff_sd
        
        
        master_df <- rbind(master_df,df)
        
        
      } # close for loop that iterates over clusters
      
    }) # Close with progress
        
        
    master_df  
        
    
  }) # close analyzed_df reactive expression
  
  
  
  ################################################################################################################################
  # Generate plotting area dynamically for individual cluster plots
  
  
  # Insert the right number of plot output objects into the web page (https://gist.github.com/wch/5436415/)
  output$plots <- renderUI({
    plot_output_list <- lapply(clusters(), function(i) {
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname, height = 500, width = 1800, brush = "brush") # optimize plotting area
    }
    ) # close lapply 
    
    # Convert the list to a tagList - this is necessary for the list of items to display properly.
    do.call(tagList, plot_output_list)
    
  }) # close output$plots renderUI
  
  
  
  ################################################################################################################################
  # Prepare individual plots
  
  
  
  # # Call renderPlot for each one. Plots are only actually generated when they
  # # are visible on the web page.
  
  observe({
    
    withProgress(message = 'Analyzing', value = 0, {
      
      for (i in clusters()) {
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/length(clusters()), detail = paste("Cluster", i))
        
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          
          df_plot <- analyzed_df() %>%
            filter(cluster == i)
          
          score_mean <- mean(df_plot$cor_coefficient)
          score_sd <- sd(df_plot$cor_coefficient)
          
          
          
          my_i <- i
          plotname <- paste("plot", my_i, sep="")
          
          output[[plotname]] <- renderPlot({
            
            df_plot_brushed <<- df_plot
            
            p <- ggdotplot(df_plot, x = "reference_id", y="cor_coefficient", 
                           fill = "reference_cell_type", xlab=F, ylab="Correlation coefficient",
                           font.y = c(14, "bold", "black"), size=1, x.text.angle=90,
                           title = paste("Cluster:",my_i), font.title = c(15, "bold.italic"),
                           font.legend = c(15, "plain", "black"))+
              theme(axis.text.x = element_text(size=10, vjust=0.5, hjust=1))+
              geom_hline(yintercept=score_mean)+
              annotate("rect", xmin = 1, xmax = length(df_plot$reference_id),
                       ymin = score_mean-score_sd, ymax = score_mean+score_sd,
                       fill = "gray50", alpha = .1)+
              annotate("rect", xmin = 1, xmax = length(df_plot$reference_id),
                       ymin = score_mean-2*score_sd, ymax = score_mean+2*score_sd,
                       fill = "gray50", alpha = .1)
            
            
            
            # Old iteration using ggdotchart function. It reorders X axis.
            # p <- ggdotchart(df_plot, x = "reference_id", y="cor_coefficient", 
            #                 group = "reference_cell_type", color = "reference_cell_type", xlab=F, ylab="Reference identity score",
            #                 font.y = c(14, "bold", "black"),
            #                 dot.size = 3, title = paste("Cluster:",my_i), font.title = c(15, "bold.italic"),
            #                 font.legend = c(15, "plain", "black"))+
            #   theme(axis.text.x = element_text(size=10, vjust=0.5, hjust=1))
            
            
            print(p)
            
          }) # close renderPlot
          
        }) # close local
        
      } # close for loop 
      
    })  # close withProgress
    
    
  }) #close observe
  
  
  ################################################################################################################################
  # Prepare top5 summary plots
  
  top_df <- reactive({
    
    
    top5_df <- analyzed_df() %>%
    group_by(cluster) %>%    #cluster
    top_n(5, wt = cor_coefficient) %>%
    arrange(as.numeric(cluster), cluster, desc(cor_coefficient))
  
  
  ordered_cluster_levels <- gtools::mixedsort(levels(as.factor(top5_df$cluster)))
  
  
  top5_df$cluster <- factor(top5_df$cluster, levels = ordered_cluster_levels)
  
  
  top5_df$index <- 1:nrow(top5_df)
  
  top5_df <- select(top5_df, cluster,
                    reference_cell_type,
                    reference_id,
                    long_name,
                    description,
                    cor_coefficient,
                    index, everything())
  
  
  
  top5_df_brush <<- top5_df
  
  top5_df
  
  })
  
  
  output$top5 <- renderPlot({
    
   top_plot <- top_df()
    
    ggdotplot(top_plot, x="index", y="cor_coefficient", 
              fill = "cluster", size=1, x.text.angle=90, 
              font.legend = c(15, "plain", "black")) +
      scale_x_discrete(labels=top_plot$reference_id)+
      theme(axis.text.x = element_text(vjust=0.5, hjust=1))
    
  })
  
  
  ################################################################################################################################
  # Download results 
  
  
  output$download_res <- downloadHandler(
    filename = "Identity_scores_all.csv",
    content = function(file) {
      write.csv(analyzed_df(), file, row.names = FALSE)
    }
  )
  
  
  output$download_top5 <- downloadHandler(
    filename = "Identity_scores_top5.csv",
    content = function(file) {
      
    write.csv(top_df(), file, row.names = FALSE)
    }
  )
  
  
} # close server function





## How to generate expression data frame summarizing gene expression per cluster:

# combined <- readRDS("combined_w_tsne_cd4_added.rds")
# 
# exprs <- df_extractor(combined, humanize = F)
# 
# saveRDS(exprs, "exprs.rds")


## Prepare data for CIPR_v2 (02/05/19)
# library(data.table)
# exprs <- readRDS("exprs.rds")
# 
# exprs <- exprs[,!colnames(exprs)%in% c("Sample", "Cell_id")]
# 
# exprs <- data.table(exprs)
# 
# exprs_mean <- exprs[,lapply(.SD, mean), by=Cluster]
# 
# exprs_mean <- data.frame(exprs_mean)
# 
# 
# cluster_names <- exprs_mean$Cluster
# gene_names <- colnames(exprs_mean)
# gene_names <- gene_names[gene_names != "Cluster"]
# 
# exprs_mean2 <- as.data.frame(t(exprs_mean[,colnames(exprs_mean) != "Cluster"]))
# colnames(exprs_mean2) <- cluster_names
# 
# exprs_mean2 <- tibble::add_column(exprs_mean2, Gene=gene_names, .after = 0)
# 
# saveRDS(exprs_mean2, "c:/OConnell Lab/R analyses/Shiny Apps/Cell_identity_predictor/deploy_v2/data/til_scseq_exprs_mean.rds")




# shinyApp(ui=ui, server=server)