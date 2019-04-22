ui <- fluidPage(
  
  shinyjs::useShinyjs(), # needed for download button to work
  
  tags$head(includeHTML(("data/google_analytics_v3.html"))), # To track user experience
  
  titlePanel(title=div(img(src="cipr_logo_small.png"), "Cluster identity predictor"), windowTitle = "CIPR"),
  
  sidebarLayout(
    
    sidebarPanel(width = 2,
                 
                 # Cluster markers (Diff exp analysis, columns must have "gene", "avg_logFC", and "cluster" data)
                 fileInput("data_file", "Upload cluster diff. exp. data",
                           multiple = F,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 # horizontal line
                 tags$hr(), 
                 
                 checkboxInput(inputId = "example_data", 
                               label = "Use example data (Ekiz HA and Huffaker TB, JCI Insight, 2019)", 
                               value = F),
                 
                 # horizontal line
                 tags$hr(), 
                 
                 radioButtons("sel_reference", 
                              label = "Select reference data", 
                              choices = c("ImmGen", "Custom"), 
                              selected = "ImmGen"), 
                 
                 
                 conditionalPanel("input.sel_reference == 'Custom'",
                                  uiOutput("ui_sel_ref")),
                 conditionalPanel("input.sel_reference == 'Custom'",
                                  uiOutput("ui_sel_ref_annot")),
                 
                 # horizontal line
                 tags$hr(), 
                 
                 radioButtons("comp_method", 
                              label = "Select method for comparisons", 
                              choices = c("logFC comparison", "Spearman", "Pearson"), 
                              selected = "logFC comparison"), 
                 
                 # horizontal line
                 tags$hr(), 
                 
                 
                 actionButton("run", label="Analyze"),
                 
                 br(),br(), br(),
                 
                 downloadButton("download_res", "Download results"),
                 
                 br(),br(),
                 
                 downloadButton("download_top5", "Download top5 hits")
                 
                 # uiOutput("downloadData_ui") #THIS ALSO WORKS TO SHOW BUTTON AFTER ANALYSIS
                 
                 
                 
                 
                 
    ),
    
    
    
    mainPanel(
      
      
      
      tabsetPanel(
        
        tabPanel("Top-5 hits",
                 
                 fluidRow(
                   plotOutput("top5", brush = "brushtop5", height="600px"),
                   tableOutput("brushtop5"),
                   p(strong("\n\nDefinition of table columns")),
                   tags$ul(
                     tags$li(strong("Cluster:"), "The name of the unknown cluster"),
                     tags$li(strong("Reference_cell_type:"), "Broad classification of the reference cell type"),
                     tags$li(strong("Reference_id:"), "Shortened unique identifier for reference cell type"),
                     tags$li(strong("Long_name:"), "Human-readable long name of the reference cell type"),
                     tags$li(strong("Description:"), "Detailed description of the reference cell type"),
                     tags$li(strong("Identity_score:"), "Identity score for the given reference cell type calculated via logFC comparison or correlation methods (see How-to-tab for details"),
                     tags$li(strong("Index:"), "Position in the dataframe (needed for ordered plotting)"),
                     tags$li(strong("Z-score:"), "How many standard deviations apart is the identity score from the average score for that cluster calculated across the reference cell subsets (see How-to tab for details)")
                     
                   )
                 )
        ),
        
        
        tabPanel("Individual Clusters",
                 
                 # This is the dynamic UI for the plots to be generated for each cluster
                 
                 fluidRow(
                   h4("Please wait while the images are being generated. This can take some time depending on the cluster numbers.\n"),
                   p("In the plots below, horizontal lines indicate the mean identity score for each experimental cell cluster, and the shaded bands indicate +/- 1 and 2 standard deviations from the mean.\n\n"),
                   uiOutput("plots")

                 )       
        ),
        
        

        tabPanel("How to use this program",
                 
                 h3("Summary"),
                 p("Understanding the biological identity of cell clusters in single cell RNA sequencing (SCseq) experiments can be challenging due to overlapping gene expression profiles. An accurate assessment of the cluster identity requires analyzing multiple differentially expressed genes as opposed to a low throughput approach of examining a few selected lineage-specific markers. Cluster Identity PRedictor (CIPR) compares user-provided SCseq cluster gene signatures with", a('Immunological Genome Project (ImmGen)', href= 'https://www.immgen.org'), "mouse immune cell datasets -- or with a user-provided reference dataset --, and calculates an", strong('identity score (IS)'), "for each SCseq cluster per reference cell subset. To calculate the IS, CIPR compares differential gene expression signatures of unknown clusters to that of known reference samples using three different methods: 1) Log fold-change comparison (dot product), 2) Spearman's correlation, and 3) Pearsons's correlation. For further information, please read the sections below"),
                 
                 br(),br(),
                 
                 h3("Quick Start"),
                 
                 p(strong(span(style="color:green", "To compare your SCseq clusters against ImmGen cell types, simply upload the cluster gene expression signatures and click 'Analyze'"))),
                 br(),
                 p(strong(span(style="color:blue", "If you'd like to compare SCseq data with a custom reference, please select the appropriate radio button and upload reference gene expression data"))),
                 br(),
                 p(strong(span(style="color:gray", "You can also run the program by using example data that our lab has generated (Ekiz HA and Huffaker TB, JCI Insight, 2019). This data were obtained by single cell transcriptomics analysis of flow-sorted murine tumor-infiltrating immune cells. Example gene expression data on this website is shortened to reduce computing times. In this example dataset, gene signatures of 4 clusters (activated T cells, natural killer cells, Langerhans dendritic cells, and plasmacytoid dendritic cells) are included, whereas the original dataset had 15 distinct cell clusters based on our analysis."))),
                 
                 br(),br(),
                 
                 
                 h3("Input data"),
                 tags$ul(
                   tags$li("Provide the expression data from SCseq cell clusters as a csv file."),
                   tags$li("This file", strong("must have at least three columns named Gene, logFC, and cluster"), ". Capitalization does not matter and extra characters before or after these column names are also allowed. The other columns are ignored."),
                   tags$li("To create plots with correct titles, avoid using spaces and special characters as column names"),
                   tags$li("Seurat's", strong("FindAllMarkers"), "function can be used to create this data frame. Please see below for example code to obtain this dataframe from your data.")
                   
                 ),
                 br(),
                 p(strong("Sample SCseq expression data (generated via Seurat's FindAllMarkers() function)")),
                 imageOutput("sample_data_file_logfc"),
                 
                 
                 h3("About reference dataset"),
                 tags$ul(
                   tags$li("This program is pre-loaded with ImmGen microarray data for easily investigating immune cell clusters in SCseq experiments."),
                   tags$li("If one would like to provide a custom reference dataset, any number of highthroughput data types can be used including microarray, RNAseq, and proteomics data, as long as data in the reference file are normalized together."),
                   tags$li("Reference data should be in linear space."),
                   tags$li("Custom reference datasets can be uploaded as a csv file and should contain gene expression data from known cell types."),
                   tags$li("Reference dataset should have genes in rows and cell types in columns."),
                   tags$li("The first column of the reference data frame must have gene names (e.g. Pdcd1, Tnfa) and must be named as", strong("Gene."), "Capitalization does not matter and characters before and after 'gene' are allowed."),
                   tags$li("This file can contain biological and technical replicates. In this case, IS is calculated and plotted for each replicate separately.")
                 ),
                 br(),
                 p(strong("Sample reference gene expression data")),
                 imageOutput("sample_reference_file"),
                 
                 
                 h3("About (optional) custom reference annotation data"),
                 tags$ul(
                   tags$li("If a custom reference dataset is used, although not necessary, an annotation file (in csv format) can be uploaded to obtain more informative plots."),
                   tags$li(strong("Annotation file must contain the columns named as 'short_name', 'long_name', 'description', and 'reference_cell_type'.")),
                   tags$li(strong("Data under 'short_name' column must EXACTLY match the column names of the reference gene expression data.")),
                   tags$li("Annotation file can have other columns as well, which will be ignored.")
                 ),
                 br(),br(),
                 p(strong("Sample reference annotation data")),
                 imageOutput("sample_annotation_file"),
                 
                 
                 h3("How is cell identity score calculated?"),
                 h4("If logFC comparison method is used:"),
                 
                 p("The algorithm works first by calculating gene expression signatures in the reference file, and then comparing the experimental cluster signatures with these reference signatures. Specifically the following processes are performed:"),
                 h5(strong("Pre-processing of reference dataset")),
                 tags$ul(
                   tags$li("The algorithm uses pre-loaded ImmGen signatures or a user-uploaded reference expression file (see below for the details of reference file format)."),
                   tags$li("The reference expression file should contain normalized gene expression values in linear space  in rows and cell types in columns. Reference data can be derived from high-throughput experiments such as microarray or RNAseq."),
                   tags$li("For each gene found in the reference file, algorithm calculates the ratio of gene expression in each cell type (i.e. individual columns) to the average gene expression value of the whole data set (i.e. all columns)."),
                   tags$li("Ratio is log transformed to obtain positive values indicating upregulation and negative values indicating downregulation in the specific cell type. Therefore, pre-processing of the reference file results in a data frame that features logFC values of each gene for each cell type.")
                 ),
                 h5(strong("Analysis of experimental cluster signatures against reference file")),
                 tags$ul(
                   tags$li("Uploaded cluster signature file should contain information about genes, logFC in expression levels and the cluster information (see below for the correct file format)."),
                   tags$li("Algorithm finds the common genes between experimental data and the reference dataset."),
                   tags$li("For each shared gene, logFC values of differentially expressed in clusters and are multiplied with the logFC values in the reference dataset. This way if a gene is upregulated or downregulated in both a given cluster and a cell type in the reference dataset (i.e. positive correlation), the multiplication will result in a positive number (i.e. multiplication of two positive or two negative numbers resulting in a positive number). Alternatively, if a gene is differentially regulated in opposite directions in cluster and reference cell type (i.e. negative correlation), multiplication of logFC values will result in a negative number."),
                   tags$li("Multiplied logFC values per each gene are added up, resulting in an aggregate IS across the dataset for each cluster in the experimental data. This way, genes that have similar expression patterns in the experimental cell cluster and the reference cell type contribute to a higher IS, whereas genes with opposite expression patterns will result in a lower IS."),
                   tags$li(strong("This way, each cluster in the SCseq experiment is analyzed against each known cell type in the reference file and scored for its overall similarity. A higher identity score corresponds to a positive correlation in gene expression profiles and hence higher chances of a cluster belonging to the particular cell type found in the reference file")),
                   tags$li("For each cell cluster in the experiment, aggregate IS of reference cell types are plotted in dot plots which shows reference cell types in the x-axis, and aggregate IS in the -axis."),
                   tags$li("A summary plot showing the reference cell types correspoint to the highest 5 IS score is also shown.")
                 ),
                 br(),
                 
                 
                 h4("If correlation methods are used:"),
                 p("The algorithm works first by subsetting the experimental and the reference datasets by using common genes found in both datasets. The second step involves a one-to-many (cluster-to-references) calculation of correlation coefficients. This way a single correlation coefficient is calculated for each reference cell type and for each cluster."),
                 
                 br(),
                 
                 p("This way, each cluster in the SCseq experiment is analyzed against each known cell type in the reference file and scored for its overall similarity. A positive and high  correlation coefficient corresponds to a positive correlation in gene expression profiles and hence higher chances of a cluster belonging to the particular cell type found in the reference file."), 
                 br(),
                 
        tags$li("For each cell cluster in the experiment, correlation coefficients of reference cell types are plotted in dot plots which shows reference cell types in the x-axis, and correlation coefficients in the y-axis."),
        tags$li("A summary plot showing the reference cell types corresponding to the highest 5 correlation coefficients is also shown."),
      
                 
                 h3("How confident is the prediction?"),
                 p("Since our algorithm compares unknown cluster gene signatures with signatures of reference samples, the confidence in the identity prediction can only be as high as the biological overlap between the experimental sample and the reference datasets. In an ideal scenario where gene expression data are available from all possible cell types under various experimental conditions, we envision that this algorithm can describe the identity of the unknown cell clusters in a highly accurate manner. However, in the absence of such data, our algorithm is still useful in characterizing unknown cluster identities by using a multiparametric approach. To help assess the confidence of the predictions, we utilize two different metrics:"),
                 tags$ul(
                   tags$li(strong("Deviation from mean:"), "For each unknown cluster, this metric is calculated by", strong(em("i) ")), "Averaging the correlation coefficients from the whole reference dataset", strong(em("ii) ")), "Subtracting the correlation coefficients of individual reference from this average", strong(em("iii) ")), "Dividing the difference by the standard deviation of the identity score across the whole reference dataset. This way, the deviation from mean is reported in the units of ", em("standard deviation."), "The higher deviation indicates a more pronounced distinction from the rest of the reference cell types, hence a higher confidence in prediction.")
                   
                 ), # close tags$ul
                 p("Although these metrics can be helpful in determining the confidence of prediction, they are not suitable to serve as a benchmark for the performance of our algorithm. A 'low-confidence' prediction can be explained by various factors including"),
                 tags$em(
                   tags$ul(
                     tags$li("A lack of comparable cell types in the reference dataset"),
                     tags$li("Using suboptimal number of metagene dimensions during clustering step of data analysis (which can result in  impure clusters composed of multiple cell types"),
                     tags$li("Presence of outliers in experimental data and clustering artifacts.")
                   )),
                 
                 p("Alternatively, a 'low-confidence' prediction can point to the discovery of new and interesting cell populations, and identify contaminating cell populations in the experiment. For instance, our experience with tumor-infiltrating lymphocyte single cell RNA-sequencing in murine B16F10 melanoma model suggests that, contaminating melanoma cells are loosely identified as ImmGen blood stem cells due to their proliferative capacity and stem-like gene expression profiles (Ekiz HA, manuscript in preparation)."),
                 
                 h3("Frequently asked questions"),
                 tags$ul(
                   tags$li(strong("Why am I getting 'Served Disconnected' error?"), "This can happen if the input and/or the custom reference data are not formatted properly. Please see guidelines above to make sure data are formatted correctly. If the problem persists, please contact us."),
                   tags$li(strong("Is CIPR limited to analyzing mouse data only?."), "Pre-loaded ImmGen data in CIPR is derived from sorted mouse immune cells. However, the pipeline converts the gene names in input and reference to lower case letters enabling the comparison of homologous genes between human and mouse. Most of the homologous gene names differ only in capitalization between human and mice, and this conversion is sufficient to compare the majority of genes found in the datasets. However, the user should use caution to derive conclusions while comparing datasets across different species. Alternatively, the user can prepare a custom reference dataset to have a more accurate cellular identity approximation.")
                 ),
        
        h3("Example R code to generate CIPR input"),
        h4("For logFC comparison method"),
        tags$pre(
"
library(Seurat)

# Generate a Seurat object by following Seurat vignettes 
# (https://satijalab.org/seurat/)
cluster_markers <- FindAllMarkers(seurat_obj)


head(cluster_markers)
# >  p_val avg_logFC pct.1 pct.2 p_val_adj   cluster  gene
# >1     0  1.794907 0.815 0.151         0   Clust_1 Cd8b1
# >2     0  1.760971 0.994 0.253         0   Clust_1  Cd3g
# >3     0  1.732805 0.832 0.116         0   Clust_1 Cxcr6
# >4     0  1.686574 0.807 0.133         0   Clust_1  Cd8a
# >5     0  1.487584 0.994 0.265         0   Clust_1  Cd3d
# >6     0  1.429320 0.763 0.119         0   Clust_1  Lag3



#  Save the cluster_markers object as .csv file
write.csv(cluster_markers, 'clustermarkers.csv', row.names=F)"),
br(),

h4("For correlation methods"),

tags$pre("

library(Seurat)
library(tibble)
library(data.table)


# Extract normalized gene expression counts
# Columns represent individual cells, rows represent genes
exprs <- as.data.frame(as.matrix(seurat_obj@data))

# Transpose the data frame
exprs <- as.data.frame(t(exprs))

# Add a column to indicate cluster belonging of each cell
exprs <- add_column(exprs, Cluster = seurat_obj@ident, .after = 0)

# Convert data frame to data.table for efficient processing
exprs <- data.table(exprs)

# Calculate the average expression of genes per cluster
exprs <- exprs[, lapply(.SD, mean), by=Cluster]

# Store cluster information as an object
cluster_names <- exprs$Cluster

# Store gene name information as an object
gene_names <- colnames(exprs)
gene_names <- gene_names[gene_names != 'Cluster']

# Discard 'Cluster' column and transpose the dataframe
exprs <- as.data.frame(t(exprs[, colnames(exprs) != 'Cluster']))

# Give cluster names to column names for CIPR to detect clusters
colnames(exprs) <- cluster_names

# Create a column named 'Gene' for CIPR to detect genes
exprs <- add_column(exprs, Gene=gene_names, .after = 0)

# Data frame in which rows represent individual cells
exprs[1:5,1:5]

# >Gene    Cluster_1  Cluster_2  Cluster_3  Cluster_4
# >Mrpl15  0.5939144  0.4352030  0.2814939  0.4863504
# >Lypla1  0.3457666  0.1579630  0.1724877  0.2503396
# >Tcea1   0.3821044  0.2566057  0.2998205  0.3587843
# >Atp6v1h 0.1757608  0.1582002  0.1243377  0.4227891
# >Rb1cc1  0.1338744  0.1592661  0.1462728  0.1576287
         

# Save the exprs object as a .csv file
write.csv(exprs, 'avg_exprs_per_cluster.csv')
         
         "),
     

                 h3("CIPR Release v4"),
        tags$ul(
          tags$li("Three different computational methods comparing differential expression data (logFC comparison, Spearman's and Pearson's correlation can now be used in the same user interface"),
          tags$li("ImmGen v1 and v2 data are renormalized from raw data files (.CEL), combined and pre-loaded as the reference. These two datasets contain partially overlapping cell types, and were generated using the same microarray platform (Affymetrix MoGene ST1.0 array). For further information please see ImmGen website and associated publications.")
        ),
                 
                 
                 h3("Contact us"),
                 p("For questions and comments:"),
                 p("atakan-dot-ekiz-at-hci.utah.edu"),
                 p("Ryan O'Connell Lab"),
                 p("Department of Pathology"),
                 p("University of Utah"),
                 p("Salt Lake City, Utah")
                 
                 
        )
        
        
        
      )
    ) 
  )
)
    

