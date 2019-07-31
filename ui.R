library(ggplot2)
library(markdown)
nki <- readRDS("breastCancerNKI_200genes.rds")
source("cluster_method.R")

navbarPage("ShinyGene",
#------tab1----------------------------                  
        tabPanel("Background Information", 
                 div(id = "description", 
                     
          h3("About the dataset"), 
          p("The data is a subset (1000 genes) of the nki data from the R package 'breastCancerNKI'. Briefly, this dataset 
            comprises the gene expression profiling of 337 breast cancer patients. There is extensive clinical data available, 
            but we will only focus on the ER (Estrogen Receptor) status, grade and overall survival of each patient. The two
           variables (ER: either positive or negative; grade: 1-3) are known to predict survival."),
          p("The original data could be found in the package below:"),
          tags$a(href="http://bioconductor.org/packages/release/data/experiment/html/breastCancerNKI.html", " breastCancerNKI"),
           
           h3("About the APP"), 
           p("In this project, we sought to produce a tool that will explore the normalized RNA expression data and their
              correlations with ER status by using two machine learning methods: tree and boosting. Finally, the user may 
              view a table of the genes selected by the threshold they have set and these may be downloaded as a CSV file 
              with a button at the bottom of the sidebar.")

          )),
#------tab2-------------------------------------------
        tabPanel("Data exploration",
             sidebarLayout(
                     sidebarPanel(
                         h3("Summary analysis"),
                         selectizeInput("gene", "Gene symbol", selected = "GREM2", choices = levels(as.factor(names(nki)))),
 #                        br(),
  #                       selectInput("group","Group:", choices = c("er","os","grade")),
                         
                         br(),
                         selectInput("plot.type","Plot Type:",
                                     list(boxplot = "boxplot", histogram = "histogram", density = "density", bar = "bar")
                         )
                         # Only show this panel if the box checked
#                        conditionalPanel(condition = "input.os == true",
#                                          checkboxInput("grade", "Also change symbol based on grade?")
                         
                         
                     ),
                     mainPanel(
                         plotOutput("boxplot1"),
                         downloadButton(outputId = "downloadPlot", label = "Download the plot"),
                         textOutput("title1"),
                         verbatimTextOutput("summaryT")
                     )
        )),
#----tab3-------------------------------------
      tabPanel("Cluster analysis",
              pageWithSidebar(
                     headerPanel('Gene clustering'),
                     sidebarPanel(
                       helpText(paste("Select distance:" )),
                       fluidRow(
                       selectInput("dmeth", NULL, choices=dmeths,
                                     selected=dmeths[1])),
                       helpText(paste("Select clustering method:" )),
                       fluidRow(
                         selectInput("meth", NULL, choices=cmeths,
                                     selected=cmeths[1])),
                       helpText(paste("Select height for cut:" )),
                       fluidRow(
                         numericInput("cutval", NULL, value=40, min=0, max=Inf, step=1)),
                       helpText(paste("Select variables for clustering from", substitute(df), ":" )),
                       fluidRow(
                         checkboxGroupInput("vars", NULL, choices=names(nki),
                                            selected=names(nki)[1:10]))
                       
                     ),
                     mainPanel(
                       tabsetPanel(
                         tabPanel("tree", 
                                  plotOutput("plot1"),
                                  downloadButton(outputId = "downloadPlot2", label = "Download the plot")),
                      
                         tabPanel("pairs", 
                                  plotOutput("pairsplot"),
                         downloadButton(outputId = "downloadPlot3", label = "Download the plot"))
#                         tabPanel("silh", 
#                                  plotOutput("silplot"))
                      )
                       )
      )),

#----tab4-------------------------------------
tabPanel("Modeling",
                 headerPanel("ML analysis"),
                 sidebarPanel(
                           h3('Gene selection for ML analysis'),
                           h3(''),                           
                           uiOutput("choose_genes"),
                           h3(''),
                           h3('Machine Learning'),
                           h3(''),
                           sliderInput('size', 'Choose sample size for training data ', min=0.2, max=0.9,
                                       value=min(0.6, 0.9), step=0.1),
                           selectInput('z','Choose ML-Method',c("Random Forest" = "rf","Decision Tree" = "rpart"),
                                       selected = "rf"),
                           conditionalPanel(condition = "input.z == 'rf'",
                                            sliderInput('m','Select mtry',min=1, max=20,
                                                        value=min(1, 20), step=1)),
                           h3(''),
                           h3('Tweek the Plot'),                           
                           uiOutput("gene_xaxis"),
                           uiOutput("gene_yaxis"),
                           downloadButton(outputId = "downloadPlot4", label = "Download the plot"),
                           h3('Math formula'),
                           uiOutput('ex1')
                         ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("ML analysis",
#                           h2('Introductory Exploratory Analysis'),
#                           plotOutput("mlplot"),
                           h4('Let\'s develop a machine learning algorithm to predict ER status. In this example, we will use all
                           200 features to develop our model. Choose the method from the sidebar.'),  
                           
                           h4('Training plot'),
                           
                           h4(''),
                           plotOutput("mlplot",click = "plot_click",  
                                      hover = hoverOpts(id = "plot_hover", delayType = "throttle")),
                           column(width = 4,  verbatimTextOutput("plot_clickinfo"),
                          verbatimTextOutput("plot_hoverinfo")),

                           h4(''),
                           h3('Prediction confusion matrix'),                   
                           verbatimTextOutput('confusionmatrix'),
                           h4('We can also see which data points are tuely/falsely predicted'),
                           plotOutput("mlplot2")  
                      ),
                     tabPanel("Data download",
                           div(style = 'height:300px; width:600px;overflow-y;overflow-y: scroll', 
                               tableOutput("table3")),
                           downloadLink("downloadDataML", "Download data in ML analysis")
                     )))
),


#-----tab5-------------------------------------      
      tabPanel("Acknowledgement",
               div(id = "description",
                   p("The code for cluster  were revised from the code in this online book: "),
               tags$a(href="https://genomicsclass.github.io/book/pages/bioc2_shiny.html", "Biomedical Data Science")
    )
)
)
      


