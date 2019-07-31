library(shiny)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RCurl)
library(psych)
library(cluster)
library(rpart)
library(visTree)
library(randomForest)
library(caret)


source("cluster_method.R")
nki <- readRDS("breastCancerNKI_200genes.rds") 
#nki <-  nki[,!duplicated(colnames(nki))]


shinyServer(function(input, output, session) {
  #---------------------subset data---------
  getData <- reactive({
            nki %>% subset(select=c(input$gene,"er","os","grade"))
	})  
  getData2 <- reactive({
  	        nki %>% subset(select=c("id",input$gene2,"er","os","grade"))
   	})


##--------- tab2 data summary characteristics--------------------
#create plot  
  boxplot1 <- function(){
    newData <- getData()
        attach(newData)
    gene1<- input$gene
 #   status <- input$group
    df1 <- newData[,c("er",gene1)]
    #df1$er <- factor(df1$er)
    colnames(df1)<-c("er","gene1")
#    plot.type <- switch(input$plot.type,
#                      "boxplot" 	= geom_boxplot(y=gene1),
#                      "histogram" =	geom_histogram(x=gene1,alpha=0.5,position="identity"),
#                      "density" 	=	geom_density(alpha=.75),
#                      "bar" 		=	geom_bar(position="dodge"))    
    df1$er <- factor(df1$er)
    colnames(df1)<-c("er","gene1")
    g <- ggplot(df1, aes(er,gene1)) + geom_boxplot(outlier.colour = "orange", outlier.shape = 1) + # plot.type +
      labs(title = paste("Plot of of the mRNA expression", input$gene, "between ER status"),y=input$gene,x="ER status")
  }
  
  output$boxplot1 <- renderPlot({
  	g <- boxplot1()
  	g
  })
  # Generate a summary of the data
  output$summaryT <- renderPrint({
    newData <- getData()
    gene1<- input$gene
    df1 <- newData[,c(gene1)]
    summary(df1)
  })

##----------------Tab 3 cluster--------------------------
#  selectedData <- reactive({
#    nki[, c(input$xcol, input$ycol)]
#  })
 df <- nki 
 nms = names(nki)
 clusterplot1 <- function(){
   xv = df[,input$vars]
    docs <- dist(as.matrix(xv), method = input$dmeth)
    hclust_dist<- as.dist(docs)
    hclust_dist[is.na(hclust_dist)] <- 0
    hclust_dist[is.nan(hclust_dist)] <- 0
    
    plot(hclust(hclust_dist, method=input$meth), labels = FALSE,
         xlab=paste(input$dmeth, "distance;", input$meth, "clustering"))
    abline(h=input$cutval, lty=2, col="gray")
   
 }
 
 output$plot1 <- renderPlot({
   clusterplot1()
  })
 pairplot <- function(){
       xv = df[,input$vars]
    pairs(data.matrix(xv))
   
 }
  output$pairsplot <- renderPlot({
    pairplot()
  })

##---------tab 4- machine learning analysis----------------

# RenderUI for gene selection
output$choose_genes <- renderUI({
   selectInput('bygene', 'Gene symbols', names(nki), multiple=TRUE, selectize=TRUE,
               selected=c("GREM2","SUHW2"))
  }) 
  
getData3 <- reactive({
    nki %>% subset(select=c("id",input$bygene,"er","os","grade"))
  })


  ###create output tablee of observations    
  output$table3 <- renderTable({
       getData3()
  })
#------------------Download sutset data used for machine learning---------------
  output$downloadDataML <- downloadHandler(
    filename = function() {
      paste("dataset", ".csv", sep ="")
    },
    content = function(file) {
      write.csv(getData3(), file, row.names = FALSE)
    },
    contentType = "csv"
  )
 
##Subset data for ML analysis 
#renderUI for gene selection for scatter plot
output$gene_xaxis <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(input$bygene))
      return()
    
    # Get the data set with the appropriate name
    dat <- getData3()  %>% select(-c(os,id,grade))#get(input$bygene)
    colnames <- names(dat)
    
    # Create the checkboxes and select them all by default
    selectInput("x", "Choose gene for x-axis", 
                       choices  = colnames,
                       selected = "er")
  })
  output$gene_yaxis <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(input$bygene))
      return()
    
    # Get the data set with the appropriate name
    dat <- getData3() %>% select(-c(os,id,grade))#get(input$bygene)
    colnames <- names(dat)
    
    # Create the checkboxes and select them all by default
    selectInput("y", "Choose gene for y-axis", 
                choices  = colnames,
                selected = "SUHW2")
  })    
  
#Split data
  df <- nki 
  inTrain <- reactive({
    set.seed(123)
    inTrain <- createDataPartition(y = df$er, p = input$size, list = FALSE) 
    inTrain
  })
  
  getFit <- reactive({
    df <- nki %>% subset(select=c(input$bygene,"er")) %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% na.omit()
    df$er <- factor(df$er)
    training <- df[inTrain(),]
    testing <- df[-inTrain(),] 
    if(input$z == 'rf'){
              modelFit <- randomForest(er ~., data = training, ntree=500,mtry=input$m)
    } else {modelFit <- rpart(er ~ ., data = training, method="class",control=rpart.control(minsplit=5,cp=0.001))
   modelFit 
   }
  })

  mlplot <-  function() {
  plot(getFit())
  }  
  output$mlplot <- renderPlot({
    mlplot()
  })

  output$plot_clickinfo <- renderPrint({
    cat("Click:\n")
    str(input$plot_click)
  })
  output$plot_hoverinfo <- renderPrint({
    cat("Hover (throttled):\n")
    str(input$plot_hover)
  })
  
##Prediction    
  output$confusionmatrix <- renderPrint({
    df <- nki %>% subset(select=c(input$bygene,"er")) %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% na.omit()
#    df_pred <- nki[,c(1:200,202)]
#    df_pred %>%  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% na.omit() -> df
     df$er <- factor(df$er)
    testing <- df[-inTrain(),]
    predictions <- predict(getFit(), newdata = testing, type="class")
    
    confusionMatrix(predictions, as.factor(testing$er))
  })
  
  mlplot2 <-  function() {
#    df_pred <- nki[,c(1:200,202)]
#    df_pred %>%  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% na.omit() -> df
    df <- nki %>% subset(select=c(input$bygene,"er")) %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% na.omit()
    df$er <- factor(df$er)
    testing <- df[-inTrain(),]
    predictions <- predict(getFit(), newdata = testing)
    testing$predRight <- ifelse(predictions == as.factor(testing$er),1,0)
    ggplot(data = testing, aes_string(x = input$x , y = input$y)) + aes(color = predRight) +
    geom_point(lwd = 4) +theme(text = element_text(size=20))

  }
  
  output$mlplot2 <- renderPlot({
    mlplot2()
  }) 
output$ex1 <- renderUI({
      withMathJax(
      helpText('The busy Cauchy distribution
               $$\\frac{1}{\\pi\\gamma\\,\\left[1 +
               \\left(\\frac{x-x_0}{\\gamma}\\right)^2\\right]}\\!$$'))
})

 
  
  #------------------Download boxplot ---------------  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("boxplot", '.png', sep='')
    },
    content = function(file) {
      ggsave(file, plot = boxplot1(), device = "png")
    }
  )
  
#------------------Download cluster plot ---------------  
  output$downloadPlot2 <- downloadHandler(
    filename = function() {
      paste("cluster", '.png', sep='')
    },
    content = function(file) {
      ggsave(file, plot = clusterplot1(), device = "png")
    }
  )  
  
  output$downloadPlot3 <- downloadHandler(
    filename = function() {
      paste("cluster", '.png', sep='')
    },
    content = function(file) {
      ggsave(file, plot = pariplot(), device = "png")
    }
  )  
 
#------------Download ML plot--------------  

  
  output$downloadPlot4 <- downloadHandler(
    filename = function() {
      paste("mlplot2", '.png', sep='')
    },
    content = function(file) {
      ggsave(file, plot=mlplot2(), h=6, w=6, units="in",dpi=300)#not work for mlplot
    }
  )  
  
  }
 
)
