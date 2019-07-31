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
    #    attach(newData)
    gene1<- input$gene
    df1 <- newData[,c("er",gene1)]
    df1$er <- factor(df1$er)
    colnames(df1)<-c("er","gene1")
    g <- ggplot(df1, aes(er, gene1)) + geom_boxplot(outlier.colour = "orange", outlier.shape = 1) +
      labs(title = paste("Boxplot of of the mRNA expression", input$gene, "between ER status"),y=input$gene,x="ER status")
    
  }
  
  output$boxplot1 <- renderPlot({
  	g <- boxplot1()
  	g
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

#source("ml_nki.R")  

  df <- nki[,c(1:200,202)] %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% na.omit() 
#  df_pred$er <- factor(df_pred$er)

    
  
  inTrain <- reactive({
    set.seed(123)
    inTrain <- createDataPartition(y = df$er, p = input$size, list = FALSE) 
    inTrain
  })
  
 
  
  getFit <- reactive({
    df$er <- factor(df$er)
    training <- df[inTrain(),]
    #  training$er <- factor(training$er, levels=c("0","1"))
    #  df$er <- factor(df$er)
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

    
  output$confusionmatrix <- renderPrint({
    df_pred <- nki[,c(1:200,202)]
    df_pred %>%  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% na.omit() -> df
     df$er <- factor(df$er)
     testing <- df[-inTrain(),]
    predictions <- predict(getFit(), newdata = testing, type="class")
    
    confusionMatrix(predictions, as.factor(testing$er))
  })
  
  mlplot2 <-  function() {
    df_pred <- nki[,c(1:200,202)]
    df_pred %>%  mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% na.omit() -> df
    df$er <- factor(df$er)
    testing <- df[-inTrain(),]
    predictions <- predict(getFit(), newdata = testing)
    testing$predRight <- ifelse(predictions == as.factor(testing$er),1,0)
    ggplot(data = testing, aes_string(x = input$x , y = input$y)) + aes(color = predRight) + geom_point(lwd = 4) +theme(text = element_text(size=20))

  }  
  output$mlplot2 <- renderPlot({
    mlplot2()
  }) 

  
##--------- tab5 data subset export--------------------  
  #update title info
  output$title1 <- renderText({
    paste0("Summary of the mRNA expression of", input$gene)
      })  
###create output for summary data    
  output$summaryT <- renderPrint({
    dt1 <- getData()
    g1 <- input$gene
    dt1 %>% subset(select=c(input$gene)) %>% summary()
    #tapply(dt$g1, dt$er, summary())
  })
  
#  observe({
#    if(input$os){
#    updateSliderInput(session,"size",max=10, min = 3)
#    } else {
#    updateSliderInput(session, "size", min = 1, max = 10)
#  }
#    })
###create output tablee of observations    
  output$table <- renderTable({
    dt2 <- getData2()

  })
  
  
#------------------Download data---------------
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("dataset", ".csv", sep ="")
    },
    content = function(file) {
      write.csv(getData2(), file, row.names = FALSE)
    },
    contentType = "csv"
  )

#------------------Download boxplot ---------------  
   output$downloadPlot <- downloadHandler(
          filename = function() {
              paste("boxplot", '.png', sep='')
          },
          content = function(file) {
              ggsave(file, plot = boxplot1(), device = "png")
          }
      )
#------------------Download plot 2---------------  
  output$downloadPlot2 <- downloadHandler(
    filename = function() {
      paste("cluster", '.png', sep='')
    },
    content = function(file) {
      ggsave(file, plot =   clusterplot1(), device = "png")
    }
  )  
  
  output$downloadPlot3 <- downloadHandler(
    filename = function() {
      paste("cluster", '.png', sep='')
    },
    content = function(file) {
      ggsave(file, plot =   pariplot(), device = "png")
    }
  )  
  
  }
 
)
