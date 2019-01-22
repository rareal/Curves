library(shiny);library(dplyr);library(magrittr);library(ggplot2);library(purrr);library(tidyr);library(tibble)
# Functions for 4PL logistic curve fitting  ========================================
# y = ((min-max)/(1+(x/EC50)^slope))+max
# x = EC50*(((min-max)/(y-max))-1)^(1/slope)

# get y
fx<-function(parm,x){
  (parm[1]-parm[4])/(1+(x/parm[3])^parm[2])+parm[4]}

# get x
xf<-function(parm,y){
  parm[3]*(((parm[1]-parm[4])/(y-parm[4]))-1)^(1/parm[2])}

# get slope
sl<-function(minY,maxY,e,x,y) {
  log((minY-maxY)/(y-maxY),x/e)  }

# sum of squares
ssq<-function(parm,x,y){
  sum((y-fx(parm,x))^2)}


pred2<-function(parms,x){
  data.frame(x=x,yp=fx(parms,x))
}


nlm2<-function(OD){
  start<-c(max(OD),3,4.5,min(OD)) # what I thought it was a good start
  #start<-c(max(OD),1,ifelse(OD[2]<0.5,max(OD)/2,4),min(OD)) # excel 4PL start
  res<-nlm(ssq,start,1:8,OD)
  c(res$estimate,res$minimum)
}

#===============================================================

# setwd("C:/Users/RBA/OneDrive/FL/R/LEARN/Shiny/curves2")
#runApp()

ui <- fluidPage(br(),
tabsetPanel(
  tabPanel("Setup",
           sidebarLayout(
            sidebarPanel(width=3,
                    radioButtons("type", "Type of Analysis",c("Concentration"="conc","Titer"="titer")),
                    fileInput("file1","Upload file"),
                    actionButton("exampledata","or use an example dataset", #"./example_data/Plate1.txt"
                                 style="background-color:#b2d7ea"),
                    hr(),
                    numericInput("blankpos", "Blank Position", 11, min = 1, max = 12),
                    numericInput("STpos", "Standard position", 12, min = 1, max = 12),
                    uiOutput("inidilv"),
                    numericInput("dilFactor", "Dilution Factor", 3),
                    conditionalPanel(
                      condition = "input.type == 'conc'",
                      textInput("STconc","Standard Concentration",2),
                      textInput("STunit","unit","ug/ml")),
                    conditionalPanel(
                      condition = "input.type == 'titer'"
                      ,numericInput("mincut", "Titer cutoff", 0.1)
                      ,numericInput("refnorm", "Reference Normalization", 1)
                      )
                    ),
            mainPanel(
              fluidRow(column(10,h1("Instructions"),
                              h4("Important!! Upload only the data, without header or row names",
                                 style="color:#7f0c18;font-style:italic")),
                       column(2,img(src = "curves_logo.jpg",height = 100, width = 269,align="right"))),
                      
                      hr(),
                      tags$b("The data will be read only if the following conditions are met:"),
                      tags$ul(
                        tags$li("	It has 8 or 16 rows"),
                        tags$li("	It has 12 or 24 columns"),
                        tags$li("	rows X columns = 96 or 384")),
                      tags$b("This app will analyse ELISA data of the following format:"),
                br(),br(),
                      img(src = "read_format.jpg",height = 145.6, width = 500)
                ,br(),br()
                      ,tags$b("It can supports data in 96 well or 384 well format")
                ,br(),tags$b("In case of 384 well format, the data will be broken into 4 x 96 well plates, like this:")
                ,br(),br(),img(src = "384.jpg",height = 209, width = 500)
              ,br(),br(),tags$b("The program will then:"),
              tags$ul(
                tags$li("	remove the average of the blanks from each value in the dataset;"),
                tags$li("	normalize all the values to the first value of the Standard (ST) per plate;"),
                tags$li("	fit a 4 parameter logistic model to each sample."))
              ,tags$b("To evaluate the fit, use the following example:")
              ,br(),img(src = "example.png",height = 285.7, width = 300)
              ,tags$div(
p("After the fit, the program will decide a Y cutoff (closest to 0.5 as possible) for each sample.")
,p("The dilution at Y = cutoff will be used to quantify the ratio sample/ST")
,p("Concentration:")
# ,tags$ul(tags$li("The final concentration of the sample will be:
# initial dilution X ratio(sample/ST) X ST concentration",style="font-weight:normal"))
,tags$ul(tags$li(withMathJax("The final concentration of the sample will be:
$$\\text{initial dilution} * 
               \\text{ratio(sample/ST)} * \\text{ST concentration}$$"),style="font-weight:normal"))
,p("Titer:")
# ,tags$ul(tags$li("The final titer of the sample will be: 
# initial dilution X ratio(sample/ST) X ST titer (defined using Titer cutoff)",style="font-weight:normal")
,tags$ul(tags$li(withMathJax("The final titer of the sample will be: 
$$log_{10}(\\text{initial dilution} * 
          \\text{ratio(sample/ST)} * \\text{ST titer})$$"),style="font-weight:normal"))
,div("ST titer is the ST dilution at Titer cutoff",style="font-weight:normal;font-style: italic",align="center")
#runApp()
,tags$ul(
tags$li("If sample maximum OD < Titer cutoff, sample Titer = 0",style="font-weight:normal")  
,tags$li("Reference Normalization will be used to normalize ODs to another reference. 
All ODs (including the current reference) will be multiplied by this number"),style="font-weight:normal")
,style="font-weight:bold")
                      )
           )
           ),
  tabPanel("Curves",plotOutput("diagnostic")),
  tabPanel("Results",
           fluidRow(
             column(2,tableOutput("resTable")),
             column(2,
                    br(),
                    actionButton("copyButton", "Copy!"))
           ))
)
  
)

### SERVER =========================================================================================

server <- function(input, output,session) {
session$onSessionEnded(stopApp)

# dinamic UI
output$inidilv<-renderUI({
  numericInput("iniDil", "Initial Dilution", ifelse(input$type=="conc",10,1))
})

#data input
  # plate1<-reactive({
  #   inFile <- input$file1
  #   
  #   if (is.null(inFile))
  #     return(NULL)
  #   read.table(inFile$datapath, header = F) 
  # })

#if (is.null(input$file1)&input$button==0) plate1<-NULL  # test!!!!!!!!!!

observeEvent(input$file1,{
  plate1 <- read.table(input$file1$datapath,header=F) })

observeEvent(input$button,{
  plate1 <- read.table("./example_data/Plate1.txt",head=F) })

####verbatimTextOutput


blankpos<-reactive(input$blankpos)
STpos<-reactive(input$STpos)
iniDil<-reactive(input$iniDil)
dilFactor<-reactive(input$dilFactor)
STconc<-reactive(input$STconc)
STunit<-reactive(input$STunit)

sep_plate_unnested <- reactive({
  if (is.null(plate1())) return(NULL)
  # separate plates, remove blank per row and normalize to ST
  sep_plate<-as.data.frame(matrix(nrow=32,ncol=11))
  sep_plate<-cbind(plate=rep(paste0("plate1.",1:4),each=8),sep_plate)
  STpos2=ifelse(STpos()<blankpos(),STpos(),STpos()-1)
  #i=1
  for(i in 1:4){
    srows=seq(ifelse(i<=2,1,9),ifelse(i<=2,8,16))
    scols=seq(ifelse(i%%2==0,2,1),24,by=2)
    #remove blanks
    int<-plate1()[srows,scols][,-blankpos()]-mean(plate1()[srows,scols][,blankpos()])
    int[int<0]=0
    #normalize to ST
    int<-int/int[1,STpos2]
    sep_plate[(i-1)*8+(1:8),2:12]<-int
  }
  
  # gather samples and nest around plates
  sep_plate_nested<-sep_plate %>% gather(col,OD,-plate) %>% group_by(plate) %>%  nest(.key=by_plate)
  
  # nest again around samples (the operator %<>% assigns the result of the right on the left)
  sep_plate_nested %<>%
    mutate(by_plate = map(by_plate, ~.x %>% group_by(col) %>% nest(.key = ODg)))
  
  #model and predictions, get cut, xcut, ratioToST, and concentration (conc)
  for(i in 1:4){
    sep_plate_nested$by_plate[[i]]$col[-STpos2]<-paste0('s',(i-1)*10+(1:10)) #name samples
    sep_plate_nested$by_plate[[i]]$col[STpos2]<-paste0('ST.',i)  #name STs
    sep_plate_nested$by_plate[[i]] %<>%
      mutate(parms=map(ODg,nlm2),   #map nlm to ODs
             pred=map(parms,pred2,seq(1,8,length.out=32)), #map predictions using parms
             cut=map(ODg,~{ifelse(min(.)>0.5,min(.),ifelse(max(.)<0.5,max(.),0.5))})) #get cut (y value to compare to ST)
    ratioToST=rep(NA,11);ratioToST[STpos2]<-1  # vector to fill inside the loop
    conc=rep(NA,11);conc[STpos2]<-STconc()       # vector to fill inside the loop
    xcut=rep(NA,11)                            # vector to fill inside the loop
    for(j in c(1:11)[-STpos2]){
      cut.=sep_plate_nested$by_plate[[i]][j,][['cut']] %>% unlist
      sDilCut=sep_plate_nested$by_plate[[i]][j,] %$% unlist(parms) %>% xf(cut.) #sample dil at cut
      xcut[j]<-sDilCut
      StDilCut=sep_plate_nested$by_plate[[i]][STpos2,] %$% unlist(parms) %>% xf(cut.) #ST dil at cut
      ratioToST[j]<-dilFactor()^(sDilCut-StDilCut)  #e.g. c^a / c^b = c^(a-b)
      conc[j]<-iniDil()*ratioToST[j]*STconc()
    }
    sep_plate_nested$by_plate[[i]] %<>% 
      mutate(ratioToST=ratioToST,conc=conc,xcut=xcut)
  }
  
  #unnest by_plate and order col factor
  sep_plate_nested %>% unnest -> sep_plate_unnested
  sep_plate_unnested$col<-factor(sep_plate_unnested$col,levels=sep_plate_unnested$col,ordered=T)
  
  sep_plate_unnested  
  
})

output$diagnostic <- renderPlot({
  if (is.null(plate1())) return(NULL)  
  #get the predictions 
  sep_plate_unnested() %>% unnest(pred) -> t1d 
  
  #get the optimization minimum for each sample
  sep_plate_unnested() %>% unnest(parms) %>% 
    slice(seq(5,nrow(.),by=5)) -> t1m 
  
  #selected cut values
  sep_plate_unnested() %>% select(plate,col,cut,xcut) %>% unnest(cut,xcut) ->selV
  
  #plot
  sep_plate_unnested() %>% unnest(ODg) %>% cbind(x=rep(1:8,nrow(.))) %>% 
    #mutate(col=factor(col,levels=unique(col),ordered=T)) %>%
    ggplot(data=.,aes(x,OD,group=col))+geom_point()+
    facet_wrap(~col,ncol=11)+
    geom_line(data=t1d,aes(x,yp,group=col))+
    geom_label(data=t1m,aes(x=Inf,y=Inf,label=format(parms,digits=1,scientific=T)),
               fill=ifelse(t1m$parms>0.01,"red","white"),hjust=1,vjust=1,size=3)+
    scale_x_continuous("Dilution",breaks=1:8)+
    labs(y="normalized OD")+theme_bw()+
    theme(panel.grid.minor = element_blank(),
          strip.background =element_rect(fill="#deebf7"),
          strip.text=element_text(face='bold',size=12))+
    #  geom_point(data=selV,aes(xcut,cut,group=col),pch=23,fill="orange",size=3)
    geom_hline(data=selV,aes(yintercept=cut),size=rep(c(rep(0.5,10),NA),4),colour="#feb24c")
  
})

resTable1<-reactive({
  if (is.null(plate1())) return(NULL)
  STpos2=ifelse(STpos()<blankpos(),STpos(),STpos()-1)
  sep_plate_unnested() %>% select(col,conc) %>% 
    slice(-seq(STpos2,44,by=11)) %>% #print(n=nrow(.))
    as.data.frame %>% mutate(conc=round(conc,1)) -> tableRes
  colnames(tableRes)<-c('samples',STunit())
  tableRes
})

output$resTable <- renderTable({resTable1()}) 

observeEvent(input$copyButton, {
  write_clip(resTable1())
})

    
}


shinyApp(ui, server)

