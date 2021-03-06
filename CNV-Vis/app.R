library(tidyverse)
library(shiny)
library(RJDBC)
library(DT)
#library(org.Hs.eg.db)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
options(scipen = 999)
source("/Users/campbellim/OneDrive - Children's Hospital of Philadelphia/Penn/BMIN502/DBPass.R")

#inner_join(select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = c("ENTREZID","SYMBOL")),
#           select(TxDb.Hsapiens.UCSC.hg19.knownGene,keys = keys(TxDb.Hsapiens.UCSC.hg19.knownGene),columns=c('GENEID', 'TXCHROM', 'TXSTART', 'TXEND', 'TXID'), keytype="GENEID"),
#           by = c("ENTREZID" = "GENEID")) %>%
#    mutate(TXCHROM = str_remove(TXCHROM,"chr")) %>%
#    group_by(SYMBOL) %>%
#    summarise(Chr = TXCHROM[1], 
#              Start = min(TXSTART,na.rm = TRUE),
#              Stop = max(TXEND,na.rm = TRUE)) -> Genes
load("/Users/campbellim/OneDrive - Children's Hospital of Philadelphia/Penn/BMIN502/GeneSymbols.hg19.RData") #Genes Table

jdbcDriver <- JDBC("oracle.jdbc.OracleDriver", classPath = "/Users/campbellim/Drivers/OJDBC-Full (v12_1_0_2)/ojdbc6.jar")
jdbcConnection <- dbConnect(jdbcDriver, "jdbc:oracle:thin:@whdb.research.chop.edu:1521:WHDB","cnv",dbpass)

Arrays <- dbGetQuery(jdbcConnection, "SELECT ARRAYID, ARRAYNAME FROM ARRAY")

CNVs <- dbGetQuery(jdbcConnection, "SELECT * FROM CNV")
CNVs %>%
    mutate(CHROMOSOME = str_remove(CHROMOSOME,"chr")) -> CNVs

ui <- fluidPage(

    titlePanel("CNV Visualization"),

    sidebarLayout(
        sidebarPanel(
            selectInput(inputId="Array", label = "Array", choices = setNames(Arrays$ARRAYID,Arrays$ARRAYNAME), selected = 1),
            selectInput(inputId="Chr", label = "Chromosome",choices = c(1:21,"X","Y","MT"), selected = 1),
            numericInput(inputId = "Start", label = "Start", value = 2700000),
            numericInput(inputId = "Stop", label = "Stop", value = 3000000),
            selectizeInput(inputId = "Gene", label = "Jump to Gene", choices = c("",Genes$SYMBOL[1:10])),
            actionButton("Left","<<"),actionButton("Out3","Out 3x"),actionButton("In3","In 3x"),actionButton("Right",">>")
        ),

        mainPanel(
           plotOutput("Plot"),
           h3("CNV Calls"),
           DT::dataTableOutput("Plot2")
        )
    )
)

server <- function(input, output, session) {


    
    query <- "SELECT A.ARRAYID, A.ARRAYNAME, PI.CHROMOSOME, PI.\"POSITION\", P.BAF, P.LRR 
    FROM ARRAY A
    LEFT JOIN PROBE_INFO PI
    ON A.VERSIONID = PI.VERSIONID 
    LEFT JOIN PROBE P
    ON PI.PROBEID = P.PROBEID 
    AND A.ARRAYID = P.ARRAYID 
    WHERE A.ARRAYID = %s
    AND PI.CHROMOSOME = '%s'
    AND PI.POSITION > %s
    AND PI.POSITION < %s"
    
    data <- reactiveVal() 
    CNVTable <- reactiveVal() 
    observe({
        complete.query <- sprintf(query,input$Array,input$Chr,input$Start,input$Stop)
        query.result <- dbGetQuery(jdbcConnection,complete.query)
        data(query.result %>%
            pivot_longer(BAF:LRR,names_to = "Variable", values_to = "Value") %>%
            mutate(Variable = factor(Variable, levels = c("LRR","BAF"))))
        
    })
    
    observeEvent(input$Gene,{
        if(input$Gene != ""){
        updateSelectInput(session, "Chr", selected = Genes[Genes$SYMBOL == input$Gene,][["Chr"]])
        updateNumericInput(session, "Start", value = Genes[Genes$SYMBOL == input$Gene,][["Start"]])
        updateNumericInput(session, "Stop", value = Genes[Genes$SYMBOL == input$Gene,][["Stop"]])
        updateSelectizeInput(session, "Gene", selected = "")
        }
    })
    
    observeEvent(input$Out3,{
            updateNumericInput(session, "Start", value = max(0,input$Start - (input$Stop-input$Start)))
            updateNumericInput(session, "Stop", value = input$Stop + (input$Stop-input$Start))
        })
    
    observeEvent(input$In3,{
        updateNumericInput(session, "Start", value = max(0,input$Start + floor((input$Stop-input$Start) * 0.33)))
        updateNumericInput(session, "Stop", value = input$Stop - floor((input$Stop-input$Start) * 0.33))
    })
    
    observeEvent(input$Left,{
        updateNumericInput(session, "Start", value = max(0,input$Start - floor((input$Stop-input$Start) * 0.5)))
        updateNumericInput(session, "Stop", value = input$Stop - floor((input$Stop-input$Start) * 0.5))
        if(CNVTable$CNVStart > input$Stop | CNVTable$CNVStart < input$Start) CNVTable$CNVStart <- NULL
        if(CNVTable$CNVStop > input$Stop | CNVTable$CNVStop < input$Start) CNVTable$CNVStop <- NULL
        })
    
    observeEvent(input$Right,{
        updateNumericInput(session, "Start", value = max(0,input$Start + floor((input$Stop-input$Start) * 0.5)))
        updateNumericInput(session, "Stop", value = input$Stop + floor((input$Stop-input$Start) * 0.5))
        if(CNVTable$CNVStart > input$Stop | CNVTable$CNVStart < input$Start) CNVTable$CNVStart <- NULL
        if(CNVTable$CNVStop > input$Stop | CNVTable$CNVStop < input$Start) CNVTable$CNVStop <- NULL
    })    
    
    observeEvent(input$select_button,{
        CNVTable$CNVStart <- CNVs[CNVs$CNVID == str_extract(input$select_button,"[0-9]+$"),"START"]
        CNVTable$CNVStop <- CNVs[CNVs$CNVID == str_extract(input$select_button,"[0-9]+$"),"STOP"]
        CNVTable$CNVLen <- CNVTable$CNVStop - CNVTable$CNVStart
        updateNumericInput(session, "Start", value = max(0, CNVTable$CNVStart - (CNVTable$CNVLen * 2)))
        updateNumericInput(session, "Stop", value = CNVTable$CNVStop + (CNVTable$CNVLen * 2))
        updateSelectInput(session, "Chr", selected = CNVs[CNVs$CNVID == str_extract(input$select_button,"[0-9]+$"),"CHROMOSOME"])
    })    
    
    output$Plot <- renderPlot({
    if(!is.null(nrow(data()))){
        p <- ggplot(data(), aes(x = POSITION, y = Value)) + 
            geom_point() + 
            facet_wrap(~Variable,ncol = 1, scale = "free_y", strip.position = "left",labeller = as_labeller(c("LRR" = "Log R Ratio","BAF" = "B Alelle Frquency"))) +
            theme_linedraw() + labs(y = NULL, x = paste0("Chromosome ",data()$CHROMOSOME[1]," Genomic Coordinate (hg19)")) +
            theme(strip.placement = "outside", strip.background = element_blank(),
                  strip.text = element_text(color = "black", size = 12))
        if(length(CNVTable$CNVStart) > 0 & length(CNVTable$CNVStop) > 0){
           p + geom_vline(aes(xintercept = CNVTable$CNVStart), lty = 3) + geom_vline(aes(xintercept = CNVTable$CNVStop), lty = 3)
        } else {p}
    }
    })
    
    CNVTable <- reactiveValues()
    observe({CNVTable$table = data.frame(
        CNVs[CNVs$ARRAYID == input$Array,c("CHROMOSOME","START","STOP","LENGTH","DEEPCNVPROB","CN")],
        Actions = sapply(CNVs[CNVs$ARRAYID == input$Array,"CNVID"],function(x)as.character(actionButton(inputId = paste0("button_",x), label = "View", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' )))
#        CNVTable$table = CNVTable$table[order(CNVTable$table$CHROMOSOME),]
        )})
    
    output$Plot2 <- DT::renderDataTable({CNVTable$table},
                                        rownames = FALSE, escape = FALSE, selection = 'none')
    
}

# Run the application 
shinyApp(ui = ui, server = server)
