sleft_join = function(...) suppressMessages(left_join(...))

exByPrimSite = function(sym) {
 stopifnot(sym %in% ccle$exgenes)
 tal = ccle$src %>% tbl("ccle_exprs_tall") %>% filter(Description == sym)
 talr = rename(tal, CCLE_name=Tumor_Sample_Barcode) 
 grp = ccle$src %>% tbl("ccle_cell_line_info") %>% select(CCLE_name, Site_Primary) 
 tmp = sleft_join(talr, grp) %>% as.data.frame()
 sps = split(tmp$Signal, tmp$Site_Primary)
 meds = sapply(sps, median)
 omeds = order(meds, decreasing=TRUE)
 sps[omeds]
}

#exboxByPrimSite = function(sym) {
# dat = exByPrimSite(sym)
# opar = par(no.readonly=TRUE)
# par(mar=c(12,4,2,2), las=2)
# on.exit(par(opar))
# boxplot(dat, ylab=paste0(sym, " expr."))
#}

exboxByPrimSite = function(sym) {
 dat = exByPrimSite(sym)
 opar = par(no.readonly=TRUE)
 par(mar=c(12,4,2,2), las=2, srt=45, font=2)
 on.exit(par(opar))
 sn = names(dat)
 names(dat) = NULL
 boxplot(dat, ylab=paste0(sym, " expr."), labels=NULL, axes=FALSE)
axis(1, labels=FALSE, at=1:length(dat))
text(1:length(dat), par("usr")[3] - 0.25, srt = 45, adj = 1,
     labels = sn, xpd = TRUE)
axis(2)
}


ic50byLineWMut = function(sym, compound) {
 stopifnot(sym %in% ccle$hcgenes)
 comp = ccle$src %>% filter_compound(compound) %>% select(CCLE_Cell_Line_Name, IC50__uM_)
 hyb = ccle$src %>% tbl("ccle_hybrid_capture") %>% filter(Hugo_Symbol==sym)
 hyb = rename(hyb, CCLE_Cell_Line_Name = Tumor_Sample_Barcode)
 sleft_join(hyb, comp)
}
 
pic50byLineWMut = function(sym, compound) {
 tab = na.omit(ic50byLineWMut(sym, compound) %>% as.data.frame())
 o = order(tab$IC)
 val = tab$IC[o]
 nam = tab$CCLE_Cell[o]
 opar = par(no.readonly=TRUE)
 par(mar=c(14,4,2,2), las=2, font=2)
 on.exit(par(opar))
 plot(val, axes=FALSE, xlab=" ", ylab=paste0(compound, " IC50"), log="y")
 axis(2)
 axis(1, at=1:length(val), labels=nam)
 invisible(tab)
}

pic50Serve = function(input, output, session) {
     output$comptitle = renderUI({HTML("<h3>IC50 for selected compound, measured on cell lines with mutation in selected gene</h3>")})
     INIpicsym = reactive({if (is.null(input$insymic)|nchar(input$insymic)==0)
                return("ABL1")
             else input$insymic})
     INIpiccomp = reactive({if (is.null(input$incomp)|nchar(input$incomp)==0)
                return("Irinotecan")
             else input$incomp})
     output$icplot = renderPlot( pic50byLineWMut( INIpicsym(), INIpiccomp() ) )
     updateSelectizeInput(session, "insymic", 
                   choices = ccle$hcgenes, server=TRUE )
     updateSelectizeInput(session, "incomp", 
                   choices = ccle$compounds, server=TRUE )
     }

exbServe = function(input, output, session) {
     output$exsititle = renderUI({HTML("<h3>RMA expression for selected gene, summarized on primary tumor sites</h3>")})
     INI = reactive({if (is.null(input$insym)|nchar(input$insym)==0)
                return("SLFN11")
             else input$insym})
     output$bplot = renderPlot( exboxByPrimSite( INI() ) )
     updateSelectizeInput(session, "insym", 
                   choices = ccle$exgenes, server=TRUE )
     }

exBoxServer = function(input, output, session) {
     exbServe(input, output, session)
     pic50Serve(input, output, session)
     }

multiWidget = function(ccle) {
  shinyApp(ui = fluidPage(
   fluidRow( column(9, htmlOutput("comptitle"))),
   fluidRow(
      column(3, selectizeInput("insymic", label="Mutated Gene",
       choices=NULL, options=list(placeholder="ABL1"),
       selected="ABL1", multiple=FALSE)),
      column(3, selectizeInput("incomp", label="Compound",
       choices=NULL, options=list(placeholder="Irinotecan"),
       selected="Irinotecan", multiple=FALSE))),
   fluidRow(column(9, plotOutput('icplot'))),
   fluidRow( column(9, htmlOutput("exsititle"))),
   fluidRow(column(3, selectizeInput("insym", label="Gene symbol",
       choices=NULL, options=list(placeholder="SLFN11"),
       selected="SLFN11", multiple=FALSE))),
   fluidRow(column(9, plotOutput('bplot')))
   )
  ,
           server = exBoxServer
  )
}

