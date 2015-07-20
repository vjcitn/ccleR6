get_affy <- function(ccle) {
  
  #get affy data for genes and cell lines of interest
  affydata <- ccle$src %>% tbl("ccle_exprs_tall") %>% 
    filter(Description %in% ccle$query_symbols & Tumor_Sample_Barcode %in% ccle$query_cell_lines ) %>%
    collect()
  
  #select and rename columns
  affydata <- affydata %>% transmute(CCLE_name = Tumor_Sample_Barcode,
                                          ID = Description,
                                          Original = as.character(Signal),
                                          Value = Signal,
                                          Type = 'affy')
  
  #check for duplicates
  dup_check <- affydata %>% group_by(ID) %>% summarise(N=n()) %>% arrange(desc(N)) %>% filter(N == median(N)) %>% as.data.frame(n=-1)
  
  #get rid of duplicates and scale values
  affydata <- affydata %>% filter(ID %in% dup_check$ID) %>% group_by(ID) %>% 
    mutate(zscore=scale(Value)) %>% ungroup() 
  return(affydata)

}


get_hybcap <- function (ccle) {
  
  #get hybrid capture sequencing data for genes and cell lines of interest
  hybcapdata <- ccle$src %>% tbl("ccle_hybrid_capture") %>% 
    filter(Hugo_symbol %in% ccle$query_symbols & Tumor_Sample_Barcode %in% ccle$query_cell_lines ) %>%
    collect()
  
  #filter rows that have a protein change and select and rename columns
  hybcapdata <- hybcapdata %>% dplyr::filter(nchar(Protein_Change) > 1) %>%
    transmute(CCLE_name = Tumor_Sample_Barcode,
              ID = Hugo_Symbol,
              Original = Protein_Change)
  
  #ensure that only one row per gene/cell line and add value column
  hybcapdata <- hybcapdata %>% group_by(CCLE_name, ID) %>% summarise(Original=paste(Original, collapse=',')) %>% mutate(Value=1)
  
  #which samples were actually sequenced?
  hybcaptested <- ccle$src %>% tbl("ccle_hybrid_capture") %>% select(Tumor_Sample_Barcode) %>%
    distinct() %>% collect() %>% filter(Tumor_Sample_Barcode %in% ccle$query_cell_lines)
  hybcaptested <- hybcaptested %>% filter(Tumor_Sample_Barcode %in% ccle$query_cell_lines)
  sequenced_ids <- hybcaptested$Tumor_Sample_Barcode
  notsequenced_ids <- setdiff(ccle$query_cell_lines, sequenced_ids)
  
  #generate dataframes for sequenced and not sequenced cell lines
  hybcapsequenced <- data.frame ( CCLE_name=rep(sequenced_ids, length(ccle$query_symbols)) , 
                                   ID=rep(ccle$query_symbols, each= length(sequenced_ids) ),
                                   Original='-',
                                   Value=0, stringsAsFactors=FALSE)
  
  if (length(notsequenced_ids) > 0 ) {
    
    hybcapnotsequenced <- data.frame ( CCLE_name=rep(notsequenced_ids, length(ccle$query_symbols)) , 
                                        ID=rep(ccle$query_symbols, each= length(notsequenced_ids) ),
                                        Original=NA,
                                        Value=NA, stringsAsFactors=FALSE )
    
  } else {
    hybcapnotsequenced <- data.frame ()
  }
  
  #get rid of rows in hybcap.sequenced which are duplicated in hybcap.data.hm
  hybcapsequenced <- hybcapsequenced %>% filter(!( paste(CCLE_name, ID) %in% paste(hybcapdata$CCLE_name, hybcapdata$ID) ))
  
  #combine hybcap dataframes and add additional standard columns
  hybcapdata <- bind_rows(hybcapdata, hybcapsequenced, hybcapnotsequenced) %>%
    mutate(Type='hybcap', zscore=Value)
  
  return ( hybcapdata )
  
}

#get compound response data
get_response <- function (ccle) {
  
  respdata <- ccle$src %>% tbl("ccle_drug_data") %>% 
    filter(Compound %in% ccle$query_compounds & CCLE_Cell_Line_Name %in% ccle$query_cell_lines ) %>%
    collect()
    
  #just get fields of interest
  respdata <- respdata %>% transmute(CCLE_name=CCLE_Cell_Line_Name, 
                                     ID=Compound, 
                                     Original=as.character(ActArea), 
                                     Value=ActArea,
                                     Type='resp')
  
  #scale
  respdata <- respdata %>% group_by(ID) %>% 
    mutate(zscore=scale(Value)) %>% ungroup() 
  
  return ( respdata )
  
}

make_df <- function (ccle, datatype=c('affy', 'hybcap'))  {
  
  require(tidyr)
  
  #make a data container list
  all_data.list <- list()
  if ('affy' %in% datatype) { all_data.list[['affy_data']]  <- get_affy (ccle) } 
  if ('hybcap' %in% datatype) { all_data.list[['hybcap_data']] <- get_hybcap (ccle) }
  all_data.list[['resp_data']] <- get_response(ccle)
  
  #combine
  all_data <- bind_rows(all_data.list)
  
  #get rid of cell lines with no response data
  present_cls <- all_data.list[['resp_data']] %>% select(CCLE_name) %>% distinct()
  missing_cls <- setdiff(ccle$query_cell_lines, present_cls$CCLE_name)
  all_data <- all_data %>% filter(CCLE_name %in% present_cls$CCLE_name)
  if (length(missing_cls) != 0) {
    warning(sprintf('No response data for following cell lines: %s', paste(missing_cls, collapse=', ')))
    
  }
  
  #create matrix
  output.df <- all_data %>% transmute(CCLE_name, fn=paste(ID,Type,sep='_'), Value) %>% spread(fn, Value)
  
  #reorder columns
  output.df <- output.df %>% dplyr::select(CCLE_name, ends_with('_resp'), everything())
  
  #make mutation fields into factors
  output.df <- output.df %>% mutate_each(funs(as.factor), ends_with('_hybcap|_cosmic'))
  
  return(output.df)
  
}


make_heatmap <- function (ccle, datatype=c('affy', 'hybcap'), compound=1) {
  
  require(reshape2)
  require(ggplot2)
  require(scales)
  
  #make a data container list
  all_data.list <- list()
  if ('affy' %in% datatype) { all_data.list[['affy_data']]  <- get_affy (ccle) } 
  if ('hybcap' %in% datatype) { all_data.list[['hybcap_data']] <- get_hybcap (ccle) }
  all_data.list[['resp_data']] <- get_response(ccle)
  
  #combine
  all_data <- bind_rows(all_data.list)
  
  #get compound to order by
  if (is.numeric(compound)) {
    if (compound <= length(ccle$query_compounds)) { compound_name <- ccle$query_compounds[compound] }
    else { stop(sprintf('There are only %s compounds but you selected compound no %s', length(ccle$query_compounds), compound))}
  } else {
    if (compound %in% ccle$query_compounds) {compound_name <- compound} 
    else {stop(sprintf('Compound %s needs to be specified in ccle object', compound))}
  }
  
  #get rid of unneeded response data
  plotdata <- all_data %>% filter( !(all_data$Type == 'resp' & all_data$ID != compound_name)  )
  
  #featurename order
  datatype <- c('resp', datatype)
  plotdata <- plotdata %>% mutate(FeatureName = paste(Type, ID, sep='_'))
                                  
  ordered_feature_names <- plotdata %>% transmute(FeatureName, Type = factor(Type, levels=datatype), ID) %>%
    distinct() %>% arrange(ID, Type)
  
  
  
  #cell lines order
  ordered_cell_lines <- all_data.list[['resp_data']] %>% filter(ID==compound_name) %>% arrange(desc(Value))

  #make new scaleinfo
  scaleinfo <- plotdata %>% group_by(Type) %>% 
    summarise(min=min(zscore, na.rm=TRUE), mean=mean(zscore, na.rm=TRUE), max=max(zscore, na.rm=TRUE)) %>%
    melt(id.vars='Type') %>% transmute(Type=as.character(Type), Level=as.character(variable), Value=value)
  
  #scale coloring
  scalecols <- data.frame(Type=c('resp', 'resp', 'resp', 'affy', 'affy', 'affy', 'hybcap', 'hybcap'),
                          Level=c('min', 'mean', 'max', 'min', 'mean', 'max', 'min', 'max' ),
                          Colour=c('red', 'yellow', 'green', 'blue', 'white', 'red', 'lightblue', 'darkblue'),
                          stringsAsFactors = FALSE)
  
  #make offset info
  offsetinfo <- data.frame(Type=c('resp', 'affy', 'hybcap'),
                           Offset=c(0,20,40), stringsAsFactors = FALSE)
  
  
  #merge scale info colours
  scaleinfo <- scaleinfo %>% inner_join(scalecols, by=c('Type', 'Level')) %>% inner_join(offsetinfo, by='Type') %>% 
    mutate(Value=Value+Offset) %>% arrange(Value)
  
  #apply offset info to plotdata
  plotdata <- plotdata %>% inner_join(offsetinfo, by='Type') %>% mutate(zscore_offset = zscore+Offset)
  
  #do the plot
  p <- ggplot(plotdata, aes(x=CCLE_name, y=FeatureName))  +
    geom_tile(aes(fill = zscore_offset), linetype=0 ) + 
    scale_fill_gradientn(colours = scaleinfo$Colour, values = rescale(scaleinfo$Value)) +
    scale_y_discrete(limits=ordered_feature_names$FeatureName) + #this orders the y axis as we want it
    scale_x_discrete(limits=ordered_cell_lines$CCLE_name) + #this orders the x axis as we want it
    coord_flip() +
    theme_bw(base_size = 9) + 
    labs(x = "", y = "") +
    theme(axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1), axis.text.y = element_text(size=rel(1))) +
    theme(panel.background=element_rect(fill="gray90", color='white'), legend.position = "none") 
  p
  
  return(p)  
  
}
