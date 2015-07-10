
ccledb = R6Class("ccledb",
   public = list(
    cell_lines = NA,
    compounds = NA,
    organs = NA,
    cngenes = NA,
    hcgenes = NA,
    src = NA,
    trimvec = function(x) {
       tmp = paste(x, collapse=" ")
       if (nchar(tmp)>50) tmp = paste0(substr(tmp, 1, 50), " ...")
       tmp
       },
    print = function(...) {
      cat("ccledb instance. Components:\n")
      flds = setdiff(ls(self), 
           c("clone", "initialize", "print", "trimvec", "src"))
      for (i in flds) {
         vals <- do.call("$", list(self, i))
         cat("  ", i, ": (", length(vals), ") ", self$trimvec(vals), "\n",
                  sep="")
         }
      cat("---\n")
      cat("use $src for reference to dplyr src_sqlite:\n")
      print(self$src)
      },
    initialize = function(src) {
      stopifnot(inherits(src, "src_sqlite"))
      self$src <- src
      message("creating guide vectors (need ~10 seconds...)")
      clinfo <- self$src %>% tbl("ccle_cell_line_info")
      self$cell_lines <- (clinfo %>% 
                   select(CCLE_name) %>% 
                      as.data.frame(n=-1))$CCLE_name
      self$organs <- sort(unique((clinfo %>% 
                   select(Site_Primary) %>% 
                      as.data.frame(n=-1))$Site_Primary))
      self$compounds <- (self$src %>% 
              tbl("ccle_drug_data") %>% 
                   select(Compound) 
                      %>% distinct() %>% as.data.frame())$Compound
      self$cngenes = sort(unique((self$src %>%
              tbl("ccle_copynumber_tall") %>%  # SLOW!!!
                   select(geneName) %>% as.data.frame(n=-1))$geneName))
      self$hcgenes = sort(unique((self$src %>%
              tbl("ccle_hybrid_capture") %>% 
                   select(Hugo_Symbol) %>% as.data.frame(n=-1))$Hugo_Symbol))
      message("done.")
      }
    )
)
