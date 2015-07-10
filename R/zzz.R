
.onAttach = function( libname, pkgname ) {
  dbloc = Sys.getenv("CCLEDB_PATH")
  stopifnot(file.exists(dbloc))
  .ccleCon <<- dbConnect(SQLite(), dbloc)
  .ccleSrc <<- src_sqlite(dbloc)
  message("defined .ccleCon as connection to SQLite db")
  message("defined .ccleSrc as dplyr src, for tbl() processing")
  invisible(NULL)
}
