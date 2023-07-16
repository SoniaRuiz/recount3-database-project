## Load the necessary tables from the database
con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
df_novel <- dplyr::tbl(con, "novel") %>% dplyr::collect()
df_intron <- dplyr::tbl(con, "intron") %>% dplyr::collect()

#' Measure the distance between a novel junction and their reference annotated
#' intron
#'
#' @param novel_start Start position of the novel junction.
#' @param novel_end End position of the novel junction.
#' @param ref_start Start position of the reference annotated intron.
#' @param ref_end End position of the reference annotated intron.
#' @param type The category of the novel junction. Either novel_acceptor or
#'   novel_donor.
#' @param strand The strand of the two junctions.
#'
#' @return The distance between the novel junction and its reference annotated
#'   intron.
#' @export
measureDistance <- function(novel_start, novel_end, ref_start, ref_end, type, strand){
  if(type == "novel_donor" & strand == "+"){
    return(ref_start - novel_start)
  }else if(type == "novel_donor" & strand == "-"){
    return(novel_end - ref_end)
  }else if(type == "novel_acceptor" & strand == "+"){
    return(novel_end - ref_end)
  }else if(type == "novel_acceptor" & strand == "-"){
    return(ref_start - novel_start)
  }else{
    return(NA)
  }
}
