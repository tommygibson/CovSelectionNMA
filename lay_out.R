#' Arranging ggplots
#' 
#' Provides a \code{layout}-like interface for arranging ggplots of different 
#' sizes.
#' 
#' @param ... Each argument should be of the form \code{list(plot, rows, 
#' columns)}, where \code{plot} is a ggplot (or similar), and \code{rows} and 
#' \code{columns} are consecutive sequences indicating the row and column 
#' numbers for \code{plot} to span.
#' 
#' @author Alan D. Jassby and James E. Cloern (originally from the \code{wq} 
#' package).
#' 
#' @examples
#' \dontrun{
#' gg <- ggplot(mtcars, aes(x = hp, y = mpg)) + geom_point()
#' layOut(list(gg, 1:2, 1:3),
#'        list(gg, 3, 1:2),
#'        list(gg, 3, 3))
#' }
#' 
#' @export
lay_out <- function(...) {
  
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]],
          vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                              layout.pos.col = x[[i]][[3]]))
  }
}