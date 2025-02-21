#' A polygon grid
#'
#' A grid with 100 cells, each \eqn{1^0} by \eqn{1^0} in size. Each cell is assigned a hypothetical local assemblage of varying species richness from a hypothetical regional pool of 20 species, where species may occupy up to 50 cells out of the 100 cells in the grid.
#'
#' \itemize{
#'    \item ID  a unique id field
#'    \item obs_rich  observed species richness
#'    \item trend  a variable with spatial trend obtained using a linear function of the coordinates}
#'
#'@format A data.frame with three data columns and a 'geometry' column containing simple feature geometries
#'@name grid_layer
NULL
