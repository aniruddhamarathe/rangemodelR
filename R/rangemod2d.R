#' Range Cohesion Models for Spatial Polygon Grids
#' @description The function takes observed site by species matrix and returns
#'               expected species richness values of each site based on user
#'               defined neighbour relationships.
#' @param spmat a site by species matrix or data frame with species in columns
#' @param shp simple features file of sites where species occurrences are recorded
#' @param reps number of replicates
#' @param nb a neighbour object similar to one generated from
#'        \code{\link[spdep]{poly2nb}} of \pkg{\link[spdep]{spdep}}
#'        If 'nb' is NULL then result is range scatter
#' @param field a number or character vector indicating which column is the unique id
#' @param var an optional vector containing explanatory variable for
#'        constraining the randomization
#' @param first If true, 'var' is used while choosing the first occurrence as
#'        well.if 'var' is null, first is always set 'FALSE'
#' @param degen If true, each randomized site by species matrix is saved and
#'        provided in output
#' @param rsize which range sizes to use for simulation, can be an integer vector
#'        of same length as number of species(columns)
#'        or either 'observed' or 'unif'. See details for explanations
#' @details Implements simulations used by Rahbeck et.al. (2007) to
#'          species distribution data on a continuous grid. In 'spmat' the sites
#'          (rows) represent each cell in the grid.The species occurrences across
#'          sites are randomly spread maintaining strict range cohesion.
#'          A neighbour object is used to limit the choice
#'          of cells during random selections to immediate neighbours. Options
#'          for creating four cell (rook) or eight cell (queen) neighbours can
#'          be accessed while creating the 'nb' object, (typically from package
#'          \pkg{\link[spdep]{spdep}}). The randomization proceeds by selecting a
#'          single site,(weighted by 'var' if provided, and first is TRUE), and
#'          then continues selecting one site at a time from a vector of
#'          available neighbours taken from 'nb' and weighted by 'var' if
#'          provided. The vector of available sites is updated after each site
#'          is selected.
#' @return  A list containing following elements:
#' \describe{
#'  \item{"out.df"}{ a data frame with four columns for mean, standard
#'                   deviation, lower and upper limits of confidence intervals
#'                   of predicted species richness }
#'  \item{"out.shp"}{ same as the input shp with additional four
#'                    columns of 'out.df'}
#'  \item{"degenerate.matrices"}{a list of all the randomized matrices - only
#'                                present if 'degen' is TRUE}
#'  }
#' @references Rahbek, C., Gotelli, N., Colwell, R., Entsminger, G., Rangel,
#'              T. & Graves, G. (2007) Predicting continental-scale patterns of
#'              bird species richness with spatially explicit
#'              models. Proceedings of the Royal Society B: Biological
#'              Sciences, 274, 165.
#'
#'              Gotelli, N.J., Anderson, M.J., Arita, H.T., Chao, A., Colwell,
#'              R.K., Connolly, S.R., Currie, D.J., Dunn, R.R., Graves, G.R. &
#'              Green, J.L. (2009) Patterns and causes of species richness:
#'              a general simulation model for macroecology. Ecology Letters,
#'              12, 873-886.
#' @examples
#'
#' if(requireNamespace('sf')&requireNamespace('spdep')){
#'   data(grid_layer)
#'   data(spmat)
#'
#'   grid_layer <- sf::st_as_sf(grid_layer,sf_column_name = 'geometry')
#'   grid_nb <- spdep::poly2nb(grid_layer)
#'   mod.out <- rangemod2d(spmat,grid_layer,"ID",nb = grid_nb,rsize = "observed", var = NULL,reps = 5)
#'   }
#'\dontrun{
#'  tmap::tm_shape(mod.out$out.shp)+
#'   tm_polygons(fill = mod.rich)
#'  }
#'
#'
#' @export
rangemod2d <- function(spmat,shp,field,nb,rsize = c("observed","unif"),
                       var = NULL, reps,degen = FALSE,first = FALSE){

  ####sanity check of arguments####
  if(any(!(shp$ID%in%rownames(spmat)))){
    stop("all unique ids in 'shp' should appear in 'spmat'")
  }
  #
  # if(!is.null(nb)&!length(nb) == nrow(spmat)){
  #   stop("length of 'nb' should be same as number of sites: ",
  #        length(nb)," and ", nrow(spmat))
  # }

  if(!is.null(var)& !length(var) == nrow(spmat)){
    stop("'var' should be of same length as number of sites: ",
         length(var),"and",nrow(spmat),".")
  }

  if(is.vector(rsize,mode = "numeric")& !length(rsize) == ncol(spmat)){
    stop("rsize should be of same length as number of species: ",
         length(rsize)," and ",ncol(spmat))
  }

  if(is.null(rownames(spmat))){
    warning("No rownames for 'spmat', setting rownames as 1:nrow(spmat)")
    rownames(spmat) <- 1:nrow(spmat)
  }

  ####chunk1 - prepare objects for input and output####
  uid <- unlist(sf::st_drop_geometry(shp[,field]))
  spmat[spmat>0] <- 1
  spmat <- as.matrix(spmat)
  keep <- which(colSums(spmat) > 0)
  spmat <- as.matrix(spmat[,keep])
  if(is.vector(rsize,mode = "numeric")){
    range.size <- rsize[keep]
  }else{
    rsize <- match.arg(rsize,c("observed","unif"))
    range.size <- switch(rsize,observed = {colSums(spmat)},
                         unif = {sample(1:nrow(spmat),ncol(spmat),replace = T)})
  }
  if(!is.null(nb)){
    #add site names to the nb object
    nb <- lapply(nb,function(x)uid[x])
    names(nb) <- uid
  }

  mat.temp <- as.matrix(spmat)
  mat.out <- matrix(nrow = nrow(spmat),ncol = reps,
                    dimnames = list(rownames(spmat),1:reps))

  degen.mats <- list()

  ####chunk2 - use 'random.range' to spread ranges on the matrix####
  for(i in 1:reps){
    mat.temp[which(mat.temp > 0)] <- 0
    for(j in 1:length(range.size)){
      temp.vec1 <- random.range(uid = uid,nb=nb,
                                range.size = range.size[j],var = var,
                                first = first)
      mat.temp[which(rownames(mat.temp)%in%as.character(temp.vec1)),j] <- 1
    }
    mat.out[,i] <- rowSums(mat.temp)
    if(degen == TRUE){degen.mats[[i]] <- mat.temp}
  }
  ####output####
  out.df <- data.frame(mod.rich = apply(mat.out,1,mean),
                       mod.sd =  apply(mat.out,1,stats::sd),
                       q2.5 = apply(mat.out,1,stats::quantile,probs = 0.025),
                       q95.5 = apply(mat.out,1,stats::quantile,probs = 0.975))

  ####chunk3 - add simulated data to shp and prepare for plotting in ggplot2
  shp <- cbind(shp,out.df)

  if(degen==TRUE){
    outlist <- list(out.df = out.df,out.shp = shp,
                    degenerate.matrices=degen.mats)
  }else{
    outlist <- list(out.df = out.df,out.shp = shp)
  }
  outlist

}
