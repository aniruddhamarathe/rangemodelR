#' random.range - supporting function for other rangemodel fucntions
#'
#' @description It is used within other functions in this package to
#'              randomly place given number of species occurrences
#' @param uid a vector of unique ids for selection
#' @param nb a neighbour object similar to generated from 'poly2nb'
#' @param range.size a vector of number of sites occupied by each species
#' @param var an optional vector of variables for constraining the randomization
#' @param first If true, var is used while choosing the first occurrence as well.
#'        if var is null, first is always set FALSE
#' @details This function is not intended for any direct use but is called
#'           within other functions of this package.
#' @return a numeric vector specifying selected values in 'uid'
#' @export

random.range <- function(uid,nb,range.size,var,first){
  if(is.null(nb)){
    sel.vec <- sample(uid,range.size,prob = var)
  }else{
    sel.vec <- NULL
    sel.nb <- NULL # objects to store selected cells and neighbours

    if(first == T){
      sel.vec <- c(sel.vec,sample(uid,1,prob = var)) # select first cell

    }else{sel.vec <- c(sel.vec,sample(uid,1))} # select first cell)

    if(range.size == 1){
      return(sel.vec)
    }else{
      sel.nb <- nb[names(nb)%in%sel.vec[1]] # query the nb object
      # with name of the selected cell in the uid vector
      sel.nb.vec <- unlist(sel.nb)
      if(is.null(var)){
        for(i in 2:range.size){
          sel.vec[i] <- sample(unique(sel.nb.vec),1)
          sel.nb <- nb[names(nb)%in%sel.vec]
          sel.nb.vec <- unlist(sel.nb)
          sel.nb.vec <- sel.nb.vec[!sel.nb.vec%in%sel.vec]
        }  #extend the selction to desired length
      }else{
        for(i in 2:range.size){
          sel.vec[i] <- sample(unique(sel.nb.vec),1,
                               prob = var[uid %in% unique(sel.nb.vec)] )
          sel.nb <- nb[names(nb)%in%sel.vec]
          sel.nb.vec <- unlist(sel.nb)
          sel.nb.vec <- sel.nb.vec[names(nb)%in%sel.vec]
        }  #extend the selction to desired length
      }
    }
  }

  return(sel.vec)
}
