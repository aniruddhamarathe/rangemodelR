#' Range shuffle models for for range extents reccorded along gradients
#' @description This function was used by Wang et. al.(2012) to test geometric
#'              constraints, on elevational gradients. The function randomizes
#'              the range extents and range location for each species and
#'              returns expected pattern in species richness under geometric
#'              constraints.
#' @param x     Input data for elevational extents of species.data.dataframe
#'              with names 'genus_species','min','max','range','mid',
#'              'num_zones'. See Wang et.al (2012) for details. Species names
#'              column is optional.
#' @param boundary nature of boundaries at the extremes of the gradien ie.
#'                 either 'hard boundaries' that species cannot cross, or
#'                 'soft boundaries' that species can move across. Can be
#'                one of the following choises "hh", "sh", "hs", or "ss"
#' @param var Predictor variable for constraining randomizations.
#'            dataframe with columns 'mid' and 'weights'. Where 'wights' provide
#'            relative chance of selecting a range location, typically based on
#'            environmental predictors
#' @param interval Numeric. Interval between
#' @param sites Numeric. Locations on domain for calculating species richness
#' @param reps  number of iterations
#' @param degen logical. If TRUE save each randomized distribution.
#' @param lowest Numeric. If analysis is only for a subset of the sampled
#'               gradient then the minimum point within the subset
#' @param highest Numeric. If analysis is only for a subset of the sampled
#'               gradient then the maximum point within the subset
#' @details range_shuffle impliments simulations described by Wang et.al (2012)
#'          to estimate effect of geometric  and environmental constraints on
#'          pattern in species richness across spatial graidients. It calculates
#'          a vector of all possible range locations for each given range extent
#'          based on the conditions for geometric constraints given by
#'          'boundary'. Range locations are randomized by sampling from this
#'          vector.
#' @return If degen is FALSE, a data frame with four colums for mean, SD and
#'        confidence intervals of expected richness
#'
#' \itemize{
#'  \item{"mod.rich"}{ mean richness of each site}
#'  \item{"mod.sd"}{ standard deviation of species richness}
#'  \item{"q2.5"}{ lower limit of the confidence interval}
#'  \item{"q97.5"}{ upper limit of the confidence interval}
#' }
#'        If degen is TRUE, then a list containing above data frame and a list
#'        of all the randomized matrices
#' @references Wang, X., and J. Fang. 2012. Constraining null models with
#'             environmental gradients: a new method for evaluating the effects
#'             of environmental factors and geometric constraints on geographic
#'             diversity patterns. Ecography 35:1147-1159.
#' @examples
#' data(rangedata)
#' range_shuffle_rnd <- range_shuffle(x=rangedata,boundary = 'hh',
#'                                    interval = 200, var = NULL,
#'                                    sites = seq(600,2400,200),
#'                                    reps =10,degen = FALSE)
#' range_shuffle_rnd
#' plot(range_shuffle_rnd$mean,ylab = "Species Richness",pch = 19)

#' @export
range_shuffle <- function(x,boundary,var,interval,sites,reps,degen,
                          lowest = NA,highest = NA){

  if(is.na(lowest)){
    min_pt <- min(x$min)
  }else{min_pt <- lowest}

  if(is.na(highest)){
    max_pt <- max(x$max)
  }else(max_pt <- highest)


  rsfd <- table(x$range)


  #make a list of feasible midpoints for each range size
  midpoints_list <- lapply(seq(min(x$range),(max_pt - min_pt),interval),
                           function(x,bnd = boundary){
                             switch(bnd,
                                    hh = {seq((x/2)+min_pt,max_pt - (x/2),
                                              by = interval)},
                                    sh = {rev(seq(max_pt - (x/2),min_pt,
                                                  by = -interval))},
                                    hs = {seq((x/2)+min_pt,max_pt,
                                              by = interval)},
                                    ss = {seq(min_pt,max_pt,by = interval)}
                             )
                           })

  midpoints_list <- stats::setNames(midpoints_list,seq(min(x$range),(max_pt - min_pt),
                                                interval))



  #sample the apropriate vector in midpoints_list with replacement using range
  #frequencies in rsfd. Then construct range min and max based on range size.
  shuffle_sprich <- list()
  shuffle_dist <- list()
  for(i in 1:reps){
    shuffle <- lapply(1:length(rsfd),
                      function(y){
                        r <- names(rsfd)[y]
                        freq <- rsfd[y]
                        if(is.null(var)){
                          weights <- NULL
                        }else{

                          weights <-var$weights[var$mid%in%midpoints_list[[r]]]
                        }
                        mid <- if(length(midpoints_list[[r]]) == 1){
                          midpoints_list[[r]]}else{
                            sample(midpoints_list[[r]],
                                      size = freq,
                                      replace = T,
                                      prob = weights )}
                        switch(boundary,
                               hh = {mid
                               },
                               sh = {
                                 mid[sapply(mid,function(a){
                                   a < min(x$min) + as.numeric(r)/2})]<-
                                   min(x$min) + as.numeric(r)/2
                               },
                               hs = {
                                 mid[sapply(mid,function(a){
                                   a > max(x$max) - as.numeric(r)/2})]<-
                                   max(x$max) - as.numeric(r)/2
                               },
                               ss = {
                                 mid[sapply(mid,function(a){
                                   a < min(x$min) + as.numeric(r)/2})]<-
                                   min(x$min) + as.numeric(r)/2
                                 mid[sapply(mid,function(a){
                                   a > max(x$max) - as.numeric(r)/2})]<-
                                   max(x$max) - as.numeric(r)/2
                               }
                        )

                        data.frame(mid = mid,
                                   min = mid - (as.numeric(r)/2),
                                   max = mid + (as.numeric(r)/2)
                        )

                      } )
    #combine all ranges
    shuffle_dat <- do.call(rbind,shuffle)
    shuffle_dat$range <- shuffle_dat$max - shuffle_dat$min

    #calculate species richness
    sites1 <- apply(shuffle_dat,1,
                    function(z){
                      sites[sites>=z[2]&sites<=z[3]]
                    })
    sites1 <- do.call(c,as.list(sites1))
    shuffle_sprich[[i]] <- data.frame(table(sites1))
    if(degen == T){
      shuffle_dist[[i]] <- shuffle_dat
    }
  }

  #bring all the species richness values together
  shuffle_sprich_dat <- do.call(rbind,shuffle_sprich)

  shuffle_sprich_out <- do.call(rbind,
                                with(shuffle_sprich_dat,
                                     tapply(Freq,sites1,
                                            function(x){
                                              data.frame(mean = mean(x),
                                                         sd = sd(x),
                                                         q025 = quantile(x,0.025),
                                                         q975 = quantile(x,0.975))
                                            }))
  )



  if(degen == T){
    out <- list(out_dist = shuffle_dist,out_dat = shuffle_sprich_out)
  }else(shuffle_sprich_out)

}
