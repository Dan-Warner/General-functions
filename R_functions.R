## Cumulative frequency ----
# does not work with adorn_totals() 
adorn_cumulative <- function(dat, colname, dir = "down"){

  if(!missing(colname)){
    colname <- rlang::enquo(colname)
  } else if("valid_percent" %in% names(dat)) {
  colname <- rlang::sym("valid_percent")
  } else if("percent" %in% names(dat)){
    colname <- rlang::sym("percent")
  } else {
    stop("\"colname\" not specified and default columns valid_percent and percent are not present in data.frame dat")
  }

  target <- dplyr::pull(dat, !! colname)

  if(dir == "up"){
    target <- rev(target)
  }
  dat$cumulative <- cumsum(ifelse(is.na(target), 0, target)) + target*0 # an na.rm version of cumsum, from https://stackoverflow.com/a/25576972
  if(dir == "up"){
    dat$cumulative <- rev(dat$cumulative)
    names(dat)[names(dat) %in% "cumulative"] <- "cumulative_up"
  }
  dat
}




## Summarizes data between subjects ----
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



# Norm the data within specified groups in a data frame; it normalizes each ---- 
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}



## Summarizes data within subjects ----
## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}





# Upper error bar
geom_uperrorbar <- function(mapping = NULL, data = NULL,
                            stat = "identity", position = "identity",
                            ...,
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomUperrorbar,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}



# Specify decimals ----

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# specify_decimal(mean(na.omit(df$a)),1)



# Percentile rank ----
# Takes into account NA's (rank()/length() and dplyr::percent_rank() do not work with NA's and NA's will be accorded a percentile rank

prct_rank <-function(p){
  r<-rank(p)/sum(!is.na(p))
  r[is.na(p)]<-NA
  r
} 

# prct_rank(df, a)



# Arc-sine transformation ----
asinTransform <- function(p) { asin(sqrt(p/100)) }  

# asinTransform(x)



# Legend on the bottom of the graph (alternatively, use Cowplot package) ----

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1, 
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}



# Color format ----

colFmt = function(x,color){
  outputFormat = opts_knit$get("rmarkdown.pandoc.to")
  if(outputFormat == 'latex')
    paste("\\textcolor{",color,"}{",x,"}",sep="")
  else if(outputFormat == 'html')
    paste("<font color='",color,"'>",x,"</font>",sep="")
  else
    x
}

# `r colFmt("write here...",'red')`



# Multi-Class Summary Function ---- 
  # Based on caret:::twoClassSummary
  require(compiler) # now comes bundled with any new R release, since R 2.13
  multiClassSummary <- cmpfun(function (data, lev = NULL, model = NULL){
    
    #Load Libraries
    require(Metrics)
    require(caret)
    
    #Check data
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
      stop("levels of observed and predicted data do not match")
    
    #Calculate custom one-vs-all stats for each class
    prob_stats <- lapply(levels(data[, "pred"]), function(class){
      
      #Grab one-vs-all data for the class
      pred <- ifelse(data[, "pred"] == class, 1, 0)
      obs  <- ifelse(data[,  "obs"] == class, 1, 0)
      prob <- data[,class]
      
      #Calculate one-vs-all AUC and logLoss and return
      cap_prob <- pmin(pmax(prob, .000001), .999999)
      prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
      names(prob_stats) <- c('ROC', 'logLoss')
      return(prob_stats) 
    })
    prob_stats <- do.call(rbind, prob_stats)
    rownames(prob_stats) <- paste('Class:', levels(data[, "pred"]))
    
    #Calculate confusion matrix-based statistics
    CM <- confusionMatrix(data[, "pred"], data[, "obs"])
    
    #Aggregate and average class-wise stats
    #Todo: add weights
    class_stats <- cbind(CM$byClass, prob_stats)
    class_stats <- colMeans(class_stats)
    
    #Aggregate overall stats
    overall_stats <- c(CM$overall)
    
    #Combine overall with class-wise stats and remove some stats we don't want 
    stats <- c(overall_stats, class_stats)
    stats <- stats[! names(stats) %in% c('AccuracyNull', 'Prevalence', 'Detection Prevalence')]
    
    #Clean names and return
    names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
    return(stats)
  })
  
  
  
# empty as NA ----
  empty_as_na <- function(x){
    if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
    ifelse(as.character(x)!="", x, NA)
  }
  
  # transform all columns: iris %>% mutate_each(funs(empty_as_na)) 
  
  
  
# favourite theme topic for plots
  my_theme <- function() {
    theme_light(base_size = 14) +
      theme(# text
            #axis.text = element_text(colour = "grey30"),
            text = element_text(colour = "grey30"), 
            # panel
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "grey50", fill = NA, size = 0.7),
            # facet appearance
            strip.background = element_rect(fill = 'transparent', colour = 'transparent', linetype = "solid"),
            strip.text = element_text(face = "plain", size = 14, colour = "grey30"), 
            # legend
            legend.background = element_rect(fill = "transparent"),
            legend.text = element_text(colour = "grey30")
            # margin
            #plot.margin = unit(c(t = 0.5, r = 0.5, b = 0.5, l = 0.5), "cm")
            ) 
  }
  # ggplot() + my_theme()
  

  
  
# mutate multiple variables based on one condition ---- 
  
  mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
    condition <- eval(substitute(condition), .data, envir)
    .data[condition, ] <- .data[condition, ] %>% mutate(...)
    .data
  }
  
  # df %>%
  #   mutate_cond(cnt_tests_5_35_dim < 4, 
  #               sd_fa_de_novo_tfa = NA, sd_fa_mixed_tfa = NA, sd_fa_preformed_tfa = NA, sd_fa_c18_1_tfa = NA,
  #               sd_fa_de_novo_day = NA, sd_fa_mixed_day = NA, sd_fa_preformed_day = NA, sd_fa_c18_1_day = NA,
  #               sd_fa_de_novo = NA, sd_fa_mixed = NA, sd_fa_preformed = NA, sd_fa_c18_1 = NA,
  #               sd_ratio_181_14 = NA, sd_ratio_181_16 = NA
  #   )
