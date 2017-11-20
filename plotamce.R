###############################################################################

## Getting conditional estimates given (1) AMCE object
## (2) LIST containing levels at which conditional effects will be calculated
## (3) NAME of current conditional variable
## (4) VALUE of current level of current conditioning variable
## (5) NAME of profile attribute being modified by conditioning variable

#############################################################################

get.conditional.effects <- function(object, conditional.levels, current.effect, current.level, mod.var) {
  
  #amce object
  amce_obj <- object
  #make dummy data and set in base levels
  cond.data <- amce_obj$data
  #set factor vars to baselines
  for (var in names(amce_obj$baselines)) {
    cond.data[[var]] <- amce_obj$baselines[[var]]
  }
  #set in conditional vars to their first given level
  for (var in names(conditional.levels)) {
    cond.data[[var]] <- conditional.levels[[var]][1]
  }        
  #set current effect to current level
  cond.data[[current.effect]] <- current.level
  
  #original levels of conditional vars
  orig.levels <- sapply(all.vars(amce_obj$cond.formula)[-1] [all.vars(amce_obj$cond.formula)[-1] %in% names(amce_obj$baselines)], function(x) levels(amce_obj$data[[x]]), simplify=F)
  #coefficients associated with conditional estimates base term
  cond.base <- unlist(sapply(amce_obj$respondent.varying,USE.NAMES = F, function(x) colnames(amce_obj$cond.estimates[[x]])))
  #estimated coefficients, adding 0 for intercept
  cond.beta <- c(0,do.call(cbind,amce_obj$cond.estimates)[1,])
  
  #blank output
  estimates.vector <- c()
  error.vector <- c()
  names.vector <- c()
  
  #quick covariance getting function
  cov.ij <- function(var1,var2) {
    out <- pred_mat[var1]*pred_mat[var2]*amce_obj$vcov.resp[var1,var2]
    return(out)
  }
  cov.ij <- Vectorize(cov.ij,vectorize.args = c("var1","var2"))
  
  #function for NA multiplication to be used in special cases
  na.multiply <- function(x,y) {
    vec <- c(x,y)
    #If either is NA and other is 0, return 0
    if (any(is.na(vec)) && vec[!is.na(vec)] == 0) {
      out <- 0
    } else {
      #otherwise normal (so 1*NA = NA)
      out <- x*y
    }
  }
  na.multiply <- Vectorize(na.multiply,vectorize.args = c("x","y"))
  
  # split up modified variable
  mod_vars <- strsplit(mod.var,":")[[1]]
  # loop over levels
  for (mod_coef in colnames(amce_obj$cond.estimates[[mod.var]])) {
    ## (1) Edit conditional data
    # split level coefficient into components
    mod_coefs <- strsplit(mod_coef,":")[[1]]
    # edit cond data to fit this level
    for (x in 1:length(mod_coefs)) {
      mod_lev <- sub(mod_vars[x],"",mod_coefs[x])            
      cond.data[[mod_vars[x]]] <- mod_lev
    }
    ## (2) Make model matrix
    pred_mat <- model.matrix(amce_obj$cond.formula,cond.data,xlev = orig.levels)
    ## (3) Turn off base term for this conditional var
    turn_off <- rep(1,ncol(pred_mat))
    names(turn_off) <- colnames(pred_mat)         
    turn_off[cond.base] <- 0
    # Use to turn off terms in pred_mat that only contain respondent varying items
    pred_mat <- pred_mat[1,]*turn_off
    ## (4) Calculate coefficient and SE
    # (a) Coefficient
    if (!any(is.na(cond.beta))) {
      pred_val <- sum(pred_mat*cond.beta)
    } else {
      #otherwise use special function
      pred_val <- sum(na.multiply(pred_mat,cond.beta))
    }
    # (b) SE
    if (!is.na(pred_val)) {     
      #variable names sans intercept
      vars <- colnames(amce_obj$vcov.resp)[2:ncol(amce_obj$vcov.resp)]
      #all other covariance combinations
      all_cov <- outer(vars,vars,FUN= function(x,y) cov.ij(x,y))
      pred_se <- sqrt(sum(all_cov))
    } else {
      pred_se <- NA
    }
    ## And print out
    estimates.vector <- c(estimates.vector,pred_val)
    error.vector <- c(error.vector,pred_se)
    names.vector <- c(names.vector,mod_coef)
  }
  
  #return list of modified coefficient estimates and se's
  names(estimates.vector) <- names(error.vector) <- names.vector
  out <- rbind(estimates.vector,error.vector)
  return(out)
  
}

########################################################
# plot amce function
#######################################################

# default will return single plot with point estimates and CI's
# to facet plots, give "facet.name" the name of the variable to facet by...
# ... by default variable facets will be ALL level combinations...
# ... or if a continuous variable is given, will use quantiles...
# ... to customize, directly give "facet.levels" a named LIST ...
# ... with desired values of each variable entered as DATA FRAME with desired names
# "display" takes one of: "all", "unconditional", "conditional"
# with no given facet and no respondent vars, this choice is irrelevant (unconditional only)
# with no given facet and respondent vars, default is all but can choose just unconditional or interaction
# with a facet given (respondent or otherwise), similarly can choose
#' @method plot amce
#' @export 
plot.amce2 <- function(x, main="", xlab="Change in E[Y]", ci=.95, colors=NULL, xlim=NULL, breaks=NULL, labels=NULL, attribute.names = NULL, level.names = NULL, label.baseline = TRUE, text.size=11, text.color = "black", point.size = .5, dodge.size=0.9, plot.theme = NULL, plot.display = "all", facet.names = NULL, facet.levels = NULL, ...) {
  
  # You need ggplot2
  amce_obj <- x
  ylim <- xlim
  
  # Make R CMD check happy
  pe <- NULL
  se <- NULL
  group <- NULL
  lower <- NULL
  upper <- NULL
  var <- NULL
  printvar <- NULL
  facet <- NULL
  
  ############################## basic set-up: get attributes and levels
  
  # Extract raw attribute names from the amce_obj$estimates object
  raw_attributes <- names(amce_obj$estimates)
  # Extract raw levels (coefficient names)
  raw_levels <- lapply(amce_obj$estimates,colnames)
  
  # Determine baseline level for each effect estimate in raw_levels and append to beginning of each vector in raw_levels
  for (effect in names(raw_levels)) {
    effect_elements <- strsplit(effect, ":")[[1]]
    baseline_interactions <- c()
    for (elem in effect_elements) {
      # get baseline, as if coefficient name
      base_coef <- paste(c(elem,amce_obj$baselines[[elem]]),collapse="")
      baseline_interactions <- c(baseline_interactions,base_coef)
    }
    interaction_str <- paste(baseline_interactions,sep="",collapse=":")
    raw_levels[[effect]] <- c(interaction_str, raw_levels[[effect]]) 
  }
  
  ################################### Incorporate and adjust user-input: general
  
  # Convert ci to z-score
  if (ci < 1 & ci > 0) {
    zscr <- qnorm(1- ((1-ci)/2))
  } else {
    cat("Invalid confidence interval -- Defaulting to 95%")
    zscr <- qnorm(1- ((1-.95)/2))
  }
  
  ################################### Incorporate and adjust user-input: naming    
  
  # Sanity check user-provided attribute.names against AMCE objects
  if (!is.null(attribute.names)) {
    attribute.names <- unique(attribute.names)
    if (length(attribute.names) != length(raw_attributes)) {
      cat(paste("Error: The number of unique elements in attribute.names ", length(attribute.names), " does not match the attributes in amce object for which estimates were obtained: ", paste(raw_attributes,collapse=", "), "\n", sep=""))
      cat("Defaulting attribute.names to attribute names in AMCE object\n")
      attribute.names <- NULL
    }
  }
  
  # Sanity check user-provided level.names against AMCE object
  if (!is.null(level.names)) {
    names(level.names) <- clean.names(names(level.names))
    for (name in names(level.names)) {
      if (name %in% names(raw_levels)) {
        if (length(level.names[[name]]) != length(raw_levels[[name]])) {
          cat(paste("Error: level.names lengths do not match levels for attribute ", name, "\n",sep=""))
          cat(paste("Defaulting level.names for attribute ", name, " to level names in AMCE object", "\n",sep=""))
          level.names[[name]] <- NULL
        }   
      } else {
        cat(paste("Error: level.names entry ",name," not in AMCE object. Removing level.names for attribute.","\n",sep=""))
        level.names[[name]] <- NULL
      }
    }
  }
  
  # If no attribute name or changed to NULL, use initial user supplied names as attribute names
  if (is.null(attribute.names)) {
    attribute.names <- c()
    for (attr in names(amce_obj$estimates)) {
      attr_split <- strsplit(attr,":")[[1]]
      attr_lookup <- paste(unlist(sapply(attr_split,function(x) amce_obj$user.names[x])), collapse=":")
      attribute.names <- c(attribute.names,attr_lookup)
    }
  }
  
  # If no level names make blank list
  if (is.null(level.names)) level.names <- list()
  # fill in blank list or missing levels, if any
  if (any(!names(raw_levels) %in% names(level.names))) {
    for (attr in names(raw_levels)[!names(raw_levels) %in% names(level.names)]) {
      attr_split <- strsplit(raw_levels[[attr]],":")
      level.names[[attr]] <- unlist(lapply(attr_split,function(x) paste(sapply(x,function(y) amce_obj$user.levels[y]),collapse=":")))
    }
  }
  
  ################################### Incorporate and adjust user-input: facetting
  
  # valid plot.display option?
  plot.display.opts <- c("all","unconditional","interaction")
  if (!is.element(plot.display,plot.display.opts)) {
    stop(paste(c("Error-- plot.display must be once of: ",paste(plot.display.opts,collapse=", ")),collapse=" "))
  }
  
  #clean facet names; if levels but no names? level names are facets
  if (!is.null(facet.names)) {
    facet.names <- clean.names(facet.names)
  } else if (!is.null(facet.levels)) {
    facet.names <- clean.names(names(facet.levels))
  }
  
  #check that they are in AMCE object
  if (!is.null(facet.names)) {
    facet.names.check <- c()
    for (facet.name in facet.names) {
      if (grepl(":",facet.name)) stop("Error-- cannot facet by interaction in current version.")
      if (!facet.name %in% names(amce_obj$estimates) & !facet.name %in% names(amce_obj$cond.estimates)) {
        stop(paste(c("Error-- cannot find facet name",facet.name,"in AMCE object output."),collapse=" "))
      } else {
        facet.names.check <- c(facet.names.check,facet.name)
      }
    }
    facet.names <- facet.names.check
  }
  
  #if no facets but there are respondent varying characteristics, use those
  if ((is.null(facet.names)) & (length(amce_obj$respondent.varying) > 0) & (plot.display != "unconditional")) {
    facet.names <- amce_obj$respondent.varying
  }
  
  #no facet name or resp var, must be unconditional
  if (is.null(facet.names) & plot.display == "interaction") {
    warning("Warning: no facet name or respondent varying characteristic provided to calculate conditional estimates. Will display unconditional only")
    plot.display <- "unconditional"
  }
  
  #unconditional but facet names given? remove facet names
  if(plot.display == "unconditional" & !is.null(facet.names)) {
    warning("Warning-- plot display is set to unconditional, facet names will be ignored")
    facet.names <- NULL
    facet.levels <- NULL
  }
  
  #check and clean facet levels if provided
  if (!is.null(facet.levels)) {
    #clean names of facet levels
    names(facet.levels) <- clean.names(names(facet.levels))
    #clean actual levels
    for (facet.name in names(facet.levels)) {
      #if it's a factor, clean up level names
      if (facet.name %in% names(amce_obj$baselines)) {
        facet.levels[[facet.name]] <- clean.names(facet.levels[[facet.name]])
        #make sure that if it's profile-varying, there's more than base
        if (facet.name %in% names(amce_obj$estimates) && is.element(amce_obj$baselines[[facet.name]],facet.levels[[facet.name]])) {
          stop (paste(c("Error: Facet level \"",as.character(amce_obj$baselines[[facet.name]]), "\" is the baseline level of a profile varying attribute. Please provide alternative facet level or use defaults."), collapse=""))               
        }
        #names from user input if none provided
        if (is.null(names(facet.levels[[facet.name]]))) {
          fac.levs <- sapply(facet.levels[[facet.name]],function(x) paste(facet.name,x,sep=""))
          names(facet.levels[[facet.name]]) <- sapply(fac.levs, USE.NAMES = F, function(x) amce_obj$user.levels[[x]])
        }
      } else if (is.null(names(facet.levels[[facet.name]]))) {
        #not a factor and no names, just take level values
        names(facet.levels[[facet.name]]) <- as.character(facet.levels[[facet.name]])
      }
    }
  }
  
  #if user didn't give any levels, make blank list
  if (is.null(facet.levels)) facet.levels <- list()
  #input missing levels if any
  if (any(!facet.names %in% names(facet.levels))) {
    for (facet.name in facet.names[!facet.names %in% names(facet.levels)]) {
      #if it's a factor, default facet levels are all levels
      if (facet.name %in%  names(amce_obj$baselines)) {
        if (facet.name %in% names(amce_obj$estimates)) {
          #if NOT respondent varying get levels and names from ESTIMATES
          fac.levs <- colnames(amce_obj$estimates[[facet.name]])
        } else {
          #get levels and names from COND.ESTIMATES
          fac.levs <- colnames(amce_obj$cond.estimates[[facet.name]])
        }
        # get pure levels
        facet.levels[[facet.name]] <- sub(facet.name,"",fac.levs)
        #add in baseline 
        facet.levels[[facet.name]] <- c(amce_obj$baselines[[facet.name]], facet.levels[[facet.name]])
        #names from user input
        fac.levs <- c(paste(facet.name,amce_obj$baselines[[facet.name]],sep=""),fac.levs)
        names(facet.levels[[facet.name]])  <- sapply(fac.levs, USE.NAMES = F,function(x) amce_obj$user.levels[[x]])
      } else if (facet.name %in% names(amce_obj$continuous)) {
        #if it's continuous, default is quantiles
        facet.levels[[facet.name]] <-  amce_obj$continuous[[facet.name]]
      }             
    }
  }
  
  #the equivalent of summary's "covariate values" are respondent-varying entries
  #so get just those
  covariate.values <- list()
  for (var in names(facet.levels)) {
    if (var %in% amce_obj$respondent.varying) {
      covariate.values[[var]] <- facet.levels[[var]]
    }
  }
  
  ################################### Compile estimates into plottable objects
  
  #blank data frame for plot data
  d <- data.frame(pe=c(), se=c(), upper=c(), lower=c(), var=c(), printvar = c(), group=c(), facet=c())
  
  ############# Unconditional estimates
  
  #only display if plot.display == all or unconditional
  if (plot.display != "interaction") {
    #if plot.display == all, add unconditional facet name (not needed for unconditional only)
    if (plot.display == "all") {
      uncond.facet.name <- "Unconditional"
    } else {
      uncond.facet.name <- NA
    }
    #if plot = all and there are non-respondent varying facet names
    #remove them from raw attributes
    if (plot.display == "all" && !is.null(facet.names)) {
      attr_remove <- c()
      for (facet.name in facet.names[!is.element(facet.names, amce_obj$respondent.varying)]) {
        attr_remove1 <- raw_attributes[grepl(":",raw_attributes)]
        attr_remove1 <- attr_remove1[grepl(facet.name,attr_remove1)]
        attr_remove <- c(attr_remove,attr_remove1)
      }
      raw_attributes <- raw_attributes[!is.element(raw_attributes,attr_remove)]
    }       
    #loop over raw attribute names
    for (i in 1:length(raw_attributes)) {
      #get raw attribute name
      attr_name <- raw_attributes[i]
      #get attribute name to print
      print_attr_name <- attribute.names[which(names(amce_obj$estimates) == raw_attributes[i])]
      #set up basic group header and add to plot
      d_head <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var= attr_name, printvar=paste(print_attr_name, ":", sep=""), group="<NA>",facet=uncond.facet.name)
      d <- rbind(d,d_head)    
      #iterate over levels
      for (j in 1:length(raw_levels[[attr_name]])) {
        #raw level name
        level_name <- raw_levels[[attr_name]][j]
        #get level name to print
        print_level_name <- level.names[[attr_name]][j]              
        #if on the first level
        if (j == 1) {
          if (label.baseline) {
            print_level_name <- paste("(Baseline = ",print_level_name,")",sep="")
          }
          #get the baseline and print a blank line
          d_lev <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=level_name, printvar=paste("   ", print_level_name,sep=""), group=print_attr_name, facet=uncond.facet.name) 
        } else {
          #retrieve estimate and SE
          val_pe <- amce_obj$estimates[[attr_name]][1,level_name]
          val_se <- amce_obj$estimates[[attr_name]][2,level_name]       
          #calculate bounds
          upper_bnd <- val_pe + zscr*val_se
          lower_bnd <- val_pe - zscr*val_se 
          #make line to add to plot data
          d_lev <- data.frame(pe=val_pe, se=val_se, upper=upper_bnd, lower=lower_bnd, var=level_name, printvar=paste("   ", print_level_name,sep=""), group=print_attr_name, facet=uncond.facet.name)       
        } #end if a baseline
        #add to plot
        d <- rbind(d,d_lev)
      } #end loop over levels
    } #end loop over non-facet related attribute names
  } #end if plot.display == all or plot.display == conditional
  
  ############# Conditional estimates
  
  #Only if plot.display is all or conditional and we got a facet name from somehere
  if (plot.display != "unconditional" & !is.null(facet.names)) {
    
    #loop over facets
    for (facet.name in facet.names) {    
      
      #how to print it
      print_facet_name <- amce_obj$user.names[[facet.name]]
      #### identify all REQUESTED terms involving facet name
      all_req_vars <- attr(terms(amce_obj$formula),"term.labels")
      all_mod <- unlist(sapply(all_req_vars,function(x) {
        y <- strsplit(x,":")[[1]]
        if (any(y == facet.name)) x
      }))
      #figure out profile attributes these refer to
      all_mod <- unlist(sapply(all_mod,function(x) {
        subs <- strsplit(x,":")[[1]]
        subs <- subs[is.element(subs,names(amce_obj$estimates))]
        subs <- subs[subs != facet.name]
        if (length(subs) > 0) paste(subs,collapse=":")
      }))
      #make sure there are some
      if (length(all_mod) == 0) {
        stop(paste(c("Error: Facet variable",facet.name,"not interacted with profile attributes"),collapse=" "))
      }
      #just unique ones
      all_mod <- unique(all_mod)
      #Temp Bug Fixing
      if (is.element(facet.name,names(amce_obj$estimates))) {
        #if ACIE, then remove baseline of facet in the d dataset filling process
        facet.start <- 2
      } else {
        #if conditional on respondent varying
        facet.start <- 1
      }
      
      #for each actual facet level make new set of plot data
      for (k in facet.start:length(facet.levels[[facet.name]])) {
        
        # set level
        facet_lev <- facet.levels[[facet.name]][k]
        #how to print facet level
        if (is.element(facet.name,names(amce_obj$estimates))) {
          #if ACIE
          print_facet_level <- paste(c("ACIE", paste(c(print_facet_name, names(facet.levels[[facet.name]])[k]), collapse = " = ")), collapse = "\n")
        } else {
          #if conditional on respondent varying
          print_facet_level <- paste(c("Conditional on",paste(c(print_facet_name, names(facet.levels[[facet.name]])[k]), collapse = " = ")), collapse = "\n")
        }
        
        #loop over variables to be modified
        for (mod_var in all_mod) {
          #how to print modified attribute
          print_attr_name <- attribute.names[which(names(amce_obj$estimates) == mod_var)]
          #set up header to reflect base (non-facet) category
          d_head <- data.frame(pe=NA, se=NA, upper=NA, lower=NA,var=mod_var, printvar=paste(print_attr_name, ":", sep=""), group="<NA>", facet=print_facet_level)
          #add new header
          d <- rbind(d, d_head)
          #Get estimates
          if (facet.name %in% names(amce_obj$estimates)) {
            #figure out interaction name
            inter_coef <- paste(sort(c(mod_var,facet.name)),collapse = ":")
            #get from unconditional estimates
            estimate.source <- amce_obj$estimates[[inter_coef]]
            estimate.source <- estimate.source[,grep(paste0(facet.name, facet_lev), colnames(estimate.source))]
          } else {
            #calculate from function if conditional
            estimate.source <- get.conditional.effects(amce_obj, covariate.values, facet.name, facet_lev, mod_var)
          }                 
          
          #split into components
          mod_vars <- strsplit(mod_var,":")[[1]]
          #iterate over levels of modified variable
          for (p in 1:length(raw_levels[[mod_var]])) {          
            #raw level name is original coefficient name
            mod_coef <- raw_levels[[mod_var]][p]
            #split it up
            mod_coefs <- strsplit(mod_coef,":")[[1]]
            #modify data dummy
            for (lev in 1:length(mod_coefs)) {
              #get level name from coefficient name 
              mod_lev <- sub(mod_vars[lev],"",mod_coefs[lev])
            }                       
            #get level name to print
            print_level_name <- level.names[[mod_var]][p] 
            #get the baseline of modified var and make a blank line
            if (p == 1) {
              if (label.baseline) {
                print_level_name <- paste("(Baseline = ",print_level_name,")",sep="")
              }
              d_lev <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var = mod_coef, printvar=paste("   ",  print_level_name,sep=""), group=print_attr_name, facet= print_facet_level)        
            } else {
              #retrieve estimate and SE
              val_pe <- estimate.source[1,p-1]
              if (!is.na(val_pe)) {
                val_se <- estimate.source[2,p-1]
                #calculate bounds
                upper_bnd <- val_pe + zscr*val_se
                lower_bnd <- val_pe - zscr*val_se 
              } else {
                val_se <- upper_bnd <- lower_bnd <- NA
              }
              #make line to add to plot data
              d_lev <- data.frame(pe=val_pe, se=val_se, upper=upper_bnd, lower=lower_bnd, var = mod_coef, printvar=paste("   ", print_level_name,sep=""), group=print_attr_name, facet=print_facet_level)    
            }
            #add level data to plot data
            d  <- rbind(d, d_lev)
          } #end loop over levels of modified var
          
        } #end loop over all modified vars
      } #end loop over level of facetted variable
    } #end loop over facets
  } else {
    #if there are no facets or plot.display is unconditional, remove that column
    d <- d[,-which(colnames(d) == "facet")]
  }
  
  #################    format "d" dataframe
  
  # Set Y bounds
  if (is.null(ylim)) {
    max_upper <- max(d$upper, na.rm=T) + .05
    min_lower <- min(d$lower, na.rm=T) - .05
    ylim <- c(min_lower, max_upper)
    d[is.na(d)] <- max_upper + 100
  } else {
    d[is.na(d)] <- max(ylim) + 100
  }
  
  # Make group factors <NA> actually NA
  d$group[d$group == "<NA>"] <- NA
  #same with facet
  if(!is.null(facet.names)) d$facet[d$facet == "<NA>"] <- NA
  
  # Reverse factor ordering
  d$var <- factor(d$var,levels=unique(d$var)[length(d$var):1])
  #make facet into factor, if it exists
  if (!is.null(facet.names)) {
    d$facet <- factor(d$facet,levels=unique(d$facet))
  }
  
  ########## plot output
  
  p = ggplot(d, aes(y=pe, x=var, colour=group))
  p = p + coord_flip(ylim = ylim) 
  p = p + geom_hline(yintercept = 0,size=.5,colour="black",linetype="dotted") 
  p = p + geom_pointrange(aes(ymin=lower,ymax=upper),position=position_dodge(width=dodge.size),size=point.size)
  
  #add facetting
  if (!is.null(facet.names)) {
    p = p + facet_wrap(~ facet)
  } 
  
  # If breaks and labels Null, use default
  if (is.null(breaks) & is.null(labels)) {
    p = p + scale_y_continuous(name=xlab)
  } else if (is.null(breaks) & !is.null(labels)) {
    p = p + scale_y_continuous(name=xlab, labels=labels)
  } else if (!is.null(breaks) & is.null(labels)) {
    p = p + scale_y_continuous(name=xlab, breaks=breaks)
  } else if (!is.null(breaks) & !is.null(labels)) {
    p = p + scale_y_continuous(name=xlab, breaks=breaks, labels=labels)
  }
  fix.xlabs <- as.character(d$printvar)[!duplicated(d$printvar)]
  p = p + scale_x_discrete(name="",labels=fix.xlabs[length(fix.xlabs):1])
  
  # If there's a title,add it
  if (!is.null(main)) {
    if (main != "") {
      p = p + ggtitle(main)
    }
  }
  # If no colors specified, use default
  if (is.null(colors)) {
    p = p + scale_colour_discrete(" ") 
  } else if (is.vector(colors)) {
    # Make manual palette
    cPal <- rep(colors, ceiling(length(unique(d$group))/length(colors)))
    # Use manual palette
    p = p + scale_colour_manual(values=cPal)
  } else {
    cat("Error: 'colors' must be a vector. Using default colors\n")
    p = p + scale_colour_discrete(" ") 
  }
  
  # colour scheme
  # if no theme specified, use default
  if (is.null(plot.theme)){
    
    theme_bw1 <- function(base_size = text.size, base_family = "") {
      theme_grey(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.text.x = element_text(size = base_size*.9, colour = text.color,  hjust = .5 , vjust=1),axis.text.y = element_text(size = base_size, colour = text.color, hjust = 0 , vjust=.5 ), axis.ticks = element_line(colour = "grey50"),axis.title.y =  element_text(size = base_size,angle=90,vjust=.01,hjust=.1),plot.title = element_text(face = "bold"),legend.position = "none")
    }
    
    p = p + theme_bw1()
    
  } else if (is.null(class(plot.theme)))  {
    
    cat("Error: 'plot.theme' is not a valid ggplot theme object. Using default theme\n")
    theme_bw1 <- function(base_size = text.size, base_family = "") {
      theme_grey(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.text.x = element_text(size = base_size*.9, colour = text.color,  hjust = .5 , vjust=1),axis.text.y = element_text(size = base_size, colour = text.color, hjust = 0 , vjust=.5 ), axis.ticks = element_line(colour = "grey50"),axis.title.y =  element_text(size = base_size,angle=90,vjust=.01,hjust=.1),plot.title = element_text(face = "bold"),legend.position = "none")
    }
    
    p = p + theme_bw1()
    
  } else if (class(plot.theme)[1] != "theme") {
    
    cat("Error: 'plot.theme' is not a valid ggplot theme object. Using default theme\n")
    theme_bw1 <- function(base_size = text.size, base_family = "") {
      theme_grey(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.text.x = element_text(size = base_size*.9, colour = text.color,  hjust = .5 , vjust=1),axis.text.y = element_text(size = base_size, colour = text.color, hjust = 0 , vjust=.5 ), axis.ticks = element_line(colour = "grey50"),axis.title.y =  element_text(size = base_size,angle=90,vjust=.01,hjust=.1),plot.title = element_text(face = "bold"),legend.position = "none")
    }
    
    p = p + theme_bw1()
    
    # otherwise use the user-passed theme   
  } else {
    p = p + plot.theme
  }
  
  #console message with level to hold resp vars as
  if (length(covariate.values) > 1) {
    resp.message <- c("Note:")
    for (this.var in names(covariate.values)) {
      resp.message <- paste(c(resp.message," For AMCE and ACIE conditional on ",this.var,", "),collapse="")
      other.vars <- names(covariate.values)[names(covariate.values) != this.var]
      
      other.levels <- c()
      for (var in other.vars) {
        other.levels <- c(other.levels,paste(c(var," will be held at level \"",names(covariate.values[[var]])[1],"\""),collapse = ""))                
      }
      other.levels <- paste(other.levels,collapse = ", and ")
      resp.message <- c(resp.message,other.levels,".")
      resp.message <- paste(resp.message,collapse = "")
    }
    cat(resp.message,"\n")
  }
  print(p)
}
