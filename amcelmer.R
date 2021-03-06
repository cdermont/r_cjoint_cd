#########################################################################
## function for removing ALL punctuation, symbols, and spaces from string
## from elements of vector "vec"
########################################################################

clean.names <- function(str) {
  #split components of interactions
  x <- strsplit(str,":")[[1]]
  #and apply cleaning separately, removing any punctuation (P), symbols (S), and separators (Z)
  x <- gsub("[\\p{P}\\p{S}\\p{Z}]","",x,perl=T)
  #re-attach
  paste(x,collapse=":")
}
clean.names <- Vectorize(clean.names,vectorize.args=("str"),USE.NAMES = F)


#####################
# dependencies
######################

compute_dependencies <- function(J, tol=1e-14){
  # Get attribute names
  attribute_names <- names(dimnames(J))
  # If only one attribute, no dependence
  if (length(attribute_names) == 1){
    dependency_list <- list()
    dependency_list[[attribute_names[1]]] <- c()
    return(dependency_list)
  }else{
    # Create list for each attribute_name
    dependency_list <- list()
    for (k in attribute_names){
      dependency_list[[k]] <- c()
    }
    # Loop over each pair of attributes - figure out if they're independent
    for (i in 1:(length(attribute_names)-1)){
      for(j in (i+1):length(attribute_names)){
        attr1 <- attribute_names[i]
        attr2 <- attribute_names[j]
        cross_tab <- apply(J, c(attr1,attr2), sum)
        # Standardize
        sums <- apply(cross_tab, 1, sum)
        cross_tab_std <- cross_tab/sums
        # Compute similarities
        is_equal = TRUE
        r_1 <- cross_tab_std[1,]
        if (nrow(cross_tab_std) > 1){
          for (m in 2:nrow(cross_tab_std)){
            if (any(as.vector(r_1) - as.vector(cross_tab_std[m,]) > tol)){
              is_equal <- FALSE 
            }
          }
        }
        
        # If not the same, append to dependency dictionary
        if (!is_equal){
          dependency_list[[attr1]] <- c(dependency_list[[attr1]], attr2)
          dependency_list[[attr2]] <- c(dependency_list[[attr2]], attr1)
        } 
      }
    }
    return(dependency_list)
  }
}

#########################
## amce.lmer function
#########################

amce.lmer <- function(formula, hierarchy, data, design="uniform", 
                      respondent.varying = NULL, subset=NULL, 
                      respondent.id=NULL, na.ignore=FALSE, 
                      baselines = NULL) {
  
  ###### Formula and Variables
  
  #we will split the formula into separate lists
  #unique_vars = all VARIABLES in formula
  #respondent_vars = VARIABLES varying by respondent
  #profile_vars = VARIABLES varying by profile
  #orig_effects = all EFFECTS in formula 
  #profile_effects = profile varying EFFECTS in formula
  #user inputted names for above end with "_user"
  
  ##### Parse formula, clean variable and input names
  
  formula_user <- formula
  #all variables in formula
  formula_char_user <- all.vars(formula)
  #lists with original names of variables and levels
  user_names <- list()
  user_levels <- list()
  for (char in formula_char_user) {
    user_names[[clean.names(char)]] <- char
    if (class(data[[char]]) == "factor") {
      old_names <- names(user_levels)
      user_levels <- c(user_levels,levels(data[[char]]))
      new_names <- sapply(clean.names(levels(data[[char]])),function(x) paste(clean.names(char),x,sep=""))
      names(user_levels) <- c(old_names,new_names)
    }
  }
  
  #make sure no duplicates after spaces and special characters removed
  formula_char <- clean.names(formula_char_user)    
  #if this makes for non-unique names, stop
  if(length(unique(formula_char)) != length(formula_char)) {
    stop("Error: Variable names must be unique when whitespace and meta-characters are removed. Please rename.")
  }
  
  #separate dependent and independent variables and clean
  y_var <- clean.names(formula_char_user[1])
  #identify ALL original effects; will add in missing base terms automatically
  orig_effects <- clean.names(attr(terms(formula_user),"term.labels"))
  #formula sorting part I: sort non-interaction terms and put them first
  orig_effects <- c(sort(orig_effects[!grepl(":",orig_effects)]), orig_effects[grepl(":",orig_effects)])
  vars_plus <- paste(orig_effects,collapse = " + ")
  form <- formula(paste(c(y_var,vars_plus),collapse = "~"))
  orig_effects <- attr(terms(form),"term.labels")
  
  #find missing base terms
  full_terms <- attr(terms(formula(paste(y_var,paste(sapply(orig_effects,function(x) gsub(":","*",x)),collapse=" + "),sep=" ~ "))),"term.labels")
  # add in any missing base terms for interactions
  missing_terms <- full_terms[!is.element(full_terms,orig_effects)]
  if (length(missing_terms > 0)) {
    orig_effects <- c(orig_effects,missing_terms)
    warning("Missing base terms for interactions added to formula")
  }
  
  #formula sorting redux: sort non-interaction terms and put them first
  orig_effects <- c(sort(orig_effects[!grepl(":",orig_effects)]), orig_effects[grepl(":",orig_effects)])
  #combine with "+"
  vars_plus <- paste(orig_effects,collapse = " + ")
  #then remake formula 
  form <- formula(paste(c(y_var,vars_plus),collapse = "~"))
  orig_effects <- clean.names(attr(terms(form),"term.labels"))
  
  #unique variables only (no interactions)
  unique_vars <- clean.names(rownames(attr(terms(form),"factor"))[-1])
  #respondent variables
  respondent_vars <- clean.names(respondent.varying)
  #profile variables
  profile_vars <- unique_vars[!is.element(unique_vars,respondent_vars)]
  
  #identify the REQUESTED profile effects and respondent effects (if any)
  if (length(respondent_vars) > 0) {
    #identify profile only effects
    profile_effects <- unlist(sapply(orig_effects,USE.NAMES = F,function(x) {
      y <- strsplit(x,":")[[1]]
      if (!any(is.element(y,respondent_vars))) x
    }))
    #terms containing a respondent var
    resp_only <- unlist(sapply(orig_effects,USE.NAMES = F, function(x) {
      y <- strsplit(x,":")[[1]]
      if(any(is.element(y,respondent_vars))) x
    }))
    #things that respondent vary is interacted with
    resp_mod <- unlist(sapply(resp_only,USE.NAMES = F,function(x) {
      y <- strsplit(x,":")[[1]]
      vars <- y[!is.element(y,respondent_vars)]
      if (length(vars) > 0) paste(vars,collapse = ":")
    }))
    resp_effects <- c(resp_mod,resp_only)  
  } else {
    profile_effects <- orig_effects
    resp_effects <- NULL
  }
  
  ### Extra name cleaning
  #cleaning additional inputs
  if (!is.null(respondent.id)) respondent.id <- clean.names(respondent.id)
  if (!is.null(baselines)) {
    names(baselines) <- clean.names(names(baselines))
    baselines <- lapply(baselines,function(x) clean.names(x))
  }
  
  #cleaning within data 
  colnames(data) <- clean.names(colnames(data))
  ## data <- data.frame(data) #in case of dplyr etc.
  var_to_check_levels <- unique(c(formula_char_user, respondent_vars, profile_vars, profile_effects, orig_effects))
  
  for (var in var_to_check_levels) {
    if (class(data[[var]]) == "factor") {
      clean.labels <- clean.names(levels(data[[var]]))
      if (length(unique(clean.labels)) != length(clean.labels)) {
        stop (paste("Error: levels of variable", var, "are not unique when whitespace and meta-characters are removed. Please rename."))
      }
      data[[var]] <- factor(data[[var]],levels=levels(data[[var]]), labels=clean.names(levels(data[[var]])))
    }
  }
  
  #######  Sanity Checks Re: Data
  
  # Are variables in data?
  for(var in formula_char) {
    if(!(var %in% colnames(data))) {
      stop(paste("Error:", var, "not in 'data'"))
    }
  }
  
  # Make sure non-respondent varying are factors
  for (var in profile_vars) {
    if (class(data[[var]]) != "factor") {
      data[[var]] <- as.factor(data[[var]])
      warning(paste(c("Warning: ",var," changed to factor"),collapse=""))
    }
  }
  
  # Is there missing data?
  if(na.ignore == FALSE){
    for(variab in formula_char){
      if (sum(is.na(data[[variab]])) != 0 ){
        stop(paste("Error:", variab, "has missing values in 'data'"))
      }
    }
  }
  
  # Is the respondent varying characteristic even in the formula obj?
  if (!is.null(respondent_vars)) {
    for (var in respondent_vars) {
      found <- 0
      for (formulavars in formula_char) {
        if (var == formulavars) {
          found <- 1
        }
      }
      if (found == 0) {
        stop(paste("Error:", var, "is specified in respondent.varying, but is not in the formula"))
      }
    }
  }
  
  # Check whether outcome variable is a binary 0-1 or numeric
  if (!is.numeric(data[[y_var]]) & !is.integer(data[[y_var]])) {
    stop(paste("Error:", y_var, "is not numeric or integer"))
  }
  
  # Are the user-supplied desired baselines in the data?
  if (!is.null(baselines)) {
    for(var in names(baselines)) {
      if(!(baselines[[var]] %in% data[[var]])) {
        stop(paste("Error: user supplied baseline",baselines[[var]],"is not a level of",var))
      }
    }      
  }
  
  ##### Sanity Checks re: design matrix
  
  # If design is already conjointDesign object, proceed to relevant sanity checks
  if (class(design) == "conjointDesign") {
    # Remove whitespaces etc from dimension names of design array 
    names(dimnames(design$J)) <- clean.names(names(dimnames(design$J)))
    dimnames(design$J) <- lapply(dimnames(design$J),function(x) clean.names(x))
    
    #and design dependencies
    names(design$depend) <- clean.names(names(design$depend))
    design$depend <- lapply(design$depend,function(x) clean.names(x))  
    #Now check to make sure profile varying attributes are in conjointDesign
    for (eff in profile_vars) {   
      if (!(eff %in% names(dimnames(design$J)))) {
        stop(paste("Error:", eff, "not in 'design' object"))
      }
    }      
    #Check to make sure conjointDesign attributes are in data and level names match
    for (eff in names(dimnames(design$J))) {
      if (!(eff %in% colnames(data))){
        stop(paste("Error: attribute", eff, "in 'design' object is not in 'data'"))
      } else {
        # Check all level names for the attribute in dataset appear in design
        for (lev in clean.names(levels(as.factor(data[[eff]])))) {
          if (!(lev %in% dimnames(design$J)[[eff]])) {
            #print(paste0('checking level ', lev))
            stop(paste("Error: factor level", lev, "of attribute", eff, "not in 'design' object"))
          }
        }
      }
    }
    
    depend_to_check <- NULL
    for (var in names(design$depend)) {
      depend_to_check <- c(depend_to_check,design$depend[[var]])
    }
    ## only check any variables that weren't already checked
    depend_to_check <- depend_to_check[!depend_to_check %in% var_to_check_levels]
    ## Check that all dependencies have unique levels
    for (var in depend_to_check){
      if (class(data[[var]]) == "factor") {
        clean.labels <- clean.names(levels(data[[var]]))
        if (length(unique(clean.labels)) != length(clean.labels)) {
          stop (paste("Error: levels of variable", var, "used as a dependency, are not unique when whitespace and meta-characters are removed. Please rename."))
        }
        data[[var]] <- factor(data[[var]],levels=levels(data[[var]]), labels=clean.names(levels(data[[var]])))
      }
    }
  } else if (design == "uniform") {    
    # else if design == "uniform", create J-dimensional array 
    design <- list()        
    # Determine dimensions
    # And create J matrix with uniform probabilities across all vars
    design.dim <- vector(length=length(profile_vars))
    dim_list <- list()
    for (i in 1:length(profile_vars)) {
      dim_list[[i]] <- levels(factor(data[[profile_vars[i]]]))
      design.dim[i] <- length(dim_list[[i]])
    }
    names(dim_list) <- profile_vars
    design$J <- array(1/prod(design.dim), dim=design.dim, dimnames=dim_list)
    design$depend <- compute_dependencies(design$J)       
  } else {
    #if neither uniform nor conjointDesign, error
    stop("Design object must be a valid character string 'uniform' or a conjointDesign object")   
  }
  
  ####### Subsetting data    
  
  if (is.null(subset)) {
    data <- data 
  } else {
    if (class(subset) == "logical") {
      if (length(subset) == nrow(data)) {
        data <- subset(data, subset) 
      } else {
        warning("Warning: invalid argument to 'subset' - must be the same length as the number of rows in data")
      }
    } else {
      warning("Warning: invalid argument to 'subset' - must be a logical")
    }
  }
  
  ###### Adjust baselines if given    
  
  if (!is.null(baselines)) {
    for (var in names(baselines)) {
      data[[var]] <- factor(data[[var]])
      data[[var]] <- relevel(data[[var]], baselines[[var]])
    } 
  }
  
  ####### Adding relevant interaction terms to model
  
  #If there are any dependencies-- only for profile-varying!
  if(any(profile_vars %in% names(design$depend))) {
    #initialize full interaction set
    depend_vars <- c()     
    #loop over effects with dependencies
    for(eff in profile_vars[profile_vars %in% names(design$depend)]) {
      #initialize interaction set for given variable
      inter <- c()
      #identify higher order occurences of variable
      #make sure it's just that variable, not followed or begun by "_"
      eff_all <- grep(paste(c(":",eff,"(?!_)","|",eff,":(?<!_)"),collapse=""),
                      orig_effects,value=T,perl=T)
      #if you find some, break up, sort and replace ":" with "*"
      if (length(eff_all) > 0) {
        eff_all <- sapply(strsplit(eff_all,":"),function(x) paste(sort(x),collapse="*"))
      }
      #combine with lower order (main effect)
      eff_all <- c(eff,eff_all)
      #for each occurrence, create interaction
      inter <- sapply(eff_all,USE.NAMES = F,function(x) {
        #get conditioning set
        T_r <- design$depend[[eff]]
        #make factors
        for (t in T_r){
          data[[t]] <- as.factor(data[[t]])
        }
        #combine name and dependency
        T_r_d <- c(x,T_r)
        #make interaction term
        paste(T_r_d,collapse="*")
      })
      #add to list
      depend_vars <- c(depend_vars,inter)
    }      
    #drop repeats
    depend_vars <- unique(depend_vars)
    #add to formula
    form_full <- formula(paste(c(form,depend_vars),collapse = " + "))
  } else {
    form_full <- form
  }
  
  #all variables to be run
  all_run_vars <- attr(terms(form_full),"term.labels")
  #formula sorting redux: sort non-interaction terms and put them first
  all_run_vars <- c(sort(all_run_vars[!grepl(":",all_run_vars)]), all_run_vars[grepl(":",all_run_vars)])
  #combine with "+"
  vars_plus <- paste(all_run_vars,collapse = " + ")
  #then remake formula 
  form_full_nohierarchy <- paste(c(y_var,vars_plus),collapse = "~")
  form_full <- formula(paste(form_full_nohierarchy, " + ", hierarchy, sep=""))
  all_run_vars <- attr(terms(formula(form_full_nohierarchy)),"term.labels")
  
  ####### If there are respondent varying terms, split into two formulas
  ######## One contains only profile effects
  ######## Second is full formula
  
  if (length(respondent_vars) > 0) {
    ### profile only formula
    #remove those involving respondent things
    prof_only <- unlist(sapply(all_run_vars,function(x) {
      y <- clean.names(strsplit(x,":")[[1]])
      if(!any(is.element(y,respondent_vars))) x
    }))
    prof_only_plus <- paste(prof_only,collapse = " + ")
    #formula with profile only
    form_prof_nohierarchy <- paste(all.vars(form_full)[1],prof_only_plus,sep=" ~ ")
    form_prof <- formula(paste(form_prof_nohierarchy, " + ", hierarchy, sep=""))
  } else {
    #otherwise use full formula
    form_prof <- form_full
  }
  all_prof <- clean.names(attr(terms(formula(form_prof_nohierarchy)),"term.labels"))
  all_run_vars <- clean.names(attr(terms(formula(form_full_nohierarchy)),"term.labels"))
  if (any(!is.element(all_prof,all_run_vars))) {
    warning("Warning: mismatch of term names between full formula and profile formula")
  }
  
  ####### Running Model
  
  lin.mod.prof <- lme4::lmer(form_prof, data=data)
  if (length(respondent_vars) > 0) {
    lin.mod.full <- lme4::lmer(form_full, data=data)
  } else {
    lin.mod.full <- NULL
  }
  
  #Get sample size
  sample_size_prof <- length(summary(lin.mod.prof)$residuals)
  if (length(respondent.varying) > 0) {
    sample_size_full <- length(summary(lin.mod.full)$residuals)
  } else {
    sample_size_full <- NULL
  }
  
  #Get coefficients
  coeffs.prof <- fixef(lin.mod.prof)  
  if (length(respondent.varying) > 0) {
    coeffs.full <- fixef(lin.mod.full)  
  } else {
    coeffs.full <- NULL
  }
  
  #Compute vcov of LMER
  vcov_mat_prof <- vcov(lin.mod.prof)
  if (length(respondent.varying) > 0) {
    vcov_mat_full <- vcov(lin.mod.full)
  } else {
    vcov_mat_full <- NULL
  }  
  ## function for fixing variance-covariance matrices
  ## input in vcov output from lm and a matrix of varprobs
  ## draws on matrix package
  fix.vcov <- function(varprob,vcov) {
    if (!requireNamespace("Matrix", quietly = TRUE)){
      stop("Matrix package needed for this function to work. Please install it.",
           call. = FALSE)
    }
    #designate inputs as sparse
    varprob2 <- Matrix(varprob,sparse=TRUE)
    vcov2 <- Matrix(vcov,sparse=TRUE)
    #calculate and add in single sum corrections
    fix1 <- varprob2 %*% vcov2 + t(varprob2 %*% vcov2) + vcov2
    #vcov matrix as vector
    vcov_vector <- matrix(as.vector(vcov2),ncol=1,nrow=length(as.vector(vcov2)))
    #use to multiply combinations of varprobs by covariances and sum
    weighted_covs <- kronecker(varprob2,varprob2) %*% vcov_vector
    #make corrections into a matrix
    weighted_covs <- matrix(weighted_covs,nrow=nrow(vcov), ncol=ncol(vcov))
    #add to single sum corrections
    fix2 <- fix1 + weighted_covs
    #return as normal matrix
    out <- matrix(fix2,ncol=ncol(vcov),nrow=nrow(vcov))
    colnames(out) <- rownames(out) <- rownames(vcov)
    return(out)
  }
  
  ######### Extract Effects from the profile-vars only linear model
  
  # proposed nomenclature here:
  # effect = attribute in question, which has "effect levels"
  # depends = attributes it depends on, each of which has "depend levels"
  
  #Make R CMD check happy
  J_baseline <- NULL
  J_effect <- NULL
  
  #before we start, make a blank call to the design array J
  J_call <- Quote(design$J[])
  J_call <- J_call[c(1, 2, rep(3, length(dim(design$J))))]
  
  #warnings counter
  warn_i <- 0
  
  ############## loop over unique profile vars only (AMCE and ACIE); interactions etc. below
  
  #blank list for output
  estimates <- list()
  #re-sort profile effects
  profile_effects <- c(sort(profile_effects[!grepl(":",profile_effects)]), profile_effects[grepl(":",profile_effects)])
  #initialize list for weighted cross-terms
  covariance_list <- list()
  #blank matrix of var probs
  varprob_mat <- matrix(0,nrow(vcov_mat_prof),ncol(vcov_mat_prof))
  colnames(varprob_mat) <- rownames(varprob_mat) <- colnames(vcov_mat_prof)
  
  for(i in 1:length(profile_effects)) {
    
    #split into sections if it's an interaction
    substrings <- strsplit(profile_effects[i], "[:*]", perl=TRUE)[[1]]
    
    #administrative loop to find levels
    all_levels <- list()
    all_levels_coefs <- list()
    for(effect in substrings) {
      #get all level names and coefficient names-- sans baseline!!!
      all_levels[[effect]] <- levels(data[[effect]])[-1]
      all_levels_coefs[[effect]] <- sapply(all_levels[[effect]], function(x) {
        paste(c(effect,x), collapse="")
      })
    }
    
    #find all combinations of level names-- add as FIRST column
    levels <- expand.grid(all_levels, stringsAsFactors = FALSE)
    #make level combos in first column
    levels <- cbind(apply(levels,1,function(x) paste(x,collapse=":")),levels)
    colnames(levels) <- c("name",substrings)
    
    #and all combinations of actual coefficient names
    coefs <- expand.grid(all_levels_coefs, stringsAsFactors = FALSE)
    coefs <- apply(coefs,1,function(x) paste(x,collapse=":"))
    
    # Initialize the results 
    results <- matrix(nrow=2, ncol = nrow(levels))
    if (length(substrings) > 1) {
      rownames(results) <- c("ACIE", "Std. Error")
    } else {
      rownames(results) <- c("AMCE", "Std. Error")
    }
    colnames(results) <- coefs
    results[2,] <- NA
    
    #### find extra times when this effect is mentioned
    all_depends <- unlist(sapply(all_prof,USE.NAMES = F,function(x) {
      y <- strsplit(x,":")[[1]]
      if (all(is.element(substrings,y))) x
    }))
    # remove the actual term
    all_depends <- all_depends[-is.element(all_depends,profile_effects[i])]
    
    #### loop over every combination of levels of component effects
    for(j in 1:nrow(levels)) {
      
      #figure out which level of inter we're doing
      ## effect_level <- as.character(levels[j,1])
      effect_level_coef <- coefs[j]
      #get its beta 
      initial_beta <- coeffs.prof[effect_level_coef]
      
      #if interaction,make sure there is baseline support for this level combination
      if (!is.na(initial_beta) & length(substrings) > 1) {
        for (effect1 in substrings) {
          #get effect base
          effect_base1 <- levels(data[[effect1]])[1]
          #subset data to that level
          base.subset <- data[which(data[[effect1]] == effect_base1),]                       
          #loop over other profile-varying vars in interaction to subset further
          for(effect in substrings[!(substrings %in% effect1)]) {
            base.subset <- base.subset[which(base.subset[[effect]] == as.character(levels[j,effect])),]
          }                        
          #if there's no support left, change beta and var to NA
          if (nrow(base.subset) == 0) {
            initial_beta <- NA
            #and give a warning that you had to do it
            warn_i <- warn_i + 1
          }
        }
      }
      
      # If initial_beta and initial_variance are not NA (are valid level combination)
      # and there are dependent variables to incorporate
      if (!is.na(initial_beta) & length(all_depends) > 0) {
        
        #get the slice of design array J associated with baseline and inter level
        J_effect_call <- J_base_call <-  J_call
        for(effect in substrings) {
          #identify its baseline and modify baseline call accordingly
          base <- levels(data[[effect]])[1]
          effect_index <- which(names(dimnames(design$J)) == effect)
          J_base_call[effect_index + 2] <- base                            
          #identify level of each effect and modify inter call accordingly
          level <- levels[j,effect]
          J_effect_call[effect_index + 2] <- level
        }
        eval(call("<-", Quote(J_baseline), J_base_call))
        eval(call("<-", Quote(J_effect), J_effect_call))
        
        #### loop over dependencies for all components of interaction
        for(k in 1:length(all_depends)) {
          
          #attribute effect is dependent on
          depend <- all_depends[[k]]
          #figure out what levels of what variables are involved
          substrings_d <- strsplit(depend,":")[[1]]
          substrings_d <- substrings_d[!is.element(substrings_d,substrings)]
          all_depend_coefs <- list()
          for (sub in substrings_d) {
            all_depend_coefs[[sub]] <-  sapply(levels(data[[sub]]), function(x) paste(c(sub,x),collapse=""))
          }
          all_depend_levels <- expand.grid(all_depend_coefs)
          substrings_l <- strsplit(effect_level_coef,":")[[1]]
          for (l in length(substrings_l):1) {
            all_depend_levels <- cbind(rep(substrings_l[l], nrow(all_depend_levels)), all_depend_levels)
          }
          colnames(all_depend_levels)[1:length(substrings_l)] <- substrings
          all_depend_levels <- all_depend_levels[sort(colnames(all_depend_levels))]
          all_depend_level_coefs <- apply(all_depend_levels, 1, function(x) paste(x,collapse=":"))
          
          #baseline support for depend attribute level 
          if (!(is.null(dim(J_baseline)))){
            baseline_support <- apply(J_baseline,substrings_d,sum)
          } else {
            baseline_support <- J_baseline
          }
          baseline_support[baseline_support != 0] <- 1
          
          #probs for depend attribute levels WITH baseline support
          if (!is.null(dim(J_effect))) {
            joint_prob <- apply(J_effect, substrings_d, sum)*baseline_support
          } else {
            joint_prob <- J_effect*baseline_support
          }
          #make it a vector
          joint_prob <- as.vector(joint_prob)
          names(joint_prob) <- all_depend_level_coefs
          
          all_depend_level_coefs <- all_depend_level_coefs[!is.na(coeffs.prof[all_depend_level_coefs])]
          varprob_mat[effect_level_coef,all_depend_level_coefs] <- as.numeric(joint_prob[all_depend_level_coefs])/as.numeric(sum(joint_prob))
          
          ##### if all_depend_level_coefs is 1 or longer 
          if (length(all_depend_level_coefs)) {
            #calculate probabilities for this effect and depend level 
            var_prob <- joint_prob[all_depend_level_coefs]
            var_prob <- as.numeric(var_prob)/as.numeric(sum(joint_prob))
            # add weighted beta to initial_beta
            depend_betas <- fixef(lin.mod.prof)[all_depend_level_coefs]
            initial_beta <- sum(initial_beta,var_prob*depend_betas,na.rm=T)
          }
          
        } #end for loop over different dependent attributes
      } #end if has valid beta, var, dependencies
      
      # Store effect and standard error estimates
      results[1,j] <- initial_beta
      
    } #end for loop over all level combinations
    
    # combine estimates + SEs into single matrix - store in list
    estimates[[profile_effects[i]]] <- results
    
  } #end for loop over profile effects      
  
  ### fix var-cov matrix
  vcov_prof <- suppressMessages(fix.vcov(varprob_mat,vcov_mat_prof))
  #write in adjusted variances to output matrix
  for (i in 1:length(estimates)) {
    coef_names <- colnames(estimates[[i]])
    variances <- sqrt(diag(vcov_prof)[coef_names])
    estimates[[i]][2,] <- ifelse(is.na(estimates[[i]][1,]),NA,variances)
  }
  
  #determine term names for profile effects (no depends) to keep
  profile_effects_plus <- paste(profile_effects,collapse=" + ")
  profile_effects_form <- formula(paste(c(y_var,profile_effects_plus),collapse = " ~ "))   
  profile_effects_terms <- colnames(model.matrix(profile_effects_form,data))   
  profile_effects_terms <- profile_effects_terms[profile_effects_terms %in% colnames(vcov_mat_prof)]
  vcov_prof <- vcov_prof[profile_effects_terms,profile_effects_terms]
  
  ######### Extract Effects from the full model (if have respondent interactions)
  
  # proposed nomenclature here:
  # effect = attribute in question, which has "effect levels"
  # depends = attributes it depends on, each of which has "depend levels"
  # inters = attributes in interaction terms each of which has "inter levels"
  
  #if there are any respondent effects
  if (length(respondent_vars) > 0) {
    
    #blank list for output
    conditional.estimates <- list()
    #re-sort respondent effects
    resp_effects <- c(sort(resp_effects[!grepl(":",resp_effects)]), resp_effects[grepl(":",resp_effects)])
    #initialize list for weighted cross-terms
    covariance_list <- list()
    #blank matrix for var probs
    varprob_mat <- matrix(0,nrow(vcov_mat_full),ncol(vcov_mat_full))
    colnames(varprob_mat) <- rownames(varprob_mat) <- colnames(vcov_mat_full)
    
    #loop over respondent-related effects
    for (i in 1:length(resp_effects)) {
      
      #split into component effects, if interaction
      substrings <- strsplit(resp_effects[i], "[:*]", perl=TRUE)[[1]]
      
      ## start by finding levels
      #administrative loop over components of interaction
      all_levels <- list()
      all_levels_coefs <- list()
      for(effect in substrings) {
        #if it's not a factor, only has the 1 "level" and coefficient name stays
        if (class(data[[effect]]) != "factor") {
          all_levels[[effect]] <- effect
          all_levels_coefs[[effect]] <- effect
        } else {
          #if it is a factor, get all level names and coefficient names-- sans baseline!!!
          all_levels[[effect]] <- levels(data[[effect]])[-1]
          all_levels_coefs[[effect]] <- sapply(all_levels[[effect]],
                                               function(x) {
                                                 paste(c(effect,x), collapse="")
                                               })
        }
      }
      
      #find all combinations of level names-- add as FIRST column
      levels <- expand.grid(all_levels, stringsAsFactors = FALSE)
      levels <- cbind(apply(levels,1,function(x) paste(x,collapse=":")),levels)
      colnames(levels) <- c("name",substrings)
      
      #and all combinations of coefficient names
      coefs <- expand.grid(all_levels_coefs, stringsAsFactors = FALSE)
      coefs <- apply(coefs,1,function(x) paste(x,collapse=":"))
      
      # Initialize the results 
      results <- matrix(nrow=2, ncol = nrow(levels))
      rownames(results) <- c("Conditional Estimate", "Std. Error")
      colnames(results) <- coefs
      #write NA to SE row to start
      results[2,] <- NA
      
      #### find extra times when this effect is mentioned in full formula
      # only if anything related to profile var is involved
      if (any(substrings %in% profile_vars)) {
        all_depends <- unlist(sapply(all_run_vars,USE.NAMES = F,function(x) {
          y <- strsplit(x,":")[[1]]
          if (all(is.element(substrings,y))) x
        }))
        # remove the actual term
        all_depends <- all_depends[!is.element(all_depends,resp_effects[i])]
        # remove any that involve any other respondent varying terms
        resp.other <- respondent_vars[!is.element(respondent_vars,substrings)]
        all_depends <- unlist(sapply(all_depends,function(x) {   
          sub_depends <- strsplit(x,":")[[1]]
          if (all(!is.element(sub_depends,resp.other))) x
        }))
      } else {
        #no profile vars, no depends
        all_depends <- c()
      }
      
      #### loop over every combination of levels of component effects
      for(j in 1:nrow(levels)) {
        
        #figure out which level of inter we're doing
        effect_level_coef <- coefs[j]
        #get its beta 
        initial_beta <-coeffs.full[effect_level_coef]
        
        #make sure there is baseline support for this level combination
        if (!is.na(initial_beta)) {                   
          for (effect1 in substrings[substrings %in% profile_vars]) {
            #get effect base
            effect_base1 <- levels(data[[effect1]])[1]
            #subset data to that level
            base.subset <- data[which(data[[effect1]] == effect_base1),]                       
            #loop over other profile-varying vars in interaction to subset further
            for(effect in substrings[substrings %in% profile_vars][!(substrings[substrings %in% profile_vars] %in% effect1)]) {
              base.subset <- base.subset[which(base.subset[[effect]] == as.character(levels[j,effect])),]
            }                        
            #if there's no support left, change beta and var to NA
            if (nrow(base.subset) == 0) {
              initial_beta <- NA
              #and give a warning that you had to do it
              warn_i <- warn_i + 1
            }
          }
        }
        
        # If initial_beta and initial_variance are not NA and there are depends
        # proceed to add to beta and var
        if (!is.na(initial_beta) & length(all_depends) > 0) {
          
          #get the slice of design array J associated with baseline and inter level
          #profile variables only!
          J_effect_call <- J_base_call <-  J_call
          for(effect in substrings[substrings %in% profile_vars]) {
            #identify its baseline and modify baseline call accordingly
            base <- levels(data[[effect]])[1]
            effect_index <- which(names(dimnames(design$J)) == effect)
            J_base_call[effect_index + 2] <- base                           
            #identify level of each effect and modify inter call accordingly
            level <- levels[j,effect]
            J_effect_call[effect_index + 2] <- level
          }
          eval(call("<-", Quote(J_baseline), J_base_call))
          eval(call("<-", Quote(J_effect), J_effect_call))
          
          #### loop over dependencies for all components of effect
          for(k in 1:length(all_depends)) {
            
            #attribute effect is dependent on
            depend <- all_depends[[k]]
            #figure out what levels of what variables are involved
            substrings_d <- strsplit(depend,":")[[1]]
            substrings_d <- substrings_d[!is.element(substrings_d,substrings)]
            all_depend_coefs <- list()
            for (sub in substrings_d) {
              all_depend_coefs[[sub]] <-  sapply(levels(data[[sub]]), function(x) paste(c(sub,x),collapse=""))
            }
            all_depend_levels <- expand.grid(all_depend_coefs)
            substrings_l <- strsplit(effect_level_coef,":")[[1]]
            for (l in length(substrings_l):1) {
              all_depend_levels <- cbind(rep(substrings_l[l], nrow(all_depend_levels)), all_depend_levels)
            }
            colnames(all_depend_levels)[1:length(substrings_l)] <- substrings
            ####put terms together in proper order
            all_depend_levels <- all_depend_levels[sort(colnames(all_depend_levels))]
            all_depend_level_coefs <- apply(all_depend_levels, 1, function(x) paste(x,collapse=":"))
            
            #baseline support for depend attribute level in inter
            if (!(is.null(dim(J_baseline)))) {
              baseline_support <- apply(J_baseline,substrings_d,sum)
            } else {
              baseline_support <- J_baseline
            }
            baseline_support[baseline_support != 0] <- 1
            
            #support for depend attribute levels WITH baseline support
            if (!is.null(dim(J_effect))) {
              joint_prob <- apply(J_effect, substrings_d, sum)*baseline_support
            } else {
              joint_prob <- J_effect*baseline_support
            }
            #make it a vector
            joint_prob <- as.vector(joint_prob)
            names(joint_prob) <- all_depend_level_coefs
            
            all_depend_level_coefs <- all_depend_level_coefs[!is.na(fixef(lin.mod.full)[all_depend_level_coefs])]
            varprob_mat[effect_level_coef,all_depend_level_coefs] <- as.numeric(joint_prob[all_depend_level_coefs])/as.numeric(sum(joint_prob))
            
            ## If there are non-null # of depend-level-coefs
            if (length(all_depend_level_coefs)) {
              #calculate probabilities for this effect and depend level 
              var_prob <- joint_prob[all_depend_level_coefs]
              var_prob <- as.numeric(var_prob)/as.numeric(sum(joint_prob))
              # add weighted beta to initial_beta
              depend_betas <- fixef(lin.mod.full)[all_depend_level_coefs]
              initial_beta <- sum(initial_beta,var_prob*depend_betas, na.rm=T)
            }        
            
          } #end for loop over different dependent attributes
          
        } #end if initial beta and var are NA, has depends
        
        # Store effect estimates
        results[1,j] <- initial_beta
        
      } #end for loop over all level combinations
      
      # combine estimates + SEs into single matrix - store in list
      conditional.estimates[[resp_effects[i]]] <- results
      
    } #end for loop over respondent related effects      
    
    ##fix variance-covariance matrix
    vcov_resp <- suppressMessages(fix.vcov(varprob_mat,vcov_mat_full))
    #write in adjusted variances
    for (i in 1:length(conditional.estimates)) {
      coef_names <- colnames(conditional.estimates[[i]])
      variances <- sqrt(diag(vcov_resp)[coef_names])
      conditional.estimates[[i]][2,] <- ifelse(is.na(conditional.estimates[[i]][1,]), NA, variances)
    }
    
    ###terms to keep (no depends)
    resp_effects_plus <- paste(resp_effects,collapse=" + ")
    resp_effects_form <- formula(paste(c(y_var,resp_effects_plus),collapse = " ~ "))
    resp_effects_terms <- colnames(model.matrix(resp_effects_form,data))
    resp_effects_terms <- resp_effects_terms[resp_effects_terms %in% colnames(vcov_mat_full)]
    vcov_resp <- vcov_resp[resp_effects_terms,resp_effects_terms]
    
  } #end if there are any respondent related effects
  
  ############  create conjoint object for output
  
  output <- list()
  class(output) <- c("amce")
  
  #saving things for unconditional estimates
  output$estimates <- estimates
  #saving profile attributes
  output$attributes <- dimnames(design$J)
  #save original profile-only vcov matrix
  output$vcov.prof <- vcov_prof
  #save sample size used for unconditional estimates
  output$samplesize_prof <- sample_size_prof
  #save style edited formula (no depends)
  output$formula <- form
  
  # save full model
  output$lmer.prof <- lin.mod.prof
  output$lmer.full <- lin.mod.full
  
  #final warning tally
  if (warn_i > 0) {
    warning(paste("Warning: ",warn_i," interaction levels lacked support at baseline, effects undefined unless alternative baseline is provided."))
  }
  
  #saving things for conditional estimates
  if (length(respondent.varying) > 0) {     
    output$cond.estimates <- conditional.estimates
    output$vcov.resp <- vcov_resp
    output$samplesize_full <- sample_size_full
    #save style edited formula (no depends), only resp-related
    output$cond.formula <- resp_effects_form
  }
  
  # Save baselines of unique (main) effects (if factor) to "baselines"
  # If continuous save summary information to "continuous"
  output$baselines <- list()
  output$continuous <- list()
  for (k in unique_vars) {
    if (class(data[[k]]) == "factor") {
      output$baselines[[k]] <- levels(data[[k]])[1]
    } else if (class(data[[k]]) == "numeric") {
      output$continuous[[k]] <- quantile(model.matrix(form,data)[,k], probs=c(0.25,0.5,0.75), na.rm=T)
    }
  }
  
  #save number of respondents if ID given
  #save number of respondents if ID given
  if (!is.null(respondent.id)) {
    output$numrespondents <- length(unique(data[[respondent.id]]))
  } else {
    output$numrespondents <- NULL
  }
  
  #save respondent variables if given
  if (!is.null(respondent.varying)) {
    output$respondent.varying <- respondent_vars
  } else {
    output$respondent.varying <- NULL
  }
  
  #save original names
  output$user.names <- user_names
  output$user.levels <- user_levels
  #save the original data
  output$data <- data
  return(output)
}


