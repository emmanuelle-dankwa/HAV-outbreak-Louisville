######### ALL FUNCTIONS USED IN MODELLING ANALYSIS #############
# Author: EMMANUELLE A. DANKWA (email:dankwa@stats.ox.ac.uk)
# IN THIS SCRIPT:
#   * HAV TRANSMISSION MODEL
#   * FUNCTIONS FOR PARAMETER ESTIMATION
#   *FUNCTION TO COMPUTE OUTBREAK TRAJECTORY
#   * MISCELLANEOUS
#   * CORRESPONDING FUNCTIONS FOR MODEL WITH SINGLE EFFECTIVE CONTRACT RATE (BETA)
###############################################################


##############################################################################
####  HAV TRANSMISSION MODEL (Equation S3 IN SUPPLEMENTARY MATERIAL)    ####
##############################################################################

#' @t times
#' @states Compartmental model states (S, L, Z, I & R)
#' @parameters model parameters
#' @fitlength The model is calibrated using case counts up to 'fitlength' week
#' @effecvanumperwk Effective number of vaccinations per week 
#' @vacc Vaccination at time, t
#' @beta Parameter: Force of infection, dynamic
#' @durL  Parameter: Duration of latency period
#' @durRel Parameter: Duration of relapse period
#' @durI Parameter: Duration of infectiousness
#' @propRelapse Parameter: Proportion of relapse period


HAVODELV2 <- function(t, states, parameters, fitlength, effecvacnumperwk){
  with(as.list(c(states,parameters)),
       {

         
         beta=a + (b - a)/(1+exp(-c*(t-d)))
        # print(beta)
         
         states[1] -> S
         states[2] -> L
         states[3] -> I
         states[4] -> Z
         states[5] -> R
         
         # sum to N
         N    <- S + L + Z + I + R
         
         
         # vaccination value
         if (t < fitlength){            
           vacc <- effecvacnumperwk[floor(t+1)]
         }else{                                   # after calibration assume same vaccination as previous week
           vacc <- effecvacnumperwk[fitlength]}
         if(vacc > S){vacc = S}
         
        # print(S/N)
         
         dS <- -(beta) * S * (I/N) - vacc*(S/N) # S/N is the variable part of theta 
         dL <-  beta * S * (I/N) - (1/durL)*L
         dI <-  (1/durL)*L +(1/durRem)*R - (1/durI)*I             
         dZ <- (1-propRelapse)*(1/durI)*I + vacc*(S/N)
         dR <- propRelapse*(1/durI)*I - (1/durRem)*R  # propRelapse = 1-eta (eta = permanent recoveries)
         
         return(list(c(dS, dL, dI, dZ, dR)))
         
       })  
}





##### FUNCTIONS FOR PARAMETER ESTIMATION #######

###  LIKELIHOOD FUNCTION  ### 

# Function to compute the log likelihood. Returns the log likelihood for a given 
# parameter set. 
#
# Inputs:
#' @x Vector of parameters to be estimated
#' @N Population size at each iteration
#' @ydata Bootstrap sample at each iteration
#' @vacnumperwk Vector of weekly vaccination counts
#' @fitlength The model is calibrated using case counts up to 'fitlength' week
#' @tau Vaccine efficacy rate


calib.function <- function(x, 
                           N, 
                           ydata,
                           vacnumperwk = vaccdata, 
                           fitlength = 69, 
                           tau=0.9,
                           fp = fixed_params){
  # As required by optimx, first argument of the optimization function is a vector of parameters of interest.
  x <- exp(x)
 
  # Parameters to be estimated
  a <- x[1]            # Effective contact rate at the start of the outbreak (β_s, in manuscript)
  b <- x[2]            # Effective contact rate later in the outbreak (β_l, in manuscript)
  c <- x[3]            # Speed of transition from a to b
  d <- x[4]            # Transition midpoint time from β_s to β_l (in weeks)
  I0 <- x[5]           # Number of infectious individuals in week 0                                             
  
  
  omega <- fp["omega", 1]     # Fraction of vaccine doses given to at-risk individuals. $\omega$ is assumed to be 1 (Table 1)

  
  # Initial conditions 
  L. = (fp["durI", 1]/fp["durL", 1])*I0     # set L_0 to be at equilibrium with I_0
  I. = I0                                   # Initial number of infectious individuals
  Z. = fp["fracImmune", 1]*N                # Number of immune individuals at the start of the outbreak
  R. = 0                                    # Number of individuals initially the remission state
  S. = N - L. - I. - Z. - R.                # Number of susceptible individuals
  
  # States
  states = c(S=S.,
             L=L.,
             I=I.,
             Z=Z.,
             R=R.)
  
  #print(paste("states", states))
  
  
  # Parameters
  parameters = c(list(durL = fp["durL", 1],
                      durI = fp["durI", 1],
                      durRem =fp["durRem", 1],
                      propRelapse =fp["propRelapse", 1],
                      a= a,
                      b = b,
                      c = c,
                      d = d))
  
  #print(parameters)
  # Calibrate up to:
  ydata <- ydata[1:fitlength]
  
  # Define time
  tspan <- seq(0, length(ydata), by=1)
  
  # Run the ODE solver
  theta_part <- vacnumperwk*tau*omega  # Constant part of theta, the effective number of vaccinations ($\theta(t)$) 
  
  
  #print(paste0("theta", theta_part))
  
  result_df <- data.frame(ode(y = states,
                              times = tspan,
                              func = HAVODELV2,
                              parms = parameters,
                              method = "bdf_d",
                              effecvacnumperwk = theta_part,
                              fitlength = fitlength
  ))
  
  # Calculate weekly **detected** incidence
  IncDetec <- ((1/fp["durL", 1])*result_df[, 3] + (1/fp["durRem", 1])*result_df[, 6])[-1]
  #print(IncDetec)
  # Calculate Poisson likelihood (**negative** log likelihood here, as optimx performs minimization)
  NLLK<- -sum((dpois(x = ydata , lambda = IncDetec, log = T))) 
  print(paste("Likelihood value:", NLLK))
  
  if(NLLK==-Inf |is.nan(NLLK)){# Return a large number if value is improbable
    NLLK = -9999
  }
  if(NLLK==Inf){
    NLLK = 9999
  }
  return(NLLK)
}




## COMPUTE THE PROFILE LIKELIHOOD CURVE FOR A PARAMETER ##

# Inputs:
#' @p.vec Parameter for which profile likelihood is computed. May be multiple parameters.
#' @pvl Parameter values over which to profile; a list, each component corresponding to a parameter. 
#' @... Other parameters to be passed to `estimate.parameters.pl` (see below)


pl_function <- function(p.vec,
                        pvl=profile_vals_list,
                        ...){
  
  V = length(pvl[[p.vec]]) # Number of values for which the profile likelihood are to be computed
  
  # Create empty list to save outputs
  master_list <- vector("list", length(p.vec))
  names(master_list) <- p.vec
  
  for(parameter in p.vec){
    # Define interval over which profile is to computed
    profile_param <- parameter
    
    profile_vals <- pvl[[profile_param]]
    
    print(profile_vals)
    
    # Create list to store values
    res <- list()                   # Parameter estimates
    
    for (iter in c(1:V)){
      print(paste0("iter: ", iter))
      
      if(iter==1){
        prev_ests_val <- NA
      }else{
        prev_ests_val <- res
        # print(paste0("prev_ests: ", prev_ests_val))
      }
      
      
      res[[iter]] <- estimate.parameters.pl(profile_est = profile_param,
                                            # pul2 =  prior_uncertainty_limits,
                                            # N.vector = N_vec,
                                            # bootstr = bootdata,
                                            # vacnumperwk = vaccdata,
                                            #fitlength = 69,
                                            #omega = 0.9,
                                            profile_lik = profile_param,
                                            profileint = profile_vals,
                                            roundx = iter,
                                            V_value = iter,
                                            prev_ests = prev_ests_val
                                            # pul1 = prior_uncertainty_limits
      )
      
      #print(res[[iter]]) # Monitor
      
    }
    
    # Save results to master list
    master_list[[parameter]]$res <- res
    master_list[[parameter]]$profile_vals <- profile_vals
    master_list[[parameter]]$name <- profile_param
    
  }
  return(master_list)
  
}






## PERFORM OPTIMIZATION AT EACH STEP OF PROFILE LIKELIHOOD COMPUTATION ##

# Inputs:
#' @profile_est Parameter for which profile likelihood is computed. May be multiple parameters.
#' @pul2 Data frame of values for other parameters.
#' @obs Observed case counts
#' @V_value Iteration number
#' @prev_ests Parameter estimates at previous optimization steps (These are used as starting points for the current optimization step)
#' @chain_algo Use the previous parameter estimates as starting points for the current optimization step


estimate.parameters.pl <- function(profile_est = profile_param,
                                   pul2 =  pul,
                                   obs = casecountperwk,
                                   V_value = 1,
                                   prev_ests,
                                   chain_algo = chain
                                   ,...){
  
  # Set initial, lower and upper values (See Table 1 in main text).
  # Initial values are sampled from a uniform distribution over the prior uncertainty intervals #
  
  
  # Exclude profiled parameter
  pul2 <- pul2[!(pul2$parameter == profile_est | pul2$status == "fixed"),]
  
  # Allocate space in memory to save optimized values
  pl <- vector("list", nrow(pul2)+4)
  # names(pl) <- c(pul2$parameter,       # for optimx
  #                "value",
  #                "fevals", 
  #                "convcode",
  #                "xtime")
  # 
  
  names(pl) <- c(pul2$parameter,  # for nl.optr
                 "value",
                 "niter", 
                 "convergence",
                 "messsage") 
  
  
  # Choose starting values
  
  # if(V_value== 1|!(chain_algo)){ # If on the first iteration, or if the chain algorithm is not being used, starting values are the optimized values 
  #   index <- 0
  # }else{
  #   convcodes <- lapply(prev_ests, function(x) x$convcode %in% c(9999, 52)) # for optimx
  #   index <- V_value - 1
  #   #print(paste0("index: ", index))
  #   #print(paste0("convcodes: ", convcodes))
  #   while(index > 0 & convcodes[[index]]){
  #     index <- index - 1
  #     # print(paste0("index: ", index))
  #   }
  # }
  
  if(V_value== 1|!(chain_algo)){ # If on the first iteration, or if the chain algorithm is not being used, starting values are the optimized values 
    index <- 0
  }else{
    convcodes <- lapply(prev_ests, function(x) x$convergence < 0) # for nloptr
    index <- V_value - 1
    #print(paste0("index: ", index))
    #print(paste0("convcodes: ", convcodes))
    while(index > 0 & convcodes[[index]]){
      index <- index - 1
      # print(paste0("index: ", index))
    }
  }
  
  # print(paste0("index: ", index))
  
  # 
  if(index==0){
    x0 <- log(pul2$opt)
  }else{
    x0 <- unlist(prev_ests[[index]][1:nrow(pul2)], use.names=FALSE)
  }
  # 
  # # print(paste0("x0: ", x0))
  # # print(paste0("lower: ", log(pul2$lower.bound) ))
  # # Quick check for suitability of values 
  # 
  # 
  # 
  # lb_logical <- x0 < log(pul2$lower.bound) # pick outlaws, those that are true do not satisfy the required  conditions
  # ub_logical <- x0 > log(pul2$upper.bound) # pick outlaws
  # 
  # while(sum(lb_logical) > 0 | sum(ub_logical) > 0){
  #   print("in while fix loop")
  #   for (i in 1:sum(lb_logical)){
  #     ind <- which(lb_logical==TRUE)     # identify indices of outlaws
  #     differences <- log(pul2$lower.bound)[ind] - x0[ind]
  #     x0[ind] <-   x0[ind] + differences 
  #   }
  #   
  #   for (i in 1:sum(ub_logical)){
  #     ind <- which(ub_logical==TRUE)     # identify indices of outlaws
  #     differences <- x0[ind] - log(pul2$upper.bound)[ind]
  #     x0[ind] <-   x0[ind] - differences 
  #     x0[ind] <- max(0, x0[ind])# Ensure non-negative values
  #     
  #   }
  #   
  #   lb_logical <- x0 < log(pul2$lower.bound) # pick outlaws, those that are true do not satisfy the required  conditions
  #   ub_logical <- x0 > log(pul2$upper.bound) # pick outlaws
  # }
  # 
  # 
  # print("lower=")
  # print(pul2$lower.bound)
  # print("x0=")
  # print(exp(x0))
  # print("upper=")
  # print(pul2$upper.bound)
  
  
calibration_results <- nloptr::bobyqa(x0 = x0,
                                        lower = log(pul2$lower.bound),
                                        upper = log(pul2$upper.bound),
                                        fn = calib.function.pl,
                                        nl.info = F,
                                        control = list(xtol_rel = 1e-15, 
                                                       maxeval = 2000,
                                                       ftol_rel=1e-15,
                                                       ftol_abs = 1e-15),
                                        ydata = obs,
                                        N = mean(N_vec),
                                        ...)
  
  print("NLLK=")
  print(calibration_results$value)
  
  print("best set:")
  print(exp(calibration_results$par))
  
  
  for(i in 1:nrow(pul2)){
    pl[[i]] <- calibration_results$par[[i]]
  }
  pl[[nrow(pul2) + 1]] <- calibration_results$value
  pl[[nrow(pul2) + 2]] <- calibration_results$iter
  pl[[nrow(pul2) + 3]] <- calibration_results$convergence
  pl[[nrow(pul2) + 4]] <- calibration_results$message
  
  
  return(pl)
  
}



## LIKELIHOOD FUNCTION FOR PROFILE LIKELIHOODS ###

calib.function.pl <- function(x,
                              N,
                              ydata,
                              vacnumperwk = vaccdata,
                              fitlength = 69,
                              tau = 0.9,
                              profile_lik = profile_param,
                              profileint = profile_vals,
                              roundx = iter,
                              pul1 = pul){


  # Parameters being optimized
  being_optimized <- pul1$parameter[!(pul1$parameter==profile_lik | pul1$status == "fixed")]

  # Satisfy requirement for optimx, where the first parameter has to be the vector of parameter values to be optimized
  for (i in 1:length(being_optimized)){
    pul1$opt[pul1$parameter==being_optimized[i]] <- exp(x[i])
  }

  # Pull index of variable being profiled
  index.opt <- pul1$index[pul1$parameter==profile_lik]

  # Set fixed value
  pul1$opt[pul1$index==index.opt] <- profileint[roundx]

  L. = (pul1$opt[pul1$parameter == "durI"]/pul1$opt[pul1$parameter == "durL"])*pul1$opt[pul1$parameter == "I0"]
  I. = pul1$opt[pul1$parameter == "I0"]
  Z. = pul1$opt[pul1$parameter == "fracImmune"]*N
  R. = 0
  S. = N - L. - I. - Z. - R.

  states = c(S=S.,
             L=L.,
             I=I.,
             Z=Z.,
             R=R.)

  # Parameters
  parameters = c(list(durL = pul1$opt[pul1$parameter == "durL"],
                      durI =pul1$opt[pul1$parameter == "durI"],
                      durRem = pul1$opt[pul1$parameter == "durRem"],
                      propRelapse = pul1$opt[pul1$parameter == "propRelapse"],
                      a= pul1$opt[pul1$parameter == "a"],
                      b =  pul1$opt[pul1$parameter == "b"],
                      c =  pul1$opt[pul1$parameter == "c"],
                      d = pul1$opt[pul1$parameter == "d"]))


  tspan = seq(0, length(ydata), by = 1)

  # Run the ODE solvers
  theta_part <-pul1$opt[pul1$parameter == "omega"]*vacnumperwk*tau    # Effective number of vaccinations ($\theta(t)$)

  result_df <- data.frame(ode(y = states,
                              times = tspan,
                              func = HAVODELV2,
                              parms = parameters,
                              method = "bdf_d",
                              effecvacnumperwk = theta_part,
                              fitlength = fitlength
  ))

  IncDetec <-((1/pul1$opt[pul1$parameter == "durL"])*result_df[, 3] + (1/pul1$opt[pul1$parameter == "durRem"])*result_df[, 6])[-1]

  # Calculate Poisson likelihood (**negative** log likelihood here, as optimx performs minimization)
  LLK<- -sum((dpois(x = ydata , lambda = IncDetec, log = T)))

  if(LLK==-Inf |is.nan(LLK)){
    LLK = -9999
  }
  if(LLK==Inf){
    LLK = 9999
  }
  return(LLK)
}



## FUNCTION TO COMPUTE OUTBREAK TRAJECTORY (MODEL OUTPUT GIVEN PARAMETER SET) ##

## GIVEN A PARAMETER SET, COMPUTE OUTBREAK TRAJECTORY AND THE LIKELIHOOD OF THE TRAJECTORY ###

# Inputs:
#' @x Vector of parameter values
#' @vacnumperwk Vector of weekly vaccination counts
#' @fitlength The model is calibrated using case counts up to `fitlength` week
#' @tau First-dose vaccine efficacy
#' @fp   Data frame of fixed parameter names and values
#' @ydata Observed case counts
#' @compute.likelihood Should likelihood be computed?
#' @sensitivity.analysis Is the function being called in the context of a sensitivity analysis?
#' @nsim Number of simulations


# Load functions 
calc.mod.estims <- function(x, 
                            vacnumperwk = vaccdata,
                            fitlength = 69, 
                            tau= 0.90,
                            fp = fixed_params, 
                            ydata = casecountperwk,
                            compute.likelihood = FALSE,
                            sensitivity.analysis = FALSE,
                            nsim = 10000){
  
  FC <- matrix(NA, nrow = fitlength, ncol = nsim)    # For saving trajectories
  NLLK <- c()                                        # For saving (negative) log likelihoods
  
  for(j in 1:nsim){
    
   # paste("Iteration:", print(j))
    
    # x <- exp(x)      # Convert to original scale 
    a <- x[j,"a"]      # Effective contact rate at the start of the outbreak (β_s, in manuscript)
    b <- x[j,"b"]      # Effective contact rate later in the outbreak (β_l, in manuscript)
    c <- x[j,"c"]      # Speed of transition from a to b
    d <- x[j,"d"]      # Transition midpoint time from β_s to β_l (in weeks)
    I0 <- x[j,"I0"]    # Number of infectious individuals in week 0
    N <- x[j, "N"]     # Population size
    
    ind <- ifelse(sensitivity.analysis, j, 1)
    
    # Initial conditions 
    L. = ((fp["durI", ind]/fp["durL", ind]))*I0    
    I. = I0                                  
    Z. = fp["fracImmune", ind]*N
    R. = 0
    S. = N - L. - I. - Z. - R.
    
    # States
    states = c(S=S.,
               L=L.,
               I=I.,
               Z=Z.,
               R=R.)
    
    #  print(paste("states", states))
    
    
    # Parameters
    parameters = c(list(durL = fp["durL",   ind],
                        durI = fp["durI",   ind],
                        durRem =fp["durRem",   ind],
                        propRelapse =fp["propRelapse",   ind],
                        a= a,
                        b = b,
                        c = c,
                        d = d))
    # print(parameters)
    
    # Define time
    tspan <- seq(0, fitlength, by=1)
    
    # Run the ODE solver
    omega <- 1
    theta_part <- tau*vacnumperwk*omega # Effective number of vaccinations ($\theta(t)$) 
    #  print(paste0("theta", sum(theta_part)))
    
    result_df <- data.frame(deSolve::ode(y = states,
                                         times = tspan,
                                         func = HAVODELV2,
                                         parms = parameters,
                                         method = "bdf_d",
                                         effecvacnumperwk = theta_part, 
                                         fitlength = fitlength
    ))
    
    # Calculate weekly **detected** incidence
    FC[,j] <- ((1/fp["durL",   ind])*result_df[, 3] + (1/fp["durRem",   ind])*result_df[, 6])[-1]
    #print(sum(FC[,j]))
    if(compute.likelihood){
      NLLK[j]  <- -sum((dpois(x = ydata , lambda = FC[,j], log = T))) 
     # print(paste("Likelihood value: ", NLLK[j]))
    }
  }  
  if(compute.likelihood){
    return(list(FC=FC, NLLK=NLLK, LHS = x))
  }else{
    return(list(FC=FC))
  }
  
}

###### FUNCTIONS FOR PLOTS #########





### FUNCTION TO PLOT PROFILE LIKELIHOODS ###
# Inputs:
#'@list_item: list of outputs from pl_function, each component of the list corresponding to a different parameter
#'@log: should x-axis be on the log scale?
#'@plot: should profile likelihood be plotted?
#'@cutoff_value: minimum cutoff for likelihood (according to profile likelihood confidence interval)


profile_function <- function(list_item = res5, log = T, plot = FALSE, cutoff_value){
  
  # Load required packages
  #print(list_item[[1]]$name)
  
  #  Obtain local likelihood values
  # Extract likelihood values from parameter under consideration
  llk_values <- c()
  
  for(j in 1:length(list_item[[1]]$res)){
    llk_values <- c(llk_values, list_item[[1]]$res[[j]]$value)
  }
  
  # Data frame of plotting values 
  df <- data.frame(profile_values = list_item[[1]]$profile_vals[1:length(list_item[[1]]$res)],
                   llk_values =  llk_values, 
                   diffs = diff(c(NA,llk_values))) # For the detection of non-monotonicities
  
  opt.y <- min(df$llk_values)
  opt.x <- df$profile_values[which(df$llk_values==opt.y)]
  
  # Exclude improbable llk values and points
  df.frac <- df[!(c(is.infinite(df$profile_values) | df$llk_values > 2e3 | is.na(df$llk_values) | abs(df$llk_values) == 9999)),  ]
  
  
  if(log){
    df$profile_values <- log(df$profile_values)
    xlabel <- xlabel.log.function(list_item[[1]]$name)
    opt.x <- log(opt.x)
    
  }else{
    xlabel <- xlabel.function(list_item[[1]]$name)
  }
  
  if(plot){
    if (!require(ggplot2)) install.packages("ggplot2")
    library(ggplot2)
    
    p <- ggplot2::ggplot(data = df.frac,
                         mapping = aes(x = profile_values,
                                       y = llk_values)) +
      scale_y_continuous(limits = c(128, 140),  breaks = seq(130, 138, 2) ) + #n.breaks = 6 # opt.y + qchisq(0.95, 9)/2
      # , breaks = seq(128, 138, 2)
      # revert to no breaks 
      # scale_x_continuous(limits = c(0.50, 1.2)) +
      geom_hline(yintercept = cutoff_value, linetype = "dashed", colour = "red")+ # Max threshold Raue(2009); 95% threshold - no theoretical justification
      # geom_hline(yintercept = opt.y + qchisq(0.95, 1)/2 , linetype = "dashed", colour = "red")+ # Min threshold
      geom_line() +
      theme_bw()+
      theme(legend.position = "none")+
      # annotate("point", x = opt.x, y = opt.y, colour = "blue") +
      labs(x = xlabel,
           y = "-log(PL)")
    return(list(optimum.info = list(minx=opt.x, miny = opt.y, name=list_item[[1]]$name), df.full = df, df = df.frac,  plt = p))
    
  }else{
    return(list(optimum.info = list(minx=opt.x, miny = opt.y, name=list_item[[1]]$name), df.full = df, df = df.frac,  plt = NULL))
    
  }
  
}




### FUNCTION TO DETERMINE X-AXIS LABELS IN PROFILE LIKELIHOOD PLOTS ###
# Inputs: 
#'@p: parameter for which axis label is to be printed

xlabel.function <- function(p = "propRelapse"){
  # Set x-axis label
  if(p == "a"){
    xlabel <- expression(beta[s])
  }else{
    if(p == "b"){
      xlabel <- expression(beta[l])
    }else{
      if(p == "c"){
        xlabel <- expression(c)
      }else{
        if(p == "d"){
          xlabel <- expression(t^"*")
        }else{
          if(p=="I0"){
            xlabel <- expression(I[0])
          }else{
            if(p == "omega"){
              xlabel <- expression(omega)
            }else{
              if(p == "fracImmune"){
                xlabel <- expression(epsilon)
              }else{
                if(p =="propRelapse"){
                  xlabel <- expression(1~-~eta)
                }else{
                  if(p == "kappa"){
                    xlabel <- expression(kappa)
                  }else{
                    if(p == "t_offset"){
                      xlabel <- expression(t["offset"])
                    }else{
                      xlabel <- "NA"
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(xlabel)
}




### FUNCTION FOR SCATTER PLOT (FIGURE S5, SUPPLEMENTARY) ###
# Note: modification of function `plot_scatter` in sensobol package (Puy et al., 2021)
# For input definition, see ?sensobol::plot_scatter

plot_scatter2 <- function (data, N, Y, params, method = "point", size = 0.7,
                           alpha = 0.2)
{
  value <- y <- NULL
  dt <- data.table::data.table(cbind(data, Y))[1:N]
  colnames(dt)[length(colnames(dt))] <- "y"
  out <- data.table::melt(dt, measure.vars = params)
  levels(out$variable) <- c(expression(beta[s]), expression(epsilon), expression(1/gamma))
  gg <- ggplot2::ggplot(out, ggplot2::aes(value, y)) + ggplot2::facet_wrap(~variable,
                                                                           scales = "free_x", labeller = label_parsed) + ggplot2::labs(x = "Value",
                                                                                                                                       y = "y") + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                                                                                                                                                                       panel.grid.minor = ggplot2::element_blank(), legend.background = ggplot2::element_rect(fill = "transparent",
                                                                                                                                                                                                                                                                              color = NA), legend.key = ggplot2::element_rect(fill = "transparent",                                                                                                                                                                                                                                                                                                                           color = NA), strip.background = ggplot2::element_rect(fill = "white"),
                                                                                                                                                                                       legend.position = "top")
  if (method == "point") {
    gg <- gg + ggplot2::geom_point(size = size, alpha = alpha) +
      ggplot2::stat_summary_bin(fun = "mean", geom = "point",
                                colour = "red", size = 0.7)
  }
  else if (method == "bin") {
    gg <- gg + ggplot2::geom_hex() + ggplot2::stat_summary_bin(fun = "mean",
                                                               geom = "point", colour = "red", size = 0.7)
  }
  else {
    stop("Method should be either point or bin")
  }
  return(gg)
}




#### MISCELLANEOUS ####


## FUNCTION TO ASSIGN "FIXED" OR "FREE" STATUS TO PARAMETERS ##
# Inputs:
#'@pul: Data frame of parameters
#'@to_fix: vector of parameters to fix
#'@fixed_values: corresponding values for parameters in `to_fix`

fix.parameters <- function(pul_arg = pul,
                           to_fix = fixed_parameters, # defined above
                           fixed_values = fixed_parameter_values){
  
  # Create two additional columns for pul_arg
  pul_arg$status <- rep("free", nrow(pul_arg)) # Indicate status of parameter: whether "free" or "fixed". You are free until you are fixed!
  pul_arg$fixed_val <- rep(NA, nrow(pul_arg)) # For fixed parameters, indicate fixed value 
  
  # Complete modification of pul: specify fixed parameter and their values 
  for(i in 1:length(to_fix)){
    pul_arg$status[pul$parameter == to_fix[i]] <-"fixed"
    pul_arg$fixed_val[pul$parameter == to_fix[i]] <- fixed_values[i]
    
  }
  return(pul_arg)
}




## FUNCTION TO COMPUTE A SEQUENCE FOR BETA, GIVEN VALUES FOR THE SIGMOIDAL PARAMETERS ##
# Inputs:
#'@x: vector of estimates of the sigmoidal parameters
#'@obs: Observed case counts

beta_function <- function(x, obs = casecountperwk){ # vector x holds the estimated sigmoidal parameters
  beta <-  x[1] + (x[2] - x[1])/(1+exp(-x[3]*(1:length(obs)-x[4])))
  return(beta)
}






### FUNCTION TO MULTIPLY VACCINATION COUNTS UP TO THE OBSERVED COUNT ####
#'@v Observed vaccination data
#'@multiplier Multiplier. To double vaccination counts, multiplier=2


multiply.vacc.counts <- function(v, multiplier=2){
  
  scenario <- c() # Vector to save updated counts 
  
  for (i in 1:length(v)){
    
    if(sum(scenario)<=sum(v)){
      scenario[i] <- round(multiplier*v[i])
    }else{
      scenario[i] <-  0
    }
  }
  scenario
}



### FUNCTION TO VARY PARAMETER VALUES BY A PERCENTAGE ####
# Inputs: 
#'@val Parameter value
#'@perc Percentage to vary 'val' by

vary <- function(val, perc=0.2){
  c(val - val*perc, val + val*perc)
}



### CORRESPONDING FUNCTIONS FOR SINGLE BETA MODEL #####


# Corresponds to "HAVODELV2" above 

HAVODELV2_old <- function(t, states, parameters, fitlength, effecvacnumperwk){
  with(as.list(c(states,parameters)),
       {


         states[1] -> S
         states[2] -> L
         states[3] -> I
         states[4] -> Z
         states[5] -> R

         # sum to N
         N    <- S + L + Z + I + R


         # vaccination value
         if (t < fitlength){
           vacc <- effecvacnumperwk[floor(t+1)]
         }else{                                   # after calibration assume same vaccination as previous week
           vacc <- effecvacnumperwk[fitlength]}
         if(vacc > S){vacc = S}

         # print(S/N)

         dS <- -(beta) * S * (I/N) - vacc*(S/N) # S/N is the variable part of theta
         dL <-  beta * S * (I/N) - (1/durL)*L
         dI <-  (1/durL)*L +(1/durRem)*R - (1/durI)*I
         dZ <- (1-propRelapse)*(1/durI)*I + vacc*(S/N)
         dR <- propRelapse*(1/durI)*I - (1/durRem)*R

         return(list(c(dS, dL, dI, dZ, dR)))

       })
}




# Corresponds to "calib.function" above 

calib.function_old <- function(x,
                           N,
                           ydata,
                           vacnumperwk = vaccdata,
                           fitlength = 69,
                           tau=0.9,
                           fp = fixed_params){
  # As required by optimx, first argument of the optimization function is a vector of parameters of interest.
  x <- exp(x)
  beta <- x[1]           # reduction in force of infection (weeks 34-69)
  I0 <- x[2]           # Initial number of infected individuals


  omega <- fp["omega", 1]     # omega is assumed to be 1
  # Initial conditions
  L. = (fp["durI", 1]/fp["durL", 1])*I0     # set L_0 to be at equilibrium with I_0
  I. = I0                                   # I_0*k
  Z. = fp["fracImmune", 1]*N
  R. = 0
  S. = N - L. - I. - Z. - R.


  states = c(S=S.,
             L=L.,
             I=I.,
             Z=Z.,
             R=R.)

  #print(paste("states", states))


  # Parameters
  parameters = c(list(durL = fp["durL", 1],
                      durI = fp["durI", 1],
                      durRem =fp["durRem", 1],
                      propRelapse =fp["propRelapse", 1],
                      beta = beta
                    ))
 # print(parameters)

  # Calibrate up to:
  ydata <- ydata[1:fitlength]

  # Define time
  tspan <- seq(0, length(ydata), by=1)

  # Run the ODE solver
  theta_part <- vacnumperwk*tau*omega  # Constant part of theta, the effective number of vaccinations ($\theta(t)$)


  #print(paste0("theta", theta_part))

  result_df <- data.frame(ode(y = states,
                              times = tspan,
                              func = HAVODELV2_old,
                              parms = parameters,
                              method = "bdf_d",
                              effecvacnumperwk = theta_part,
                              fitlength = fitlength
  ))

  # Calculate weekly **detected** incidence
  IncDetec <- ((1/fp["durL", 1])*result_df[, 3] + (1/fp["durRem", 1])*result_df[, 6])[-1]
  #print(IncDetec)
  # Calculate Poisson likelihood (**negative** log likelihood here, as optimx performs minimization)
  NLLK<- -sum((dpois(x = ydata , lambda = IncDetec, log = T)))
  print(paste("Likelihood value", NLLK))

  if(NLLK==-Inf |is.nan(NLLK)){# Return a large number if value is improbable
    NLLK = -9999
  }
  if(NLLK==Inf){
    NLLK = 9999
  }
  return(NLLK)
}







# Corresponds to "calc.mod.estims" above 

calc.mod.estims_old <- function(x,
                            vacnumperwk = vaccdata,
                            fitlength = 69,
                            tau=0.9,
                            fp = fixed_params,
                            ydata = casecountperwk,
                            compute.likelihood = FALSE,
                            sensitivity.analysis = FALSE,
                            nsim = 10000){

  FC <- matrix(NA, nrow = fitlength, ncol = nsim) # For saving trajectories
  NLLK <- c()                                        # For saving (negative) log likelihoods

  for(j in 1:ncol(FC)){

   # print(j)
    # x <- exp(x)         # Convert to original scale
    beta <- x[j,"beta"]      # reduction in force of infection (weeks 34-69)
    I0 <-   x[j,"I0"]        # Initial number of infected individuals


    N <- ifelse(sensitivity.analysis, fp["N", 1] ,  x[j, "N"]) # If sensitivity analysis, do not use fixed value in x vector
    #print(N)
    # Initial conditions
    L. = (fp["durI", 1]/fp["durL", 1])*I0     # set L_0 to be at equilibrium with I_0
    I. = I0                                   # I_0*k
    Z. = fp["fracImmune", 1]*N
    R. = 0
    S. = N - L. - I. - Z. - R.


    states = c(S=S.,
               L=L.,
               I=I.,
               Z=Z.,
               R=R.)

    #print(paste("states", states))


    # Parameters
    parameters = c(list(durL = fp["durL", 1],
                        durI = fp["durI", 1],
                        durRem =fp["durRem", 1],
                        propRelapse =fp["propRelapse", 1],
                        beta = beta))
    #print(parameters)

    # Define time
    tspan <- seq(0, fitlength, by=1)

    # Run the ODE solver
    theta_part <- tau*vacnumperwk*fp["omega", 1]  # Effective number of vaccinations ($\theta(t)$)
    #print(paste0("theta", sum(theta)))

    result_df <- data.frame(deSolve::ode(y = states,
                                         times = tspan,
                                         func = HAVODELV2_old,
                                         parms = parameters,
                                         method = "bdf_d",
                                         effecvacnumperwk = theta_part,
                                         fitlength = fitlength
    ))

    # Calculate weekly **detected** incidence
    FC[,j] <- ((1/fp["durL", 1])*result_df[, 3] + (1/fp["durRem", 1])*result_df[, 6])[-1]
    #print(FC)
    if(compute.likelihood){
      NLLK[j]  <- -sum((dpois(x = ydata , lambda = FC[,j], log = T)))
      print(NLLK[j])
    }
  }
  if(compute.likelihood){
    return(list(FC=FC, NLLK=NLLK, LHS = x))
  }else{
    return(list(FC=FC))
  }

}





# Corresponds to "estimate.parameters.pl" above 


estimate.parameters.pl_old <- function(profile_est = profile_param,
                                   pul2 =  pul,
                                   obs = casecountperwk,
                                   V_value = 1,
                                   prev_ests,
                                   chain_algo = chain
                                   ,...){
  
  # Set initial, lower and upper values (See Table 1 in main text).
  # Initial values are sampled from a uniform distribution over the prior uncertainty intervals #
  
  
  # Exclude profiled parameter
  pul2 <- pul2[!(pul2$parameter == profile_est | pul2$status == "fixed"),]
  
  # Allocate space in memory to save optimized values
  pl <- vector("list", nrow(pul2)+4)
  # names(pl) <- c(pul2$parameter,       # for optimx
  #                "value",
  #                "fevals", 
  #                "convcode",
  #                "xtime")
  # 
  
  names(pl) <- c(pul2$parameter,  # for nl.optr
                 "value",
                 "niter", 
                 "convergence",
                 "messsage") 

  if(V_value== 1|!(chain_algo)){ # If on the first iteration, or if the chain algorithm is not being used, starting values are the optimized values 
    index <- 0
  }else{
    convcodes <- lapply(prev_ests, function(x) x$convergence < 0) # for nloptr
    index <- V_value - 1
    #print(paste0("index: ", index))
    #print(paste0("convcodes: ", convcodes))
    while(index > 0 & convcodes[[index]]){
      index <- index - 1
      # print(paste0("index: ", index))
    }
  }
  
  if(index==0){
    x0 <- log(pul2$opt)
  }else{
    x0 <- unlist(prev_ests[[index]][1:nrow(pul2)], use.names=FALSE)
  }

  
  calibration_results <- nloptr::bobyqa(x0 = x0,
                                        lower = log(pul2$lower.bound),
                                        upper = log(pul2$upper.bound),
                                        fn = calib.function.pl_old,
                                        nl.info = F,
                                        control = list(xtol_rel = 1e-15, 
                                                       maxeval = 2000,
                                                       ftol_rel=1e-15,
                                                       ftol_abs = 1e-15),
                                        ydata = obs,
                                        N = mean(N_vec),
                                        ...)
  
  print("NLLK=")
  print(calibration_results$value)
  
  print("best set:")
  print(exp(calibration_results$par))
  
  
  for(i in 1:nrow(pul2)){
    pl[[i]] <- calibration_results$par[[i]]
  }
  pl[[nrow(pul2) + 1]] <- calibration_results$value
  pl[[nrow(pul2) + 2]] <- calibration_results$iter
  pl[[nrow(pul2) + 3]] <- calibration_results$convergence
  pl[[nrow(pul2) + 4]] <- calibration_results$message
  
  
  return(pl)
  
}




calib.function.pl_old <- function(x,
                              N,
                              ydata,
                              vacnumperwk = vaccdata,
                              fitlength = 69,
                              tau = 0.9,
                              profile_lik = profile_param,
                              profileint = profile_vals,
                              roundx = iter,
                              pul1 = pul){
  
  
  # Parameters being optimized
  being_optimized <- pul1$parameter[!(pul1$parameter==profile_lik | pul1$status == "fixed")]
  
  # Satisfy requirement for optimx, where the first parameter has to be the vector of parameter values to be optimized
  for (i in 1:length(being_optimized)){
    pul1$opt[pul1$parameter==being_optimized[i]] <- exp(x[i])
  }
  
  # Pull index of variable being profiled
  index.opt <- pul1$index[pul1$parameter==profile_lik]
  
  # Set fixed value
  pul1$opt[pul1$index==index.opt] <- profileint[roundx]
  
  L. = (pul1$opt[pul1$parameter == "durI"]/pul1$opt[pul1$parameter == "durL"])*pul1$opt[pul1$parameter == "I0"]
  I. = pul1$opt[pul1$parameter == "I0"]
  Z. = pul1$opt[pul1$parameter == "fracImmune"]*N
  R. = 0
  S. = N - L. - I. - Z. - R.
  
  states = c(S=S.,
             L=L.,
             I=I.,
             Z=Z.,
             R=R.)
  
  # Parameters
  parameters = c(list(durL = pul1$opt[pul1$parameter == "durL"],
                      durI =pul1$opt[pul1$parameter == "durI"],
                      durRem = pul1$opt[pul1$parameter == "durRem"],
                      propRelapse = pul1$opt[pul1$parameter == "propRelapse"],
                      beta= pul1$opt[pul1$parameter == "beta"]
                     ))
  
  
  tspan = seq(0, length(ydata), by = 1)
  
  # Run the ODE solvers
  theta_part <-pul1$opt[pul1$parameter == "omega"]*vacnumperwk*tau    # Effective number of vaccinations ($\theta(t)$)
  
  result_df <- data.frame(ode(y = states,
                              times = tspan,
                              func = HAVODELV2_old,
                              parms = parameters,
                              method = "bdf_d",
                              effecvacnumperwk = theta_part,
                              fitlength = fitlength
  ))
  
  IncDetec <-((1/pul1$opt[pul1$parameter == "durL"])*result_df[, 3] + (1/pul1$opt[pul1$parameter == "durRem"])*result_df[, 6])[-1]
  
  # Calculate Poisson likelihood (**negative** log likelihood here, as optimx performs minimization)
  LLK<- -sum((dpois(x = ydata , lambda = IncDetec, log = T)))
  
  if(LLK==-Inf |is.nan(LLK)){
    LLK = -9999
  }
  if(LLK==Inf){
    LLK = 9999
  }
  return(LLK)
}



# Corresponds to "pl_function" above 

pl_function_old <- function(p.vec,
                        pvl=profile_vals_list,
                        ...){

  V = length(pvl[[p.vec]]) # Number of values for which the profile likelihood are to be computed

  # Create empty list to save outputs
  master_list <- vector("list", length(p.vec))
  names(master_list) <- p.vec

  for(parameter in p.vec){
    # Define interval over which profile is to computed
    profile_param <- parameter

    profile_vals <- pvl[[profile_param]]

    print(profile_vals)

    # Create list to store values
    res <- list()                   # Parameter estimates

    for (iter in c(1:V)){
      print(paste0("iter: ", iter))

      if(iter==1){
        prev_ests_val <- NA
      }else{
        prev_ests_val <- res
        # print(paste0("prev_ests: ", prev_ests_val))
      }


      res[[iter]] <- estimate.parameters.pl_old(profile_est = profile_param,
                                            # pul2 =  prior_uncertainty_limits,
                                            # N.vector = N_vec,
                                            # bootstr = bootdata,
                                            # vacnumperwk = vaccdata,
                                            #fitlength = 69,
                                            #omega = 0.9,
                                            profile_lik = profile_param,
                                            profileint = profile_vals,
                                            roundx = iter,
                                            V_value = iter,
                                            prev_ests = prev_ests_val
                                            # pul1 = prior_uncertainty_limits
      )

      #print(res[[iter]]) # Monitor

    }

    # Save results to master list
    master_list[[parameter]]$res <- res
    master_list[[parameter]]$profile_vals <- profile_vals
    master_list[[parameter]]$name <- profile_param

  }
  return(master_list)

}




#Reference
# 1. Puy A, Piano SL, Satelli A, Levin SA. sensobol: an R package to compute variance-based sensitivity indices 2021. arXiv:2101.10103v2
