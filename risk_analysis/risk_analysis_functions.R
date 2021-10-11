######### MISCELLANEOUS FUNCTIONS  FOR RELATIVE RISK ANALYSIS ###

# Author: EMMANUELLE A. DANKWA (email:dankwa@stats.ox.ac.uk)

###############################################################



# 1. Functions to convert dates to year-week format ####
# Arguments:
# df is the data frame to be converted
# yc - year column

func1 <- function(df, yc = 4){
  
  df$Year_Week <- NA
  
  for(i in 1: dim(df)[1]){
    
    df$Year_Week[i] = paste0(df[i,yc] , '-', df[i,2] )
    
  }
  return(df)
}



# Arguments:
# dat is the standard (with all weeks represented)
# dat_10 is the data frame to be made standard

aggregate <- function(dat_10, dat = dat2_yrwk){
  
  dat_10 <- dat_10 [, c('MMWRYear','MMWRWeek')]
  
  # Aggregating over weeks
  
  dat_10 <- as.data.frame(table(dat_10$MMWRYear, dat_10$MMWRWeek ))
  
  # Assign column names
  
  colnames(dat_10) <- c( 'Year', 'Week', 'numcases')
  
  # Order by year
  
  dat_10 <- dat_10[order(dat_10$Year), ]
  rownames(dat_10) <- NULL
  
  n.cases <- sum(dat_10$numcases)
  
  ############################################################
  # Convert to year-week format
  
  ###########################################################  
  dat_10 <- func1(df=dat_10, yc=1)
  dat <- func1(df=dat, yc=1)
  
  cases_10 <- c()   # Store case counts
  
  for(i in 1:length(dat$Year_Week)){ # Compare to (full) dat_m and put 0 in week positions where no case was recorded.
    if(dat$Year_Week[i]%in%dat_10$Year_Week){
      # Pick up index
      index <- which(dat_10$Year_Week == dat$Year_Week[i])
      cases_10[i] <- dat_10$numcases[index]
    }else{
      cases_10[i] <- 0       
    }
    n.cases2 <- sum(cases_10) # this should be equal to n.cases
  }
  
  dat_10 <- data.frame(Year = dat$Year, Week= dat$Week, numcases = cases_10, Year_Week = dat$Year_Week)
  
  return(list(final_dat=dat_10, n.cases = n.cases, n.cases2 = n.cases2, cases = cases_10))
}






# 3. Function to calculate the number of cases in each period.
# Arguments:
# vec = vector of cases 

cases.by.period <- function(vec){
  # Pre-peak (up to week 14 of 2018),
  # peak (weeks 15-17 of 2018) 
  # Post-peak (week 18 of 2018 to week 22 of 2019) 
  
  prp <- sum(vec[1:31])
  peak <- sum(vec[32:34])
  pop <- sum(vec[35:91])
  total <- sum(prp, peak, pop)
  
  return(list(prepeak=prp, peak = peak, postpeak = pop, total = total))
  
}



### RISK ANALYSIS ###



# 4. Function to compute relative risk ratios and odd ratios by age and sex ###

relative.risk.results <- function(population = 'risk', totalprepeak, totalpostpeak, var = 'age'){
  
  if(population=='general'){
    # set reference data frame based on population group
    ref_df <- dat2
  }else{
    if(population=='risk'){
      ref_df <- target
    }else{
      stop('population must be either `general` or `risk`')
    }
  }
  
  if(var == 'age'){ # If Age
    # Create data frames for all age levels in year-week-numcases format. 
    var_dfs <- automate.age.dfs(ref = ref_df)
    num_levels <- length(levels(ref_df$Age_group))
    levs <- levels(ref_df$Age_group)
  }else{
    if(var=='sex'){ # If sex
      # Create data frames for males and females
      
      dat_m1 <- aggregate(ref_df[ref_df$Sex=='Male', ])$final_dat
      dat_f1 <- aggregate(ref_df[ref_df$Sex=='Female',])$final_dat
      var_dfs <- list(dat_m1, dat_f1)
      num_levels <- 2 # Male and Female (No NAs in data)
      levs <- ref_df$Sex
    }else{
      stop('var must be either `age` or `sex`')
    }
  }
  
  # Dataframe to store results
  res <- data.frame(category = numeric(num_levels),
                    prepeak = numeric(num_levels),
                    peak = numeric(num_levels),
                    postpeak = numeric(num_levels),
                    total = numeric(num_levels),
                    relative_risk = numeric(num_levels),
                    relative_risk_CI = numeric(num_levels))
  
  for(i in 1:num_levels){
    temp <- cases.by.period(var_dfs[[i]]$numcases)
    res$category[i] <-  paste0(var, '_', levs[i])
    res$prepeak[i] <- temp$prepeak
    res$peak[i] <- temp$peak
    res$postpeak[i] <- temp$postpeak
    res$total[i] <- temp$total
    
    if(res$total[i] > 10){
      rr <- riskratio(res$prepeak[i], res$postpeak[i], totalprepeak, totalpostpeak)
      res$relative_risk[i] <- round(rr$estimate, 2)
      res$relative_risk_CI[i] <- paste0(round(rr$conf.int[1],2),  '-',  round(rr$conf.int[2], 2))
    }else{
      res$relative_risk[i] <- NA
      res$relative_risk_CI[i] <- NA
    }
  }
  
  # Compute odds ratios
  
  if(var=='sex'){
    odds_ratio <- oddsratio(a=matrix(data=c(res$prepeak[1], res$postpeak[1], res$prepeak[2], res$postpeak[2]), nrow=2, ncol = 2))
    odds_ratio_estimate <- round(odds_ratio$estimate, 2)
    odds_ratio_CI <-  paste0(round(odds_ratio$conf.int[1],2),  '-',  round(odds_ratio$conf.int[2], 2))
    results <- list(relative_risk_results = res, odds_ratio = odds_ratio_estimate, odds_ratio_CI = odds_ratio_CI)
    
  }else{ 
    # when var='age', do:
    # Create list to save odds ratios
    odds_ratios <- list()
    length(odds_ratios) <- 4
    for(i in 3:(num_levels-1)){
      temp_name <- paste0('odds_ratio_', levs[i], 'vs')
      temp_list <- list() # list to save odds ratios for each age level
      
      for( j in 2:(i-1)){
        sub <- paste0(levs[j], ':')
        temp_res <- oddsratio(a=matrix(data=c(res$prepeak[i], res$postpeak[i], res$prepeak[j], res$postpeak[j]), nrow=2, ncol = 2))
        odds_ratio_estimate <- round(temp_res$estimate, 2)
        odds_ratio_CI <- paste(round(temp_res$conf.int[1], 2), '-', round(temp_res$conf.int[2], 2))
        temp_res <- paste(paste0(temp_name, sub), 'estimate = ', paste0(odds_ratio_estimate,','), 'CI = ', odds_ratio_CI)
        odds_ratios[[i-2]][j-1] <- temp_res
      }}
    results <- list(relative_risk_results = res, odds_ratios = odds_ratios)
    
  }
  results
}



# 5. Function to compute moving averages

ma = function(x, k) {
  if (k == 1)
    return(x)
  n = length(x)
  out = x[-(1:(k - 1))]/k
  
  for (i in 2:k) {
    out = out + x[seq(from = k + 1 - i, to = n + 1 - i)]/k
  }
  out
  
}


# 6. Function to clean age-focused data frames and produce results in the year-week-number of cases format. 



automate.age.dfs <- function(ref, dat = dat2_yrwk){
  
  # Create empty list for storing data frames for age
  
  results <- list()
  for(i in 1:length(levels(ref$Age_group))){
    
    # Subset ages for that age group
    dat_10 <- ref[ref$Age_group ==levels(ref$Age_group)[i],]
    results[[i]] <- aggregate(dat_10, dat)$final_dat 
    
  }
  results
}





### PLOTTING FUNCTIONS ###

# Function to plot epidemic curves by risk group (PEH/PWUD or not)
# Argument
# dat: weekly detected cases by risk group 

main.risk.plot <- function(dat){
  
  windowsFonts(A = windowsFont("Times New Roman")) # Set font
  par(family = "serif", cex= 1.1)
  plot(x = ma(1:91, 2), y = ma(as.numeric(dat$numcases), 2), type = 'l', xlab = 'Week', ylab = 'Number of cases', lwd = 2, panel.first = rect(c(1, 35), -1e6, c(31, 91), col='light gray', border=NA, ytop = 35 ), family = "A") # Total
  lines(x = ma(1:91, 2), y = ma(dat$targetcases, 2), type = 'l', col = 'red', lwd=2)  # Risk group 
  lines(x = ma(1:91, 2), y = ma(dat$non_target_cases, 2), type = 'l', col = 'blue', lwd=2) # Other
  legend('topright', legend = c('Total', 'PEH/PWUD', 'Other'), col = c('black', 'red', 'blue') , title = 'All cases' , lwd = 2, cex = 0.8)
  title(main = 'A)', adj = 0, cex = 2.0 )
}






# Function to plot epidemic curves by housing status among PEH/PWUD 
# Argument
# dat: weekly detected cases among PEH/PWUD with housing and drug use data

housing.plot <- function(dat){
  
  # Select cases who report homelessness
  hml <- target[which(target$Homelessness == 'Homeless - Shelter/Streets' | target$Homelessness== 'Unstable Housing'), ]
  # Not homeless
  nhml <-  target[-which(target$Homelessness == 'Homeless - Shelter/Streets' | target$Homelessness== 'Unstable Housing'), ]
  
  # Convert to year-week-numcases format 
  hml <- aggregate(hml)$final_dat
  nhml <- aggregate(nhml)$final_dat
  dat <- aggregate(dat)$final_dat
  
  # Plot
  plot(x = ma(c(1:91), 2), y = ma(as.numeric(dat$numcases), 2), type = 'l', xlab = 'Week', ylab = 'Number of cases', lwd = 2,
       panel.first = rect(c(1, 35), -1e6, c(31, 91), col='light gray', border=NA, ytop = 35) ) # Total
  lines(x = ma(c(1:91), 2), y = ma(hml$numcases, 2), type = 'l', col = 'red', lwd=2) # Homeless
  lines(x = ma(c(1:91), 2), y = ma(nhml$numcases, 2), type = 'l', col = 'blue', lwd = 2) # Not homeless
  legend('topright', legend = c('Total',  'Homeless' , 'Not homeless'), col = c('black', 'red', 'blue') , title = 'PEH/PWUD', lwd = 2, cex = 0.8)
  title(main = 'B)', adj = 0, cex = 2.0 )
  
}




# Function to plot epidemic curves by drug use status among PEH/PWUD 
# Argument
# dat: weekly detected cases among PEH/PWUD with housing and drug use data

druguse.plot <- function(dat){
  
  # Drug using
  dru <- target[which(target$IVDU=='Yes' | target$Non.IV.DU=='Yes' ), ]
  # Non-drug using
  ndru <- target[-which(target$IVDU=='Yes' | target$Non.IV.DU=='Yes' ), ]
  
  # Convert to year-week-numcases format 
  dru <- aggregate(dru)$final_dat
  ndru <- aggregate(ndru)$final_dat
  dat <- aggregate(dat)$final_dat
  
  # Plot
  plot(x = ma(c(1:91), 2), y = ma(as.numeric(dat$numcases), 2), type = 'l', xlab = 'Week', ylab = 'Number of cases', lwd = 2,
       panel.first = rect(c(1, 35), -1e6, c(31, 91), col='light gray', border=NA, ytop = 35)) # Total
  lines(x = ma(c(1:91), 2), y = ma(dru$numcases, 2), type = 'l', col = 'red', lwd=2) # Drug using 
  lines(x = ma(c(1:91), 2), y = ma(ndru$numcases, 2), type = 'l', col = 'blue', lwd = 2) # Non-drug using
  legend('topright', legend = c('Total',  'Drug using' , 'Non-drug using'), col = c('black', 'red', 'blue') , title = 'PEH/PWUD', lwd = 2, cex = 0.8)
  title(main = 'C)', adj = 0, cex = 2.0 )
}