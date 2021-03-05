library(data.table)
library(stringr)
library(effectsize)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggeasy)
library(ggpubr)
library(matrixStats)
library(metaviz)

# Save data folder
savePath <- 'your-path-to-data'

## Global Variables
Png <- '.png'; Csv <- '.csv'
Z_crit <- 1.96 # for two-tailed; 1.65 for one-tailed

## Functions
# Check if all elec in SumElec correspond to the ones in d
checkElec <- function(d, SumElec){
  ## Check if all electrodes in SumElec have been used in d
  AllElec <- NULL;
  for(di in d ){ El <- NULL;
  for(e in di$Elec){ if( e!= 'NP'){ El <- c( El, str_split(e, ", ")[[1]]) }}
  AllElec <- c(AllElec, El ) }
  AllElec <- unique( AllElec[order(AllElec)])
  SElec <- NULL; for(e in SumElec){ SElec <- c( SElec, str_split(e, ", ")[[1]])  }
  SElec <- unique(SElec[order(SElec)])

  if( length(setdiff(AllElec, SElec)) + length(setdiff(SElec, AllElec)) != 0 ){
    if( length(setdiff(AllElec, SElec)) > 0 ){ print('Problem with Elec from d: '); print(setdiff(AllElec, SElec)) }
    if( length(setdiff(SElec, AllElec)) > 0 ){ print('No data for: '); print(setdiff(SElec, AllElec)) }
  }else{
    print('OK')
  }
}

# Aggregate data in single matrix
getAgg <- function(d, SumElec){
  Agg <- matrix(data=NA, nrow = length(SumElec), ncol = 3*length(names(d)))

  for(nr in 1:length(SumElec) ){
    e <- str_split(SumElec[nr], ", ")[[1]]
    for(nd in names(d) ){
      de <- d[[nd]]$Elec
      nb <- NULL; for(ei in e){
        for(del in de){
          if(any( ei == str_split(del, ', ')[[1]]) ){
            nb <- c(nb, which(del==de)) } } }
      nb <- unique(nb)
      nc <- seq(1,3) + (which(nd==names(d))-1)*3

      if( length(nb) > 1 ){
        Agg[nr, nc ] <- colMeans(as.matrix(d[[nd]][nb,c('d','CI_low','CI_high'),with=F]))
      }else if( length(nb) >  0){
        Agg[nr, nc ] <- as.matrix(d[[nd]][nb,c('d','CI_low','CI_high'),with=F])
      }
    }
  }
  Agg <- data.table(Agg)
  Names <- NULL; for( n in names(d)){ Names <- c( Names, paste(n, 'd') , paste(n, 'CI_low') , paste(n, 'CI_high') ) }
  names(Agg) <- Names
  Agg <- cbind(data.table(Elec = SumElec), Agg)
}

# Calculate Cohen's d for Motor Imagery
getEffectSizeMI <- function(SumElec, Z_crit){
  d <- NULL
  ## From "Partially Overlapping Neural Networks for Real and Imagined Hand Movements"
  # From Table 1
  # two-tailed Z-score (standard normal deviate)
  # Imagined movement with rest compared
  Z2 <- c( 3.76, 4.58, 4.75, 4.31, 5.69, 5.35, 6.29, 5.81, 4.61,
          4.3, 5.18, 3.6, 4.63, 4.53, 4.73, 4.65, 4.81,
          4.3, 5.47, 5.47 )
  Elec <- c('F5', 'FP2', 'F8', 'F7', 'FC2, FC4', 'FC3', 'FCZ', 'FCZ', 'FC1',
            'C4', 'C3, C5', 'CP2, P2', 'P1', 'CP4', 'CP3', 'TP7, T7', 'NP',
            'NP', 'NP', 'NP' )

  # eight right-handed healthy volunteers
  n <- 8
  ## For right-handers only
  d[['Partially_Both']] <- data.table( Subjects = 8, Trials = 96, Elec = Elec, z_to_d(Z2, n) ) # z_to_d(Z2, 2*n-2)

  
  ## From "Localization of grasp representations in humans by positron emission tomography"
  # From Table 2
  Control <- c( 59.9, 56.3, 50.3, 67.4, 55.6, 47.6, 53.9, 61.3, 51.4, 51.8, 58.6, 61.6 )
  MI <- c( 62.7, 58.2, 51.6, 70.5, 57.4, 49.3, 55.6, 64.2, 52.6, 54.1, 60.1, 63.6 )
  ttest <- c( 6.47, 5.87, 4.18, 12.21, 4.36, 4.52, 4.76, 4.78, 5.1, 3.9, 4.13,5.06 )
  # Seven subjects
  n <- 7
  Elec <- c( 'FCZ', 'C1', 'C2', 'FCZ', 'CP3', 'PO3', 'P1', 'FC3', 'F3', 'NP', 'C3', 'NP' )

  d[['Localization_Right_Dom']] <- data.table( Subjects = 7, Trials = 130, Elec = Elec, t_to_d(ttest, 2*n-2) )
  # # From compute.es, equivalent results to effectsize
  # d <- tes(t=ttest, n.1=n, n.2=n)
  # print(d$d)

  
  ## From "Somatotopic mapping of the human primary sensorimotor cortex during motor imagery and motor execution by functional magnetic resonance imaging"
  # From Table 1
  r <- c( 0.54 )
  r_sd <- c(0.13 )
  dS <- c( 0.73 )
  dS_sd <- c( 0.48 )

  Cohen_d <- r_to_d(r)
  CI <- rep(0.95, length(r))
  # Standard deviation: With probability about 95%  (r - 2 * SD; r + 2 * SD)
  CI_low <- r_to_d(r - 2* r_sd)
  CI_high <- r_to_d(r + 2* r_sd)
  Elec <- c( 'C2, C4' )
  d[['Somatotopic_Left_NonDom']] <- data.table( Subjects = 14, Trials = 90, Elec = Elec, d = Cohen_d, CI = CI, CI_low = CI_low, CI_high = CI_high)

  # From Table 1
  r <- c( 0.59 )
  r_sd <- c( 0.09 )
  dS <- c( 1.47 )
  dS_sd <- c( 1.27 )

  Cohen_d <- r_to_d(r)
  CI <- rep(0.95, length(r))
  # Standard deviation: With probability about 95%  (r - 2 * SD; r + 2 * SD)
  CI_low <- r_to_d(r - 2* r_sd)
  CI_high <- r_to_d(r + 2* r_sd)
  Elec <- c( 'C3' )
  d[['Somatotopic_Right_Dom']] <- data.table( Subjects = 14, Trials = 90, Elec = Elec, d = Cohen_d, CI = CI, CI_low = CI_low, CI_high = CI_high)
  

  ## From "Body-specific motor imagery of hand actions - neural evidence from right- and left-handers"
  # number of words = 96
  # half MAN half NONMAN: 96/2 
  # 32 subjects: 16 left-handed, 16-right handed
  # 4 groups: left-handers left hemisphere, left-handers right hem, 
  # right-handers left hem, right-handers right hem
  n <- 32/2
  # From Figure 3 (approximate values)
  # Right Hem Left-handers MAN > NONMAN and Left Hem Right-handers MAN > NONMAN
  Perc_signal_change <- c( 0.105, 0.17, 0.085  )
  # Standard Error
  Error_bar <- c( 0.01, 0.015, 0.005 )
  nrVoxels <- c( 619, 1052, 364 )

  # SD = SE * sqrt(N)
  Perc_signal_change_sd <- Error_bar*sqrt(n)

  # From "How to Analyze Change from Baseline: Absolute or Percentage Change?" equation 4
  ttest <- Perc_signal_change/Perc_signal_change_sd
  Elec <- c( 'FC2, FC4', 'CP4, C4', 'P8' )
  d[['Body_Left_Dom']] <- data.table( Subjects = 32, Trials = 96, Elec = Elec, t_to_d(ttest, 2*n-2))

  
  # From Figure 3 (approximate values)
  # Right Hem Left-handers MAN > NONMAN and Left Hem Right-handers MAN > NONMAN
  Perc_signal_change <- c( 0.12, 0.12, 0.12  )
  Error_bar <- c( 0.008, 0.015, 0.015 )
  nrVoxels <- c( 1022, 404, 1213 )

  # SD = SE * sqrt(N)
  Perc_signal_change_sd <- Error_bar*sqrt(n)

  # From "How to Analyze Change from Baseline: Absolute or Percentage Change?" equation 4
  ttest <- Perc_signal_change/Perc_signal_change_sd
  Elec <- c( 'FC1, F1', 'CP3', 'P7' )
  d[['Body_Right_Dom']] <- data.table( Subjects = 32, Trials = 96, Elec = Elec, t_to_d(ttest, 2*n-2))

  
  ## From "Motor imagery-based brain activity parallels that of motor execution - Evidence from magnetic source imaging of cortical oscillations"
  # By definition: ERS_ERD = (P_event - P_rest)/P_rest * 100%
  # From Figure 2
  ERS_ERD <- c( 1.5, -1.5, -2.5, 2, -3, -3, -1.5, -2.5, -2.5 )
  # Hypotheses from Figure 4
  Sigma <- 0.5
  P_rest <- 3
  # Calculations
  P_event <- abs(ERS_ERD*P_rest + P_rest)
  if( any(P_event==P_rest) ){print('Wrong hypothesis for P_rest')}
  Cohen_d <- abs(sqrt(P_event) - sqrt(P_rest))/Sigma
  # 16 subjects
  N <- 16
  n1 <- N/2
  n2 <- N/2
  # From "The Handbook of Research synthesis" p238
  d_var <- ((n1+n2)/(n1*n2) + (Cohen_d^2)/(2*(n1+n2-2)))*((n1+n2)/(n1+n2-2))
  d_sd <- sqrt(d_var)
  # z-score at 95% confidence interval
  CI <- rep(0.95, length(Cohen_d))
  # d +/- d_sd*Z_crit
  CI_low <- Cohen_d - d_sd*Z_crit
  CI_high <- Cohen_d + d_sd*Z_crit
  Elec <- c( 'O2', 'FC4', 'C4', 'O1', 'CP1, CP3', 'CP1, CP3', 'F1', 'CP3', 'FC1' )
  d[['Motor_Left_NonDom']] <- data.table( Subjects = 16, Trials = 21, Elec = Elec, d = Cohen_d, CI = CI, CI_low = CI_low, CI_high = CI_high)

  
  ## From "Multimodal functional imaging of motor imagery using a novel paradigm"
  # From Table 3
  q <- c( -42, -54 )
  q <- abs(q)/100 # in percentage
  # From Figure 6 (extremely approximate)
  q_sd <- 0.03

  # from https://www2.psych.ubc.ca/~schaller/528Readings/RosnowRosenthal2003.pdf
  # Fisher zr-transformed r
  r <- (exp(2*q)-1) / (exp(2*q)+1)
  r_sd <- (exp(2*q_sd)-1) / (exp(2*q_sd)+1)
  r_low  <- (exp(2*(q-2*q_sd))-1) / (exp(2*(q-2*q_sd))+1)
  r_high  <- (exp(2*(q+2*q_sd))-1) / (exp(2*(q+2*q_sd))+1)
  Cohen_d <- r_to_d(r)
  CI_low <- r_to_d(r_low)
  CI_high <- r_to_d(r_high)
  CI <- rep(0.95, length(r))
  Elec <- c( 'FC1', 'FC4' )
  d[['Multimodal_Left_NonDom']] <- data.table( Subjects = 14, Trials = 10, Elec = Elec, d = Cohen_d, CI = CI, CI_low = CI_low, CI_high = CI_high)

  # From Table 3
  q <- c( -69, -43 )
  q <- abs(q)/100 # in percentage
  # From Figure 6 (extremely approximate)
  q_sd <- 0.03

  # from https://www2.psych.ubc.ca/~schaller/528Readings/RosnowRosenthal2003.pdf
  # Fisher zr-transformed r
  r <- (exp(2*q)-1) / (exp(2*q)+1)
  r_sd <- (exp(2*q_sd)-1) / (exp(2*q_sd)+1)
  r_low  <- (exp(2*(q-2*q_sd))-1) / (exp(2*(q-2*q_sd))+1)
  r_high  <- (exp(2*(q+2*q_sd))-1) / (exp(2*(q+2*q_sd))+1)
  Cohen_d <- r_to_d(r)
  CI_low <- r_to_d(r_low)
  CI_high <- r_to_d(r_high)
  CI <- rep(0.95, length(r))
  Elec <- c(  'FC1', 'FC4' )
  d[['Multimodal_Right_Dom']] <- data.table( Subjects = 14, Trials = 10, Elec = Elec, d = Cohen_d, CI = CI, CI_low = CI_low, CI_high = CI_high)

  
  ## From "Increased motor cortex excitability during motor imagery in brain-computer interface trained subjects"
  # Educated guess from summary
  Elec <- c( 'C4, C6', 'C2', 'FC4, FC2', 'FCZ', 'F1', 'CP4' )
  d[['Increased_Right_Dom']] <- data.table( Subjects = 11, Trials = NA, Elec = Elec, d = rep(NA,length(Elec)),
                             CI = rep(NA,length(Elec)), CI_low = rep(NA,length(Elec)),
                             CI_high = rep(NA,length(Elec)))
  # Educated guess from summary
  Elec <- c( 'C5, C3', 'C1', 'FC3', 'FCZ', 'F5, F3', 'CP3, CP1' )
  d[['Increased_Left_NonDom']] <- data.table( Subjects = 11, Trials = NA, Elec = Elec, d = rep(NA,length(Elec)),
                                   CI = rep(NA,length(Elec)), CI_low = rep(NA,length(Elec)),
                                   CI_high = rep(NA,length(Elec)))


  checkElec(d, SumElec)

  Agg <- getAgg(d, SumElec)

  return(list(Agg, d))
}

# Calculate Cohen's d for Motor Execution
getEffectSizeME <- function(SumElec, Z_crit){
  d <- NULL
  ## From "Partially Overlapping Neural Networks for Real and Imagined Hand Movements"
  # From Table 1
  # two-tailed Z-score (standard normal deviate)
  # Executed movement with rest compared
  Z2 <- c( 5.62, 7.14, 6.55, 5.67, 5.19, 5.11, 4.37, 4.59, 6.01, 5.03, 4.94, 4,
           3.66, 4.95, 3.95, 3.93, 4.7, 5.73, 3.25, 4.2, 4.83, 4.67, 5.35, 6.23 )
  Elec <- c('F8', 'C2', 'C1', 'FC2, FC4', 'FC1', 'FCZ', 'FC3', 'FC2, FCZ',
            'FC1, FCZ', 'C4', 'C3', 'CP2', 'CP1', 'CP3', 'CP6', 'CP5', 'NP',
            'NP', 'NP', 'NP', 'NP', 'NP', 'NP', 'NP' )

  # eight right-handed healthy volunteers
  n <- 8
  ## For right-handers only
  d[['Partially_Both']] <- data.table( Subjects = 8, Trials = 96, Elec = Elec, z_to_d(Z2, n) ) # z_to_d(Z2, 2*n-2)

  
  ## From "Somatotopic mapping of the human primary sensorimotor cortex during motor imagery and motor execution by functional magnetic resonance imaging"
  # From Table 1
  r <- c( 0.83 )
  r_sd <- c( 0.05 )
  dS <- c(  3.13 )
  dS_sd <- c( 1.22 )

  Cohen_d <- r_to_d(r)
  CI <- rep(0.95, length(r))
  # Standard deviation: With probability about 95%  (r - 2 * SD; r + 2 * SD)
  CI_low <- r_to_d(r - 2* r_sd)
  CI_high <- r_to_d(r + 2* r_sd)
  Elec <- c(  'C2, C4' )
  d[['Somatotopic_Left_NonDom']] <- data.table( Subjects = 14, Trials = 90, Elec = Elec, d = Cohen_d, CI = CI, CI_low = CI_low, CI_high = CI_high)

  # From Table 1
  r <- c( 0.83 )
  r_sd <- c( 0.06 )
  dS <- c( 3.18 )
  dS_sd <- c( 1.74 )

  Cohen_d <- r_to_d(r)
  CI <- rep(0.95, length(r))
  # Standard deviation: With probability about 95%  (r - 2 * SD; r + 2 * SD)
  CI_low <- r_to_d(r - 2* r_sd)
  CI_high <- r_to_d(r + 2* r_sd)
  Elec <- c( 'C3, C1' )
  d[['Somatotopic_Right_Dom']] <- data.table( Subjects = 14, Trials = 90, Elec = Elec, d = Cohen_d, CI = CI, CI_low = CI_low, CI_high = CI_high)
  
  
  ## From "Motor imagery-based brain activity parallels that of motor execution - Evidence from magnetic source imaging of cortical oscillations"
  # By definition: ERS_ERD = (P_event - P_rest)/P_rest * 100%
  # From Figure 1
  ERS_ERD <- c(-2.5, -3, -3, -2.5)
  # Hypotheses from Figure 4
  Sigma <- 0.5
  P_rest <- 3
  # Calculations
  P_event <- abs(ERS_ERD*P_rest + P_rest)
  if( any(P_event==P_rest) ){print('Wrong hypothesis for P_rest')}
  Cohen_d <- abs(sqrt(P_event) - sqrt(P_rest))/Sigma
  # Variance of Cohen's d
  # 16 subjects, N=480-144=336 trials
  N <- 16  
  n1 <- N/2
  n2 <- N/2
  # From "The Handbook of Research synthesis" p238
  d_var <- ((n1+n2)/(n1*n2) + (Cohen_d^2)/(2*(n1+n2-2)))*((n1+n2)/(n1+n2-2))
  d_sd <- sqrt(d_var)
  # z-score at 95% confidence interval
  CI <- rep(0.95, length(Cohen_d))
  # d +/- d_sd*Z_crit
  CI_low <- Cohen_d - d_sd*Z_crit
  CI_high <- Cohen_d + d_sd*Z_crit
  Elec <- c( 'CPZ, CZ', 'C2', 'C2', 'C3, C1' )
  d[['Motor_Left_NonDom']] <- data.table( Subjects = 16, Trials = 21, Elec = Elec, d = Cohen_d, CI = CI, CI_low = CI_low, CI_high = CI_high)

  
  ## From "Multimodal functional imaging of motor imagery using a novel paradigm"
  # From Table 3
  # Left ME
  q <- c(-90, -97 )
  q <- abs(q)/100 # in percentage
  # From Figure 6 (extremely approximate)
  q_sd <- 0.03

  # From https://www2.psych.ubc.ca/~schaller/528Readings/RosnowRosenthal2003.pdf
  # Fisher zr-transformed r
  r <- (exp(2*q)-1) / (exp(2*q)+1)
  r_sd <- (exp(2*q_sd)-1) / (exp(2*q_sd)+1)
  r_low  <- (exp(2*(q-2*q_sd))-1) / (exp(2*(q-2*q_sd))+1)
  r_high  <- (exp(2*(q+2*q_sd))-1) / (exp(2*(q+2*q_sd))+1)
  Cohen_d <- r_to_d(r)
  CI_low <- r_to_d(r_low)
  CI_high <- r_to_d(r_high)
  CI <- rep(0.95, length(r))
  Elec <- c( 'FC1',  'FC4, FC2' )
  d[['Multimodal_Left_NonDom']] <- data.table( Subjects = 14, Trials = 10, Elec = Elec, d = Cohen_d, CI = CI, CI_low = CI_low, CI_high = CI_high)

  # From Table 3
  # Right ME
  q <- c( -87, -86 )
  q <- abs(q)/100 # in percentage
  # From Figure 6 (extremely approximate)
  q_sd <- 0.03

  # From https://www2.psych.ubc.ca/~schaller/528Readings/RosnowRosenthal2003.pdf
  # Fisher zr-transformed r
  r <- (exp(2*q)-1) / (exp(2*q)+1)
  r_sd <- (exp(2*q_sd)-1) / (exp(2*q_sd)+1)
  r_low  <- (exp(2*(q-2*q_sd))-1) / (exp(2*(q-2*q_sd))+1)
  r_high  <- (exp(2*(q+2*q_sd))-1) / (exp(2*(q+2*q_sd))+1)
  Cohen_d <- r_to_d(r)
  CI_low <- r_to_d(r_low)
  CI_high <- r_to_d(r_high)
  CI <- rep(0.95, length(r))
  Elec <- c( 'C1', 'FC4, FC2' )
  d[['Multimodal_Right_Dom']] <- data.table( Subjects = 14, Trials = 10, Elec = Elec, d = Cohen_d, CI = CI, CI_low = CI_low, CI_high = CI_high)

  
  ## From "Cortical Activations in Humans Grasp-Related Areas Depend on Hand Used and Handedness"
  # 34 subjects
  n <- 34
  # From Table 1
  # Main effect of performing hand
  # (Right-handers/LH+Left-handers/LH) > (Right-handers/RH+Left-handers/RH)
  ttest <- c( 8.16, 7.7 )
  Z <- c( 6.58, 6.32 )
  Elec <- c( 'C4', 'C2' )
  d[['Cortical_Left_Both']] <- data.table( Subjects = n, Trials = 140, Elec = Elec, t_to_d(ttest, 2*n-2))

  # From Table 1
  # Main effect of performing hand
  # (Right-handers/RH+Left-handers/RH) > (Right-handers/LH+Left-handers/LH)
  ttest <- c( 9.66, 8.81 )
  Z <- c( 7.36, 6.93 )
  Elec <- c( 'C1', 'C3' )
  d[['Cortical_Right_Both']] <- data.table( Subjects = 34, Trials = 140, Elec = Elec, t_to_d(ttest, 2*n-2))

  
  ## From "The Neural Correlates of Grasping in Left-Handers: When Handedness Does Not Matter"
  # 16 subjects
  n <- 16
  # From Table 2
  ttest <- c(10.39, 8.36, 6.66, 6.65, 7.78, 7.26, 5.67, 7.08, 6.77, 6.5 )
  Z <- c(6.72, 5.96, 5.18, 5.12, 5.71, 5.47, 4.64, 5.39, 5.23, 5.1 )
  Elec <- c( 'C2', 'C4', 'CP4', 'FC4', 'C3, C1', 'FC3', 'FC1', 'F8', 'FC5', 'CP3' )
  ## For left-handers only
  d[['Neural_Both']] <- data.table( Subjects = 16, Trials = 60, Elec = Elec, t_to_d(ttest, 2*n-2))

  
  checkElec(d, SumElec)

  Agg <- getAgg(d, SumElec)

  return(list(Agg, d))

}

# Calculate Cohen's for P300
getEffectSizeP300 <- function(SumElec, Z_crit){
  d <- NULL
  ## From "Effect of the green/blue flicker matrix for P300-based brain–computer interface: an EEG–fMRI study"
  # From Table 1
  # white/gray
  Z2 <- c( 4.71, 4.27, 3.74, 3.67, 3.66, 3.29 )
  Elec <- c( 'P7, PO7', 'PO8', 'PO4', 'P2', 'C2', 'PO3, O1' )
  
  # the residual artifacts were so severe that we had to discard the data
  # of three of the 12 subjects
  n <- 9
  ## For right-handers only
  d[['Effect']] <- data.table( Subjects = 9, Trials = 63, Elec = Elec, z_to_d(Z2, n) ) # z_to_d(Z2, 2*n-2)
  
  
  ## From "Brain-Computer Interface System Based on P300 Processing with Convolutional Neural Network, Novel Speller, and Low Number of Electrodes"
  # # From Figure 5
  # # Amplitude (microVolts)
  # AmpP300 <- c( 0.21492, 0.26217, 0.40276, 0.39004, 0.5779 )
  # AmpNoP300 <- c( 0.042623, 0.0085862, 0.066889, -0.014932, 0.053107 )
  # # Standard deviation (microVolts): calculated on no P300, and assumed to be equal to P300
  # AmpSD <- 0.06
  # 
  # Elec <- c( 'C3', 'P3', 'P4', 'O1', 'O2' )
  # 
  # d[['BCI']] <- data.table( Subjects = 8, Trials = 34, Elec = Elec, d = (AmpP300-AmpNoP300)/mean(AmpSD), 
  #                            CI = 0.95, 
  #                            CI_low = ((AmpP300-AmpNoP300)-Z_crit*mean(AmpSD))/mean(AmpSD),
  #                            CI_high = ((AmpP300-AmpNoP300)+Z_crit*mean(AmpSD))/mean(AmpSD) ) 
  # From data given by author
  Elec <- c( 'O1', 'OZ', 'O2' )
  d[['BCI']] <- data.table( Subjects = 8, Trials = 34, Elec = Elec, d = c(6.68453649495561,
                                                                     5.17625013672877,
                                                                     5.25275144271726 ),
                             CI = 0.95,
                             CI_low = c(4.03148567186625,
                                        3.02247199291775,
                                        3.07424483565272 ),
                             CI_high = c(9.30466614738737,
                                         7.29118769060386,
                                         7.39274476936806 ) )
  
  
  # ## From "Functional magnetic resonance imaging of brain activity in the visual oddball task"
  # # Not enough data found in paper
  
  
  ## From "Combined event-related fMRI and EEG evidence for temporal—parietal cortex activation during target detection"
  # From Table 1
  Z2 <- c(6.1, 4.79, 5.8 )
  Elec <- c( 'CP5, C3, P3', 'CP5, CP3, C3, P3', 'CP6, CP4, C4, P4' )
  # Number of subjects
  n <- 11

  ## For right-handers only
  d[['Combined']] <- data.table( Subjects = 11, Trials = 20, Elec = Elec, z_to_d(Z2, n) ) # z_to_d(Z2, 2*n-2)
  
  
  ## From "Responses to Rare Visual Target and Distractor Stimuli Using Event-Related fMRI"
  # From Figure 1B
  Z2 <- c(3, 3, 4.5, 3.5, 3.5, 5.5, 3.5, 5, 4, 4 )
  # From Table 1: Target Stimulus
  Elec <- c('OZ', 'P3', 'P2', 'CP5', 'C3, CP3', 'CP4, P4', 'FC4',
            'FCZ', 'F5, FC5', 'F8' )
  # number of subjects
  n <- 6
  ## For right-handers only
  # 8 runs per subject. each run 90 stimuli
  # rare target, the letter “X”: P = 0.09
  d[['Responses']] <- data.table( Subjects = 6, Trials = round((90*8)*0.09), Elec = Elec, z_to_d(Z2, n) ) # z_to_d(Z2, 2*n-2)
  
   
  
  # ## From "The neural basis of the P300 potential"
  ## Not enough data to calculate standard deviation
  # # From Table 3: Current source density
  # CSD <- c(2.06, 1.77, 1.47, 1.14, 1.13,
  #          2.02, 1.77, 1.69,
  #          1.5, 1.41, 1.4, 1.33,
  #          2, 1.11, 1.09, 1.09, 1.07,
  #          2.24, 1.79, 1.16, 1.06, 1.06 )
  # Elec <- c( 'FCZ', 'CP3', 'T7, TP7', 'CP6', 'FC1, FC3',
  #            'FT7, TP7', 'P6', 'FZ, FPZ',
  #            'CP5', 'PZ, POZ', 'CP6', 'FT7',
  #            'FCZ, FZ', 'TP7, T7', 'P7', 'CP3', 'CP6',
  #            'FCZ', 'CP3', 'FC1', 'CP6', 'CP4' )
  # 
  # 
  # d[['Neural']] <- data.table( Subjects = 50, Trials = 100, Elec = Elec, d = CSD/CSD_SD, CI = 0.95, 
  #                              CI_low = CSD-Z_crit*CSD_SD,
  #                              CI_high = CSD+Z_crit*CSD_SD ) 
  
  
  ## From "Localizing P300 Generators in Visual Target and Distractor Processing - A Combined Event-Related Potential and Functional Magnetic Resonance Imaging Study"
  # From Table 2: Target Amplitude (nanoamperes)
  TA <- c( 8.6, 11.5, 13.7, 16.1, 25.5, 44.2,
           7.2, 6.4, -10.4, -9.1 )
  # “Min” and “max” denote the borders of the 95% confidence interval
  TA_Min <- c( 4.9, 8.9, 10.3, 11, 19.5, 36,
               3.4, 3.6, -4.8, -6.7 )
  TA_Max <- c( 13.3, 15.1, 18.6, 21.8, 34.6, 53.9,
               11.5, 9.3, -18.1, -12.7 )
  
  Elec <- c( 'CP3', 'CP4', 'P1, P3', 'P2, P4', 'P7', 'P8',
             'F3, AF3', 'F4, AF4', 'FCZ', 'CP4, P6' )
  
  TA_SD <- (TA_Max-TA_Min)/2/Z_crit
  
  d[['Localizing']] <- data.table( Subjects = 10, Trials = 35, Elec = Elec, d = TA/TA_SD, CI = 0.95, 
                               CI_low = TA_Min/TA_SD, CI_high = TA_Max/TA_SD ) 
  
  
  
  checkElec(d, SumElec)
  
  Agg <- getAgg(d, SumElec)
  
  return(list(Agg, d))
  
}

# Create bold legend
bold_legend <- function(value) {
  x <- as.character(value)
  do.call(expression, lapply(x, function(y) if(str_detect(y, 'Mean')){
    diffNchar <- first(nchar(str_split(y, '\n')[[1]])) - last(nchar(str_split(y, '\n')[[1]]))
    if( abs(diffNchar)-2 >= 0 ){ Spaces <- paste0(rep(" ", abs(diffNchar)-2 ), collapse = ""); Spaces2 <- ""
    }else{Spaces <- ""; Spaces2 <- paste0(rep(" ", abs(diffNchar) ), collapse = ""); }
    bquote(bold(.(paste( Spaces , first(str_split(y, '\n')[[1]]), '\n', Spaces2, last(str_split(y, '\n')[[1]]) ))) )}else{y} ))
}

# Forest plot
ForestPlot <- function( toPlot , El, boolPlot = T, boolSave = F, nameFile = "" ){ 
  # red, blue, grey
  TextColor <- "#353238"
  myColors <- c("#92140C","#247BA0", "#353238")
  myColorsSig <- c( 'Dominant', 'Non Dominant', 'Both' )
  myLabels <- c( 'Left', 'Right', 'Dominant', 'Non Dominant', 'Mixed Dominant \n          and \n Non Dominant')
  myLabs <- c("ᐊ","ᐅ")
  
  # Add number of subjects under papers names
  xLegends <- NULL
  for( xl in 1:nrow(toPlot) ){ # xl <- 3
    if( !str_detect(toPlot$Study[xl], 'Mean') ){
      nbSpaces <- abs( nchar( toPlot$Study[xl] ) - nchar( paste0( 'N = ', toPlot$Subjects[xl] )  ) )
      if( round(nbSpaces/2) == 0 ){ nbSpaces <- 1 }else{ nbSpaces <- round(nbSpaces/2) }
      nbSpaces <- paste0( rep( ' ' , nbSpaces ) , collapse = "")
      if( nchar( toPlot$Study[xl] ) > nchar( paste0( 'N = ', toPlot$Subjects[xl] )  ) ){
        xLegends <- c( xLegends , paste0( toPlot$Study[xl], '\n', nbSpaces, 'N = ', toPlot$Subjects[xl]  )  )
      }else{
        xLegends <- c( xLegends , paste0( nbSpaces, toPlot$Study[xl],'\n',    'N = ', toPlot$Subjects[xl] )  )
      }
    }else{
      xLegends <- c( xLegends , toPlot$Study[xl]  )
    }
  }
  toPlot$Study <- xLegends
  colS <- NULL; for( xl in 1:nrow(toPlot) ){
    colS <- c(colS, myColors[which(myColorsSig == toPlot$Dom[xl])] ) }
  toPlot <- cbind(toPlot, Col = colS)
  
  # size texts: x axis, y axis labels and ticks
  sizeText <- 20
  # size text under triangles
  sizeText2 <- 5
  # size legend right
  sizeText3 <- 6
  # size symbol legend right
  sizeSymbol <- 7
  # size triangles
  sizeTriangle <- 10
  # size lines
  sizeLines <- 2
  
  Length <- length(myLabels)
  pos.x <- rep(0, Length )
  # To change position legend, change from = XX
  pos.y <- rev(seq(from = Length - 1, by = 0.7, length.out = Length ))
  pText <- ggplot( data.table(x = pos.x-0.1, y = pos.y) , aes(x,y)) + theme_void() +
    xlab(NULL) + ylab(NULL) + xlim(-0.5, 0.5) + ylim(0,ceiling(max(pos.y))) +
    geom_text(aes( pos.x, pos.y, label = myLabels ), colour = TextColor, size = sizeText3, hjust = 'left' ) +
    geom_point(col = c(rep(TextColor , 2), rep('NA' , 3) ),
               shape = c(myLabs, rep(NA,3)) ,
               size = sizeSymbol) + 
    geom_point(col = c(rep('NA' , 2), myColors ), 
               shape = c(rep(NA,2), rep(19,3)) ,
               size = sizeSymbol) + 
    theme(legend.position = 'none') 
  if( min(toPlot$CI_low) > 0 ){ RangeD <- 0 }else{ RangeD <- -ceiling(abs(min(toPlot$CI_low))/0.5)*0.5 }
  RangeD <- c( RangeD, ceiling(max(toPlot$CI_high)/0.5)*0.5 )
  
  # Vertical adjustment depending on whether there are two cohen's d to write on same line
  if( any(as.double(table(toPlot$Study)) != 1) ){
    repName <- names(table(toPlot$Study))[which(as.double(table(toPlot$Study)) != 1)]
    VerticalAdg <- rep(3, nrow(toPlot) )
    VerticalAdg[toPlot$Study == repName][2] <- 5
  }else{
    VerticalAdg <- 3
  }
  pPlot <- toPlot %>% ggplot( aes( x = factor( Study, 
                                               # Ordering so that y-axis is alphabeticall ordered
                                               # And Means are all at the end
                                               level = c( Study[str_detect(Study, 'Mean')] ,
                                                          rev(unique(Study[-which(str_detect(Study, 'Mean'))])[order(unique(Study[-which(str_detect(Study, 'Mean'))]))])
                                               ) ),
                                   y=d, ymin=CI_low, ymax=CI_high )) +
    geom_linerange(aes(colour = Dom), size = sizeLines ) +
    scale_color_manual(values = toPlot$Col ) +
    geom_text(aes(x=Study, y = d, colour = Dom, label = Symbols),
              vjust = 0.3, hjust = 0.5,
              size = sizeTriangle, family = "HiraKakuPro-W3") + 
    geom_text(aes(x=Study, y = d, label = Text, vjust = VerticalAdg ), size = sizeText2, color = TextColor ) +  
    geom_hline(yintercept = 0, linetype=2) +
    coord_flip() +
    xlab('Study') +
    ylab(paste0("Cohen's d for ", El)) + 
    # Mean in bold
    scale_x_discrete( labels =  bold_legend ) +
    scale_y_continuous( limits = RangeD ) + 
    theme_light(base_size = sizeText, base_family = "Roboto Condensed") +
    theme(legend.position = 'none', panel.border = element_blank(), axis.text.y = element_text(hjust = 0.5) ) # ,
  
  if(any(toPlot$CI_low<0)){ print('Negative values!')}
  
  g <-  ggarrange( pPlot , pText, ncol = 2, nrow = 1, widths = c(3,1) )
  if(boolPlot){ dev.new(width = 40, height = 20, units = "cm" ); print(g) }
  if( boolSave && !file.exists(nameFile) ){ ggsave(file = nameFile, g , width = 40, height = 20, units = "cm" , dpi = 300  ) }
  
}

# Meta-analysis
metaAnalysis <- function(nameTask, Agg, d, Z_crit, Distinction = F, boolPlot = F, boolSave = F, savePath){ 
  savePath <- file.path(savePath, paste0( 'Cohen_' , nameTask))
  if(boolSave && !dir.exists(savePath)){ dir.create(savePath) }
  
  for(nr in 1:nrow(Agg)){ # nr <- 1
    El <- Agg$Elec[nr]
    nA <- matrix( unlist(Agg[nr,2:ncol(Agg),with=F]) , ncol = 3, byrow = TRUE)
    Study <- Hand <- Dom <- Sub <- Trials <- NULL; 
    for(nd in names(d)){ # nd <- names(d)[1]
      if( length(str_split(nd, '_')[[1]]) > 1 ){
        Study <- c(Study, str_split(nd, '_')[[1]][1] )
        Hand <- c(Hand, str_split(nd, '_')[[1]][2] )
        Dom <- c(Dom, str_split(nd, '_')[[1]][3] )
      }else{
        Study <- c(Study, nd)
        Hand <- c(Hand, 'NA')
        Dom <- c(Dom, 'NA')   
      }
      Sub <- c(Sub, unique(d[[nd]]$Subjects) )
      Trials <- c(Trials, unique(d[[nd]]$Trials) )
    }
    # Changing Dom so the labels are c("Both", "Dominant", "Non Dominant")
    if(any(is.na(Dom))){ Dom[which(is.na(Dom))] <- 'Both' }
    Dom[which(Dom=='Dom')] <- "Dominant"
    Dom[which(Dom=='NonDom')] <- "Non Dominant"
    if( nrow(data.table(nA)) == 0 ){ next }
    nA <- cbind(Study = Study, Subjects = Sub, Trials = Trials, Hand = Hand, Dom = Dom, data.table(nA))
    names(nA)[(ncol(nA)-2):ncol(nA)] <- c('d', 'CI_low', 'CI_high')
    nA <- as.data.table(nA[which(!is.na(nA$d))])
    Symbols <- NULL; for(nh in nA$Hand){ 
      if(nh == "Left"){ Symbols <- c(Symbols , "◀" )
      }else if(nh == "Right"){ Symbols <- c(Symbols , "▶" )
      }else{ Symbols <- c(Symbols , "◀▶" ) } }
    # the 95% confidence interval is 3.92 standard errors wide (3.92 = 2 × 1.96).
    # CI = Cohen_d +/- SD*Z_crit
    # Standard error = Standard Dev / sqrt(N)
    # SD = sqrt(N) * (upper limit - lower limit) / (Z_crit*2)
    # SE = (upper limit - lower limit) / (Z_crit*2)
    myText <- NULL; for(nr in 1:nrow(nA)){ myText <- c( myText, 
                                                        paste0(round(nA$d[nr], digits = 2), ' [',
                                                               round(nA$CI_low[nr], digits = 2), ' ; ', 
                                                               round(nA$CI_high[nr], digits = 2), '] ' ) ) }
    
    nA <- cbind(nA, data.table(Text = myText, SE = ( nA$CI_high - nA$CI_low )/(Z_crit*2) , Symbols = Symbols ) )
    # Number of studies providing effect size for this particular electrode
    nbStudies <- nrow( nA[!is.na(nA$d)] )
    
    if( nbStudies > 0 ){
      print(El)
      SimpleMean <- WeightedMean <- NULL
      
      if( Distinction ){
        if( str_detect(nameTask, 'M') ){
          for(ha in unique(nA$Hand)[order(unique(nA$Hand))] ){ # ha <- unique(nA$Hand)[order(unique(nA$Hand))][1]
            if( ha == 'Both' ){
              nASub <- nA
            }else{
              nASub <- subset(nA, Hand == ha )
            }
            if( nrow(nASub) > 1 ){
              if( length(unique(nASub$Hand))==1 ){ H <- unique(nASub$Hand) }else{ H <- 'Both'}
              if( H == 'Left'){ S <- "◀" }else if( H == 'Right'){ S <- "▶" }else{ S <- "◀▶" }
              if( length(unique(nASub$Dom))==1 ){ D <- unique(nASub$Dom) }else{ D <- 'Both'}
              # Simple mean effect size
              SimpleMean <- rbind(SimpleMean, 
                                  data.table( Study = paste0('Simple Mean \nEffect Size ', ha), 
                                              Subjects = sum(nASub$Sub),
                                              Hand = H, Dom = D ,
                                              d = mean(nASub$d), CI_low = mean(nASub$CI_low), 
                                              CI_high = mean(nASub$CI_high),
                                              Text = paste0(round(mean(nASub$d),digits = 2), ' [', 
                                                            round(mean(nASub$CI_low),digits = 2), ' ; ',
                                                            round(mean(nASub$CI_high),digits = 2), ']' ), 
                                              NbStudies = nbStudies,
                                              SE = mean(nASub$SE), Symbols = S ) )
              
              # Weighted mean effect size
              W <- nASub$Sub
              WeightedMean <- rbind( WeightedMean ,
                                     data.table( Study = paste0('Weighted Mean \nEffect Size ', ha), 
                                                 Subjects = sum(W),
                                                 Hand = H, Dom = D ,
                                                 d = weighted.mean( nASub$d, W ), 
                                                 CI_low = weighted.mean( nASub$CI_low, W ), 
                                                 CI_high = weighted.mean( nASub$CI_high, W ),
                                                 Text = paste0(round(weighted.mean( nASub$d, W ),digits = 2), ' [', 
                                                               round(weighted.mean( nASub$CI_low, W ),digits = 2), ' ; ',
                                                               round(weighted.mean( nASub$CI_high, W ),digits = 2), ']' ),
                                                 NbStudies = nbStudies,
                                                 SE = weighted.mean( nASub$SE, W ), Symbols = S ) )
            }else{
              SimpleMean <- rbind(SimpleMean, 
                                  data.table( Study = paste0('Simple Mean \nEffect Size ', ha), nASub[,c(2,4:(ncol(nASub)-2)), with=F],
                                              NbStudies = nrow(nASub), nASub[,(ncol(nASub)-1):ncol(nASub),with=F] ) )
              WeightedMean <- rbind( WeightedMean ,
                                     data.table( Study = paste0('Weighted Mean \nEffect Size ', ha),nASub[,c(2,4:(ncol(nASub)-2)), with=F],
                                                 NbStudies = nrow(nASub), nASub[,(ncol(nASub)-1):ncol(nASub),with=F] ) )
            }
          }
        }else{
          # Simple mean effect size
          SimpleMean <- data.table( Study = 'Simple Mean \nEffect Size All', 
                                    Subjects = nrow(nA),
                                    Hand = 'Both', Dom = 'Both' ,
                                    d = mean(nA$d), CI_low = mean(nA$CI_low), CI_high = mean(nA$CI_high),
                                    Text = paste0(round(mean(nA$d),digits = 2), ' [', 
                                                  round(mean(nA$CI_low),digits = 2), ' ; ',
                                                  round(mean(nA$CI_high),digits = 2), ']' ), 
                                    NbStudies = nbStudies,
                                    SE = mean(nA$SE), Symbols = "◀▶" )
          # Weighted mean effect size
          # if( !is.null(nA$Sub) && !is.null(nA$Trials) ){ W <- nA$Sub*nA$Trials
          if( !is.null(nA$Sub) ){ W <- nA$Sub
          }else{ W <- 1 }
          WeightedMean <- data.table( Study = 'Weighted Mean \nEffect Size All', 
                                      Subjects = sum(W),
                                      Hand = 'Both', Dom = 'Both' ,
                                      d = weighted.mean( nA$d, W ), 
                                      CI_low = weighted.mean( nA$CI_low, W ), 
                                      CI_high = weighted.mean( nA$CI_high, W ),
                                      Text = paste0(round(weighted.mean( nA$d, W ),digits = 2), ' [', 
                                                    round(weighted.mean( nA$CI_low, W ),digits = 2), ' ; ',
                                                    round(weighted.mean( nA$CI_high, W ),digits = 2), ']' ), 
                                      NbStudies = nbStudies,
                                      SE = weighted.mean( nA$SE, W ), Symbols = "◀▶" )
        }
      }else{
        # Simple mean effect size
        SimpleMean <- data.table( Study = 'Simple Mean \nEffect Size All', 
                                  Subjects = nrow(nA),
                                  Hand = 'Both', Dom = 'Both' ,
                                  d = mean(nA$d), CI_low = mean(nA$CI_low), CI_high = mean(nA$CI_high),
                                  Text = paste0(round(mean(nA$d),digits = 2), ' [', 
                                                round(mean(nA$CI_low),digits = 2), ' ; ',
                                                round(mean(nA$CI_high),digits = 2), ']' ), 
                                  NbStudies = nbStudies,
                                  SE = mean(nA$SE), Symbols = "◀▶" )
        # Weighted mean effect size
        # if( !is.null(nA$Sub) && !is.null(nA$Trials) ){ W <- nA$Sub*nA$Trials
        if( !is.null(nA$Sub) ){ W <- nA$Sub
        }else{ W <- 1 }
        WeightedMean <- data.table( Study = 'Weighted Mean \nEffect Size All', 
                                    Subjects = sum(W),
                                    Hand = 'Both', Dom = 'Both' ,
                                    d = weighted.mean( nA$d, W ), 
                                    CI_low = weighted.mean( nA$CI_low, W ), 
                                    CI_high = weighted.mean( nA$CI_high, W ),
                                    Text = paste0(round(weighted.mean( nA$d, W ),digits = 2), ' [', 
                                                  round(weighted.mean( nA$CI_low, W ),digits = 2), ' ; ',
                                                  round(weighted.mean( nA$CI_high, W ),digits = 2), ']' ), 
                                    NbStudies = nbStudies,
                                    SE = weighted.mean( nA$SE, W ), Symbols = "◀▶" )
      }
      # If no distinction is made
      nameFile <- file.path(savePath, paste0('Simple Mean ', nameTask, ' ', El, Png))
      ForestPlot(rbind( nA, SimpleMean, fill = TRUE ) , El , boolPlot, boolSave , nameFile = nameFile)
      
      nameFile <- file.path(savePath, paste0('Weighted Mean ', nameTask, ' ', El, Png))
      ForestPlot( rbind( nA, WeightedMean, fill = TRUE ) , El , boolPlot , boolSave , nameFile = nameFile)
      
      if( !exists('AllSimpleMean') ){
        AllSimpleMean <- cbind(Elec = El, SimpleMean )
        AllWeightedMean <- cbind(Elec = El, WeightedMean )
      }else{
        AllSimpleMean <- rbind( AllSimpleMean, cbind(Elec = El, SimpleMean ) )
        AllWeightedMean <- rbind( AllWeightedMean,cbind(Elec = El, WeightedMean ) )
      }
    } # End loop if(nbStudies > 1)
  } # End loop nr
  if( exists('AllSimpleMean')){
    nameFile <- file.path(savePath, paste0('Simple Mean ', nameTask, Csv))
    if( boolSave && !file.exists(nameFile) ){ write.csv(AllSimpleMean, file = nameFile)}
    nameFile <- file.path(savePath, paste0('Weighted Mean ', nameTask, Csv))
    if( boolSave && !file.exists(nameFile)){ write.csv(AllWeightedMean, file = nameFile)}
  }else{
    # print(d[[1]])
    AllSimpleMean <- AllWeightedMean <- NULL
    nameFile <- file.path(savePath, paste0("Cohen's d ", nameTask, Csv))
    # print(nameFile)
    if( boolSave && !file.exists(nameFile) ){ write.csv(d[[1]], file = nameFile)}
  }
  return(list(AllSimpleMean, AllWeightedMean))
}

# For motor-related task, distinction between hand used
Distinction <- T

### Motor Imagery
SumElec <- c( 'FP2', 'F1', 'F5, F3', 'F8', 'F7', 'FC4, FC2',
              'FC3', 'FCZ', 'FCZ', 'FC2', 'FC1', 'C2', 'C1',
              'C4, C6', 'C5, C3', 'P2, CP2', 'P1',
              'CP4', 'CP3, CP1', 'TP7, T7',
              'P8', 'P7', 'PO3', 'O2', 'O1' )
Res <- getEffectSizeMI(SumElec, Z_crit)
MI <- Res[[1]]
d <- Res[[2]]
Means_MI <- metaAnalysis('MI', MI, d, Z_crit, Distinction, boolPlot = F, boolSave = T, savePath)


### Motor Execution
SumElec <- c( 'F8', 'FC5', 'FC4, FC2', 'FC1, FC3', 'FCZ', 'FCZ',
              'FCZ', 'FCZ', 'C2', 'C1', 'C4', 'C3', 'CPZ, CZ', 'CP2',
              'CP1', 'CP4', 'CP3', 'CP6', 'CP5' )
Res <- getEffectSizeME(SumElec, Z_crit)
ME <- Res[[1]]
d <- Res[[2]]
Means_ME <- metaAnalysis('ME', ME, d, Z_crit, Distinction, boolPlot = F, boolSave = T, savePath)

## P300
SumElec <- c( 'F5, FC5', 'F3, AF3', 'F4, AF4', 'F8', 'FC1, FC3',
              'FC4', 'FZ, FPZ', 'FCZ', 'FCZ', 'C2',
              'C4','C3','CP3', 'CP4, P4', 'P1, P3', 'P2, P4', 
              'P7', 'P8', 'T7, TP7', 'P7, PO7', 'PO8', 'CP4, P6',
              'FT7, TP7', 'CP5', 'CP6, P6', 'PO4, O2',
              'PO3, O1', 'OZ' )
Res <- getEffectSizeP300(SumElec, Z_crit)
P300 <- Res[[1]]
d <- Res[[2]]
Means_P300 <- metaAnalysis('P300', P300, d, Z_crit, Distinction = F, boolPlot = F, boolSave = T, savePath)
