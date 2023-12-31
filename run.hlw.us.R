## PREAMBLE: clear screen/workspace -----
cat("\014"); rm(list = ls()); gc();
# SET DEFAULTS: display options, font and y axis label rotation
options(digits = 14); options(scipen = 999);  options(max.print=10000)
# INSTALL PACMAN: if not installed (note: may neeD:\_teaching\_current.teaching\_SU.EFF\code-EFF\helper_functions\print_results.Rd to disable windows firewall for packages to install)
if(!"pacman" %in% installed.packages()){install.packages("pacman")}
# LOAD/INSTALL: other required packages
pacman::p_load(tis,mFilter,nloptr,openxlsx)
# LOAD HELPER FUNCTIONS STORED IN functions_path
functions_path  = c("./helper_functions/")
invisible(lapply( paste0(functions_path, list.files(functions_path, "*.R")), source ))

# NOTE: the sample dates MUST correspond to data in input file

# Set the start and end dates of the estimation sample (format is c(year,quarter))
sample.start <- c(1961,1)
sample.end   <- c(2019,4)

# The estimation process uses data beginning 4 quarters prior to the sample start
data.start    <- shiftQuarter(sample.start,-4)

# Initialization of state vector and covariance matrix
# Set as NA to follow procedure in HLW paper
# Or can input values manually
xi.00.stage1 <- NA
xi.00.stage2 <- NA
xi.00.stage3 <- NA

P.00.stage1 <- NA
P.00.stage2 <- NA
P.00.stage3 <- NA

# Upper bound on a_3 parameter (slope of the IS curve)
a.r.constraint <- -0.0025

# Lower bound on b_2 parameter (slope of the Phillips curve)
b.y.constraint <- 0.025

# Set start index for g.pot series; used in state vector initialization
g.pot.start.index <- 1 + ti(shiftQuarter(sample.start,-3),'quarterly')-ti(data.start,'quarterly')

# Because the MC standard error procedure is time consuming, we include a run switch
# Set run.se to TRUE to run the procedure
run.se <- FALSE

# Set number of iterations for Monte Carlo standard error procedure
niter <- 5

# =================
# COVID-ADJUSTED MODEL SETTINGS
# =================

# Set to TRUE if using time-varying volatility; FALSE if not
# Must specify kappa.inputs if TRUE
use.kappa <- TRUE

# fix.phi must be set at NA or a numeric value
# Set as NA to estimate the COVID indicator coefficient
# Set at a numeric value to fix phi at that value (e.g. 0)
fix.phi <- NA


# =================
# VARIANCE SCALE PARAMETERS
# =================

# SETTINGS:
# kappa.inputs DESCRIPTIONS:
# name: used as label in param.num
# year: assumes kappa applies to full year, unless manually adjusted
# T.start: time series index start; will be set in subsequent loop for YYYY:Q1
# T.end: time series index end; will be set in subsequent loop for YYYY:Q4
# init: value to initialize kappa in parameter estimation; default of 1
# lower.bound : lower bound for kappa in maximum likelihood estimation; default 1
# upper.bound : upper bound for kappa in maximum likelihood estimation; default Inf (no bound)
# theta.index: leave as NA; will be filled in within each stage

# NOTE: fix kappa at value by setting lower.bound=upper.bound=value

kappa.inputs <- data.frame('name'=c('kappa2020Q2-Q4','kappa2021','kappa2022'),
                           'year'=c(2020,2021,2022),
                           'T.start'=c(NA,NA,NA),
                           'T.end'=c(NA,NA,NA),
                           'init'=c(1,1,1),
                           'lower.bound'=c(1,1,1),
                           'upper.bound'=c(Inf,Inf,Inf),
                           'theta.index'=c(NA,NA,NA))

# NOTE: Sets Q1-Q4 of years provided
if (use.kappa) {

  # Number of kappas introduced
  n.kappa <- dim(kappa.inputs)[1]
  for (k in 1:n.kappa) {
    # Indexing to start of y_t vector
    covid.variance.start.yq <- c(kappa.inputs$year[k],1) - sample.start

    kappa.inputs$T.start[k] <- max(covid.variance.start.yq[1]*4 + covid.variance.start.yq[2] +1,0)

    covid.variance.end.yq <- c(kappa.inputs$year[k],4) - sample.start

    kappa.inputs$T.end[k] <- max(covid.variance.end.yq[1]*4 + covid.variance.end.yq[2] +1,0)

    rm(covid.variance.start.yq, covid.variance.end.yq)

    # Manual adjustment to start Kappa_2020 in second quarter
    # Comment out under alternative specifications
    if (kappa.inputs$year[k]==2020) {
      kappa.inputs$T.start[k] <- kappa.inputs$T.start[k] + 1
    }
  }
}


# =================
# INPUT DATA
# =================

# Read input data from FRBNY website
us.data <- read.xlsx("inputData/Holston_Laubach_Williams_current_estimates.xlsx", sheet="US input data",
                      na.strings = ".", colNames=TRUE, rowNames=FALSE, detectDates = TRUE)

us.log.output             <- us.data$gdp.log
us.inflation              <- us.data$inflation
us.inflation.expectations <- us.data$inflation.expectations
us.nominal.interest.rate  <- us.data$interest
us.real.interest.rate     <- us.nominal.interest.rate - us.inflation.expectations
us.covid.indicator        <- us.data$covid.ind


# =================
# ESTIMATION
# =================

us.estimation <- run.hlw.estimation(log.output=us.log.output,
                                    inflation=us.inflation,
                                    real.interest.rate=us.real.interest.rate,
                                    nominal.interest.rate=us.nominal.interest.rate,
                                    covid.indicator=us.covid.indicator,
                                    a.r.constraint=a.r.constraint,
                                    b.y.constraint=b.y.constraint,
                                    g.pot.start.index=g.pot.start.index,
                                    use.kappa=use.kappa,
                                    kappa.inputs=kappa.inputs,
                                    fix.phi=fix.phi,
                                    xi.00.stage1=xi.00.stage1,
                                    xi.00.stage2=xi.00.stage2,
                                    xi.00.stage3=xi.00.stage3,
                                    P.00.stage1=P.00.stage1,
                                    P.00.stage2=P.00.stage2,
                                    P.00.stage3=P.00.stage3,
                                    run.se=run.se,
                                    sample.end=sample.end)

# One-sided (filtered) estimates
one.sided.est.us <- cbind(us.estimation$out.stage3$rstar.filtered,
                          us.estimation$out.stage3$trend.filtered,
                          us.estimation$out.stage3$z.filtered,
                          us.estimation$out.stage3$output.gap.filtered)


# =================
# OUTPUT
# =================

# Set up output for export
output.us <- format.output(country.estimation=us.estimation,
                           one.sided.est.country=one.sided.est.us,
                           real.rate.country=us.real.interest.rate,
                           start=sample.start,
                           end=sample.end,
                           run.se=run.se)

# Save output to CSV
# write.table(output.us, 'output/output.us.csv', col.names=TRUE, quote=FALSE, row.names=FALSE, sep = ',', na = '')
# 
# 
# #### WRITE FULL KF STATES TO CSV to check if KF of g(t) and g(t-1) as well as z(t) and z(t-1) are identical -----
# KFstates = us.estimation$out.stage3$states            # make KF data.frame to print to csv with headers and dates
# KF.df = as.data.frame(KFstates$filtered$xi.tt)        # remove xi rownames
# rownames(KF.df) = NULL                                # add proper column names for the states
# colnames(KF.df) = c("y*(t)","y*(t-1)","y*(t-2)","g(t)","g(t-1)","g(t-2)","z(t)","z(t-1)","z(t-2)")
# KF = data.frame(Date = us.data$date[5:nrow(us.data)]) # make df and join
# KF = cbind(KF,KF.df)                                  # write to CSV
# write.csv(KF,"output/KF.csv",row.names=FALSE)
# 
# # repeat for Smoothed states
# KS.df = as.data.frame(KFstates$smoothed$xi.tT)
# rownames(KS.df) = NULL
# colnames(KS.df) = c("y*(t)","y*(t-1)","y*(t-2)","g(t)","g(t-1)","g(t-2)","z(t)","z(t-1)","z(t-2)")
# KS = data.frame(Date = us.data$date[5:nrow(us.data)])
# KS = cbind(KS,KS.df)
# write.csv(KS,"output/KS.csv",row.names=FALSE)


# store extra output
# # Signal-to-noise Ratios
lambda.g	= 0.072769547
lambda.z	= 0.020609208
# write.table(out.stage3$xi.00,"xi00.csv", sep = ",", row.names = FALSE, col.names = FALSE)
# write.table(out.stage3$P.00,"P00.csv",   sep = ",", row.names = FALSE, col.names = FALSE)
xi.00 = as.matrix(read.csv("./init.Vals/xi00.csv", header=FALSE))
P.00  = as.matrix(read.csv("./init.Vals/P00.csv" , header=FALSE))
# P00  = 0.2*diag(9)
