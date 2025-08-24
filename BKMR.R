#### 1. Load Required Packages ####
# If 'bkmr' is not installed, uncomment the next line to install it
# devtools::install_github("jenfb/bkmr")

setwd("Your Own") # Set working directory

# Load essential libraries
library(bkmr)    # BKMR modeling
library(fields)  # For cover.design function
library(ggplot2) # For data visualization

#### 2. Read and Explore Data ####
dat <- read.csv("ASD-Co exposure.csv")   # Import dataset
head(dat)                         # Preview first few rows
str(dat)                          # Check data structure

#### 3. Data Preparation for BKMR ####

# 3.1. Exposure Variables (First 5 columns: PM10, PM25, CO, SO2, NO2)
mixture <- as.matrix(dat[, 1:5])
colnames(mixture) <- c("PM10", "PM25", "CO", "SO2", "NO2")
mixture <- scale(mixture) # Standardize exposure variables

# 3.2. Outcome Variable (CARS column, continuous outcome)
outcome <- dat$CARS

# 3.3. Covariates (Columns 6-10: Age, Sex, BMI, Race, Edu)
covariates <- as.matrix(dat[, 6:10])
colnames(covariates) <- c("Age", "Sex", "BMI", "Race", "Edu")

# 3.4. Check Data Dimensions
cat("Exposure matrix dimensions:", dim(mixture), "\n")
cat("Outcome variable length:", length(outcome), "\n") 
cat("Covariate matrix dimensions:", dim(covariates), "\n")

# 3.5. Check for Missing Values
cat("Missing values in exposure:", sum(is.na(mixture)), "\n")
cat("Missing values in outcome:", sum(is.na(outcome)), "\n")
cat("Missing values in covariates:", sum(is.na(covariates)), "\n")

#### 4. Set BKMR Knots to Speed Up Model ####
set.seed(123) # For reproducibility
n_samples <- nrow(mixture)
cat("Sample size:", n_samples, "\n")
nd <- max(10, round(n_samples / 10)) # Number of knots: at least 10 or 1/10 of samples
knots <- cover.design(mixture, nd = nd)$design

#### 5. Fit BKMR Model ####
fit_bkmr <- kmbayes(
  y = outcome, 
  Z = mixture, 
  X = covariates, 
  iter = 5000,            # Number of iterations (increase for real analysis)
  family = "gaussian",    # For continuous outcomes; use "binomial" for binary
  verbose = FALSE, 
  varsel = TRUE,          # Enable variable selection
  knots = knots           # Use specified knots
)

#### 6. Extract Variable Importance (Posterior Inclusion Probabilities, PIP) ####
pips <- ExtractPIPs(fit_bkmr)
print(pips)

#### 7. Plot Univariate Exposure-Response Relationships ####
pred.resp.univar <- PredictorResponseUnivar(
  fit = fit_bkmr,
  q.fixed = 0.5           # Fix other exposures at their median
)

ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  geom_smooth(stat = "identity") +
  ylab("h(z)") +
  facet_wrap(~variable)

#### 8. Bivariate Exposure-Response and Quantile Level Plots ####
# Bivariate exposure-response
pred.resp.bivar <- PredictorResponseBivar(
  fit = fit_bkmr,
  min.plot.dist = 1,
  q.fixed = 0.5
)

# Plot at specific quantile levels
pred.resp.bivar.levels <- PredictorResponseBivarLevels(
  pred.resp.df = pred.resp.bivar,
  Z = mixture,
  qs = c(0.10, 0.5, 0.90)
)

ggplot(pred.resp.bivar.levels, aes(z1, est)) +
  geom_smooth(aes(col = quantile), stat = "identity") +
  facet_grid(variable2 ~ variable1) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1") + 
  theme_minimal()

#### 9. Summarize Overall Mixture Risk Effects ####
risks.overall <- OverallRiskSummaries(
  fit = fit_bkmr,
  qs = seq(0.25, 0.75, by = 0.05),
  q.fixed = 0.5
)

p1 <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96 * sd, ymax = est + 1.96 * sd)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_pointrange()
print(p1)

#### 10. Summarize Single Variable Risk Effects ####
risks.singvar <- SingVarRiskSummaries(
  fit = fit_bkmr,
  y = outcome, 
  Z = mixture, 
  X = covariates, 
  qs.diff = c(0.25, 0.75), 
  q.fixed = c(0.25, 0.50, 0.75), 
  method = "exact"
)

ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96 * sd, 
                          ymax = est + 1.96 * sd, 
                          col = q.fixed)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_pointrange(position = position_dodge(width = 0.75)) +
  coord_flip()
