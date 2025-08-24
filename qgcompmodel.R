#### 1. Load Required Packages ####
library(qgcomp)      # For Quantile G-Computation mixture analysis
library(ggplot2)     # For visualization
library(dplyr)       # For data manipulation
library(tibble)      # For tidy data frames

#### 2. Read and Explore Data ####
data <- read.csv("ASD-Co exposure.csv") # Import dataset
head(data)                        # Preview first few rows
str(data)                         # Check data structure

#### 3. Define Exposure Variables ####
exposure <- names(data)[1:5]      # Select first 5 columns as exposures

#### 4. Function to Extract Weights from qgcomp Model ####
extract_weights <- function(fit) {
  # Compatible with qgcomp.noboot and qgcomp.boot
  if (!is.null(fit$weights)) {
    w <- fit$weights
    tibble(
      exposure = names(w),
      value = as.numeric(w)
    )
  } else if (!is.null(fit$pos.weights) && !is.null(fit$neg.weights)) {
    expnms <- fit$expnms
    npos <- length(fit$pos.weights)
    nneg <- length(fit$neg.weights)
    tibble(
      exposure = c(expnms[1:npos], expnms[(npos+1):(npos+nneg)]),
      value = c(fit$pos.weights, -fit$neg.weights)
    )
  } else {
    stop("Unrecognized weight structure")
  }
}

#### 5. Fit qgcomp Model for Continuous Outcome (with Covariates) ####
fit1 <- qgcomp.noboot(
  CARS ~ PM10 + PM25 + CO + SO2 + NO2 + Age + Sex + BMI + Race + Edu,
  expnms = exposure,
  data = data,
  family = gaussian(),
  q = 4
)
print(fit1)    # Print model summary
plot(fit1)     # Default plot

#### 6. Beautify the Weights Plot ####
wdat1 <- extract_weights(fit1) %>% arrange(desc(abs(value)))
ggplot(wdat1, aes(x = reorder(exposure, value), y = value, fill = value > 0)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c("#E57373", "#64B5F6")) + # Red for negative, blue for positive
  labs(
    title = "Beautified: Exposure Weights for Outcome (qgcomp)",
    x = "Exposure",
    y = "Weight"
  ) +
  theme_minimal(base_size = 14)
