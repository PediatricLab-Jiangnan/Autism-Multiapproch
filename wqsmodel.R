#### 1. Load Required Packages ####
library(gWQS)      # For Weighted Quantile Sum regression
library(ggplot2)   # For visualization

#### 2. Read and Explore Data ####
data <- read.csv("ASD-Co exposure.csv") # Import dataset
head(data)                        # Preview first few rows
str(data)                         # Check data structure

#### 3. Define Exposure Mixture Variables ####
# Select the first 5 columns as mixture exposures (modify as needed)
mixed_exposure <- names(data)[1:5]

#### 4. Fit the WQS Regression Model ####
wqsmodel1 <- gwqs(
  formula = CARS ~ wqs + Age + Sex + BMI + Race + Edu, # The formula must include 'wqs'
  mix_name = mixed_exposure,    # Specify mixture exposures
  data = data,                  # Input data
  q = 4,                        # Number of quantiles (quartiles)
  validation = 0.6,             # 60% data for training, 40% for validation
  b = 100,                      # Number of bootstrap samples
  b1_pos = TRUE,                # Assume positive association
  b_constr = FALSE,             # No direction constraint on weights
  family = gaussian,            # For continuous outcome; use 'binomial' for binary
  seed = 1234                   # Set random seed for reproducibility
)
# Note: For binary outcome, replace 'CARS' with your binary variable and set family=binomial

#### 5. Model Summary and Results ####
summary(wqsmodel1)          # Show model summary
confint(wqsmodel1)          # Confidence intervals for coefficients
wqsmodel1$final_weights     # Show component weights

#### 6. Extract and Sort Weights ####
weights <- wqsmodel1$final_weights
weights <- weights[order(-weights$mean_weight), ] # Sort by descending weight

#### 7. Calculate Weight Threshold ####
# The threshold is the reciprocal of the number of mixture components
tau <- 1 / length(weights$mix_name)

#### 8. Plot Component Weights (Barplot) ####
ggplot(weights, aes(x = reorder(mix_name, mean_weight), y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity", color = "black") +    # Draw barplot with black border
  theme_bw() +                                      # Black and white theme
  theme(
    axis.ticks = element_blank(),                   # Remove axis ticks
    axis.text.x = element_text(color = 'black'),    # X-axis text color
    axis.title.x = element_text(size = 14, margin = margin(t = 10)), # X-axis title style
    axis.title.y = element_text(size = 14, margin = margin(r = 10)), # Y-axis title style
    legend.position = "none"                        # Hide legend
  ) +
  coord_flip() +                                    # Flip coordinates for horizontal barplot
  labs(
    title = "WQS Component Weight Distribution",     # Plot title
    x = "Pollutant",                                # X-axis label
    y = "Estimated Weight"                          # Y-axis label
  ) +
  # Add red dashed threshold line
  geom_hline(yintercept = tau, linetype = "dashed", color = "red", alpha = 0.7)
