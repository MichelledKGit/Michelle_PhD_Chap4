library(tinytex)
library(sp)
library(spData)
library(sf)
library(spdep)
library(spatialreg)
library(gstat)
library(readxl)

data_full<-filepathtodata

data_full$x <- as.numeric(data_full$x)
data_full$y <- as.numeric(data_full$y)

# Convert to factors in the training dataset (data)
data_full$combination_number <- as.factor(data_full$combination_number)

data <- data_full[!is.na(data_full$average_sales), ]
data <- data_full[!is.na(data_full$Hours_Open), ]

sps<-SpatialPoints(data[,c("x","y")],proj4string = CRS("+proj=utm +zone=32")) 
spst <- spTransform(sps, CRS("+proj=longlat +datum=WGS84")) 
data[, c("long", "lat")] <- coordinates(spst) 
cords=cbind(data$long,data$lat)

# Distance calculations and weight matrices
Eucl.dist <- as.matrix(dist(cords))

Inv_w<-1/Eucl.dist
diag(Inv_w) <- 0 
rs_inv <- rowSums(Inv_w); 
Inv_w <- apply(Inv_w, 2, function(q) q/rs_inv)


# Assuming `combination_number` is a vector of categorical values for each location
combination_number <- data$combination_number  # or replace with your actual combination_number vector

# Create binary difference matrix for the categorical covariate
N <- length(combination_number)
pairwise_diff_matrix <- matrix(0, N, N)

for (i in 1:N) {
  for (j in 1:N) {
    pairwise_diff_matrix[i, j] <- ifelse(combination_number[i] == combination_number[j], 0, 1)
  }
}

rs_inv <- rowSums(Inv_w); 
Inv_w <- apply(Inv_w, 2, function(q) q/rs_inv)

pairwise_diff_matrix_inv=1/pairwise_diff_matrix
pairwise_diff_matrix_inv[which(!is.finite(pairwise_diff_matrix_inv))] <- 0

# Row-standardize Eucl.dist
row_sums_Eucl <- rowSums(Eucl.dist)  # Calculate the sum of each row for Eucl.dist
standardized_Eucl <- apply(Eucl.dist, 2, function(q) q/row_sums_Eucl)  # Divide each row by its sum
standardized_Eucl[is.na(standardized_Eucl)] <- 0  # Handle NaNs if any rows sum to zero

# Row-standardize pairwise_diff_matrix
row_sums_pairwise <- rowSums(pairwise_diff_matrix)  # Calculate the sum of each row for pairwise_diff_matrix
standardized_pairwise <- apply(pairwise_diff_matrix, 2, function(q) q/row_sums_pairwise)  # Divide each row by its sum
standardized_pairwise[is.na(standardized_pairwise)] <- 0  # Handle NaNs if any rows sum to zero

# Display the standardized matrices
standardized_Eucl<-standardized_Eucl*1000
standardized_pairwise<-standardized_pairwise*1000
Inv_w<-Inv_w*1000

##################################################################################################
########################################     Moran's I    ########################################
##################################################################################################

coordinates(data) <- ~x+y
nb <- knn2nb(knearneigh(coordinates(data), k = 5))  # 4-nearest neighbours as an example
lw <- nb2listw(nb, style = "W")

glm_model <- lm(log(data$average_sales)~combination_number+Hours_Open, data = data)  # Replace other_covariates with your actual covariates

residuals_glm <- residuals(glm_model)
data$residuals <- residuals_glm

moran_test <- moran.test(residuals_glm, lw)

print(moran_test)

spplot(data, "residuals", 
       col.regions = colorRampPalette(c("blue", "white", "red"))(100),
       at = seq(min(data$residuals), max(data$residuals), length.out = 100),
       main = "OLS model residual",
       sp.layout = list("sp.points", data, pch = 16, cex = abs(data$residuals) * 0.5))


###################################################################################################
############ Calculating the R square at every value of alpha ####################################
###################################################################################################

# Initialize vectors to store alpha values and corresponding R-squared values
alpha_values <- seq(0, 1, by = 0.05)
r_squared_values <- numeric(length(alpha_values))

# Loop over alpha values
for (i in seq_along(alpha_values)) {
  alpha <- alpha_values[i]
  
  # Compute weight matrix W1
  W1 <- exp(- (alpha * standardized_pairwise + (1 - alpha) * Inv_w))
  diag(W1) <- 0
  rs_w1 <- rowSums(W1)
  W1 <- apply(W1, 2, function(x) ifelse(rs_w1 == 0, 0, x / rs_w1))  # Avoid division by zero
  
  # Fit null and covariate models using spatial autoregressive model
  pw1_null <- spautolm(log(average_sales) ~ 1, family = "SAR", listw = mat2listw(W1, style = 'W'), data = data)
  pw1_cov <- spautolm(log(average_sales) ~ combination_number + Hours_Open, family = "SAR", listw = mat2listw(W1, style = 'W'), data = data)
  
  # Calculate and store R-squared value
  r_squared <- (1 - exp(-(2 / nrow(data)) * (logLik(pw1_cov) - logLik(pw1_null)))) * 100
  r_squared_values[i] <- r_squared
  
  # Print output for each alpha
  cat(paste("Alpha:", alpha, "\n"))
  cat(paste("PW1-R2:", round(r_squared, 3), "\n"))
  cat(paste("PW1-AIC:", round(AIC(pw1_cov), 3), "\n"))
  cat(paste("PW1-RMSE:", round(sqrt(pw1_cov$fit$SSE / nrow(data)), 3), "\n"))
}


# Find the maximum R-squared value and corresponding alpha
max_r_squared <- max(r_squared_values, na.rm = TRUE)
max_r_squared_alpha <- alpha_values[which.max(r_squared_values)]

# Plot R-squared values against alpha
plot(alpha_values, r_squared_values, type = "b", pch = 19, col = "blue",
     xlab = "Alpha", ylab = "R-Squared Value", main = "R-Squared vs Alpha",
     cex = 1.2, lwd = 2)
grid()  # Add a grid for better readability

# Add a horizontal line at the maximum R-squared value
abline(h = max_r_squared, col = "red", lty = 2, lwd = 2)

# Add a vertical line at the alpha value corresponding to the maximum R-squared
abline(v = max_r_squared_alpha, col = "red", lty = 2, lwd = 2)

# Annotate the maximum R-squared value and corresponding alpha on the plot
text(0.5, max_r_squared, labels = paste("Max R² =", round(max_r_squared, 3)),
     pos = 3, col = "red")
text(max_r_squared_alpha, min(r_squared_values), 
     labels = paste("Alpha =", round(max_r_squared_alpha, 2)), pos = 4, col = "red")


###############################################################################
########################### Overlaying the p-values with the R square ########
###############################################################################

# Disable scientific notation globally
options(scipen = 999)

# Initialize vectors to store alpha values, corresponding R-squared values and p-values
p_values <- numeric(length(alpha_values))  # To store p-values

# Loop over alpha values
for (i in seq_along(alpha_values)) {
  alpha <- alpha_values[i]
  
  # Compute weight matrix W1
  W1 <- exp(- (alpha * standardized_pairwise + (1 - alpha) * Inv_w))
  diag(W1) <- 0
  rs_w1 <- rowSums(W1)
  W1 <- apply(W1, 2, function(x) ifelse(rs_w1 == 0, 0, x / rs_w1))  # Avoid division by zero
  
  # Fit null and covariate models using spatial autoregressive model
  pw1_null <- spautolm(log(average_sales) ~ 1, family = "SAR", listw = mat2listw(W1, style = 'W'), data = data)
  pw1_cov <- spautolm(log(average_sales) ~ combination_number + Hours_Open, family = "SAR", listw = mat2listw(W1, style = 'W'), data = data)
  
  # Calculate and store R-squared value
  r_squared <- (1 - exp(-(2 / nrow(data)) * (logLik(pw1_cov) - logLik(pw1_null)))) * 100
  r_squared_values[i] <- r_squared
  
  # Calculate and store the p-value for the Likelihood Ratio Test (LRT)
  lrt_statistic <- 2 * (logLik(pw1_cov) - logLik(pw1_null))
  k <- length(levels(data$combination_number))  # Number of levels in the categorical variable
  df <- k - 1  # Degrees of freedom for the likelihood ratio test
  p_value <- 1 - pchisq(lrt_statistic, df)  # p-value for the LRT
  p_values[i] <- p_value
  
  # Print output for each alpha
  cat(paste("Alpha:", alpha, "\n"))
  cat(paste("PW1-R2:", round(r_squared, 3), "\n"))
  cat(paste("PW1-p-value:", format(p_value, scientific = FALSE), "\n"))  # Ensure p-value is not in scientific notation
}

# Find the maximum R-squared value and corresponding alpha
max_r_squared <- max(r_squared_values, na.rm = TRUE)
max_r_squared_alpha <- alpha_values[which.max(r_squared_values)]

# Adjust the plot margins to ensure labels don't overlap
par(mar = c(5, 5, 4, 4) + 0.1)  # Adjust margins (left and right)

# Plot R-squared values against alpha
plot(alpha_values, r_squared_values, type = "b", pch = 19, col = "blue",
     xlab = "Alpha", ylab = "R-Squared", 
     cex = 1.2, lwd = 2, xaxt = "n")  # Disable default x-axis for custom ticks

# Customize the x-axis to show every 0.1 mark
axis(1, at = seq(0, 1, by = 0.1))

# Add a grid for better readability
grid()

# Add a dotted horizontal line at p-value = 0.05 to show the significance threshold
abline(h = 40.85, col = "red", lty = 2, lwd = 2)

# Add a vertical line at the alpha value corresponding to the maximum R-squared
abline(v = max_r_squared_alpha, col = "dark green", lty = 2, lwd = 2)

# Annotate the maximum R-squared value and corresponding alpha on the plot
text(x = 0.5, y = 40.85, labels = "p-value = 0.05", col = "red", pos = 3, cex = 0.8)
#text(max_r_squared_alpha, min(r_squared_values), 
#     labels = paste("Alpha =", round(max_r_squared_alpha, 2)), pos = 4, col = "red")

# Add p-value to a secondary y-axis
par(new = TRUE)
plot(alpha_values, p_values, type = "b", pch = 19, col = "black", 
     xlab = "", ylab = "", yaxt = "n", lwd = 2)  # Changed to black

# Add a legend to the plot
legend("topright", legend = c("R-Squared", "p-value(LRT)"), col = c("blue", "black"), pch = 19, lty = 1, cex = 0.8)


#############################################################################################
################ Creating output where the refernece category changes ######################
#############################################################################################

# Add a new column 'attseg' based on 'combination_number' using base R
data$attseg <- NA  # Create an empty column first

data$attseg[data$combination_number == 2] <- "_PNH"
data$attseg[data$combination_number == 5] <- "_PHH"
data$attseg[data$combination_number == 7] <- "_PNM"
data$attseg[data$combination_number == 8] <- "_PHM"

data$OpHours <- data$Hours_Open

data$attseg <- as.factor(data$attseg)

############################################################################################
############### Change reference category to get the 4 different graphs ####################
############################################################################################

data$attseg <- relevel(data$attseg, ref = "_PNM")

alpha=0.65
W1 <- exp(- (alpha * standardized_pairwise + (1 - alpha) * Inv_w))
diag(W1) <- 0
rs_w1 <- rowSums(W1) 
W1 <- apply(W1, 2, function(x) x/rs_w1)

pw1_null<-spautolm(log(average_sales) ~1,family="SAR", listw =mat2listw(W1, style ='W'), data=data)
pw1_cov<-spautolm(log(average_sales) ~ OpHours + attseg,family="SAR", listw =mat2listw(W1, style ='W'), data=data)


# Print the summary of the model to check if coefficients are being extracted
summary(pw1_cov)

# Extract the coefficients using the coef() function
model_coef <- coef(pw1_cov)

# Check the extracted coefficients
print(model_coef)

# Extract coefficients for 'combination_number' categories (excluding the intercept)
coef_values <- model_coef[-1]  # Excluding the intercept

# Get the names of the coefficients (i.e., the variable names)
coef_names <- names(coef_values)

# Store the results in a data frame
coef_results <- data.frame(Covariates = coef_names, Estimate = coef_values)

coef_results <- coef_results[coef_results$Covariates != "lambda", ]

# Creating the bar plot to visualize the coefficients for each reference group
library(ggplot2)

ggplot(coef_results, aes(x = Covariates, y = Estimate, fill = Covariates)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  # Bar plot without error bars
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +  # Darker line at y = 0.0 (using linewidth instead of size)
  theme_minimal() +
  labs(
       x = "Covariates", y = "Coefficient Estimate") +
  scale_fill_manual(values = c("red", "green", "blue", "purple")) +  # Color for different reference categories
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


