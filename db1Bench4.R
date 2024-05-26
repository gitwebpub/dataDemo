
p.model.d <- function(model) {
  # Extract data used in the model
  data <- model.frame(model)
  oldpar <- par(mfrow = c(1, 2))
  
  # Standardized residuals plot
  plot(rstandard(model), main = 'standardRes by Index')
  text(x = 1:length(rstandard(model)), 
       y = rstandard(model), 
       labels = rownames(data), 
       pos = 3, cex = 0.6)
  
  # Hat values plot
  plot(hatvalues(model), main = 'hatVals by Index')
  text(x = 1:length(hatvalues(model)), 
       y = hatvalues(model), 
       labels = rownames(data), 
       pos = 3, cex = 0.6)
  
  # Calculate the mean of hat values
  mean_hatvalues <- mean(hatvalues(model))
  # Draw horizontal lines at mean and 2 * mean of hat values
  abline(h = mean_hatvalues, col = "red", lty = 2)
  abline(h = 2 * mean_hatvalues, col = "blue", lty = 2)
  
  # Reset to old parameters
  par(oldpar)
}

# Define the function to calculate influence measures
table.influence <- function(model) {
  # Calculate required statistics
  standard_res <- rstandard(model)
  hat_values <- hatvalues(model)
  cooks_distance <- cooks.distance(model)
  residuals <- rstudent(model)
  dfits_values <- residuals * sqrt(hat_values / (1 - hat_values))
  # Calculate Hadi's influence measure using olsrr package
  hadi_values <- olsrr::ols_hadi(model)$hadi
  # Combine the statistics into a data frame
  influence_measures <- data.frame(
    Std_Res = standard_res,
    p_ii = hat_values,
    C_i = cooks_distance,
    DFITS_i = dfits_values,
    Hadi = hadi_values
  )
  return(influence_measures)
}


# Define the function to plot influence measures
plot.inf <- function(model) {
  # Calculate required statistics
  standard_res <- rstandard(model)
  hat_values <- hatvalues(model)
  cooks_distance <- cooks.distance(model)
  dfits_values <- standard_res * sqrt(hat_values / (1 - hat_values))
  hadi_values <- olsrr::ols_hadi(model)$hadi
  
  # Combine the statistics into a data frame
  cal <- data.frame(
    Std_Res = standard_res,
    p_ii = hat_values,
    C_i = cooks_distance,
    DFITS_i = dfits_values,
    Hadi = hadi_values
  )
  
  par(mfrow = c(2, 2))
  n <- nrow(cal)
  Param.Model <- length(coef(model))
  dfits_threshold <- 2 * sqrt(Param.Model / (n - Param.Model))
  
  # p : number of parameters
  # Plot Cook's Distance
  plot(cal$C_i,
       xlab = "Index", ylab = "Cook's Distance", 
       main = "Cook's Distance")
  text(x = 1:n, 
       y = cal$C_i, 
       labels = 1:n, 
       pos = 1, cex = 0.6)
  
  # Plot DFITS
  plot(cal$DFITS_i,
       xlab = "Index", ylab = "DFITS",
       main = "DFITS")
  abline(h = c(dfits_threshold, -dfits_threshold), col = "red", lty = 2)
  # abline(h = 2*sqrt(length(cal$DFITS_i)/nrow(t.diag)), col = "red", lty = 2)
  # abline(h = -2*sqrt(length(cal$DFITS_i)/nrow(t.diag)), col = "red", lty = 2)
  text(x = 1:nrow(t.diag),
       y = cal$DFITS_i,
       labels = 1:n,
       pos = 1, cex = 0.6)
  # p : number of parameters
  DFITSth <- 2*sqrt((Param.Model)/(n-Param.Model))
  abline(h = DFITSth, col = "red", lty = 2)
  abline(h = -DFITSth, col = "red", lty = 2)
  
  # Plot Hadi's influence measure
  plot(cal$Hadi,
       xlab = "Index", ylab = "Hadi's Measure",
       main = "Hadi's Influence Measure")
  #  abline(h = 1, col = "red", lty = 2)
  #  abline(h = 0.5, col = "blue", lty = 2)
  text(x = 1:n,
       y = cal$Hadi,
       labels = 1:n,
       pos = 1, cex = 0.6)
  
  # Plot PR plot
  Pii <- cal$p_ii
  HX <- sqrt(Pii / (1 - Pii)) # Potential
  HY <- cal$Hadi - HX # Residual
  plot(HY, HX, 
       xlab = "Residual (Y)", 
       ylab = "Potential (X)", 
       main = "Potential-Residual plot")
  text(HY, HX, labels = 1:n, pos = 1, cex = 0.6)
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
}
