library(DBI)
library(odbc)
library(rlang)
library(geoR) # For variogram analysis
library(gstat) # For kriging
library(sp)
library(maps)
library(maptools)
library(spatstat)
library(sf)
library(tmap)
library(tmaptool)
library(readxl)
library(dplyr)
library(spdep)
library(gstat)
library(plyr)
library(gridExtra)
library(ggplot2)
library(spaMM)
library(RSpectra)

######################### Max temp ##########################################

data  <- filepathtodata

coordinates(data) <- ~ X+Y

data_new_pop<-data[which(is.na(data$y_max) == FALSE),]
data_new_notpop<-data[which(is.na(data$y_max) == TRUE),]

data.var <- variogram(log(y_max)~1, data=data_new_pop)
plot(data.var)
data.fit <- fit.variogram(data.var, vgm(psill = 0.004, model = "Gau", range = 0.5, nugget = 0.0005)) #Matter
plot(data.var, model=data.fit)

# Compute the variogram
data.var <- variogram(log(y_max) ~ 1, data = data_new_pop)

# Fit the variogram model
data.fit <- fit.variogram(data.var, vgm(psill = 0.004, model = "Gau", range = 0.5, nugget = 0.0005)) 

# Plot the empirical and fitted variogram
plot(data.var, model = data.fit, main = "Variogram for Maximum Temperatures")

# Add a label to the plot using annotation
par(new = TRUE)  # Ensure we're adding to the existing plot
text(x = 0.1, y = 0.003, 
     labels = "Gaussian Model\nPartial Sill: 0.004\nRange: 0.5\nNugget: 0.0005", 
     adj = c(0, 0), cex = 0.8, col = "blue")


lzn.kriged <- krige(log(y_max)~1,data_new_pop,data_new_notpop,model=data.fit)
print(lzn.kriged)

bt <- exp( lzn.kriged@data$var1.pred + (lzn.kriged@data$var1.var / 2) )
mu_bt <- mean(bt)
mu_original <- mean(data_new_pop$y_max)

print(mu_bt)
print(mu_original)

# correction to remove kriging bias in sample mean
btt <- bt * (mu_original/mu_bt)             # correct backtransfomed vals
lzn.kriged@data$var1.pred <- btt                    # overwrite w/ correct vals 
names(lzn.kriged) <- c("y_max", "variance")
spplot(lzn.kriged, "y_max", main = "HH estimates")


data_new_notpop <- data_new_notpop[,!names(data_new_notpop) %in% c("y_max")]

krig.output <- cbind(data_new_notpop,lzn.kriged[1])


Full_y_max <- rbind(krig.output,data_new_pop)


library(tmap)

# Convert Full_y_max to sf object
Full_y_max_sf <- st_as_sf(Full_y_max, coords = c("X", "Y"), crs = 4326)

# Set tmap mode to interactive
tmap_mode("view")

# Create custom colour palette and breaks
custom_palette <- c("blue", "dodgerblue","yellow", "#FF8C00", "red")
breaks <- c(33, 34, 35, 36, 37, 38)  # Define custom breaks for y_max values

# Create the heatmap with custom colours
tm_shape(Full_y_max_sf) +
  tm_dots(
    col = "y_max",                     # Column to use for coloring
    palette = custom_palette,          # Custom palette
    breaks = breaks,                   # Define value ranges for the palette
    size = 0.1,                        # Dot size
    alpha = 0.5,                       # Transparency for the dots
    title = "Maximum temperature in °C"             # Legend title
  ) +
  tm_layout(
    title = "Observed and predicted maximum temperatures",
    legend.outside = TRUE              # Place the legend outside the map
  )


######################### Min temp ##########################################
data  <- filepathtodata

coordinates(data) <- ~ X+Y

data_new_pop<-data[which(is.na(data$y_min) == FALSE),]
data_new_notpop<-data[which(is.na(data$y_min) == TRUE),]



data.var <- variogram(log(y_min)~1, data=data_new_pop)
plot(data.var)
data.fit <- fit.variogram(data.var, vgm(psill = 0.004, model = "Exp", range = 0.25, nugget = 0.0005)) #Matter
plot(data.var, model=data.fit)



# Compute the variogram for minimum temperatures
data.var <- variogram(log(y_min) ~ 1, data = data_new_pop)

# Fit the variogram model using the Exponential structure
data.fit <- fit.variogram(data.var, vgm(psill = 0.004, model = "Exp", range = 0.25, nugget = 0.0005))

# Plot the empirical and fitted variogram
plot(data.var, model = data.fit, main = "Variogram for Minimum Temperatures")

# Add a label to the plot using annotation
text(x = 0.1, y = 0.003, 
     labels = "Exponential Model\nPartial Sill: 0.004\nRange: 0.25\nNugget: 0.0005", 
     adj = c(0, 0), cex = 0.8, col = "darkgreen")



lzn.kriged <- krige(log(y_min)~1,data_new_pop,data_new_notpop,model=data.fit)
print(lzn.kriged)

bt <- exp( lzn.kriged@data$var1.pred + (lzn.kriged@data$var1.var / 2) )
mu_bt <- mean(bt)
mu_original <- mean(data_new_pop$y_min)

print(mu_bt)
print(mu_original)

# correction to remove kriging bias in sample mean
btt <- bt * (mu_original/mu_bt)             # correct backtransfomed vals
lzn.kriged@data$var1.pred <- btt                    # overwrite w/ correct vals 
names(lzn.kriged) <- c("y_min", "variance")
spplot(lzn.kriged, "y_min", main = "HH estimates")


data_new_notpop <- data_new_notpop[,!names(data_new_notpop) %in% c("y_min")]

krig.output <- cbind(data_new_notpop,lzn.kriged[1])


Full_y_min <- rbind(krig.output,data_new_pop)


library(tmap)

# Convert Full_y_min to sf object
Full_y_min_sf <- st_as_sf(Full_y_min, coords = c("X", "Y"), crs = 4326)

# Set tmap mode to interactive
tmap_mode("view")

# Create custom colour palette and breaks
custom_palette <- c("blue", "dodgerblue","yellow", "orange")
breaks <- c(23, 24,25, 26, 27)  # Define custom breaks for y_min values

# Create the heatmap with custom colours
tm_shape(Full_y_min_sf) +
  tm_dots(
    col = "y_min",                     # Column to use for coloring
    palette = custom_palette,          # Custom palette
    breaks = breaks,                   # Define value ranges for the palette
    size = 0.1,                        # Dot size
    alpha = 0.5,                       # Transparency for the dots
    title = "Minimum temperature in °C"             # Legend title
  ) +
  tm_layout(
    title = "Observed and predicted minimum temperatures",
    legend.outside = TRUE              # Place the legend outside the map
  )





