# ---
# title: "Gibbs Regression Test"
# author: "mcl"
# date: "April 22, 2016"
# output: html_document
# ---
#   
#   With Gibbs sampling, we draw from the joint conditional distribution of unknowns. In our case, the unknowns include the model parameters as well as the data. Here we test the effectiveness of using Gibbs sampling in assorted levels of spatial uncertainty and model error.
# 
# Import libraries
# ```{r}
library(RCurl)
library(ggplot2)
library(sp)
library(rgeos)
library(maptools)
library(mvtnorm)
library(MCMCpack)
#```

# Probability model
# 
# For each project location, we construct an area of coverage where the project location may in fact be. We assume that for any region of interest, the probability that the location is in the roi is equal to the area of overlap between the region of interest and the area of coverage divided by the area of coverage. The parameters function describes the joint probability distribution across regions of interest for each project location. dollar_expected_values and dollar_covariances calculate the expected value and the variance-covariance matrix for the joint probability distribution of aid.
# 
# ```{r}
parameters <- function(Aid.Data,P1_shp,P2_shp,P3_shp,P45_shp,P68_shp){ 
  
  Aid.Data <- Aid.Data[order(Aid.Data$project_location_id),] 
  
  # Fill in missing precision codes
  Aid.Data$latitude[Aid.Data$precision_code == 6] <- 0.30000
  Aid.Data$longitude[Aid.Data$precision_code == 6] <- 32.77500
  Aid.Data$latitude[Aid.Data$precision_code == 8] <- 0.30000
  Aid.Data$longitude[Aid.Data$precision_code == 8] <- 32.77500
  Aid.Data <- Aid.Data[(!is.na(Aid.Data$precision_code)),]
  
  Aid.PrjLoc <- Aid.Data
  
  # number of projects and rois
  num_rois <- length(P1_shp$NAME_3)
  num_locations <- nrow(Aid.PrjLoc)
  
  # dollar amounts
  dollar_set <- Aid.PrjLoc$even_split_commitments
  
  # project precisions
  precision_set <- Aid.PrjLoc$precision_code
  
  # number locations
  num_loc_set <- sample(x=1,prob=c(1),replace=T,size=num_locations)
  
  # locus of placement
  P1_shp_locus <- P1_shp
  P1_shp_locus$Dist.ID <- row.names(P1_shp_locus)
  P1_shp_locus$Dist.ID <- as.integer(P1_shp_locus$Dist.ID)
  
  coordinates(Aid.PrjLoc) <- ~longitude+latitude
  Aid.PrjLoc.extract <- over(Aid.PrjLoc,P1_shp_locus)
  locus_set <- subset(Aid.PrjLoc.extract, select = c(Dist.ID))
  
  ##### Defining Probability of Landing in Each ROI #####
  # Calculating Reliability for Project-Based Analysis
  # Using districts (N = 75) as the unit of analysis
  
  ##### Subsetting Datasets 
  Aid.Data <- Aid.Data[!(is.na(Aid.Data$precision_code)),]
  
  # Precision Code 1
  Aid.Data.P1 <- Aid.Data[Aid.Data$precision_code == 1,]
  
  # Precision Code 2
  Aid.Data.P2 <- Aid.Data[Aid.Data$precision_code == 2,]
  
  # Precision Code 3
  Aid.Data.P3 <- Aid.Data[Aid.Data$precision_code == 3,]
  
  # Precision Code 45
  Aid.Data.P45 <- Aid.Data[Aid.Data$precision_code == 4 | Aid.Data$precision_code == 5,]
  
  # Precision Code 68
  Aid.Data.P68 <- Aid.Data[Aid.Data$precision_code == 6 | Aid.Data$precision_code == 8,]
  
  ##### Calculating Spatial Area of Polygons #  ,,,
  # Calculations are by current projection
  P68_shp$Area <- gArea(P68_shp, byid=TRUE)
  P45_shp$Area <- gArea(P45_shp, byid=TRUE) 
  P3_shp$Area <- gArea(P3_shp, byid=TRUE) 
  P2_shp$Area <- gArea(P2_shp, byid=TRUE) 
  P1_shp$Area <- gArea(P1_shp, byid=TRUE) 
  
  # Joining Location to Points #
  # NOTE: Making assuming that smaller administrative zones are
  # completely embedded in larger zones. For example, a district
  # is not split between two zones. 
  
  ##### Projecting Datasets 
  
  
  # PRECISION CODE 1
  if (nrow(Aid.Data.P1) > 0){
    coordinates(Aid.Data.P1) <- ~longitude+latitude
    #Aid.Data.P1$Area_project <- over(Aid.Data.P1,P1_shp)$Area  
    Aid.Data.P1$ADM_p1 <- over(Aid.Data.P1,P1_shp)$NAME_3
    Aid.Data.P1$ADM_p2 <- over(Aid.Data.P1,P2_shp)$NAME_2
    Aid.Data.P1$ADM_p3 <- over(Aid.Data.P1,P3_shp)$NAME_1
    Aid.Data.P1$ADM_p45 <- over(Aid.Data.P1,P45_shp)$Adm_Region
    Aid.Data.P1$ADM_p68 <- "Country"
    
    Aid.Data.P1 <- subset(Aid.Data.P1@data,select = c(project_id,
                                                      project_location_id,
                                                      precision_code,
                                                      total_commitments,
                                                      total_disbursements,
                                                      ADM_p1,
                                                      ADM_p2,
                                                      ADM_p3,
                                                      ADM_p45,
                                                      ADM_p68))
    
  }
  
  # PRECISION CODE 2
  if (nrow(Aid.Data.P2) > 0){
    coordinates(Aid.Data.P2) <- ~longitude+latitude
    #Aid.Data.P2$Area_project <- over(Aid.Data.P2,P2_shp)$Area  
    Aid.Data.P2$ADM_p1 <- ""
    Aid.Data.P2$ADM_p2 <- over(Aid.Data.P2,P2_shp)$NAME_2
    Aid.Data.P2$ADM_p3 <- over(Aid.Data.P2,P3_shp)$NAME_1
    Aid.Data.P2$ADM_p45 <- over(Aid.Data.P2,P45_shp)$Adm_Region
    Aid.Data.P2$ADM_p68 <- "Country"
    
    Aid.Data.P2 <- subset(Aid.Data.P2@data,select = c(project_id,
                                                      project_location_id,
                                                      precision_code,
                                                      total_commitments,
                                                      total_disbursements,
                                                      ADM_p1,
                                                      ADM_p2,
                                                      ADM_p3,
                                                      ADM_p45,
                                                      ADM_p68))
    
  }
  
  # PRECISION CODE 3
  if (nrow(Aid.Data.P3) > 0){
    coordinates(Aid.Data.P3) <- ~longitude+latitude
    #Aid.Data.P3$Area_project <- over(Aid.Data.P3,P3_shp)$Area  
    Aid.Data.P3$ADM_p1 <- ""
    Aid.Data.P3$ADM_p2 <- ""
    Aid.Data.P3$ADM_p3 <- over(Aid.Data.P3,P3_shp)$NAME_1
    Aid.Data.P3$ADM_p45 <- over(Aid.Data.P3,P45_shp)$Adm_Region
    Aid.Data.P3$ADM_p68 <- "Country" 
    
    Aid.Data.P3 <- subset(Aid.Data.P3@data,select = c(project_id,
                                                      project_location_id,
                                                      precision_code,
                                                      total_commitments,
                                                      total_disbursements,
                                                      ADM_p1,
                                                      ADM_p2,
                                                      ADM_p3,
                                                      ADM_p45,
                                                      ADM_p68))                                                    
  }
  
  # PRECISION CODE 4
  if (nrow(Aid.Data.P45) > 0){
    coordinates(Aid.Data.P45) <- ~longitude+latitude
    #Aid.Data.P45$Area_project <- over(Aid.Data.P45,P45_shp)$Area  
    Aid.Data.P45$ADM_p1 <- ""
    Aid.Data.P45$ADM_p2 <- ""
    Aid.Data.P45$ADM_p3 <- ""
    Aid.Data.P45$ADM_p45 <- over(Aid.Data.P45,P45_shp)$Adm_Region
    Aid.Data.P45$ADM_p68 <- "Country" 
    
    Aid.Data.P45 <- subset(Aid.Data.P45@data,select = c(project_id,
                                                        project_location_id,
                                                        precision_code,
                                                        total_commitments,
                                                        total_disbursements,
                                                        ADM_p1,
                                                        ADM_p2,
                                                        ADM_p3,
                                                        ADM_p45,
                                                        ADM_p68))                                                    
  }
  
  # PRECISION CODE 6 AND 8
  if (nrow(Aid.Data.P68) > 0){
    coordinates(Aid.Data.P68) <- ~longitude+latitude
    #Aid.Data.P68$Area_project <- over(Aid.Data.P68,P68_shp)$Area  
    Aid.Data.P68$ADM_p1 <- ""
    Aid.Data.P68$ADM_p2 <- ""
    Aid.Data.P68$ADM_p3 <- ""
    Aid.Data.P68$ADM_p45 <- ""
    Aid.Data.P68$ADM_p68 <- "Country" 
    
    Aid.Data.P68 <- subset(Aid.Data.P68@data,select = c(project_id,
                                                        project_location_id,
                                                        precision_code,
                                                        total_commitments,
                                                        total_disbursements,
                                                        ADM_p1,
                                                        ADM_p2,
                                                        ADM_p3,
                                                        ADM_p45,
                                                        ADM_p68))                                                    
  }
  
  Aid.Data.All <- rbind(Aid.Data.P1,Aid.Data.P2, Aid.Data.P3, Aid.Data.P45, Aid.Data.P68)
  Aid.Data.All <- Aid.Data.All[!(is.na(Aid.Data.All$precision_code)),]
  Aid.Data.All <- Aid.Data.All[order(Aid.Data.All$project_location_id),] 
  
  #write.csv(Aid.Data.All, file = "~/Desktop/NepalAidDataAll.csv")
  
  #### Making Prevision Code 1 Dataset to Iterate Through
  P1_shp$ADM_p1 <- as.character(P1_shp$NAME_3)
  P1_shp$ADM_p2 <- as.character(P1_shp$NAME_2)
  P1_shp$ADM_p3 <- as.character(P1_shp$NAME_1)
  P1_shp$ADM_p45 <- as.character(P1_shp$Adm_Region)
  P1_shp$ADM_p68 <- "Country"
  P1_shp$AreaROI <- P1_shp$Area
  ROI.Data <- subset(P1_shp@data,select = c(ADM_p1,
                                            ADM_p2,
                                            ADM_p3,
                                            ADM_p45,
                                            ADM_p68,
                                            AreaROI))  
  
  ROI.Data$ROI.ID <- row.names(ROI.Data)
  ROI.Data$ROI.ID <- as.integer(ROI.Data$ROI.ID)
  
  Aid.Data.All$project_location_id <- as.character(Aid.Data.All$project_location_id)
  
  roi_num=1
  roiValues <- function(i){
    
    # Make P1 Matrix    
    make_roi_params.col <- subset(P1_shp@data, select=c("NAME_3","AreaROI"))
    names(make_roi_params.col) <- c("ADM_p1","AreaROI")
    make_roi_params.col$Var <- 0
    
    # Which precision code is project
    project.pc <- Aid.Data.All[Aid.Data.All$project_location_id == i,]$precision_code 
    
    if (project.pc == 1){
      p1.name <- as.character(Aid.Data.All[Aid.Data.All$project_location_id == i,]$ADM_p1) 
      make_roi_params.col[make_roi_params.col$ADM_p1 == p1.name,]$Var <- 1
    }
    
    if (project.pc == 2){
      p2.name <- as.character(Aid.Data.All[Aid.Data.All$project_location_id == i,]$ADM_p2)
      p2.dataset.P1names <- P1_shp@data[P1_shp@data$NAME_2 == p2.name,]$ADM_p1
      make_roi_params.col[make_roi_params.col$ADM_p1 %in% p2.dataset.P1names,]$Var <- 
        make_roi_params.col[make_roi_params.col$ADM_p1 %in% p2.dataset.P1names,]$AreaROI / 
        sum(make_roi_params.col[make_roi_params.col$ADM_p1 %in% p2.dataset.P1names,]$AreaROI)
    }
    
    if (project.pc == 3){
      p3.name <- as.character(Aid.Data.All[Aid.Data.All$project_location_id == i,]$ADM_p3)
      p3.dataset.P1names <- P1_shp@data[P1_shp@data$NAME_1 == p3.name,]$ADM_p1
      make_roi_params.col[make_roi_params.col$ADM_p1 %in% p3.dataset.P1names,]$Var <- 
        make_roi_params.col[make_roi_params.col$ADM_p1 %in% p3.dataset.P1names,]$AreaROI / 
        sum(make_roi_params.col[make_roi_params.col$ADM_p1 %in% p3.dataset.P1names,]$AreaROI)
    }
    
    if (project.pc == 4 | project.pc == 5){
      p4.name <- as.character(Aid.Data.All[Aid.Data.All$project_location_id == i,]$ADM_p4)
      p4.dataset.P1names <- P1_shp@data[P1_shp@data$Adm_Region == p4.name,]$ADM_p1
      make_roi_params.col[make_roi_params.col$ADM_p1 %in% p4.dataset.P1names,]$Var <- 
        make_roi_params.col[make_roi_params.col$ADM_p1 %in% p4.dataset.P1names,]$AreaROI / 
        sum(make_roi_params.col[make_roi_params.col$ADM_p1 %in% p4.dataset.P1names,]$AreaROI)
    }
    
    if (project.pc == 6 | project.pc == 8){
      make_roi_params.col$Var <- make_roi_params.col$AreaROI / sum(make_roi_params.col$AreaROI)
    }
    
    row.names(make_roi_params.col) <- make_roi_params.col$ADM_p1
    make_roi_params.col <- subset(make_roi_params.col,select="Var")
    names(make_roi_params.col) <- i
    
    return(make_roi_params.col)
  }
  
  require(parallel)
  make_roi_params.col.all = mclapply(Aid.Data.All$project_location_id, roiValues)
  
  make_roi_params <- subset(P1_shp@data, select="NAME_3")
  row.names(make_roi_params) <- make_roi_params$NAME_3
  
  for(i in 1:length(make_roi_params.col.all)){
    make_roi_params <- cbind(make_roi_params, as.data.frame(make_roi_params.col.all[i]))
  }
  
  param_set <- make_roi_params[ , !(names(make_roi_params) == "NAME_3")]
  
  # Sort param set
  param_set <- param_set[order(row.names(param_set)),] 
  
  ###### Sum of area overlaps #####
  param_set.matrix <- param_set
  param_set.matrix[(param_set.matrix > 0)] <- 1
  numOverlap <- rowSums(param_set.matrix)
  sumAreaOverlap <- numOverlap*P1_shp$Area
  
  return(list(num_rois=num_rois,
              num_locations=num_locations,
              precision_set=precision_set,
              num_loc_set=num_loc_set,
              locus_set=locus_set,
              param_set=param_set,
              dollar_set=dollar_set,
              sumAreaOverlap=sumAreaOverlap,
              area=P1_shp$Area))
}

dollar_expected_value <- function(param_set, num_loc_set, dollar_set, num_rois){
  param_loc <- rbind(param_set, num_loc_set)
  expected_vals_by_project <- apply(param_loc, 2, function(x) x[1:num_rois] * x[(num_rois+1)])
  expected_vals_dollar_bind <- rbind(expected_vals_by_project, dollar_set)
  expected_vals_dollar <- apply(expected_vals_dollar_bind, 2, function(x) x[1:num_rois] * x[(num_rois+1)])
  expected_vals_dollar_summed <- apply(expected_vals_dollar, 1, sum)
  return(expected_vals_dollar_summed)
}


dollar_covariances <-function(param_set,dollar_set){
  dollar_variances.vector <- apply(param_set, 1, function(x) sum(x*(1-x)*dollar_set^2))
  dollar_covariances.matrix <- matrix(rep(0,nrow(param_set)^2),nrow=nrow(param_set))
  for(project in 1:ncol(param_set)){
    covariances <- matrix(as.numeric(param_set[,project]),ncol=1) %*% matrix(as.numeric(param_set[,project]),nrow=1)
    dollar_covariances <- covariances * -(dollar_set[project]^2)
    dollar_covariances.matrix <- dollar_covariances.matrix + dollar_covariances
  }  
  
  for(diag_spot in 1:nrow(param_set)){
    dollar_covariances.matrix[diag_spot,diag_spot] <- dollar_variances.vector[diag_spot]
  }
  return (dollar_covariances.matrix)
} 
#```

#The Gibbs sampler incorporates a multivariate normal model of the spread of aid, drawing a new data set from the distribution of believed aid at each run, and then drawing the parameters from their conditional distribution.

#```{r}
sim_x <- function(expected_vals,error_matrix){
  sim_vals <- rmvn(1,expected_vals,error_matrix,isChol=TRUE)
  return(matrix(sim_vals,ncol=1))
}

sim_beta <- function(x, y, sigma_epsilon, sigma_beta){
  delta <- sigma_epsilon^2/sigma_beta^2
  
  eye <- diag(rep(1,dim(x)[2]))
  
  mu <- solve(t(x) %*% x + delta*eye) %*% t(x) %*% y
  
  sigma <- sigma_epsilon^2 * solve(t(x) %*% x + delta * eye)
  
  return(rmvnorm(n=1,mean=mu,sigma=sigma))
}

sim_model_error <- function(x,y,lil_delta_1,lil_delta_2,beta){
  n <- length(y)
  a <- lil_delta_1 + n/2
  b <- lil_delta_2 + 1/2 * sum((y - x %*% beta)^2)
  return(rinvgamma(1,a,b))
}

Gibbs_Sampler <- function(known_y, expected_x, meas_errors){
  
  # hyper parameters
  sigma_beta <- 1000
  lil_delta_1 <- 1
  lil_delta_2 <- 1
  
  # Generate starting vals
  
  fit <- lm(known_y ~ expected_x)
  beta_0 <- coefficients(fit)['(Intercept)']
  beta_x <- coefficients(fit)[
    'expected_x']
  x <- cbind(1,expected_x)
  sigma_epsilon <- runif(1,.05,100)
  
  #Gibbs sampler
  
  n.samples <- 15000
  samples <- matrix(0,n.samples,3)
  colnames(samples) <- c('beta_0','beta_x','sigma_epsilon')
  
  for(i in 1:n.samples){
    # draw data
    x <- cbind(1,sim_x(expected_vals=expected_x,error_matrix=meas_errors))
    
    # draw beta coefficients
    beta <- matrix(sim_beta(x=x,y=known_y,sigma_epsilon=sigma_epsilon,sigma_beta=sigma_beta),ncol=1)
    
    # draw model error
    sigma_epsilon <- sim_model_error(x=x,y=known_y,lil_delta_1
                                     =lil_delta_1,lil_delta_2=lil_delta_2,beta=beta)
    
    # update
    samples[i,] <- c(t(beta),sigma_epsilon)
  }
  
  return(samples)
}
#```

#To test the Gibbs sampler, we first create a set of known project locations. 

#```{r}
##### Importing Shapefiles for Each Precision Code #####
UGA_p68 <- readShapePoly("UGA_shapes/UGA_adm_shp/UGA_adm0")
UGA_p45 <- readShapePoly("UGA_shapes/UG_regions/Uganda_regions_2014.shp")
UGA_p3 <- readShapePoly("UGA_shapes/UGA_adm_shp/UGA_adm1")
UGA_p2 <- readShapePoly("UGA_shapes/UGA_adm_shp/UGA_adm2")
UGA_p1 <- readShapePoly("UGA_shapes/UGA_adm_shp/UGA_adm3")



# Random Points Across Uganda (not stratified into sub-counties)
randAcrossUganda = "Yes"
if(randAcrossUganda == "Yes"){
  UGA_p68.Area <- gArea(UGA_p68, byid=TRUE) 
  UGA_p1$Area <- gArea(UGA_p1, byid=TRUE) 
  #UGA_p1$percentArea <- log(UGA_p1$Area * 10^40)
  
  subcountiesWithAid <- rmultinom(n=1,size=1000,prob=UGA_p1$Area)
  
  UGA_p1.centroids <- getSpPPolygonsLabptSlots(UGA_p1)
  UGA_p1.centroids <- as.data.frame(UGA_p1.centroids)
  UGA_p1.centroids <- cbind(UGA_p1.centroids, UGA_p1@data, subcountiesWithAid)
  
  # Confirming that large subcounties are more likely to get aid
  table(UGA_p1.centroids$subcountiesWithAid)
  #UGA_p1.centroids$subcountiesWithAid[UGA_p1.centroids$subcountiesWithAid >= 10] <- 10
  summary(lm(subcountiesWithAid ~ UGA_p1$Area))
  
  UGA_p1.centroids <- subset(UGA_p1.centroids, select = c("V1","V2","subcountiesWithAid"))
  randPoints.1 <- UGA_p1.centroids[UGA_p1.centroids$subcountiesWithAid == 1,]
  randPoints.2 <- UGA_p1.centroids[UGA_p1.centroids$subcountiesWithAid == 2,]
  randPoints.3 <- UGA_p1.centroids[UGA_p1.centroids$subcountiesWithAid == 3,]
  randPoints.4 <- UGA_p1.centroids[UGA_p1.centroids$subcountiesWithAid == 4,]
  randPoints.5 <- UGA_p1.centroids[UGA_p1.centroids$subcountiesWithAid == 5,]
  randPoints.6 <- UGA_p1.centroids[UGA_p1.centroids$subcountiesWithAid == 6,]
  randPoints.7 <- UGA_p1.centroids[UGA_p1.centroids$subcountiesWithAid == 7,]
  randPoints.8 <- UGA_p1.centroids[UGA_p1.centroids$subcountiesWithAid == 8,]
  randPoints.9 <- UGA_p1.centroids[UGA_p1.centroids$subcountiesWithAid == 9,]
  randPoints.10 <- UGA_p1.centroids[UGA_p1.centroids$subcountiesWithAid == 10,]
  
  randPoints <- rbind(randPoints.1,
                      randPoints.2,
                      randPoints.2,
                      randPoints.3,
                      randPoints.3,
                      randPoints.3,
                      randPoints.4,
                      randPoints.4,
                      randPoints.4,
                      randPoints.4,
                      randPoints.5,
                      randPoints.5,
                      randPoints.5,
                      randPoints.5,
                      randPoints.5,
                      randPoints.6,
                      randPoints.6,
                      randPoints.6,
                      randPoints.6,
                      randPoints.6,
                      randPoints.6,
                      randPoints.7,
                      randPoints.7,
                      randPoints.7,
                      randPoints.7,
                      randPoints.7,
                      randPoints.7,
                      randPoints.7,
                      randPoints.2,
                      randPoints.8,
                      randPoints.8,
                      randPoints.8,
                      randPoints.8,
                      randPoints.8,
                      randPoints.8,
                      randPoints.8,
                      randPoints.8,
                      randPoints.9,
                      randPoints.9,
                      randPoints.9,
                      randPoints.9,
                      randPoints.9,
                      randPoints.9,
                      randPoints.9,
                      randPoints.9,
                      randPoints.9,
                      randPoints.10,
                      randPoints.10,
                      randPoints.10,
                      randPoints.10,
                      randPoints.10,
                      randPoints.10,
                      randPoints.10,
                      randPoints.10,
                      randPoints.10,
                      randPoints.10)
  
  randPoints <- subset(randPoints, select=c("V1","V2"))
  names(randPoints) <- c("x","y")
}

randPoints$latitude <- randPoints$y
randPoints$longitude <- randPoints$x
Uganda.AidData <- randPoints
# Uganda.AidData$total_commitments <- 1*runif(nrow(Uganda.AidData), min=0, max=1) 
# Uganda.AidData$total_disbursements <- 1*runif(nrow(Uganda.AidData), min=0, max=1) 
# Uganda.AidData$even_split_commitments <- 1*runif(nrow(Uganda.AidData), min=0, max=1) 

Uganda.AidData$total_commitments <- rep(1,nrow(Uganda.AidData))
Uganda.AidData$total_disbursements <- rep(1,nrow(Uganda.AidData)) 
Uganda.AidData$even_split_commitments <- rep(1,nrow(Uganda.AidData))
Uganda.AidData$precision_code <- 1

# randomly give 100 projects
Uganda.AidData$project_id <- round(runif(nrow(Uganda.AidData),1,100))
Uganda.AidData$project_location_id <- 1:nrow(Uganda.AidData)

# removing points over lakes
#Uganda.AidData.GIS <- Uganda.AidData
#coordinates(Uganda.AidData.GIS) <- ~longitude+latitude
#Uganda.AidData$p1adm <- over(Uganda.AidData.GIS, UGA_p1)$NAME_1
#Uganda.AidData$p1adm_area <- over(Uganda.AidData.GIS, UGA_p1)$Area
#Uganda.AidData <- Uganda.AidData[Uganda.AidData$p1adm_area < 0.05,]
#Uganda.AidData <- Uganda.AidData[Uganda.AidData$p1adm != "Lake Albert",]
#Uganda.AidData <- Uganda.AidData[Uganda.AidData$p1adm != "Lake Victoria",]
#Uganda.AidData <- subset(Uganda.AidData, select = -c(p1adm))

# Making Sure UGA_p1 dataset has UGA_p45 names in there
UGA_p1.centroids <- getSpPPolygonsLabptSlots(UGA_p1)
UGA_p1.centroids <- cbind(UGA_p1.centroids,UGA_p1@data)
names(UGA_p1.centroids)[names(UGA_p1.centroids)=="1"] <- "latitude"
names(UGA_p1.centroids)[names(UGA_p1.centroids)=="2"] <- "longitude"
coordinates(UGA_p1.centroids) <- ~latitude+longitude
UGA_p1$Adm_Region <- over(UGA_p1.centroids,UGA_p45)$Adm_Region

# Names of Admins for P1 are not unique, so we need to make them unique.
# Merge sub-county and county name
mean(as.data.frame(table(UGA_p1$NAME_3))$Freq)
UGA_p1$NAME_3 <- paste(UGA_p1$NAME_3,UGA_p1$NAME_2,sep="")
mean(as.data.frame(table(UGA_p1$NAME_3))$Freq)

# Doing same thing with covariate data, and sort by that name. 
#Uganda.data$NAME_3 <- paste(Uganda.data$NAME_3,Uganda.data$NAME_2,sep="")
#Uganda.data <- Uganda.data[order(Uganda.data$NAME_3),] 
UGA_p1 <- UGA_p1[order(UGA_p1$NAME_3),]

# Defining True Relation
precisionCodeMatrix = 1
i = 1
#Uganda.AidData$precision_code <- sample(size=nrow(Uganda.AidData), x=p.codes.ascend, prob=p.probs.all[[precisionCodeMatrix]][,i], replace=TRUE)
Uganda.AidData$precision_code <- rep(1,nrow(Uganda.AidData))
paramVals <- parameters(Uganda.AidData, P1_shp=UGA_p1, P2_shp=UGA_p2,P3_shp=UGA_p3,P45_shp=UGA_p45,P68_shp=UGA_p68)
aid.expected.true <- dollar_expected_value(paramVals$param_set,paramVals$num_loc_set, paramVals$dollar_set, paramVals$num_rois)
beta_0 <- 0
beta_1 <- 1
poverty <- beta_0 + beta_1*aid.expected.true

plot(aid.expected.true,poverty)
#```

#Then, we create a set of scenarios involving various levels of spatial uncertainty

#```{r}
##### Increasing Precision Code #####
p.codes.ascend <- c(1,2,3,4,5,6,8)
p.probs.12 <- matrix(data = 0, nrow =7, ncol =6)
p.probs.13 <- matrix(data = 0, nrow =7, ncol =6)
p.probs.14 <- matrix(data = 0, nrow =7, ncol =6)
p.probs.16 <- matrix(data = 0, nrow =7, ncol =6)

for(i in 1:6){  
  prob_prec_1 <- 1 - (i-1)*.2
  
  p.probs.12[1,i] <- prob_prec_1
  p.probs.12[2,i] <- 1 - p.probs.12[1,i]
  
  p.probs.13[1,i] <- prob_prec_1
  p.probs.13[3,i] <- 1 - p.probs.13[1,i]
  
  p.probs.14[1,i] <- prob_prec_1
  p.probs.14[4,i] <- 1 - p.probs.14[1,i]
  
  p.probs.16[1,i] <- prob_prec_1
  p.probs.16[6,i] <- 1 - p.probs.13[1,i]
}

p.probs.all <- rep(NA,4)
p.probs.all <- as.list(p.probs.all)
p.probs.all[[1]] <- p.probs.12
p.probs.all[[2]] <- p.probs.13
p.probs.all[[3]] <- p.probs.14
p.probs.all[[4]] <- p.probs.16

# Making matrix to iterate through
p.i.matrix <- matrix(NA,ncol=3,nrow=45)
p.i.matrix <- as.data.frame(p.i.matrix)
names(p.i.matrix) <- c("precision code matrix","percent not pc1","model error")
p.i.matrix[1:15,3] <- 0.01
p.i.matrix[16:30,3] <- 0.1
p.i.matrix[31:45,3] <- 1
p.i.matrix[1:5,1] <- 1
p.i.matrix[6:10,1] <- 2
p.i.matrix[11:15,1] <- 4
p.i.matrix[16:20,1] <- 1
p.i.matrix[21:25,1] <- 2
p.i.matrix[26:30,1] <- 4
p.i.matrix[31:35,1] <- 1
p.i.matrix[36:40,1] <- 2
p.i.matrix[41:45,1] <- 4
p.i.matrix[1:5,2] <- 2:6
p.i.matrix[6:10,2] <- 2:6
p.i.matrix[11:15,2] <- 2:6
p.i.matrix[16:20,2] <- 2:6
p.i.matrix[21:25,2] <- 2:6
p.i.matrix[26:30,2] <- 2:6
p.i.matrix[31:35,2] <- 2:6
p.i.matrix[36:40,2] <- 2:6
p.i.matrix[41:45,2] <- 2:6

makeData <- function(j){
  print(paste('Generating Scenario',j))
  precisionCodeMatrix = p.i.matrix[j,1]
  i = p.i.matrix[j,2]
  model.error <- p.i.matrix[j,3]
  Uganda.AidData$precision_code <- sample(size=nrow(Uganda.AidData), x=p.codes.ascend, prob=p.probs.all[[precisionCodeMatrix]][,i], replace=TRUE)
  
  paramVals <- parameters(Uganda.AidData,P1_shp=UGA_p1,P2_shp=UGA_p2,P3_shp=UGA_p3,P45_shp=UGA_p45,P68_shp=UGA_p68)
  aid.expected <- dollar_expected_value(paramVals$param_set,paramVals$num_loc_set, paramVals$dollar_set, paramVals$num_rois)
  aid.dollar_covariances <- dollar_covariances(paramVals$param_set,paramVals$dollar_set)
  poverty <- beta_0 + aid.expected.true*beta_1 + model.error*rnorm(length(aid.expected))
  
  return(cbind(poverty,aid.expected,aid.dollar_covariances))
}

scenario_data <- mclapply(1:45,makeData)

#```

#We analyze the data using Gibbs sampling and compare with a simple naive regression on the expected values.

#```{r}
results_dir <- 'Gibbs Results'
makeGibMatrices.have.data <- function(j,poverty,aid.expected,aid.dollar_covariances){
  
  ### Scenario label ###
  print(paste('Assessing Scenario',j))
  precisionCodeMatrix = p.i.matrix[j,1]
  i = p.i.matrix[j,2]
  model.error <- p.i.matrix[j,3]
  
  ### Regression ###
  print('Performing Regression Model')
  gibbs.results <- Gibbs_Sampler(poverty,aid.expected,aid.dollar_covariances)
  naive.coef <- lm(poverty ~ aid.expected)
  
  ### Print Graphics ###
  print('Printing Results')
  
  percent_alt <- sum(p.probs.all[[precisionCodeMatrix]][,i][2:7])*100
  alt_code <- p.codes.ascend[2:7][p.probs.all[[precisionCodeMatrix]][,i][2:7] != 0]
  
  scenario_code <- paste(percent_alt,'percent','code',alt_code,'with',model.error,'modelError',sep='_')
  
  dir.create(paste(results_dir,scenario_code,sep='/'))
  
  scenario_code_pretty <-paste(percent_alt,'Percent','Code',alt_code,'with',model.error,'Model Error',sep=' ')
  scatter_plot_file <- paste(results_dir,'/',scenario_code,'/','scatter_plot','_',scenario_code,'.png',sep='')
  intercept_trace_file <- paste(results_dir,'/',scenario_code,'/','intercept_trace','_',scenario_code,'.png',sep='')
  coefficient_trace_file <- paste(results_dir,'/',scenario_code,'/','coefficient_trace','_',scenario_code,'.png',sep='')
  error_trace_file <- paste(results_dir,'/',scenario_code,'/','error_trace','_',scenario_code,'.png',sep='')
  intercept_hist_file <- paste(results_dir,'/',scenario_code,'/','intercept_hist','_',scenario_code,'.png',sep='')
  coefficient_hist_file <- paste(results_dir,'/',scenario_code,'/','coefficient_hist','_',scenario_code,'.png',sep='')
  error_hist_file <- paste(results_dir,'/',scenario_code,'/','error_hist','_',scenario_code,'.png',sep='')
  
  # Scatter plot
  aid.dollar_var <- sapply(1:dim(aid.dollar_covariances)[1], function(i) aid.dollar_covariances[i,i])
  aid.dollar_sd <- sqrt(aid.dollar_var)
  png(filename = scatter_plot_file)
  plot(aid.expected,poverty,main=scenario_code_pretty,xlab='Expected Aid',ylab='Poverty')
  segments(aid.expected-aid.dollar_sd,poverty,aid.expected+aid.dollar_sd,poverty)
  abline(naive.coef)
  dev.off()
  
  # Trace plots
  png(filename = intercept_trace_file)
  plot(burnin:n.samples,gibbs.results[burnin:n.samples,1],main=paste('Intercept, ',scenario_code_pretty,sep=''),xlab='',ylab='')
  abline(h=mean(gibbs.results[burnin:n.samples,1]),col=19)
  dev.off()
  
  png(filename = coefficient_trace_file)
  plot(burnin:n.samples,gibbs.results[burnin:n.samples,2],main=paste('Coefficient, ',scenario_code_pretty,sep=''),xlab='',ylab='')
  abline(h=mean(gibbs.results[burnin:n.samples,2]),col=19)
  dev.off()
  
  png(filename = error_trace_file)
  plot(burnin:n.samples,gibbs.results[burnin:n.samples,3],main=paste('Error, ',scenario_code_pretty,sep=''),xlab='',ylab='')
  abline(h=mean(gibbs.results[burnin:n.samples,3]),col=19)
  dev.off()
  
  png(filename = intercept_hist_file)
  hist(gibbs.results[burnin:n.samples,1],main=paste('Intercept, ',main=scenario_code_pretty,sep=''),xlab='Inercept Value',ylab='Frequency')
  abline(v=mean(gibbs.results[burnin:n.samples,1]))
  dev.off()
  
  png(filename = coefficient_hist_file)
  hist(gibbs.results[burnin:n.samples,2],main=paste('Coefficient, ',scenario_code_pretty,sep=''),xlab='Coefficient Value',ylab='Frequency')
  abline(v=mean(gibbs.results[burnin:n.samples,2]))
  dev.off()
  
  png(filename = error_hist_file)
  hist(gibbs.results[burnin:n.samples,3],main=paste('Model Error, ',scenario_code_pretty,sep=''),xlab='Error Value',ylab='Frequency')
  abline(v=mean(gibbs.results[burnin:n.samples,3]))
  dev.off()
  
  # csv of gibbs results
  results_file <- paste(results_dir,'/',scenario_code,'/',scenario_code,'.csv',sep='')
  results.df <- t(data.frame('intercept'=quantile(gibbs.results[burnin:n.samples,1],c(.025,.975)),'coef'=quantile(gibbs.results[burnin:n.samples,2],c(.025,.975),'error'=quantile(gibbs.results[burnin:n.samples,3],c(.025,.975))),'naive'=coefficients(naive.coef)))
  write.csv(results.df,file=results_file)
  
  return(list("gibbs.results" = gibbs.results,
              "naive.coef" = naive.coef))
}

burnin=5000

scenario_analyses <- mclapply(1:45, function(i) makeGibMatrices.have.data(i,poverty=scenario_data[[i]][,1],aid.expected=scenario_data[[i]][,2],aid.dollar_covariances = scenario_data[[i]][,3:dim(scenario_data[[i]])[2]]))
#```

#And we print the results

#```{r}

listScenarioCode <- function(k){
  precisionCodeMatrix = p.i.matrix[k,1]
  i = p.i.matrix[k,2]
  model.error <- p.i.matrix[k,3]
  percent_alt <- sum(p.probs.all[[precisionCodeMatrix]][,i][2:7])*100
  alt_code <- p.codes.ascend[2:7][p.probs.all[[precisionCodeMatrix]][,i][2:7] != 0]
  scenario_code_pretty <-paste(percent_alt,'Percent','Code',alt_code,'with',model.error,'Model Error',sep=' ')
  
  return(scenario_code_pretty)
}

scenarios <- lapply(1:45, listScenarioCode)
naive.intercepts <- sapply(1:45, function(i) coefficients(scenario_analyses[[i]]$naive.coef)['(Intercept)'])
naive.coefs <- sapply(1:45, function(i) coefficients(scenario_analyses[[i]]$naive.coef)['aid.expected'])
intercepts.25 <- sapply(1:45, function(i) quantile(scenario_analyses[[i]]$gibbs.results[,1],c(.025)))
intercepts.975 <- sapply(1:45, function(i) quantile(scenario_analyses[[i]]$gibbs.results[,1],c(.975)))
intercepts.mean <- sapply(1:45, function(i) mean(scenario_analyses[[i]]$gibbs.results[,1]))
coef.25 <- sapply(1:45, function(i) quantile(scenario_analyses[[i]]$gibbs.results[,2],c(.025)))
coef.975 <- sapply(1:45, function(i) quantile(scenario_analyses[[i]]$gibbs.results[,2],c(.975)))
coef.mean <- sapply(1:45, function(i) mean(scenario_analyses[[i]]$gibbs.results[,2]))

result.summary <- cbind(naive.intercepts,naive.coefs,intercepts.25,intercepts.975,intercepts.mean,coef.25,coef.975,coef.mean)
rownames(result.summary) <- scenarios


summary_result_dir <- 'Gibbs Results/Summary'
create.dir(summary_result_dir)

coef.filenames <- c('Coefficient Code 2 with .01 Model Error','Coefficient Code 3 with .01 Model Error',
                    'Coefficient Code 6 with .01 Model Error','Coefficient Code 2 with .1 Model Error',
                    'Coefficient Code 3 with .1 Model Error','Coefficient Code 6 with .1 Model Error',
                    'Coefficient Code 2 with 1 Model Error','Coefficient Code 3 with 1 Model Error',
                    'Coefficient Code 6 with 1 Model Error')

xlabels <- rep(c('Percent Code 2', 'Percent Code 3', 'Percent Code 6'),3)

for(l in 1:9){
  png(paste(summary_result_dir,'/',coef.filenames[l],'.png',sep=''))
  plot(c(20,40,60,80,100),coef.mean[l:(l+4)],main=coef.filenames[l],
       xlab=xlabels[l] ,ylab='Coefficient',ylim=c(0,2),pch=19,col='red')
  segments(c(20,40,60,80,100),coef.25[l:(l+4)],c(20,40,60,80,100),coef.975[l:(l+4)],col='red')
  abline(h=1)
  par(new=T)
  plot(c(20,40,60,80,100),naive.coefs[l:(l+4)],axes=F,xlab='',ylab='',ylim=c(0,2),pch=19,col='green')
  par(new=F)
  dev.off()
}

intercept.filenames <- c('Intercept Code 2 with .01 Model Error','Intercept Code 3 with .01 Model Error',
                         'Intercept Code 6 with .01 Model Error','Intercept Code 2 with .1 Model Error',
                         'Intercept Code 3 with .1 Model Error','Intercept Code 6 with .1 Model Error',
                         'Intercept Code 2 with 1 Model Error','Intercept Code 3 with 1 Model Error',
                         'Intercept Code 6 with 1 Model Error')

for(m in 1:9){
  png(paste(summary_result_dir,'/',intercept.filenames[m],'.png',sep=''))
  plot(c(20,40,60,80,100),coef.mean[m:(m+4)],main=intercept.filenames[m],
       xlab=xlabels[m],ylab='Coefficient',ylim=c(0,2),pch=19,col='red')
  segments(c(20,40,60,80,100),coef.25[m:(m+4)],c(20,40,60,80,100),coef.975[m:(m+4)],col='red')
  abline(h=1)
  par(new=T)
  plot(c(20,40,60,80,100),naive.coefs[m:(m+4)],axes=F,xlab='',ylab='',ylim=c(0,2),pch=19,col='green')
  par(new=F)
  dev.off()
}
#```