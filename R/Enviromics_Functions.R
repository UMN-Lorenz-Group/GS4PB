
#' Title of get_helper
#'
#' Description of what get_helper does.
#'
#' @param lon Description of lon.
#' @param lat Description of lat.
#' @param variables.names Description of variables.names.
#' @param start.day Description of start.day.
#' @param end.day Description of end.day.
#' @param env.id Description of env.id.
#' @param save Description of save.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of get_helper
#' result <- get_helper(...)
#' }

get_helper <- function(lon, lat, variables.names, start.day, 
        end.day, env.id, save) {
        CL = data.frame(nasapower::get_power(community = "ag", 
            lonlat = c(lon, lat), pars = variables.names, dates = c(start.day, 
                end.day), temporal_api = "daily"))
        cor_rain_name = which(names(CL) %in% "PRECTOTCORR")
        names(CL)[cor_rain_name] = "PRECTOT"
        CL$daysFromStart = 1:nrow(CL)
        CL$env <- env.id
        CL <- CL[, c(which(colnames(CL) == "env"), which(colnames(CL) != 
            "env"))]
        if (isTRUE(save)) {
            utils::write.csv(file = paste(env.id, ".csv", sep = ""), 
                row.names = F, x = CL)
        }
        return(CL)
    }

	# Function to initialize progress bar settings
	#' Title of progress
#'
#' Description of what progress does.
#'
#' @param min = 0 Description of min = 0.
#' @param max = 100 Description of max = 100.
#' @param leftd = "|" Description of leftd = "|".
#' @param rightd = "|" Description of rightd = "|".
#' @param char = "=" Description of char = "=".
#' @param style = 2 Description of style = 2.
#' @param width = getOption("width") Description of width = getOption("width").
#' @param time = Sys.time() Description of time = Sys.time().
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of progress
#' result <- progress(...)
#' }
progress <- function(min = 0, max = 100, leftd = "|", rightd = "|", 
						 char = "=", style = 2, width = getOption("width"), time = Sys.time()) {
	  if (max <= 0) {
		stop("max must be a positive number")
	  }
	  if (width <= 0) {
		stop("width must be a positive number")
	  }
	  return(list(min = min, max = max, leftd = leftd, rightd = rightd, 
				  char = char, style = style, width = width, time = time))
	}

	# Function to run and update the progress bar
	#' Title of run_progress
#'
#' Description of what run_progress does.
#'
#' @param pb Description of pb.
#' @param actual Description of actual.
#' @param text = "" Description of text = "".
#' @param digits = 0 Description of digits = 0.
#' @param sleep = 0 Description of sleep = 0.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of run_progress
#' result <- run_progress(...)
#' }

run_progress <- function(pb, actual, text = "", digits = 0, sleep = 0){
	  Sys.sleep(sleep)
	  elapsed <- sec_to_hms(as.numeric(difftime(Sys.time(), pb$time, units = "secs")))
	  
	  temp <- switch(pb$style,
					 list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd),
						  text = paste(text, paste(pb$leftd, "%s%s", pb$rightd, sep = ""))),
					 list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 6,
						  text = paste(text, paste(pb$leftd, "%s%s", pb$rightd, sep = ""), "%s%%")),
					 list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 9,
						  text = paste(text, paste(pb$leftd, "%s%s", pb$rightd, sep = ""), elapsed)),
					 list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 15,
						  text = paste(text, paste(pb$leftd, "%s%s", pb$rightd, sep = ""), "%s%%", elapsed))
	  )
	  
	  step <- round(actual / pb$max * (pb$width - temp$extra))
	  
	  # Ensure step and width calculation are valid
	  if (step < 0) step <- 0
	  if (pb$width - step - temp$extra < 0) {
		# Adjust pb$width to a minimum valid size
		pb$width <- step + temp$extra + 1
	  }
	  
	  # Print the progress bar
	  cat(sprintf(temp$text, strrep(pb$char, step), strrep(" ", pb$width - step - temp$extra),
				  round(actual / pb$max * 100, digits = digits)), "\r")
	  
	  if (actual >= pb$max) {
		cat("\n")
	  }
	}

# Helper function to convert seconds to HH:MM:SS format
#' Title of sec_to_hms
#'
#' Description of what sec_to_hms does.
#'
#' @param seconds Description of seconds.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of sec_to_hms
#' result <- sec_to_hms(...)
#' }
sec_to_hms <- function(seconds) {
	  h <- floor(seconds / 3600)
	  m <- floor((seconds %% 3600) / 60)
	  s <- seconds %% 60
	  sprintf("%02f:%02f:%02f", h, m, s)
	}

	
	
#' Title of split_chunk
#'
#' Description of what split_chunk does.
#'
#' @param vec Description of vec.
#' @param length Description of length.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of split_chunk
#' result <- split_chunk(...)
#' }
split_chunk <- function(vec, length) {
        split(vec, ceiling(seq_along(vec)/length))
}


#' Download and Process Weather Data from NASA POWER API
#'
#' This function retrieves weather data from the NASA POWER API based on input geographical coordinates 
#' and date ranges. It supports parallel processing to handle large datasets efficiently and includes 
#' functionality to integrate elevation data from geospatial sources.
#'
#' @param env.id A vector of environment IDs (e.g., "env1"). If `NULL`, they will be generated automatically.
#' @param lat A vector of latitudes for the locations.
#' @param lon A vector of longitudes for the locations.
#' @param start.day A vector of start dates for data retrieval in "YYYY-MM-DD" format. If `NULL`, defaults to 1000 days prior to the current date.
#' @param end.day A vector of end dates for data retrieval in "YYYY-MM-DD" format. If `NULL`, defaults to 30 days after `start.day`.
#' @param variables.names A character vector of variable names to retrieve (e.g., "T2M", "PRECTOTCORR"). Defaults to common meteorological variables.
#' @param dir.path The directory path for saving output files. Defaults to the current working directory.
#' @param save Logical. If `TRUE`, saves the retrieved data to files in the specified directory.
#' @param temporal.scale The temporal resolution of the data. Defaults to `"daily"`.
#' @param country A vector of country codes (e.g., "USA", "BRA") for retrieving elevation data. If `NULL`, elevation data is not included.
#' @param parallel Logical. If `TRUE`, enables parallel processing for faster data retrieval.
#' @param workers Number of workers/cores to use for parallel processing. Defaults to 90\% of available cores.
#' @param chunk_size Number of locations to process in each parallel chunk. Defaults to 29.
#' @param sleep Time in seconds to wait between API requests. Defaults to 60.
#'
#' @return A data frame containing the retrieved weather data, including elevation if `country` is specified.
#'
#' @details 
#' This function interacts with the NASA POWER API to retrieve weather data, such as temperature, precipitation, and wind speed, for specific locations and time periods. It optionally integrates elevation data from geospatial sources. Parallel processing is supported to improve performance when handling large datasets.
#'
#' @examples
#' \dontrun{
#' # Retrieve daily weather data for two locations
#' weather_data <- get_weather_terra(
#'   env.id = c("env1", "env2"),
#'   lat = c(40.7128, 34.0522),
#'   lon = c(-74.0060, -118.2437),
#'   start.day = c("2023-01-01", "2023-01-01"),
#'   end.day = c("2023-01-31", "2023-01-31")
#' )
#' }
#'
#' @export
 
get_weather_terra <- function(env.id = NULL, lat = NULL, lon = NULL, start.day = NULL, 
    end.day = NULL, variables.names = NULL, dir.path = NULL, 
    save = FALSE, temporal.scale = "daily", country = NULL, parallel = TRUE, 
    workers = NULL, chunk_size = 29, sleep = 60) 
{
    if (!requireNamespace("doParallel", quietly = TRUE)) {
        utils::install.packages("doParallel")
    }
    if (!requireNamespace("parallel", quietly = TRUE)) {
        utils::install.packages("doParallel")
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
        utils::install.packages("foreach")
    }
   
    cat("------------------------------------------------ \n")
    cat("ATTENTION: This function requires internet access \n")
    cat("------------------------------------------------  \n")
    cat("Connecting to the NASA POWER API Client, Sparks et al 2018 \n")
    cat("https://docs.ropensci.org/nasapower \n")
    cat("------------------------------------------------  \n")
    if (is.null(env.id)) {
        env.id <- paste0("env", seq_along(lat))
    }
    if (!(is.character(env.id) || is.factor(env.id))) {
        stop("env.id should be a vector of characters (e.g. 'env1') or factors")
    }
    if (!requireNamespace("nasapower", quietly = TRUE)) {
        utils::install.packages("nasapower")
    }
    if (!requireNamespace("plyr", quietly = TRUE)) {
        utils::install.packages("plyr")
    }
    if (is.null(dir.path)) {
        dir.path = getwd()
    }
    if (is.null(start.day)) {
        start.day <- Sys.Date() - 1000
        cat(paste0("start.day is NULL", "\n"))
        cat(paste0("matched as ", start.day, "\n"))
    }
    if (is.null(end.day)) {
        end.day <- start.day + 30
        cat(paste0("end.day is NULL", "\n"))
        cat(paste0("matched as ", end.day, "\n"))
    }
    if (is.null(variables.names)) {
        variables.names = c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOTCORR", 
            "WS2M", "RH2M", "T2MDEW", "ALLSKY_SFC_LW_DWN", "ALLSKY_SFC_SW_DWN")
    }
    variables.names[grepl(variables.names, pattern = "PRECTOT")] = "PRECTOTCORR"
    env.id = as.factor(env.id)
    if (parallel == FALSE) {
        results <- list()
        pb <- progress(max = length(env.id), style = 4)
        init_time <- Sys.time()
        iter <- 0
        for (i in 1:length(env.id)) {
            iter <- iter + 1
            query_issue <- as.numeric(difftime(Sys.time(), init_time, 
                units = "secs")) > 60
            if (iter >= 30 & query_issue) {
                message("Waiting ", sleep, "s for a new query to the API.")
                Sys.sleep(sleep)
                iter <- 0
                init_time <- Sys.time()
            }
            results[[i]] <- get_helper(lon = lon[i], lat = lat[i], 
                variables.names = variables.names, start.day = start.day[i], 
                end.day = end.day[i], env.id = env.id[i], save = save)
            msg <- paste0("Env ", env.id[i], " (", i, "/", length(env.id), 
                ") ", "downloaded")
            run_progress(pb, actual = i, text = msg)
        }
		cat("\nNASA POWER: Done!")
    }
    if (parallel == TRUE) {
        env.id_par = split_chunk(env.id, length = chunk_size)
        lat_par = split_chunk(lat, length = chunk_size)
        lon_par = split_chunk(lon, length = chunk_size)
        start.day_par = split_chunk(start.day, length = chunk_size)
        end.day_par = split_chunk(end.day, length = chunk_size)
        nworkers <- ifelse(is.null(workers), trunc(parallel::detectCores() * 
            0.9), workers)
        clust <- parallel::makeCluster(nworkers)
        on.exit(parallel::stopCluster(clust))
        results <- list()
        pb <- progress(max = length(env.id_par), style = 4)
        for (i in 1:length(env.id_par)) {
            env.id_par_tmp <- env.id_par[[i]]
            lat_par_tmp <- lat_par[[i]]
            lon_par_tmp <- lon_par[[i]]
            start.day_par_tmp <- start.day_par[[i]]
            end.day_par_tmp <- end.day_par[[i]]
            parallel::clusterExport(clust, varlist = c("get_helper", 
                "lat_par_tmp", "lon_par_tmp", "variables.names", 
                "start.day_par_tmp", "end.day_par_tmp", "env.id_par_tmp"), 
                envir = environment())
            length_chunk <- length(env.id_par[[i]])
            temp <- parallel::parLapply(clust, 1:length_chunk, 
                function(j) {
                  get_helper(lon = lon_par_tmp[j], lat = lat_par_tmp[j], 
                    variables.names = variables.names, start.day = start.day_par_tmp[j], 
                    end.day = end.day_par_tmp[j], env.id = env.id_par_tmp[j], 
                    save = save)
                })
            results[[i]] <- plyr::ldply(temp)
            if (i < length(env.id_par)) {
                message("Waiting ", sleep, "s for a new query to the API.")
                msg <- paste0("Chunk ", i, "/", length(env.id_par), 
                  " (", length_chunk, " points) downloaded")
                run_progress(pb, actual = i, text = msg)
                Sys.sleep(sleep)
            }
        }
        cat("\nNASA POWER: Done!")
    }
   
    df <- plyr::ldply(results)
    if (is.null(country)) {
        cat("\n")
        cat("country = NULL ; Not possible to download ALT data. Please provide a country ID (e.g., USA1, BRA, FRA)")
        variables.names[grepl(variables.names, pattern = "PRECTOTCORR")] = "PRECTOT"
        ids = which(names(df) %in% variables.names)
        df[, ids][df[, ids] == -999] = NA
        return(df)
    }
    if (!is.null(country)) {
        if (!requireNamespace("geodata")) 
            utils::install.packages("geodata")
        suppressMessages("geodata")
        suppressWarnings("geodata")
        cat("\n")
        cat("Connecting to https://geodata.ucdavis.edu/ using Hijmans 2021")
        cat("------------------------------------------------  \n")
        unique_country <- unique(country)
        geodata_alt <- vector(mode = "list", length = length(unique_country))
        names(geodata_alt) = unique_country
        for (n in 1:length(unique_country)) {
            if (isTRUE(grepl(x = unique_country[n], pattern = "USA"))) {
                id_region = as.numeric(gsub(x = unique_country[n], 
                  pattern = "USA", replacement = ""))
				geodata_alt[[n]] <- suppressMessages(geodata::elevation_30s("USA",getwd(),mask=TRUE)[[id_region]])
            }
            else {
                geodata_alt[[n]] <- suppressMessages(geodata::elevation_30s(country = unique_country[n], getwd(),mask = TRUE))
            }
        }
        df_no_alt <- df
        df <- c()
        for (i in 1:length(env.id)) {
		
            country_raster <- geodata_alt[[which(names(geodata_alt) %in% 
                country[i])]]
            id <- which(df_no_alt$env %in% env.id[i])
            df <- rbind(df, extract_GIS_terra(covraster = country_raster, 
                env.data = df_no_alt[id, ], name.out = "ALT"))
			
        }
		colnames(df) <- gsub("ALT.*","ALT",colnames(df))
        variables.names <- c(variables.names, "ALT")
        cat("\nSRTM: Done!")
        variables.names[grepl(variables.names, pattern = "PRECTOTCORR")] = "PRECTOT"
        ids = which(names(df) %in% variables.names)
        df[, ids][df[, ids] == -999] = NA
        return(df)
    }
}


####

#' Title of extract_GIS_terra
#'
#' Description of what extract_GIS_terra does.
#'
#' @param covraster = NULL Description of covraster = NULL.
#' @param Latitude = NULL Description of Latitude = NULL.
#' @param Longitude = NULL Description of Longitude = NULL.
#' @param env.data = NULL Description of env.data = NULL.
#' @param env.id = NULL Description of env.id = NULL.
#' @param name.out = NULL Description of name.out = NULL.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of extract_GIS_terra
#' result <- extract_GIS_terra(...)
#' }

extract_GIS_terra <- function(covraster = NULL, Latitude = NULL, Longitude = NULL, 
                              env.data = NULL, env.id = NULL, name.out = NULL) 
{
    # Setting default names if not provided
    if (is.null(name.out)) 
        name.out = "ALT"
    if (is.null(Latitude)) 
        Latitude <- "LAT"
    if (is.null(Longitude)) 
        Longitude <- "LON"
    if (is.null(env.id)) 
        env.id <- "env"
    
    # Creating a vector data frame from the environmental data for latitude and longitude
    loc <- data.frame(x = env.data[, Longitude], y = env.data[, Latitude], id = env.data[, env.id])
    points <- terra::vect(loc, geom=c("x", "y"), crs="+proj=longlat +datum=WGS84")
    
    # Extracting data from each raster in the list
    extracted_data <- lapply(covraster, function(rast) terra::extract(rast,points, fun=mean, na.rm=TRUE, df=TRUE))
    extracted_data_df <- do.call(rbind, extracted_data)
    names(extracted_data_df)[-1] <- sapply(names(covraster), function(nm) paste(name.out, nm, sep="_"))

    # Ensuring the ID column from points is included for merging
    extracted_data_df$id <- points$id
    
    # Merging extracted data with the environmental data using the unique identifier
    final_data <- merge(env.data, extracted_data_df, by.x = env.id, by.y = "id")
	
	idCol <- which(colnames(final_data) %in% "ID")
	if(length(idCol)>0) { final_data_V2 <- final_data[,-idCol]}else{final_data_V2 <- final_data} 
    
    return(final_data_V2)
}



#' Title of getEnvData
#'
#' Description of what getEnvData does.
#'
#' @param Coord_Data_Tab Description of Coord_Data_Tab.
#' @param startDate Description of startDate.
#' @param endDate Description of endDate.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getEnvData
#' result <- getEnvData(...)
#' }

getEnvData <- function(Coord_Data_Tab,startDate,endDate){
    
	Dat_Loc_Coords_Filter <- Coord_Data_Tab

	Dat_Env <- as.character(Dat_Loc_Coords_Filter[,"Location"])
    Dat_Count <- as.character(Dat_Loc_Coords_Filter[,"Country"])
	lat <- as.numeric(as.character(Dat_Loc_Coords_Filter[,"Lat"]))
	long <- as.numeric(as.character(Dat_Loc_Coords_Filter[,"Long"]))
	startDay <- rep(as.character(startDate),length(Dat_Env)) #as.vector(paste(yrNamesSpl,"-05-01",sep=""))
	endDay <- rep(as.character(endDate),length(Dat_Env)) #paste(yrNamesSpl,"-11-20",sep="")
		
       
    #EnvRtype::
    ### Get weather data for the loc coords within the time interval 

   

	env.data <- get_weather_terra(env.id = Dat_Env,lat = lat,
								  lon =long,start.day = startDay,
								  end.day = endDay, country = Dat_Count,parallel=FALSE)


    return(env.data)

}

###

#' Title of getEnvKernel
#'
#' Description of what getEnvKernel does.
#'
#' @param Env.Data Description of Env.Data.
#' @param process=FALSE Description of process=FALSE.
#' @param Gaussian=FALSE Description of Gaussian=FALSE.
#'
#' @importFrom EnvRtype processWTH W_matrix env_kernel
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of getEnvKernel
#' result <- getEnvKernel(...)
#' }

getEnvKernel <- function(Env.Data,process=FALSE,Gaussian=FALSE){
   
    env.Dat <- Env.Data
	
	if(process==FALSE){
	   
		var3 <- c("T2M","T2M_MAX","T2M_MIN","PRECTOT","WS2M","RH2M","T2MDEW","ALLSKY_SFC_LW_DWN","ALLSKY_SFC_SW_DWN")
		EC.Dat = W_matrix(env.data = env.Dat,env.id="env", var.id = var3) 
		KE <- list(W = env_kernel(env.data = EC.Dat,gaussian=Gaussian)[[2]])
    }else if(process==TRUE){
	
		prData <- processWTH(env.data = env.Dat)
		var4 <- c("T2M","T2M_MAX","T2M_MIN","PRECTOT","WS2M","RH2M","T2MDEW","ALLSKY_SFC_LW_DWN","ALLSKY_SFC_SW_DWN","RTA","VPD","SPV","ETP","PETP","GDD","FRUE","T2M_RANGE")
			
		EC.Dat.Pr = W_matrix(env.data = prData,env.id="env", var.id = var4) 
		rmECVar <- which(apply(EC.Dat.Pr,2,function(x) length(which(x %in% NaN)))!=0)
		EC.Dat.Pr.Filt <- EC.Dat.Pr[,-rmECVar]

		KE = list(W = env_kernel(env.data = EC.Dat.Pr.Filt,gaussian=Gaussian)[[2]])
	}
	
	return(KE)
   }
	
# ### Based on four variables

	# var2 = c("PRECTOT","T2M","ALLSKY_SFC_SW_DWN","RH2M")
	# EC.Dat = W_matrix(env.data = env.data.Dat,env.id="env", var.id = var2) 
	# KE.Raw <- list(W = env_kernel(env.data = EC.Dat,gaussian=FALSE)[[2]])

###

#' Title of plotEnvRel
#'
#' Description of what plotEnvRel does.
#'
#' @param KE Description of KE.
#' @param outDirPath Description of outDirPath.
#'
#' @importFrom gplots heatmap.2 bluered
#' @importFrom grDevices png dev.off
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of plotEnvRel
#' result <- plotEnvRel(...)
#' }

plotEnvRel <- function(KE,outDirPath){
    #dev.off()
	#### Heatmap of Raw Weather Data  ##key.par=list(mar=c(1,4,2,1))
	par(oma = c(1, 1, 1, 1), mar = c(5, 4, 4, 2) + 0.1) # Adjust margins to avoid the 'pin' issue
	Hmp_Plot <- heatmap.2(KE$W,dendrogram="row",Colv=NA,col=bluered(75),scale="none",margins=c(10,10),key=TRUE,keysize=0.8,key.title="Env Relationship",
			  symm=TRUE,key.par=list(mar=c(5,1,1,1)),cexRow=1.5,cexCol=1.5,trace='none')

    OFN <- paste(outDirPath,"/","Environmental_Relationship.png",sep="")
	
    png(OFN,width=1024,height=768,pointsize=20)
    print(Hmp_Plot)
	Sys.sleep(5)
	dev.off()
    return(Hmp_Plot)

}

 







#' Title of syncEnvPhenoDat
#'
#' Description of what syncEnvPhenoDat does.
#'
#' @param KE Description of KE.
#' @param LocCoord Description of LocCoord.
#' @param OtherLoc Description of OtherLoc.
#'
#' @return Description of what is returned.
#' @examples
#' \dontrun{
#' # Example usage of syncEnvPhenoDat
#' result <- syncEnvPhenoDat(...)
#' }

syncEnvPhenoDat <- function(KE,LocCoord,OtherLoc){

 ke <- KE
 # Create a named vector to map abbreviations to full names
  if(OtherLoc==TRUE){
	  level_mapping <- setNames(as.character(LocCoord[,"OtherLocName"]),as.character(LocCoord[,"Location"]))
	  colnames(ke$W) <- (level_mapping[colnames(ke$W)])
	  rownames(ke$W) <- (level_mapping[rownames(ke$W)])
  }
  
  return(ke)
 }

### 
