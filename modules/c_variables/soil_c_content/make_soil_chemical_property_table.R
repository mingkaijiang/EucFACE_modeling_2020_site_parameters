make_soil_chemical_property_table <- function(bk_density){

    ### download the data
    infile <- "FACE_P0014_ALL_BasicSoilProperties_L1_2012.csv"
    if(!file.exists(paste0("download/", infile))) {
        download_soil_p_data()
    }

    ### read in data - soil property data
    myDF2 <- read.csv(file.path(getToPath(), 
                                "FACE_P0014_ALL_BasicSoilProperties_L1_2012.csv"))
    
    myDF3 <- read.csv(file.path(getToPath(), 
                                "FACE_P0014_ALL_BasicSoilProperties_L1_2013.csv"))   
    
    myDF4 <- read.csv(file.path(getToPath(), 
                                "FACE_P0014_ALL_BasicSoilProperties_L1_2014.csv"))
    
    myDF5 <- read.csv(file.path(getToPath(), 
                                "FACE_P0014_ALL_BasicSoilProperties_L1_2015.csv"))
    
    myDF2 <- myDF2[,1:10]
    myDF3 <- myDF3[,1:10]
    myDF4 <- myDF4[,1:10]
    myDF5 <- myDF5[,1:10]
    
    colnames(myDF2) <- c("Date", "SampleNumber", "ring", "plot", "depth", "pH", "gmc", "totC", "totN", "totP_ppm")
    colnames(myDF3) <- c("Date", "SampleNumber", "ring", "plot", "depth", "pH", "gmc", "totC", "totN", "totP_ppm")
    colnames(myDF4) <- c("Date", "SampleNumber", "ring", "plot", "depth", "pH", "gmc", "totC", "totN", "totP_ppm")
    colnames(myDF5) <- c("Date", "SampleNumber", "ring", "plot", "depth", "pH", "gmc", "totC", "totN", "totP_ppm")
    
    myDF <- rbind(myDF2, myDF3, myDF4, myDF5)
    myDF$Date <- dmy(myDF$Date)
    
    ### get rid of spaces in the variable "Depth"
    myDF$depth <- as.character(myDF$depth)
    myDF$depth <- factor(gsub(" ", "", myDF$depth, fixed = TRUE)) 
    
    ### note that all data in 2015 are missing. Remove them.
    myDF <- subset(myDF,Date<as.Date("2015-01-01"))
    
    ### merge soil C with bulk density
    mydat <- merge(myDF,bk_density,by.x=c("depth", "ring"), by.y=c("Depth", "ring"))

    ### only choose the start of the experiment
    #outDF <- subset(mydat, Date=="2012-06-17")
    
    outDF2 <- summaryBy(pH+totC+totN+totP_ppm~depth, FUN=c(mean,sd),
                        data=mydat, keep.names=T, na.rm=T)
    
    return(outDF2)
    
}
