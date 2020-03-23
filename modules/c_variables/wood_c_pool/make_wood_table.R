
# Make the live wood C pool
make_wood_table <- function(ring_area){
    
    #### download the data from HIEv
    download_diameter_data()
    
    #### read in 2012-15 data sets
    f13 <- read.csv(file.path(getToPath(), "FACE_P0025_RA_TREEMEAS_2012-13_RAW-V1.csv"))
    f14 <- read.csv(file.path(getToPath(), "FACE_P0025_RA_TREEMEAS_2013-14_RAW_V1.csv"))
    f15 <- read.csv(file.path(getToPath(), "FACE_P0025_RA_TREEMEAS_2015_RAW_V1.csv"))
    f16 <- read.csv(file.path(getToPath(), "FACE_P0025_RA_TREEMEAS_2016_RAW_V1.csv"))
    # this file is not on HIEv yet!
    f12 <- read.csv("temp_files/EucFACE_dendrometers2011-12_RAW.csv")
    
    #### Read in additional files that I used when doing the data analysis
    classif <- read.csv("download/FACE_AUX_RA_TREE-DESCRIPTIONS_R_20130201.csv",stringsAsFactors = FALSE)
    classif$Active.FALSE.means.dead.[classif$Tree == 608] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 125] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 206] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 210] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 212] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 510] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 518] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 520] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 524] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 527] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 531] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 605] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 615] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 616] <- FALSE  # This tree dead
    classif$Active.FALSE.means.dead.[classif$Tree == 617] <- FALSE  # This tree dead
    #classif$Active.FALSE.means.dead.[classif$Tree == 101] <- FALSE  # This tree dead in 2018
    #classif$Active.FALSE.means.dead.[classif$Tree == 219] <- FALSE  # This tree dead in 2018
    #classif$Active.FALSE.means.dead.[classif$Tree == 220] <- FALSE  # This tree dead in 2018
    #classif$Active.FALSE.means.dead.[classif$Tree == 621] <- FALSE  # This tree dead in 2018
    
    
    #### Merge the files
    all <- merge(classif,f12,by=c("Tree","Ring","CO2.trt"))
    all <- merge(all,f13,by=c("Tree","Ring","CO2.trt")) 
    all <- merge(all,f14,by=c("Tree","Ring","CO2.trt"))  
    all <- merge(all,f15,by=c("Tree","Ring","CO2.trt"))
    all <- merge(all,f16,by=c("Tree","Ring","CO2.trt"))
    
    #### remove dead trees
    all$Active.FALSE.means.dead.[is.na(all$Active.FALSE.means.dead.)] <- "TRUE"
    all <- subset(all, Active.FALSE.means.dead.== TRUE)
    #all <- all[complete.cases(all),]
    
    #### remove "CORR" columns and dead column
    uncorr <- all[,-grep("CORR",names(all))]
    uncorr <- uncorr[,-grep("Coor",names(uncorr))]
    uncorr <- uncorr[,names(uncorr) != "Active.FALSE.means.dead."]
    
    #### make a long-form version of dataframe
    long <- reshape(uncorr,idvar="Tree",varying=list(7:58),direction="long")
    dates <- names(uncorr)[7:58]
    long$Date <- c(rep(Sys.Date(),length(long$time)))  #wasn't sure how else to make this column date type
    for (i in (1:length(long$time))) {
        long$Date[i] <- as.Date(dates[long$time[i]],format="X%d.%m.%Y")
    }
    long <- renameCol(long,c("X17.02.2011"),c("diam"))
    
    long$diam <- as.numeric(long$diam)
    
    #### add biomass to long-form dataframe
    long$biom <- allom_agb(long$diam)  # in kg DM
    
    #### The bark removal affects the diameters mid-year. 
    #### Hence, just calculate biomass once per year 
    #### Specify dates here - may update this to March in future
    dates <- c(as.Date("2012-12-20"),as.Date("2013-12-20"),
               as.Date("2014-12-23"),as.Date("2015-12-14"),
               as.Date("2016-12-21"))
    data <- long[long$Date=="2012-12-20",]
    
    ### calculate basal area
    data$basal_area <- (pi/4) * data$Diameter^2
    
    ### calculate total number of trees per ring
    outDF1 <- summaryBy(Diameter~Ring, FUN=mean, data=data, keep.names=T, na.rm=T)
    outDF2 <- summaryBy(biom+basal_area~Ring, FUN=sum, data=data, keep.names=T, na.rm=T)
    outDF3 <- summaryBy(Height~Ring, FUN=max, data=data, keep.names=T, na.rm=T)

    ### return biomass in unit of kg /m-2
    outDF1$Biomass <- outDF2$biom / ring_area
    outDF1$BA <- outDF2$basal_area / ring_area
    outDF1$Height <- outDF3$Height
    
    ### count number of trees per plot
    for (i in 1:6) {
        tmpDF <- subset(data, Ring==i & Date=="2012-12-20")
        outDF1[outDF1$Ring==i, "Trees"] <- nrow(tmpDF)
    }
    
    ### unit conversions
    # from no. tree per ring to no. tree per hecture
    outDF1$Trees <- outDF1$Trees / ring_area * 10000
    
    # from kg DM per m-2 to mg DM per hectare
    outDF1$Biomass <- outDF1$Biomass * 10000 / 1000
    
    outDF1$Trt <- "aCO2"
    outDF1$Trt[outDF1$Ring%in%c(1,4,5)] <- "eCO2"
    
    
    ### outDF
    outDF <- summaryBy(Diameter+Biomass+BA+Height+Trees~Trt, FUN=c(mean, sd), data=outDF1, keep.names=T, na.rm=T)
    

    return(outDF)
    
}