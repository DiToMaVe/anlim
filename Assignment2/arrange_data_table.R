# Collect spatial data from workshop data (lecture 6)
#####################################################

# Load workshop data
load("Data.Rdata")

# Link coordinates to postal codes
pc_long_lat <- subset(DT, select=c(PC,LONG,LAT))
pc_long_lat <- subset(pc_long_lat, !duplicated(PC))
# Convert data.table to data.frame
pc_long_lat <- as.data.frame(pc_long_lat)

# Load assignment data
DT_A2 = read.csv("Assignment2.csv")

# Add spatial coordinates corresponding with postal code
DT_A2$id <- 1:nrow(DT_A2)
DT_long_lat_merged <- merge(DT_A2,pc_long_lat,by="PC", sort=FALSE)

# the merge() function messed up the original order of the data frame
DT_long_lat <- DT_long_lat_merged[order(DT_long_lat_merged$id), ]
DT_A2$id <- NULL
DT_long_lat$id <- NULL

# Check wether the reshuffeling was done correctly
head(DT_long_lat)
head(DT_A2)
check <- match_df(DT_A2, DT_long_lat)

# Store appended data frame under correct name
DT <- DT_long_lat

# Add id column (for reasons of compatibility)
DT$ID <- seq(from=1, to = dim(DT)[1], by = 1)

# Add column with average cost of claim
DT$AVG <- DT$AMOUNT/DT$NCLAIMS
# Replace NaN (due to division through zero) with zero
DT[is.na(DT)] <- 0

# Create data frame for severity modeling by removing rows without reported claims
DT.sev <- DT[DT$NCLAIMS>0,]

# Convert PC from integer to factor for reasons of compatibility
DT$PC <- factor(DT$PC)

# Convert data frame to data table
DT <- as.data.table(DT)
DT.sev <- as.data.table(DT.sev)

