##
# script that downloads and stores the data needed to run the crossvalstudy
# Note prescified location ../data/
##

data_location = '../data/'

data.url <- "http://www.math.chalmers.se/~bodavid/data/downscaling/"
files <- c("index.norway.1x1kmAggregert.sav",
           "lat.norway.sav",
           "lon.norway.sav",
           "dataTrainingNotRandom.obs.sav",
           "dataTrainingNotRandom.bcm.sav",
           "dataTrainingNotRandom.era.sav",
           "dataTestNotRandom.bcm.sav",
           "dataTestNotRandom.obs.sav",
           "dataTestNotRandom.era.sav")

cat(paste('Downloading files to ',data_location,'. May take sometime.\n'))
for(i in 1:length(files))
  download.file(paste(data.url,files[i],sep=""), paste(data_location,files[i],sep=""))

cat('Donwloading files done\n')
