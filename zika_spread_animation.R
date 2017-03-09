# The original codes were created by Rob Harrand on Kagg website.  

library(ggplot2) # Data visualization
library(readr) 
library(animation)
library(maptools)
setwd("~/Downloads")
zika = read.csv("cdc_zika.csv", stringsAsFactors = F, header = T) #Load the data

#Load the map data,
s = map_data("world")

zika$value = as.numeric(zika$value)
zika$report_date = as.Date(zika$report_date, '%Y-%m-%d')

#Exclude rows with missing data,
exc.Numbers = (is.na(zika$report_date) | (is.na(zika$value)))
zika = zika[!exc.Numbers, ]

s$colour = 0

firstdate = zika$report_date[order(zika$report_date)][1]
lastdate = zika$report_date[order(zika$report_date)][length(zika$report_date)]
duration = lastdate - firstdate

i=0
runsum_usa = 0
runsum_mexico = 0
runsum_panama = 0
runsum_nicaragua = 0
runsum_haiti = 0
runsum_guatemala = 0
runsum_el_salvador = 0
runsum_ecuador = 0
runsum_dominican_republic = 0
runsum_colombia = 0
runsum_argentina = 0
runsum_brazil = 0

totals = data.frame(runsum_brazil, runsum_argentina, runsum_colombia, runsum_dominican_republic, runsum_ecuador, runsum_el_salvador,
                    runsum_guatemala, runsum_haiti, runsum_nicaragua, runsum_panama, runsum_usa, runsum_mexico)

usa =  zika[grep("United_States", zika$location),]
mexico = zika[grep("Mexico", zika$location),]
panama = zika[grep("Panama", zika$location),]
nicaragua = zika[grep("Nicaragua", zika$location),]
haiti = zika[grep("Haiti", zika$location),]
guatemala = zika[grep("Guatemala", zika$location),]
el_salvador = zika[grep("El_Salvador", zika$location),]
ecuador = zika[grep("Ecuador", zika$location),]
dominican_republic = zika[grep("Dominican_Republic", zika$location),]
colombia = zika[grep("Colombia", zika$location),]
argentina = zika[grep("Argentina", zika$location),]
brazil = zika[grep("Brazil", zika$location),]


#Set animation time,
ani.options(interval = 0.025)


#Loop through the rows and save the gif...

saveGIF(while (i<duration) {
  
  temp_sum_usa = sum(usa$value[usa$report_date == firstdate+i])
  runsum_usa = runsum_usa + temp_sum_usa
  s$colour[grep("USA", s$region)] = log(runsum_usa)
  
  temp_sum_mexico = sum(mexico$value[mexico$report_date == firstdate+i])
  runsum_mexico = runsum_mexico + temp_sum_mexico
  s$colour[grep("Mexico", s$region)] = log(runsum_mexico)
  
  temp_sum_panama = sum(panama$value[panama$report_date == firstdate+i])
  runsum_panama = runsum_panama + temp_sum_panama
  s$colour[grep("Panama", s$region)] = log(runsum_panama)
  
  temp_sum_nicaragua = sum(nicaragua$value[nicaragua$report_date == firstdate+i])
  runsum_nicaragua = runsum_nicaragua + temp_sum_nicaragua
  s$colour[grep("Nicaragua", s$region)] = log(runsum_nicaragua)
  
  temp_sum_haiti = sum(haiti$value[haiti$report_date == firstdate+i])
  runsum_haiti = runsum_haiti + temp_sum_haiti
  s$colour[grep("Haiti", s$region)] = log(runsum_haiti)
  
  temp_sum_guatemala = sum(guatemala$value[guatemala$report_date == firstdate+i])
  runsum_guatemala = runsum_guatemala + temp_sum_guatemala
  s$colour[grep("Guatemala", s$region)] = log(runsum_guatemala)
  
  temp_sum_el_salvador = sum(el_salvador$value[el_salvador$report_date == firstdate+i])
  runsum_el_salvador = runsum_el_salvador + temp_sum_el_salvador
  s$colour[grep("El_Salvador", s$region)] = log(runsum_el_salvador)
  
  temp_sum_ecuador = sum(ecuador$value[ecuador$report_date == firstdate+i])
  runsum_ecuador = runsum_ecuador + temp_sum_ecuador
  s$colour[grep("Ecuador", s$region)] = log(runsum_ecuador)
  
  temp_sum_dominican_republic = sum(dominican_republic$value[dominican_republic$report_date == firstdate+i])
  runsum_dominican_republic = runsum_dominican_republic + temp_sum_dominican_republic
  s$colour[grep("Dominican_Republic", s$region)] = log(runsum_dominican_republic)
  
  temp_sum_colombia = sum(colombia$value[colombia$report_date == firstdate+i])
  runsum_colombia = runsum_colombia + temp_sum_colombia
  s$colour[grep("Colombia", s$region)] = log(runsum_colombia)
  
  temp_sum_argentina = sum(argentina$value[argentina$report_date == firstdate+i])
  runsum_argentina = runsum_argentina + temp_sum_argentina
  s$colour[grep("Argentina", s$region)] = log(runsum_argentina)
  
  temp_sum_brazil = sum(brazil$value[brazil$report_date == firstdate+i])
  runsum_brazil = runsum_brazil + temp_sum_brazil
  s$colour[grep("Brazil", s$region)] = log(runsum_brazil)
  
  print(m <- ggplot(s, aes(x=long, y=lat, group=group, fill=colour)) + #Set ggplot2
          
          geom_polygon(alpha=1) + #Set transparency
          
          geom_path(data = s, aes(x=long, y=lat, group=group), colour="black") + #Plot the Earth
          
          scale_fill_gradient(low = "green", high = "red", guide = "colourbar", limits=c(0,14.52374)) + #Set the colours
          
          labs(colour = "Zika cases") + 
          
          theme(plot.title = element_text(size = rel(2), colour = "blue")) + #Change the text size,
          
          ggtitle(paste("Spread of Zika Virus: ", firstdate+i))) 
  
  
  #Wait a while between plots,
  ani.pause()
  
  i=i+1
  
}, movie.name = "Zika_spread.gif", interval = 0.025, convert = "convert", ani.width = 500, 
ani.height = 350)