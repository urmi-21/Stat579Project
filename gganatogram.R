#gganatogram
#devtools::install_github("jespermaag/gganatogram")
library(gganatogram)
library(dplyr)
library(viridis)
library(gridExtra)


#example
organPlot <- data.frame(organ = c("heart", "leukocyte", "nerve", "brain", "liver", "stomach", "colon"), 
                        type = c("circulation", "circulation",  "nervous system", "nervous system", "digestion", "digestion", "digestion"), 
                        colour = c("red", "red", "purple", "purple", "orange", "orange", "orange"), 
                        value = c(100, 5, 100, 8000, 2000, 5, 5), 
                        stringsAsFactors=F)

head(organPlot)

gganatogram(data=organPlot, fillOutline='#a6bddb', organism='human', sex='male', fill="colour")+ 
  theme_void()+  scale_fill_gradient(low = "white", high = "red")


hgMale_key$organ

gganatogram(data=hgMale_key, fillOutline='#a6bddb', organism='human', sex='male', fill="colour") +theme_void()

gganatogram(data=organPlot, fillOutline='#a6bddb', organism='human', sex='male', fill="value") + 
  theme_void() +
  scale_fill_gradient(low = "white", high = "red")
