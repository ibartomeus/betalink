#install.packages(devtools)
library(devtools)
install_github("ibartomeus/betalink", ref = "new_features") #just once

library(betalink)

#test functions
#fake data
w1 <- matrix(ncol = 6, nrow = 6, c(1,1,1,1,1,0,
                                   1,1,0,1,0,1,
                                   1,0,1,0,0,0,
                                   1,0,0,0,0,0,
                                   0,1,0,1,0,0,
                                   1,0,0,0,0,0), 
             dimnames = list(c("a","b", "c", "d", "e","f"), 
                             c("A", "B", "C","D","E","F")), byrow = TRUE)
 
w2 <- matrix(ncol = 5, nrow = 5, c(1,1,1,1,1,
                                   0,1,0,1,0, #small change in OS (2 links)
                                   1,0,0,0,0, 
                                   1,0,1,0,0, #I also add one link to turnover
                                   0,1,0,1,0),
             dimnames = list(c("a","b", "c", "g", "h"), 
                             c("A", "B", "C","G","H")), byrow = TRUE)

w3 <- matrix(ncol = 6, nrow = 6, c(2,1,2,1,2,0,
                                   2,1,0,1,0,2,
                                   1,0,2,0,0,0,
                                   1,0,0,0,0,0,
                                   0,2,0,1,0,0,
                                   2,0,0,0,0,0), 
             dimnames = list(c("a","b", "c", "d", "e","f"), 
                             c("A", "B", "C","D","E","F")), byrow = TRUE)
 
#reproduce Poissot claulction
betalink(w1,w2, bf = B15)
#quantitative version
betalink(w1,w2, quant = TRUE) #same result for links (species are quant also, so they differ)
betalink(w1,w3, quant = TRUE)

#second decomposition
betalink(w1, w2, calculate_2nd_decomposition = TRUE)

#second decomposition for quant
betalink(w1, w2, calculate_2nd_decomposition = TRUE, quant = TRUE)

#extract links
extract.links(W = list(w1,w2,w3)) 

