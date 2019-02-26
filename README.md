# R package betalink
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2577565.svg)](https://doi.org/10.5281/zenodo.2577565)

Measures of network dissimiliraty using R.

Building upon original package developed by Tim Poisot, I added the quantitative forms of calculating betadiversity, and a second decomposition of betadiversity in richness and turnover components. Those are on this branch. 

For this I am using diversity calculations based on the Colwell and Coddington measure (`B15` in Koleff et al. 2003) because it allows for paralel calculations of binary and quantitative data and it's decomposition. The measures is calculated as follow for binary data:

`B15 = (b+c)/(a+b+c)`

where `a` is the number of species common to both sites, `b` is the number of species exclusive to the first site, and `c` is the number of species exclusive to the second site (Koleff et al. 2003; Carvalho et al. 2012).   

Quantitative diversity calculations were based on the Ruzicka distance coefficient, which takes a parallel notation:

`βRuz = (B+C)/(A+B+C)`

where `A` is the sum of intersections (or minima) of species abundances at two sites, `B` is the sum at site 1 minus `A`, and `C` is the sum at site 2 minus `A` (Legendre 2014). 

This two proposed beta diversity metrics for both species and interactions can also be further partitioned to see if possible differences in beta diversity are due to real turnover of species/links, or to differences in number of species/links, following Carvalho et al. (2012).

Those are calculated as: 

`B15 = B_3 + Brich`  

where 

`B_3 = 2*min(b,c)/(a+b+c)` 

and 

`Brich = abs(b-c)/(a+b+c)`

This is analogous to the quantititive version `βRuz = B_3 + Brich` where `B_3 = 2*min(B,C)/(A+B+C)` and `Brich = abs(B-C)/(A+B+C)`


#How to use it?

```
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

w3 <- matrix(ncol = 6, nrow = 6, c(2,1,2,1,2,0, #quantitative
                                   2,1,0,1,0,2,
                                   1,0,2,0,0,0,
                                   1,0,0,0,0,0,
                                   0,2,0,1,0,0,
                                   2,0,0,0,0,0), 
             dimnames = list(c("a","b", "c", "d", "e","f"), 
                             c("A", "B", "C","D","E","F")), byrow = TRUE)
 
#reproduce Poisot origial claulction
betalink(w1,w2, bf = B15)
#quantitative version
betalink(w1,w2, quant = TRUE) #same result for links (species are quant also, and hence they differ)
betalink(w1,w3, quant = TRUE)

#second decomposition
betalink(w1, w2, calculate_2nd_decomposition = TRUE)

#second decomposition for quant
betalink(w1, w2, calculate_2nd_decomposition = TRUE, quant = TRUE)

#finally I made a fucntion to extract for each web its links and see which ones are lost or not when comparing with the metaweb.
extract.links(W = list(w1,w2,w3)) 

#Other functions like betalink.dist() work as expected.
```

