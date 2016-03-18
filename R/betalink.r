#Poisot et al. 2012 Ecol Lett. separate link beta diversity (WN) in Dissimilarity of 
    #interactions established between species common to both networks (OS)
    #and Dissimilarity of interactions due to species turnover (ST), 
#Carvalho et al 2012 GEB. proposed that species betadiversity can be decomposed 
    #in betadiversity due to turnover and due to species richness change. 
    #This same framework can be applied to betadiversity of links.
#I further separate here OS in diferences due to link replacement,
    #and diferences due to increases in conectance. That can be done applying Carvalho 
    #formulas on the common subweb.
    #However, for separating interactions change due to species turnover (ST), 
    #link dissimilarity due to links replacement and link dissimilarity due to 
    #link connectance change I use the same substracting aproach than Poisot.

#using poisot naming:
#WN = OS + ST (subweb + species turnover) #as per Poisot 2012
#CC = rich + _3 (Link richness + total link turnover) #as per Carvalho 2012
#CC and WN are equal when B15 is used

#CC = (b+c)/(a+b+c) #B15 in measures.r
#_3 = 2*min(b,c)/(a+b+c) #B_3 in measures.r
#rich = |b-c|/(a+b+c) #Brich in measures.r

#Hence, for ST I define the component of link replacement (ST_3) as:
#ST_3 = WN_3 - OS (Link turnover not including the subwew) #not works, as OS < > WN_3
#And the complementary component of link richness would be: 
#STrich = ST - ST_3

#Total decomposition:
#WN = OS + ST
#where
#OS = OS_3 + OSrich
#ST = ST_3 + STrich

#I also add the quantitative form based on of calculating a,b and c from 
  #the Ruzicka = (B+C)/(A+B+C) distance coefficients, where
# A = sum of the intersections (or minima) of species abundances at two sites,
# B = sum at site 1 minus A, C = sum at site 2 minus A.#Ružicˇka differences (Legendre et al 2014 GEB),
# which can be used with the second decomposition also.

betalink  <-  function(w1, w2, bf=B01, quant = FALSE, calculate_2nd_decomposition = FALSE){
  if(quant == FALSE){
    #need to force matrices to binary  
    w1[w1>1] <- 1
    w2[w2>1] <- 1
    if(calculate_2nd_decomposition == TRUE){ 
        bf = B15
        print("bf overwritten when calculate_2nd_decomposition = TRUE. Using bf = B15 (=(b+c)/(a+b+c))")
    }
    pmb = function(A,B) list(b=sum(!(A %in% B)),c=sum(!(B %in% A)),a=sum(B %in% A))
    sp1 = list(top=rownames(w1),bottom=colnames(w1),all=unique(c(colnames(w1),rownames(w1))))
    sp2 = list(top=rownames(w2),bottom=colnames(w2),all=unique(c(colnames(w2),rownames(w2))))
    beta_U = bf(pmb(sp1$top,sp2$top))
    beta_L = bf(pmb(sp1$bottom,sp2$bottom))
    beta_S = bf(pmb(sp1$all,sp2$all))
    #calculate species beta diversity decompositions
    if(calculate_2nd_decomposition == TRUE){ 
    beta_U_3 = B_3(pmb(sp1$top,sp2$top))
    beta_Urich = Brich(pmb(sp1$top,sp2$top))
    beta_L_3 = B_3(pmb(sp1$bottom,sp2$bottom))
    beta_Lrich = Brich(pmb(sp1$bottom,sp2$bottom))
    beta_S_3 = B_3(pmb(sp1$all,sp2$all))
    beta_Srich = Brich(pmb(sp1$all,sp2$all))
    }
    # Common species
    Csp = sp1$all[sp1$all %in% sp2$all]
    CUsp = sp1$top[sp1$top %in% sp2$top]
    CLsp = sp1$bottom[sp1$bottom %in% sp2$bottom]
    if((length(CUsp)>0) & (length(CLsp)>0))
    {
        w1Con = w1[CUsp,CLsp]
        w2Con = w2[CUsp,CLsp]
        nCon = sum((w1Con == w2Con) & (w1Con == 1))
        pmBos = list(b=sum(w1Con)-nCon,c=sum(w2Con)-nCon,a=nCon)
        pmBwn = list(b=sum(w1)-nCon,c=sum(w2)-nCon,a=nCon)
        beta_OS = bf(pmBos)
        beta_WN = bf(pmBwn)
        if(is.na(beta_OS)) beta_OS = 0
        if(is.na(beta_WN)) beta_WN = 0
        beta_ST = beta_WN - beta_OS
        #calculate link second decomposition
        if(calculate_2nd_decomposition == TRUE){ 
            beta_3 = B_3(pmBwn)
            beta_rich = Brich(pmBwn)
            beta_OS_3 = B_3(pmBos)
            beta_OSrich = Brich(pmBos)    
            #beta_ST_3 = beta_3 - beta_OS
            #beta_STrich = beta_ST - beta_ST_3
            #beta_ST_3 = beta_3 - beta_OS_3
            #beta_STrich = beta_rich - beta_OSrich
            if(is.na(beta_3)) beta_3 = 0
            if(is.na(beta_rich)) beta_rich = 0
            if(is.na(beta_OS_3)) beta_OS_3 = 0
            if(is.na(beta_OSrich)) beta_OSrich = 0
            #if(is.na(beta_ST_3)) beta_ST_3 = 0
            #if(is.na(beta_STrich)) beta_STrich = 0
        }            
        if(beta_WN > 0){
            b_contrib = beta_ST / beta_WN
            #reltive contribs may be added here
        } else {
            b_contrib = 0
        }
    } else {
        beta_WN = 1 #why can you not calculate WN when all species are different?
        #NB makes this change after thinking carefully. I think no overlap implies maximum 
        #beta, not minimum!!
        beta_OS = 0
        beta_ST = 1
        b_contrib = 1
        #calculate link second decomposition 
        if(calculate_2nd_decomposition == TRUE){ 
            beta_3 = NA #when there is no overlap can't say why, right?
            beta_rich = NA
            beta_OS_3 = 0
            beta_OSrich = 0    
            #beta_ST_3 = 0
            #beta_STrich = 0
        }            
    }
    if(calculate_2nd_decomposition == TRUE){
        return(list(U = beta_U, L = beta_L, S = beta_S, OS = beta_OS, 
                    WN = beta_WN, ST = beta_ST, contrib = b_contrib,
                    U_3 = beta_U_3, Urich = beta_Urich, L_3 = beta_L_3, 
                    Lrich = beta_Lrich, S_3 = beta_S_3, Srich = beta_Srich, 
                    WN_3 = beta_3, WNrich = beta_rich, OS_3 = beta_OS_3, 
                    OSrich = beta_OSrich)) #, ST_3 = beta_ST_3,
                    #STrich = beta_STrich))
    }else{    
        return(list(U = beta_U, L = beta_L, S = beta_S, OS = beta_OS, 
                    WN = beta_WN, ST = beta_ST, contrib = b_contrib))
    }
  } 
  if(quant == TRUE){
    bf = B15
    print("bf overwritten when quant = TRUE. Using bf = B15 (=(b+c)/(a+b+c)). I am not sure if other indexes can be used")
    if(calculate_2nd_decomposition == TRUE){ 
      bf = B15
      print("bf overwritten when calculate_2nd_decomposition = TRUE. Using bf = B15 (=(b+c)/(a+b+c))")
    }
    # a sum of minima in among-site comparisons; b and c are its sum minus a.
    pmb = function(A,B) list(b= sum(A) - sum(pmin(A, B)), c = sum(B) - sum(pmin(A, B)), a = sum(pmin(A, B)))
    #species level
    rows_ <- unique(c(rownames(w1), rownames(w2)))
    cols_ <- unique(c(colnames(w1), colnames(w2)))  
    empty <- matrix(0, nrow = length(rows_), ncol = length(cols_), dimnames = list(rows_, cols_))
    w1full <- empty
    w2full <- empty
    w1full[rownames(w1), colnames(w1)]  <- w1
    w2full[rownames(w2), colnames(w2)]  <- w2
    sp1 = list(top=rowSums(w1full),bottom=colSums(w1full),all=c(rowSums(w1full),colSums(w1full)))
    sp2 = list(top=rowSums(w2full),bottom=colSums(w2full),all=c(rowSums(w2full),colSums(w2full)))
    beta_U = bf(pmb(sp1$top,sp2$top))
    beta_L = bf(pmb(sp1$bottom,sp2$bottom))
    beta_S = bf(pmb(sp1$all,sp2$all))
    #calculate species beta diversity decompositions
    if(calculate_2nd_decomposition == TRUE){ 
      beta_U_3 = B_3(pmb(sp1$top,sp2$top))
      beta_Urich = Brich(pmb(sp1$top,sp2$top))
      beta_L_3 = B_3(pmb(sp1$bottom,sp2$bottom))
      beta_Lrich = Brich(pmb(sp1$bottom,sp2$bottom))
      beta_S_3 = B_3(pmb(sp1$all,sp2$all))
      beta_Srich = Brich(pmb(sp1$all,sp2$all))
    }
    # Common species
    sp1 = list(top=rownames(w1),bottom=colnames(w1),all=unique(c(colnames(w1),rownames(w1))))
    sp2 = list(top=rownames(w2),bottom=colnames(w2),all=unique(c(colnames(w2),rownames(w2))))
    Csp = sp1$all[sp1$all %in% sp2$all]
    CUsp = sp1$top[sp1$top %in% sp2$top]
    CLsp = sp1$bottom[sp1$bottom %in% sp2$bottom]
    if((length(CUsp)>0) & (length(CLsp)>0))
    {
      w1Con = w1[CUsp,CLsp]
      w2Con = w2[CUsp,CLsp]
      nCon = sum(pmin(w1Con, w2Con))
      nfull = sum(pmin(w1full, w2full))
      #sum((w1Con == w2Con) & (w1Con == 1)) #this wont work with freq data!!!
      pmBos = list(b=sum(w1Con)-nCon,c=sum(w2Con)-nCon,a=nCon)
      #need to use them comparable in size
      pmBwn = list(b=sum(w1full)-nfull,c=sum(w2full)-nfull,a=nfull)
      beta_OS = bf(pmBos)
      beta_WN = bf(pmBwn)
      if(is.na(beta_OS)) beta_OS = 0
      if(is.na(beta_WN)) beta_WN = 0
      beta_ST = beta_WN - beta_OS
      #calculate link second decomposition
      if(calculate_2nd_decomposition == TRUE){ 
        beta_3 = B_3(pmBwn)
        beta_rich = Brich(pmBwn)
        beta_OS_3 = B_3(pmBos)
        beta_OSrich = Brich(pmBos)    
        #beta_ST_3 = beta_3 - beta_OS
        #beta_STrich = beta_ST - beta_ST_3
        #beta_ST_3 = beta_3 - beta_OS_3
        #beta_STrich = beta_rich - beta_OSrich
        if(is.na(beta_3)) beta_3 = 0
        if(is.na(beta_rich)) beta_rich = 0
        if(is.na(beta_OS_3)) beta_OS_3 = 0
        if(is.na(beta_OSrich)) beta_OSrich = 0
        #if(is.na(beta_ST_3)) beta_ST_3 = 0
        #if(is.na(beta_STrich)) beta_STrich = 0
      }            
      if(beta_WN > 0){
        b_contrib = beta_ST / beta_WN
        #reltive contribs may be added here
      } else {
        b_contrib = 0
      }
    } else {
      beta_WN = 1 #why can you not calculate WN when all species are different?
      #idem reasoning as before, no overlap, max WN
      beta_OS = 0
      beta_ST = 1
      b_contrib = 1
      #calculate link second decomposition 
      if(calculate_2nd_decomposition == TRUE){ 
        beta_3 = NA
        beta_rich = NA
        beta_OS_3 = 0
        beta_OSrich = 0    
        #beta_ST_3 = 0
        #beta_STrich = 0
      }            
    }
    if(calculate_2nd_decomposition == TRUE){
      return(list(U = beta_U, L = beta_L, S = beta_S, OS = beta_OS, 
                  WN = beta_WN, ST = beta_ST, contrib = b_contrib,
                  U_3 = beta_U_3, Urich = beta_Urich, L_3 = beta_L_3, 
                  Lrich = beta_Lrich, S_3 = beta_S_3, Srich = beta_Srich, 
                  WN_3 = beta_3, WNrich = beta_rich, OS_3 = beta_OS_3, 
                  OSrich = beta_OSrich)) #, ST_3 = beta_ST_3,
      #STrich = beta_STrich))
    }else{    
      return(list(U = beta_U, L = beta_L, S = beta_S, OS = beta_OS, 
                  WN = beta_WN, ST = beta_ST, contrib = b_contrib))
    }
  }
}
# 
# data to test
# 
# w1 <- matrix(ncol = 6, nrow = 6, c(1,1,1,1,1,0,
#                                   1,1,0,1,0,1,
#                                   1,0,1,0,0,0,
#                                   1,0,0,0,0,0,
#                                   0,1,0,1,0,0,
#                                   1,0,0,0,0,0),
#             dimnames = list(c("a","b", "c", "d", "e","f"),
#                             c("A", "B", "C","D","E","F")), byrow = TRUE)
# 
# w2 <- matrix(ncol = 5, nrow = 5, c(1,1,1,1,1,
#                                   0,1,0,1,0, #small change in OS (2 links)
#                                   1,0,0,0,0,
#                                   1,0,1,0,0, #I also add one link to turnover
#                                   0,1,0,1,0),
#             dimnames = list(c("a","b", "c", "g", "h"),
#                             c("A", "B", "C","G","H")), byrow = TRUE)
# 
# 
# x <- betalink(w1,w2, calculate_2nd_decomposition = TRUE)
# 
# predictions:
# (poisot decomposition)
# x$WN #is equal to 
# x$OS + x$ST
# 
# 
# (Carvalho decomposition)
# x$WN #is equal to 
# x$WN_3 + x$WNrich
# 
# (Total decomposition)
# x$WN #is equal to
# x$OS_3 + x$OSrich + x$ST_3 + x$STrich
#     
# calculate contributions:
# x$OSrich/x$OS # contribution of having higher conectance in one network among the species forming the common subweb
# x$STrich/x$ST # contribution of having higher number of links in one network among the non common species
