betalink.dist <-  function(W, calculate_2nd_decomposition = TRUE, ...){
	dWN = matrix(NA,ncol=length(W),nrow=length(W))
	colnames(dWN) = names(W)
	rownames(dWN) = names(W)
	dOS = dWN
	dS = dWN
	dST = dWN
	dContrib = dWN
    #add new metrics
	dU_3 = dWN
    dUrich = dWN
    dL_3 = dWN
    dLrich = dWN
    dS_3 = dWN
    dSrich = dWN
    dWN_3 = dWN
    dWNrich = dWN
    dOS_3 = dWN
    dOSrich = dWN
    #dST_3 = dWN
    #dSTrich = dWN
	for(i in c(1:(length(W)-1))){
		for(j in c((i+1):(length(W)))){
			partition = betalink(W[[i]],W[[j]], calculate_2nd_decomposition = calculate_2nd_decomposition)
            dWN[j,i]		= partition$WN
			dOS[j,i]		= partition$OS
			dS[j,i]			= partition$S
			dST[j,i]		= partition$ST
			dContrib[j,i]    = partition$contrib
			dU_3[j,i]	    = partition$U_3
			dUrich[j,i]    = partition$Urich
			dL_3[j,i]    = partition$L_3
			dLrich[j,i]    = partition$Lrich
			dS_3[j,i]    = partition$S_3
			dSrich[j,i]    = partition$Srich
			dWN_3[j,i]    = partition$WN_3
			dWNrich[j,i]    = partition$WNrich
			dOS_3[j,i]    = partition$OS_3
			dOSrich[j,i]    = partition$OSrich
			#dST_3[j,i]    = partition$ST_3
			#dSTrich[j,i]    = partition$STrich
		}
	}
	distances = list(WN=dWN, OS=dOS, S=dS, ST=dST, contrib=dContrib,
	                 U_3=dU_3, Urich=dUrich, L_3=dL_3, Lrich=dLrich,
	                 S_3=dS_3, Srich=dSrich, WN_3=dWN_3, WNrich=dWNrich,
	                 OS_3=dOS_3, OSrich=dOSrich) #, ST_3=dST_3, STrich=dSTrich)
	distances = lapply(distances,as.dist)
	return(distances)
}
