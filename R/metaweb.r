metaweb = function(W, force_binary = TRUE,...){
	Lo = unique(unlist(lapply(W,colnames)))
	Up = unique(unlist(lapply(W,rownames)))
	meta = matrix(0,ncol=length(Lo),nrow=length(Up))
	colnames(meta) = Lo
	rownames(meta) = Up
	co.oc = meta
	for(w in W){
        if(force_binary) w[w>0] = 1
		meta[rownames(w),colnames(w)] = meta[rownames(w),colnames(w)] + w
		co.oc[rownames(w),colnames(w)] = co.oc[rownames(w),colnames(w)] + 1
	}
	null.template = meta/co.oc
	null.template[is.nan(null.template)] = 0
	if(force_binary) meta[meta>0] = 1
	return(list(web=meta,template=null.template, cooc = co.oc))
}
