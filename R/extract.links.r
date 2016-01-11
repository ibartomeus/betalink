#extract betaOS links

extract.links  <-  function (W, ...){
  require(reshape)
  i = 1
  metaweb = metaweb(W, force_binary = FALSE)$web
  Out <- data.frame(w = NA, P = NA, A = NA, lost = NA, freq = NA)
  for (j in 1:length(W)) {
    w <- W[[j]]
    meta_melted <- subset(melt(metaweb[rownames(w), colnames(w), drop=FALSE]), 
                          value > 0)
    meta_melted$link <- paste(meta_melted$X1, meta_melted$X2, 
                              sep = "_")
    w_melted <- subset(melt(as.matrix(w)), value > 0)
    w_melted$link <- paste(w_melted$X1, w_melted$X2, sep = "_")
    n <- which(!meta_melted$link %in% w_melted$link)
    meta_melted$lost <- ifelse(c(1:nrow(meta_melted)) %in% 
                                 n, 1, 0)
    w_Out <- data.frame(w = rep(i, nrow(meta_melted)), P = meta_melted$X1, 
                        A = meta_melted$X2, lost = meta_melted$lost, freq = meta_melted$value)
    Out <- rbind(Out, w_Out)
    i = i + 1
  }
  Out <- Out[-1, ]
  Out
}
