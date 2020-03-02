library(dplyr)
library(data.table)
epsilon <- 1e-10 # A very small number

getMatrix <- function(m,d1,d2){
  m = matrix(m, nrow = length(d1), ncol = length(d2))
  colnames(m) = d2
  row.names(m) = d1
  return(m)
}

# Read inputs -------------------------------------------------------------

drugs <- read.csv("./Data/external/drug links.csv")
drugs <- subset(drugs, Drug.Type == "SmallMoleculeDrug")

drugMetaData = read.csv("./Data/DB_full_meta_data.csv",header = TRUE,sep = ",", stringsAsFactors = FALSE)
drugs.ws <- drugMetaData[which(is.na(drugMetaData["Smiles"])),1] %>%
            setdiff(drugs[,1], .)

drugs.ws.app <- intersect(drugs.ws, drugMetaData[grepl("approved",as.character(drugMetaData[,3])),1])

sc.m <- as.matrix(read.csv("./Data/chem_similarity.csv", row.names = 1))
sp.m <- as.matrix(read.csv("./Data/path_similarity_v2.csv",  row.names = 1)) # improved: used BioCor package (for KEGG)
st.m <- as.matrix(read.csv("./Data/target_similarity.csv",  row.names = 1))
sgoCC.m <- as.matrix(read.csv("./Data/GO_Sim/CC/GO_Sim_CC_combined.csv", row.names = 1))
sgoMF.m <- as.matrix(read.csv("./Data/GO_Sim/MF/GO_Sim_MF_combined.csv", row.names = 1))
sgoBP.m <- as.matrix(read.csv("./Data/GO_Sim/BP/GO_Sim_BP_combined.csv", row.names = 1))

# Remove self-similarities, sub-setting drug-pairs, -----------------------------
# denoising similarity matrices
# and finally Aggregate and p-value calculation     -----------------------------

diag(sc.m) <- 0; diag(sp.m) <- 0; diag(st.m) <- 0; diag(sgoCC.m) <- 0; diag(sgoMF.m) <- 0; diag(sgoBP.m) <- 0;
sc.m <- sc.m[which(rownames(sc.m) %in% drugs.ws),]
sc.m <- sc.m[,which(colnames(sc.m) %in% drugs.ws.app)]
sp.m <- sp.m[which(rownames(sp.m) %in% drugs.ws),]
sp.m <- sp.m[,which(colnames(sp.m) %in% drugs.ws.app)]
st.m <- st.m[which(rownames(st.m) %in% drugs.ws),]
st.m <- st.m[,which(colnames(st.m) %in% drugs.ws.app)]
sgoCC.m <- sgoCC.m[which(rownames(sgoCC.m) %in% drugs.ws),]
sgoCC.m <- sgoCC.m[,which(colnames(sgoCC.m) %in% drugs.ws.app)]
sgoMF.m <- sgoMF.m[which(rownames(sgoMF.m) %in% drugs.ws),]
sgoMF.m <- sgoMF.m[,which(colnames(sgoMF.m) %in% drugs.ws.app)]
sgoBP.m <- sgoBP.m[which(rownames(sgoBP.m) %in% drugs.ws),]
sgoBP.m <- sgoBP.m[,which(colnames(sgoBP.m) %in% drugs.ws.app)]

# although it's less likely to have any NAs is n  sc.m
# since we've discared drugs without any structures.
sc.m <- ifelse(is.na(sc.m), 0, sc.m) 
# epsilon instead here just be fair for some drugs
# there may not have any pathway associated
# sp.m <- ifelse(is.na(sp.m), epsilon, sp.m) 
# sp.m <- ifelse(is.na(sp.m), 0, sp.m) 
st.m <- ifelse(is.na(st.m), 0, st.m)
sgoCC.m <- ifelse(is.na(sgoCC.m), epsilon, sgoCC.m) 
sgoMF.m <- ifelse(is.na(sgoMF.m), epsilon, sgoMF.m) 
sgoBP.m <- ifelse(is.na(sgoBP.m), epsilon, sgoBP.m) 

m <- rowMeans(cbind(c(sc.m), c(st.m), c(sp.m), c(sgoCC.m), c(sgoMF.m), c(sgoBP.m)), na.rm = TRUE)
tmp <- m
tmp[tmp<0.05] = NA
p <- pnorm(scale(tmp),lower.tail = FALSE)
p.adj <- p.adjust(p, method = "fdr")
rm(tmp)
# Build the network
d1    = drugs.ws 
d2    = drugs.ws.app
p     = getMatrix(p, d1,d2)
p.adj = getMatrix(p.adj, d1, d2)
m     = getMatrix(m,d1, d2)


# Output network ----------------------------------------------------------
# percentile = quantile(m, 0.99, na.rm = TRUE)

# library(foreach)
# library(doParallel)
# cl <- makeCluster(122)
# registerDoParallel(cl)
# 
# out <- c()
# s <- Sys.time()
# out <- foreach(i = 1:length(d1), .combine='rbind') %:% 
#   foreach(j = 1:length(d2), .combine='rbind') %dopar% {
#     di <- d1[i]
#     dj <- d2[j]
#     if(!is.na(m[di,dj]) & m[di,dj] > percentile){
#       cbind(di,dj,sc.m[di,dj],st.m[di,dj],sp.m[di,dj], m[di,dj], p.adj[di,dj], p[di,dj])
#       # r
#     }
#   }
# e <- Sys.time()
# print(e-s)

# out_perc <- c()
out_pval <- c()
s <- Sys.time()
counter <- 0
elapsed <- 0
for (i in d1) {
  for (j in d2) {
    # if(!is.na(m[i,j]) & m[i,j] > percentile){
    #   r <- c(i,j,sc.m[i,j],st.m[i,j],sp.m[i,j],sgoCC.m[i,j], sgoMF.m[i,j], sgoBP.m[i,j], m[i,j], p[i,j], p.adj[i,j])
    #   out_perc <- rbind(out_perc,r)
    # }
    if(!is.na(p[i,j]) & p[i,j] < 0.05){
      r <- c(i,j,sc.m[i,j],st.m[i,j],sp.m[i,j], sgoCC.m[i,j], sgoMF.m[i,j], sgoBP.m[i,j], m[i,j], p[i,j], p.adj[i,j])
      out_pval <- rbind(out_pval,r)
    }
    # if(!is.na(p[i,j])){
    #   r <- c(i,j,sc.m[i,j],st.m[i,j],sp.m[i,j], sgoCC.m[i,j], sgoMF.m[i,j], sgoBP.m[i,j], m[i,j], p[i,j], p.adj[i,j])
    #   out_pval <- rbind(out_pval,r)
    # }
  }
  e <- Sys.time()
  elapsed <- elapsed + (e-s)
  print(paste0((counter/length(d1))*100, "% fininsed in: ", elapsed, " seconds: number of drug-pairs: ", nrow(out_pval)))
  counter <- counter+1
}
# out1_perc <- as.data.frame(out_perc)
out1_pval <- as.data.frame(out_pval)
# colnames(out1_perc) <- c("ID1","ID2","Chem_similarity","Target_similarity","Pathway_similarity","GO_CC_Similarity", "GO_MF_Similarity","GO_BP_Similarity", "rowMeans","p_value","adjP_value")
colnames(out1_pval) <- c("ID1","ID2","Chem_similarity","Target_similarity","Pathway_similarity","GO_CC_Similarity", "GO_MF_Similarity","GO_BP_Similarity", "rowMeans","p_value","adjP_value")
# fwrite(out1_perc,"./Data/new_net_info_V4_perc.csv", na="NA", quote = F, row.names = F)
fwrite(out1_pval,"./Data/new_net_info_V7_pval.csv", na="NA", quote = F, row.names = F)
# write.csv(out, "new_net_info.csv",quote = F, row.names = F, col.names = F)
