library(msa)
library(seqinr)
library(foreach)
library(doParallel)

newFunc <- function(mat, fasta, startEpoch, endEpoch){
  require(msa)
  require(seqinr)
  require(foreach)
  require(doParallel)
  # fasta <- readAAStringSet("/short/yr31/aa7970/azData/DD/Data/protein.fasta")

protSim <- function(seq1,seq2){
  return(Biostrings::pid(Biostrings::pairwiseAlignment(seq1,seq2)))
}

foreach(i = startEpoch:endEpoch, .combine='rbind') %:% 
  foreach(j = 1:nrow(mat), .combine='cbind', .packages='Biostrings') %dopar% {
    mat[i,j] = protSim(fasta[i],fasta[j]);
  }

write.csv(mat, file = paste0("protSim_", startEpoch, "_", endEpoch, ".csv"), quote = F, row.names = F)
}

fasta <- readAAStringSet("protein.fasta")
cl <- makeCluster(64)
registerDoParallel(cl)

n = 100
jump = 10
x = matrix(0, nrow = n, ncol = n)
start_time <- Sys.time()
# foreach(i = 1:(n-jump), .combine = "rbind") %dopar% {
#   newFunc(x,i,i+jump)
#}
newFunc(x,fasta, 1, 10)

end_time <- Sys.time()
print(end_time - start_time)
stopCluster(cl)
