library(pbapply)
getWeight <- function(value, quant){
  if(is.na(value)){
    return(NA)
  }
  return(which((quant - value)>0)[1]*0.001-0.001)
}

drugs <- read.csv("./data/external/drug links.csv")
drugs <- subset(drugs, Drug.Type == "SmallMoleculeDrug")

# Step 1: Read files

sc.m <- as.matrix(read.csv("./data/chem_similarity.csv", row.names = 1))
diag(sc.m) = 0
# Step2: produce quantiles corresponding to probabilities 0 to 1 with stepsize 0.001
sc.q <- quantile(sc.m,seq(from = 0, to = 1, by = 0.001), na.rm = TRUE)

# Step 3: Estimate weights as the score corresponding quantile  
sc.w <- pbsapply(sc.m, function(x) sapply(x, getWeight, sc.q))

# Step 4: Replace zero similarity (or similarities < epsilon) with NAs

# Step 5: Compute p-values
z = (w-mean(w, na.rm = TRUE))/sd(w, na.rm = TRUE) # check for hist(z)
p = pnorm(-z)

#Adjust FDR

# for aggregation apply w = apply(x, 1, function(x) mean(x, na.rm=TRUE)) and then compute z-scores

