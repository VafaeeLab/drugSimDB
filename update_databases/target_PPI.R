library(data.table)
library(dplyr)

# read Human PPI from i2d database
i2d.db = read.delim("./data/external/i2d.2_9.Public.HUMAN.tab", header = TRUE, sep = "\t")
# length(unique(i2d.db$SwissProt2))
# [1] 17394
# length(unique(i2d.db$SwissProt1))
# [1] 14939

# read target protein lists for each Drugbank drugs
uni.db = read.delim("./data/external/uniprot links.csv", header = TRUE, sep = ",")
# length(unique(uni.db$UniProt.ID))
# [1] 4763

# check how many target proteins have PPI info in the i2d.db
# length(unique(i2d.db$SwissProt1) %in% unique(uni.db$UniProt.ID))
# [1] 14939
# length(unique(i2d.db$SwissProt2) %in% unique(uni.db$UniProt.ID))
# [1] 17394

