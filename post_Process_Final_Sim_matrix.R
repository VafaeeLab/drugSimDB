library(tidyverse)
library(dplyr)
library(data.table)

sim.df = data.table::fread("Data/new_net_info_V7_pval.csv")
drug.link = data.table::fread("Data/external/drug links.csv") %>%
            dplyr::select(c("DrugBank ID","Name"))

sim.df.drugNames = dplyr::left_join(sim.df,drug.link,by=c("ID1"="DrugBank ID")) %>%
        dplyr::left_join(drug.link, by=c("ID2" = "DrugBank ID")) %>%
        dplyr::select("ID1","Name.x","ID2","Name.y","Chem_similarity","Target_similarity","Pathway_similarity",
                      "GO_CC_Similarity","GO_MF_Similarity","GO_BP_Similarity","rowMeans","p_value","adjP_value")
  

fwrite(sim.df.drugNames, "Data/new_net_info_V7_pval_with_Names.csv")
