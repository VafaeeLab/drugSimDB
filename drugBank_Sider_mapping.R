# load Drug info
library(data.table)
library(dplyr)

drug.dat = fread("Data/external/drug links.csv", sep = "")
# meddra.all.indication.dat = fread("Data/external/Sider/meddra_all_indications.tsv", header = F)
# colnames(meddra.all.indication.dat) = c("CID", "SE.ID1", "MENTION", "SE.NAME1","LABEL","SE.ID2","SE.NAME2")
drug.name.sider = fread("Data/external/Sider/drug_names.tsv", header=F)
colnames(drug.name.sider) = c("CID","DRUG.NAME")

temp = dplyr::inner_join(drug.name.sider,meddra.all.indication.dat,"CID")
temp2 = dplyr::inner_join(temp,drug.dat, by = c("DRUG.NAME" = "Name"))
fwrite(temp2, file="Data/external/drugBank_Sider_mapping.csv")
