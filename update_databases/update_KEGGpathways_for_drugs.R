update_KEGGpathways_for_drugs <- function(){
  
  require(KEGGREST)
  DK_mappingFile = "./input/external/drug links.csv"

  # Load KEGG IDs (Drugs) ---------------------------------------------------- |
  
  drugs <- subset(read.csv("./data/external/drug links.csv"), Drug.Type == "SmallMoleculeDrug")
  DK_ids <- drugs[,c(1,6)]
  
  DK_ids = DK_ids[which(DK_ids[,2]!=""),]   # filter out empty KEGG drug IDs
  DK_ids[,2] = paste0("dr:",DK_ids[,2])   # make KEGG compatible KEGG drug IDs
  colnames(DK_ids) = c("Drug.ID","Kegg.Drug.ID")
  print(paste("total KEGG drug entries to be retrieved:",length(DK_ids[,2])))
  
  # For showing progress ----------------------------------------------------- |
  env=environment()
  counter=0
  pb <- txtProgressBar(min = 0, max = length(DK_ids[,2]), style = 3)
  
  # Retrieve info from KEGG using REST APIs ---------------------------------- |
  kegg.info <- lapply(DK_ids[,2], function(kegg_drug_id){
    ent = keggGet(kegg_drug_id)
    
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir= env)
    setTxtProgressBar(get("pb", envir= env),
                      curVal +1)
    if(!is.atomic(ent[[1]]$TARGET)){
      if(!is.null(ent[[1]]$TARGET$PATHWAY)){
        tdf = as.data.frame(do.call(rbind,strsplit(ent[[1]]$TARGET$PATHWAY, "  ")))
        tdf[,1] = gsub("\\(.*","",tdf[,1])  # remove anything after "(", e.g. hsa0001(1244+343+422)
        tdf[,1] = gsub("hsa","ko",tdf[,1])  # replace 'hsa' with 'ko', as required in downstream analyses
        ret = cbind(kegg_drug_id, tdf)
      }
    }
  })
  
  # Convert output into DataFrame and dump into file ------------------------- |
  kegg.info = as.data.frame(do.call(rbind,kegg.info))
  colnames(kegg.info) = c("Kegg.Drug.ID","Kegg.ID","Kegg.Def")
  
  
  # join 
  kegg.info[,1] = gsub("dr:", "", kegg.info[,1])
  DK_ids[,2] = gsub("dr:", "", DK_ids[,2])
  fo = merge(DK_ids,kegg.info, by="Kegg.Drug.ID")
  
  write.csv(fo[c("Drug.ID","Kegg.ID")],file = "./data/durg_KEGGpathways.csv",quote = FALSE,row.names = FALSE)
  print("DONE !")
  
  # End ---------------------------------------------------------------------- |
}