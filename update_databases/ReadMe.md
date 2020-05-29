# Updating similarity measures based on the latest version of KEGG and Drugbank databases
 
## Updating databases:
Run the following cURL commands to download latest versions of required DrugBank databases (you need to make a free account) and save them into ./data/external:

### For Structureal similarities:
Downlod it from: https://www.drugbank.ca/releases/5-1-3/downloads/all-structures 

### For Target-based protein sequence similarities:

- Step #1: find protein similarity first
	- step #1.1: "multiple_jobs_prot.sh" to invoke batches of proteins to find similarities with others
	- step #1.2: "parallel_target_Sim_Raijin.sh" is invoked by the previous step, which runs a batch, which
			eventually invokes a R code "parallel_target_Sim_Raijin.R" > the output is a chunk of protein
			similarity matrix
	- step #1.3: ( Optionally: download all those chunks into local machine and )
			Run "combine_proteinSimFiles_and_prepare_4_targetSim.R" file, which yeilds "protein_similarity.csv"

- Step #2: Find Target similarity
	- step #2.1: "multiple_jobs_target.sh" to invoke batches of drugs to find similairties with others
	- step #2.2: "newTargetSim.sh" is invoked by the previous step, which runs a batch, which
			eventually invokes a R code "newTargetSim.R" > output a chunk of drug similarity matrix
	- step #2.3:  Optionally: download all those chunks into local machine and )
			Run "combine_proteinSimFiles_and_prepare_4_targetSim.R", OR "combineFiles.R" file, 
				which would yeilds "target_similarity.csv"

### For Target-based pathway similarities:
- Steps #1: Run "update_KEGGpathways_for_drugs_cov.R" for getting KEGG pathways annotated as targetted by the drugs
- Steps #2: Run "update_pathway_similarity_v2.sh" (for "update_pathway_similarity_v2.R") in HPC that finds target pathway 
	similarities of all the drugs
