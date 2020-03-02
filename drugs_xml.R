library("xml2")
library(dplyr)
library(tidyr)
data <- read_xml("./data/external/full database.xml")
drugs <- xml_find_all(data, "//d1:drugbank/d1:drug")

name <- xml_find_first(drugs, ".//d1:name") %>%
    xml_text() %>%
    tbl_df() %>%
    rename("Generic name" = value)
formula <- xml_find_first(drugs, ".//d1:experimental-properties/d1:property/d1:kind[text()='Molecular Formula']/following-sibling::d1:value") %>%
    xml_text() %>%
    tbl_df() %>%
    rename(Formula = value)
molecular_weight <- xml_find_first(drugs, ".//d1:experimental-properties/d1:property/d1:kind[text()='Molecular Weight']/following-sibling::d1:value") %>%
    xml_text() %>%
    tbl_df() %>%
    rename("Molecular weight" = value)
logp1 <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='logP']/following-sibling::d1:value") %>%
    xml_text() %>%
    tbl_df() %>%
    rename("LogP1" = value)
logp1_source <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='logP']/following-sibling::d1:source") %>%
    xml_text() %>%
    tbl_df() %>%
    rename("LogP1_source" = value)

logp2 <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='logP']/../following-sibling::d1:property/d1:kind[text()='logP']/following-sibling::d1:value") %>%
    xml_text() %>%
    tbl_df() %>%
    rename("LogP2" = value)
logp2_source <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='logP']/../following-sibling::d1:property/d1:kind[text()='logP']/following-sibling::d1:source") %>%
    xml_text() %>%
    tbl_df() %>%
    rename("LogP2_source" = value)

h_bond_acc <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='H Bond Acceptor Count']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Hydrogen bond acceptors" = value)


h_bond_don <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='H Bond Donor Count']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Hydrogen bond donors" = value)
#chiral <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Chiral center count" = value)
ring <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='Number of Rings']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Ring count" = value)
polar <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='Polar Surface Area (PSA)']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Polar surface area" = value)
refrac <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='Refractivity']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Molecular Refractivity" = value)
polarizability <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='Polarizability']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Molecular polarizability" = value)
#electrotopological <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Electrotopological states" = value)
rotatable_bonds <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='Rotatable Bond Count']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Rotatable bonds" = value)
physiological_charge <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='Physiological Charge']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Physiological charge" = value)
pka_acid <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='pKa (strongest acidic)']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("pKa (strongest acidic)" = value)
pka_basic <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='pKa (strongest basic)']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("pKa (strongest basic)" = value)
cas <- xml_find_first(drugs, ".//d1:cas-number") %>% xml_text() %>% tbl_df() %>% rename("CAS number" = value)
smiles <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='SMILES']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("Smiles" = value)
inchi <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='InChI']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("InChI" = value)
iupac <- xml_find_first(drugs, ".//d1:calculated-properties/d1:property/d1:kind[text()='IUPAC Name']/following-sibling::d1:value") %>% xml_text() %>% tbl_df() %>% rename("IUPAC name" = value)
descript <- xml_find_first(drugs, ".//d1:description") %>% xml_text() %>% tbl_df() %>% rename("Description" = value)
indic <- xml_find_first(drugs, ".//d1:indication") %>% xml_text() %>% tbl_df() %>% rename("Indication" = value)
mech <- xml_find_first(drugs, ".//d1:mechanism-of-action") %>% xml_text() %>% tbl_df() %>% rename("Mechanism of action" = value)

target <- xml_find_first(drugs, ".//d1:targets/d1:target[not(@*)]/d1:name") %>% xml_text() %>% tbl_df() %>% rename("Target names" = value)
carrier <- xml_find_first(drugs, ".//d1:carriers/d1:carrier[not(@*)]/d1:name") %>% xml_text() %>% tbl_df() %>% rename("Carrier names" = value)
enzyme <- xml_find_first(drugs, ".//d1:enzymes/d1:enzyme[not(@*)]/d1:name") %>% xml_text() %>% tbl_df() %>% rename("Enzyme names" = value)
transporter <- xml_find_first(drugs, ".//d1:transporters/d1:transporter[not(@*)]/d1:name") %>% xml_text() %>% tbl_df() %>% rename("Transporter names" = value)
for (i in 1:400 ) {
  print(i)
    target <- xml_find_first(drugs, paste0(".//d1:targets/d1:target[@position='", i, "']/d1:name")) %>% xml_text() %>% tbl_df() %>% rename("Target names2" = value) %>% bind_cols(target) %>% unite("Target names", c("Target names", "Target names2"), sep=',')
    carrier <- xml_find_first(drugs, paste0(".//d1:carriers/d1:carrier[@position='", i, "']/d1:name")) %>% xml_text() %>% tbl_df() %>% rename("Carrier names2" = value) %>% bind_cols(carrier) %>% unite("Carrier names", c("Carrier names", "Carrier names2"), sep=',')
    enzyme <- xml_find_first(drugs, paste0(".//d1:enzymes/d1:enzyme[@position='", i, "']/d1:name")) %>% xml_text() %>% tbl_df() %>% rename("Enzyme names2" = value) %>% bind_cols(enzyme) %>% unite("Enzyme names", c("Enzyme names", "Enzyme names2"), sep=',')
    transporter <- xml_find_first(drugs, paste0(".//d1:transporters/d1:transporter[@position='", i, "']/d1:name")) %>% xml_text() %>% tbl_df() %>% rename("Transporter names2" = value) %>% bind_cols(transporter) %>% unite("Transporter names", c("Transporter names", "Transporter names2"), sep=',')
}
pharmacodynamics <- xml_find_first(drugs, ".//d1:pharmacodynamics") %>% xml_text() %>% tbl_df() %>% rename("Pharmacodynamics" = value)
toxicity <- xml_find_first(drugs, ".//d1:toxicity") %>% xml_text() %>% tbl_df() %>% rename("Toxicity" = value)
absorption <- xml_find_first(drugs, ".//d1:absorption") %>% xml_text() %>% tbl_df() %>% rename("Absorption" = value)
protein_binding <- xml_find_first(drugs, ".//d1:protein-binding") %>% xml_text() %>% tbl_df() %>% rename("Protein Binding" = value)
met <- xml_find_first(drugs, ".//d1:metabolism") %>% xml_text() %>% tbl_df() %>% rename("Metabolism" = value)
half <- xml_find_first(drugs, ".//d1:half-life") %>% xml_text() %>% tbl_df() %>% rename("Half Life" = value)
route <- xml_find_first(drugs, ".//d1:route-of-elimination") %>% xml_text() %>% tbl_df() %>% rename("Route of Elimination" = value)
vol <- xml_find_first(drugs, ".//d1:volume-of-distribution") %>% xml_text() %>% tbl_df() %>% rename("Volume of Distribution" = value)
clear <- xml_find_first(drugs, ".//d1:clearance") %>% xml_text() %>% tbl_df() %>% rename("Clearance" = value)

tbl <- name %>%
    bind_cols(formula) %>%
    bind_cols(molecular_weight) %>%
    bind_cols(logp1) %>%
    bind_cols(logp1_source) %>%
    bind_cols(logp2) %>%
    bind_cols(logp2_source) %>%
    bind_cols(h_bond_acc) %>%
    bind_cols(h_bond_don) %>%
#    bind_cols(chiral) %>%
    bind_cols(ring) %>%
    bind_cols(polar) %>%
    bind_cols(refrac) %>%
    bind_cols(polarizability) %>%
#    bind_cols(electrotopological) %>%
    bind_cols(rotatable_bonds) %>%
    bind_cols(physiological_charge) %>%
    bind_cols(pka_acid) %>%
    bind_cols(pka_basic) %>%
    bind_cols(cas) %>%
    bind_cols(smiles) %>%
    bind_cols(inchi) %>%
    bind_cols(iupac) %>%
    bind_cols(descript) %>%
    bind_cols(indic) %>%
    bind_cols(mech) %>%
    bind_cols(target) %>%
    bind_cols(carrier) %>%
    bind_cols(enzyme) %>%
    bind_cols(transporter) %>%
    bind_cols(pharmacodynamics) %>%
    bind_cols(toxicity) %>%
    bind_cols(absorption) %>%
    bind_cols(protein_binding) %>%
    bind_cols(met) %>%
    bind_cols(half) %>%
    bind_cols(route) %>%
    bind_cols(vol) %>%
    bind_cols(clear)

tbl