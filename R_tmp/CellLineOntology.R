



# Ontology Lookup Service API :
# GiHub: https://github.com/lgatto/rols/
# BioC Tutorial: http://bioconductor.org/packages/release/bioc/vignettes/rols/inst/doc/rols.html
# Cell Ontology: https://www.ebi.ac.uk/ols/ontologies/cl
# CEll Line Ontology: https://www.ebi.ac.uk/ols/ontologies/clo
library(rols)# BiocManager::install("rols")
 
# CLO <- Ontology("clo") 
cell_lines <- c("A549", "GM12878","HMEC","MCF-7","MCF10A-Er-Src") 

CLO.dat <- lapply(1:length(cell_lines), function(i){
  print(cell_lines[i])
  qry <- OlsSearch(q = cell_lines[i], ontology = "CLO") #fieldList = c("label","short_form")
  if(qry@numFound > 0 ){
    (qry <- olsSearch(qry))
    df <- as(qry, "data.frame") 
    df.summ <- data.table::data.table(query=cell_lines[i],
                                      label = df$label[1], 
                                      description = unlist(unique(df$description))[1],
                                      CLO_id = df$short_form[1],
                                      url=df$iri[1])
  }
  return(df.summ) 
}) %>% data.table::rbindlist(fill=T)

# CELLOSAURUS R API
## https://github.com/jimvine/rcellosaurus

# library(dplyr)
# download.file("ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml",
#               "~/Downloads/cellosaurus.xml")
# 

# library(rcellosaurus) # devtools::install_github("jimvine/rcellosaurus")
 
# 
# cell_lines <- xml2::read_xml("~/Downloads/cellosaurus.xml") %>% 
#   xml2::xml_find_all("./cell-line-list/*")
#  
# 
# # One at a time (slow as hell)
# cell.dat <- lapply(1:10, function(i){ #length(cell_lines)
#   # comments <- cell_lines[i] %>% xml2::xml_find_all("comment-list") %>%  xml2::xml_find_all("comment")
#   species <- cell_lines[i] %>% xml2::xml_find_all("species-list") %>%  xml2::xml_text()
#   # hla <- cell_lines[i] %>% xml2::xml_find_all("hla-lists") %>%  xml2::xml_text() 
#   # xref <- cell_lines[i] %>% xml2::xml_find_all("xref-list") %>%  xml2::xml_text()
#   disease <- cell_lines[i] %>% xml2::xml_find_all("xref-list") %>%  xml2::xml_text()
#   
#   
#   accession <- xml2::xml_find_all(cell_lines[i],"accession-list") %>% xml2::xml_text()
#   print(paste("Cell line",i,":",accession))
#   type <- xml2::xml_find_all(cell_lines[i],"accession-list") %>% 
#     xml2::xml_find_all("accession") %>%
#     xml2::xml_attr("type")
#   category <- xml2::xml_attr(cell_lines[i],"category")
#   age <- xml2::xml_attr(cell_lines[i],"age")
#   
#   synonyms <- xml2::xml_find_all(cell_lines[i], "name-list") %>% 
#     xml2::xml_find_all("name") %>% 
#     xml2::xml_text() 
#   id <- synonyms[1]
#   dat <- data.table::data.table(accession=accession,
#                          type=type,
#                          category=category,
#                          age=age,
#                          id=id, 
#                          synonyms=paste(synonyms,collapse="|"))
#   return(dat)
# }) %>% data.table::rbindlist()
# 
# # All at once (far faster)
# accession <- xml2::xml_find_all(cell_lines,"accession-list") %>% xml2::xml_text()
# category <- xml2::xml_attr(cell_lines,"category")
# type <- xml2::xml_find_all(cell_lines,"accession-list") %>% 
#   xml2::xml_find_all("accession") %>%
#   xml2::xml_attr("type")
# age <- xml2::xml_attr(cell_lines,"age")
# 
# dat <- data.table::data.table(accession=accession, 
#                               category=category,
#                               type=type, 
#                               age=age)
#     

# Human Protein Atlas
# hpa <- readxl::read_excel("./echolocatoR/tools/Annotations/CellLineOntology/HumanProteinAtlas.xlsx")
# hpa %>% dplyr::group_by(Name) %>% summarise_each(function(x){paste(x,collapse=" ")})

# Saudi study
# c1 <- readxl::read_excel("echolocatoR/tools/Annotations/CellLineOntology/Cells_Representation_Analysis.V2.xlsx", sheet = "Cell Ontology")
# c2 <- readxl::read_excel("echolocatoR/tools/Annotations/CellLineOntology/Cells_Representation_Analysis.V2.xlsx", sheet = "Cell Line Ontology")


