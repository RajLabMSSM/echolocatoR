
# PICS
top_sub <- top_SNPs[Gene=="LRRK2"][1]

lead_snp = top_sub$SNP
neg_log_pval =  -log10(top_sub$P)
superpopulation <- "EUR"


library(httr)
read_html()
url <- "https://pubs.broadinstitute.org/pubs/finemapping/pics.php"

fd <- list(command1  = lead_snp,
           command2 = neg_log_pval,
           command3 = superpopulation
)

html <- RCurl::postForm(url, 
                        .params = fd, 
                        style ="POST")

# 
# resp <- POST(url, body=fd, encode="form")
# html <- content(resp)
# html %>% html_nodes("table") %>% html_children()


library(rvest) 
session <- rvest::html_session(url = url) 
form.unfilled <- sess %>% html_node("form") %>% html_form()
form.filled <- form.unfilled %>%
  set_values("command1" = lead_snp,
             "command2" = toString(neg_log_pval) ,
             "command3" = superpopulation)
form.filled$url <- url
results <- submit_form(session, form.filled, submit = "submit") 
results %>% html_nodes("body") %>% 
  html_nodes(".content-area") %>% 
  html_nodes(".row")  


 
