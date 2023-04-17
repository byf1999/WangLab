
\name{KEGG_abundance}
\alias{KEGG_abundance}
\title{extract pathway abundance and barplot}
\usage{
https://www.kegg.jp/kegg/mapper/search.html upload your KO id and paste the web to 'kegg_txy.txt'; 'pathway_list.txt'is a table of all_kegg pathway data

KEGG_abundance(mapper='kegg_txy.txt',targetLev = 'L3',pathway_list = 'pathway_list.txt',plot = T)
}
\description{
extract pathway abundance from mapper result of KEGG web
}
\examples{
library(ggplot2)
library(readr)
#assert you have 'kegg_txy.txt' and 'pathway_list.txt' in your work path.
KEGG_abundance(mapper='kegg_txy.txt',targetLev = 'L3',pathway_list = 'pathway_list.txt',plot = T)
}

