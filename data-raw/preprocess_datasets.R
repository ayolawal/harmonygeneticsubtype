# Copy and run this code on the HARMONY platform to preprocess your datasets to be used as inputs arguments in your tools

library(tidyverse)
library(data.table)
library("lubridate")
`%ni%`<- Negate(`%in%`) # init

# ============== load Tables =======================
{
  # tables in OMOP
  tables <- (as.data.frame(SparkR::sql("show tables in ekberfall")))


  # Load all the OMOP tables and save them in a data frame
  ## Example: concept <- as.data.frame(SparkR::sql('select * from ekberfall.concept'))
  names_cols <- data.frame()
  for (i in tables$tableName){
    eval(parse(text=sprintf("%s <- as.data.frame(SparkR::sql('select * from ekberfall.%s')); x=%s",i,i,i)))
    names_cols <- rbind(names_cols, data.frame(table=i, col=names(x)))  # Creates a dataframe with 2 cols containing each OMOP table with their columns
  }

  measurement %>% unique -> measurement

} # code names and load data.frames

# Solving issue with POSIXct and Spark of empty columns
special_cols <- function(x)  all(is.na(x)) && class(x) == "POSIXct"

for (i in unique(tables$tableName)){
  # names_ <- concept %>% select(is.POSIXct) %>% select(Cols_AllMissing(.)) %>% names
  text=sprintf("names_ <- %s %%>%% select_if(special_cols)  %%>%% names",i)
  eval(parse(text=text))
  if (length(names_) > 0)
    # concept[names_] <- as.POSIXct(NA)
    text=sprintf("%s[names_]  <- as.POSIXct(NA)", i)
  eval(parse(text=text))
}

# ============== Missing concepts =======================

# In the table concept some concepts are missing. Thus, we added them manually
aux_ <- read.table(header=T, stringsAsFactors = F, text='
                           concept_id      concept_name
                       3028956   "t(9;22)(q34.1;q11)(ABL1BCR) fusion transcript in Blood or Tissue by Molecular genetics method"
                       32976 "UK Biobank"
                       4265453 "Age"') %>% filter(!duplicated(concept_id))

concept <- bind_rows(concept, aux_ ) %>% unique
