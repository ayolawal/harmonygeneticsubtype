## extract KMT2A::MLLT1 principal genetic abnormality
#'
#' Function to extract KMT2A::MLLT1 principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#'
#' @return A data frame containing the KMT2A::MLLT1 principal genetic abnormality
#' @export
#'

extract_kmt2a_mllt1 <- function(person, measurement) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## extract karyotype
  kar_refined <- extract_karyotype(measurement)

  ## extract karyotype status
  karyotype_status <- gen_kary_status(person, measurement)

  kmt2a_mllt1_kar <- kar_refined %>%
    mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    #mutate(karyotype = gsub("((\\.).{0,5}[Ii][Ss][Hh]|\\.([Aa][Rr][Rr])|\\.([Mm][Oo][Ll])|(\\]\\.)).*$", "", karyotype)) %>%
    filter(grepl("11.19|19.11", karyotype) & grepl("q23.?p13", karyotype))

  df4 <- measurement %>% filter(grepl("3002279", measurement_concept_id) | grepl("37021664", measurement_source_concept_id)) %>%
    filter(!duplicated(person_id)) %>%
    select(person_id, value_as_concept_id) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135, FALSE, TRUE)) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    mutate(KMT2A_MLLT1 = ifelse(is.na(value_as_concept_id), "Unknown", ifelse(value_as_concept_id == TRUE, "Present", "Absent"))) %>%
    mutate(KMT2A_MLLT1 = ifelse(person_id  %in% kmt2a_mllt1_kar$person_id, "Present", KMT2A_MLLT1)) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(KMT2A_MLLT1 = ifelse(KMT2A_MLLT1 == "Unknown" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", KMT2A_MLLT1)) %>%
    select(-c(kary_status, value_as_concept_id))

  return(df4)
}
