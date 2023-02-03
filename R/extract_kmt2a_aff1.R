## extract KMT2A::AFF1 principal genetic abnormality
#'
#' Function to extract KMT2A::AFF1 principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#'
#' @return A data frame containing the KMT2A::AFF1 principal genetic abnormality
#' @export
#'

extract_kmt2a_aff1 <- function(person, measurement) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## extract karyotype
  kar_refined <- extract_karyotype(measurement)

  ## extract karyotype status
  karyotype_status <- gen_kary_status(person, measurement)

  kmt2a_aff1_kar <- kar_refined %>%
    mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    #mutate(karyotype = gsub("((\\.).{0,5}[Ii][Ss][Hh]|\\.([Aa][Rr][Rr])|\\.([Mm][Oo][Ll])|(\\]\\.)).*$", "", karyotype)) %>%
    filter(grepl("4.11|11.4", karyotype) & grepl("q21.q23|q23.q21|isht\\(4.11\\)", karyotype))

  df3 <- measurement %>% filter(grepl("3012425", measurement_concept_id) | grepl("37030733|2000109034", measurement_source_concept_id)) %>%
    select(person_id, value_as_concept_id) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    arrange(person_id, desc(value_as_concept_id)) %>%             ## SELECT ANY POSITIVE INSTANCE OF THE ABNORMALITY
    filter(!duplicated(person_id)) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135, FALSE, TRUE)) %>%
    mutate(KMT2A_AFF1 = ifelse(is.na(value_as_concept_id), "Not_done", ifelse(value_as_concept_id == TRUE, "Present", "Absent"))) %>%
    mutate(KMT2A_AFF1 = ifelse(person_id  %in% kmt2a_aff1_kar$person_id, "Present", KMT2A_AFF1)) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(KMT2A_AFF1 = ifelse(KMT2A_AFF1 == "Not_done" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", KMT2A_AFF1)) %>%
    select(-c(kary_status, value_as_concept_id))

  return(df3)
}
