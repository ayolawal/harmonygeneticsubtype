## extract KMT2A::AFF1 principal genetic abnormality
#'
#' Function to extract KMT2A::AFF1 principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#' @param visit_occurrence The visit_occurrence data frame from the OMOP table
#'
#' @return A data frame containing the KMT2A::AFF1 principal genetic abnormality
#' @export
#'

extract_kmt2a_aff1 <- function(person, measurement, visit_occurrence) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame") | !is(visit_occurrence, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## extract karyotype
  kar_refined <- extract_karyotype(measurement)

  ## extract karyotype status
  karyotype_status <- gen_kary_status(person, measurement)

  Visit_Date <- visit_occurrence %>% rename(VSD = "visit_start_date", VOI = "visit_occurrence_id", VSV = "visit_source_value") %>%
    filter(grepl("2000100001", visit_source_concept_id)) %>% select(person_id, VOI, VSD) %>%  arrange(person_id) %>%
    mutate(VSD = as.Date(ymd_hms(VSD)))

  kmt2a_aff1_kar <- kar_refined %>%
    mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    filter(grepl("4.11|11.4", karyotype) & grepl("q21.q23|q23.q21|isht\\(4.11\\)", karyotype))

  df3 <- measurement %>% filter(grepl("3012425", measurement_concept_id) | grepl("37030733|2000109034", measurement_source_concept_id)) %>%
    select(person_id, value_as_concept_id, visit_occurrence_id) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    arrange(person_id, desc(value_as_concept_id)) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135, FALSE, TRUE)) %>%
    mutate(KMT2A_AFF1 = ifelse(is.na(value_as_concept_id), "Unknown", ifelse(value_as_concept_id == TRUE, "Present", "Absent"))) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(KMT2A_AFF1 = ifelse(KMT2A_AFF1 == "Unknown" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", KMT2A_AFF1)) %>%
    merge(Visit_Date, by = "person_id", all.x = T) %>%
    filter(!(KMT2A_AFF1 == "Present" & visit_occurrence_id != VOI)) %>%
    filter(!duplicated(person_id)) %>%
    mutate(KMT2A_AFF1 = ifelse(person_id  %in% kmt2a_aff1_kar$person_id, "Present", KMT2A_AFF1)) %>%
    select(person_id, KMT2A_AFF1)

  return(df3)
}
