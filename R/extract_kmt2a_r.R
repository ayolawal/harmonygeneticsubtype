## extract KMT2A_r principal genetic abnormality
#'
#' Function to extract KMT2A_r principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#' @param visit_occurrence The visit_occurrence data frame from the OMOP table
#'
#' @return  A data frame containing the KMT2A_r principal genetic abnormality
#' @export
#'

extract_kmt2a_r <- function(person, measurement, visit_occurrence) {

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

  kmt2a_r_kar <- kar_refined %>%
    mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    filter(grepl("\\(9.11\\)\\(p([1-2][1-3])?.q23|t\\(11.1[0-9]\\)\\(q23|([0-9]|1[0-9]|2[0-2]).11\\)\\(.{0,8}q23|(inv|ins)\\(.{0,3}11\\)\\(.+q23", karyotype)) %>%
    filter(!(grepl("del.11.{0,15}q2.", karyotype) | grepl("4.11|1.11\\)\\(q11", karyotype) | grepl("11.19|9.22\\)\\(q34", karyotype) | grepl("del.11q23", karyotype)))

  df8 <-  measurement %>% filter(grepl("3022443|36017933|3037009|40764964", measurement_concept_id) | grepl("2000000269|2000000263", measurement_source_concept_id)) %>%
    select(person_id, value_as_concept_id, visit_occurrence_id, measurement_date) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    mutate(MD = as.numeric(measurement_date)/86400) %>%
    arrange(person_id, MD, desc(value_as_concept_id)) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135, FALSE, TRUE)) %>%
    mutate(KMT2A_r = ifelse(is.na(value_as_concept_id), "Unknown", ifelse(value_as_concept_id == TRUE, "Present", "Absent"))) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(KMT2A_r = ifelse(KMT2A_r == "Unknown" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", KMT2A_r)) %>%
    merge(Visit_Date, by = "person_id", all.x = T) %>%
    filter(!(KMT2A_r == "Present" & visit_occurrence_id != VOI)) %>%
    filter(!duplicated(person_id)) %>%
    mutate(KMT2A_r = ifelse(person_id  %in% kmt2a_r_kar$person_id, "Present", KMT2A_r)) %>%
    merge(extract_kmt2a_aff1(person, measurement, visit_occurrence), by = "person_id", all.x = T) %>%
    merge(extract_kmt2a_mllt1(person, measurement, visit_occurrence), by = "person_id", all.x = T) %>%
    mutate(KMT2A_r = ifelse(KMT2A_AFF1 =="Present" | KMT2A_MLLT1 == "Present", "Absent", KMT2A_r)) %>%
    select(person_id, KMT2A_r)

  return(df8)
}
