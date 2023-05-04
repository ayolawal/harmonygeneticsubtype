## extract TCF3::PBX1 principal genetic abnormality => t(1;19)(q23;p13)
#'
#' Function to extract TCF3::PBX1 principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#' @param visit_occurrence The visit_occurrence data frame from the OMOP table
#'
#' @return A data frame containing the TCF3::PBX1 principal genetic abnormality
#' @export
#'

extract_tcf3_pbx1 <- function(person, measurement, visit_occurrence) {

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

  tcf3_pbx1_kar <- kar_refined %>%
    filter(grepl("t\\(1.19\\)\\(q23.p13|t\\(1.19\\)", karyotype) & !grepl("11.19|t\\(1.19\\)\\([pq](21|43|25|13|21|1|36)|t\\(1.19\\)\\(q23.(p11|q13)", karyotype))

  df5 <- measurement %>% filter(grepl("3000296", measurement_concept_id) | grepl("2000108068", measurement_source_concept_id)) %>%
    select( person_id, value_as_concept_id, visit_occurrence_id) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135, FALSE, TRUE)) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    mutate(TCF3_PBX1 = ifelse(is.na(value_as_concept_id), "Unknown", ifelse(value_as_concept_id == TRUE, "Present", "Absent"))) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(TCF3_PBX1 = ifelse(TCF3_PBX1 == "Unknown" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", TCF3_PBX1)) %>%
    merge(Visit_Date, by = "person_id", all.x = T) %>%
    filter(!(TCF3_PBX1 == "Present" & visit_occurrence_id != VOI)) %>%
    filter(!duplicated(person_id)) %>%
    mutate(TCF3_PBX1 = ifelse(person_id  %in% tcf3_pbx1_kar$person_id, "Present", TCF3_PBX1)) %>%
    select(person_id, TCF3_PBX1)

  return(df5)
}
