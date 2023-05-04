## tcf3_hlf => t(17;19)(q22;p13)
## extract TCF3::HLF principal genetic abnormality
#'
#' Function to extract TCF3::HLF principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#' @param visit_occurrence The visit_occurrence data frame from the OMOP table
#'
#' @return  A data frame containing the TCF3::HLF principal genetic abnormality
#' @export
#'

extract_tcf3_hlf <- function(person, measurement, visit_occurrence) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame") | !is(visit_occurrence, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  Visit_Date <- visit_occurrence %>% rename(VSD = "visit_start_date", VOI = "visit_occurrence_id", VSV = "visit_source_value") %>%
    filter(grepl("2000100001", visit_source_concept_id)) %>% select(person_id, VOI, VSD) %>%  arrange(person_id) %>%
    mutate(VSD = as.Date(ymd_hms(VSD)))

  ## extract karyotype status
  karyotype_status <- gen_kary_status(person, measurement)

  df6 <- measurement %>% filter(grepl("42868760", measurement_concept_id)) %>%
    select( person_id, value_as_concept_id, visit_occurrence_id) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    mutate(TCF3_HLF = ifelse(is.na(value_as_concept_id), "Unknown", ifelse(value_as_concept_id == 4181412, "Present", "Absent"))) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(TCF3_HLF = ifelse(TCF3_HLF == "Unknown" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", TCF3_HLF)) %>%
    merge(Visit_Date, by = "person_id", all.x = T) %>%
    filter(!(TCF3_HLF == "Present" & visit_occurrence_id != VOI)) %>%
    filter(!duplicated(person_id)) %>%
    select(person_id, TCF3_HLF)

  return(df6)
}
