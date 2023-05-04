## extract ETV6::RUNX1 principal genetic abnormality => t(12;21)(p13;q22.3)(ETV6,RUNX1)

#' Function to extract the ETV6::RUNX1 principal genetic abnormality
#'
#' @param person The person data frame from the OMOP table
#' @param measurement The measurement data frame from the OMOP table
#' @param visit_occurrence The visit_occurrence data frame from the OMOP table
#'
#' @return A data frame containing the ETV6::RUNX1 principal genetic abnormality
#' @export
#'

extract_etv6_runx1 <- function(person, measurement, visit_occurrence) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame") | !is(visit_occurrence, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  Visit_Date <- visit_occurrence %>% rename(VSD = "visit_start_date", VOI = "visit_occurrence_id", VSV = "visit_source_value") %>%
    filter(grepl("2000100001", visit_source_concept_id)) %>% select(person_id, VOI, VSD) %>%  arrange(person_id) %>%
    mutate(VSD = as.Date(ymd_hms(VSD)))

  df1 <- measurement %>% filter(grepl("3036903", measurement_concept_id) | grepl("35977015", measurement_source_concept_id ) ) %>%
    select(person_id, value_as_concept_id, visit_occurrence_id) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    arrange(person_id, desc(value_as_concept_id)) %>%
    mutate(ETV6_RUNX1 = ifelse(is.na(value_as_concept_id), "Unknown", ifelse(value_as_concept_id == 4181412, "Present", "Absent"))) %>%
    merge(Visit_Date, by = "person_id", all.x = T) %>%
    filter(!(ETV6_RUNX1 == "Present" & visit_occurrence_id != VOI)) %>%
    filter(!duplicated(person_id)) %>%
    select(person_id, ETV6_RUNX1)

  return(df1)
}
