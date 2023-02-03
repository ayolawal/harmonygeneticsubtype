## extract ETV6::RUNX1 principal genetic abnormality => t(12;21)(p13;q22.3)(ETV6,RUNX1)

#' Function to extract the ETV6::RUNX1 principal genetic abnormality
#'
#' @param person The person data frame from the OMOP table
#' @param measurement The measurement data frame from the OMOP table
#'
#' @return A data frame containing the ETV6::RUNX1 principal genetic abnormality
#' @export
#'

extract_etv6_runx1 <- function(person, measurement) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  df1 <- measurement %>% filter(grepl("3036903", measurement_concept_id) | grepl("35977015", measurement_source_concept_id ) ) %>%
    select(person_id, value_as_concept_id) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    arrange(person_id, desc(value_as_concept_id)) %>%             ## SELECT ANY POSITIVE INSTANCE OF THE ABNORMALITY
    filter(!duplicated(person_id)) %>%
    mutate(ETV6_RUNX1 = ifelse(is.na(value_as_concept_id), "Not_done", ifelse(value_as_concept_id == 4181412, "Present", "Absent"))) %>%
    select(-c(value_as_concept_id))

  return(df1)
}
