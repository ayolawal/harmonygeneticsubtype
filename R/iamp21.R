## extract IAMP21 principal genetic abnormality
#'
#' Function to extract IAMP21 principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#'
#' @return  A data frame containing the IAMP21 principal genetic abnormality
#' @export
#'

iamp21 <- function(person, measurement) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## extract karyotype status
  karyotype_status <- kary_status(person, measurement)

  df8 <- measurement %>% filter(grepl("35977099", measurement_concept_id) | grepl("2000000466", measurement_source_concept_id)) %>%
    filter(!duplicated(person_id)) %>%
    select( person_id, value_as_concept_id) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    mutate("iAMP21" = ifelse(is.na(value_as_concept_id), "Not_done", ifelse(value_as_concept_id == 4181412, "Present", "Absent"))) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(iAMP21 = ifelse(iAMP21 == "Not_done" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", iAMP21)) %>%
    select(-c(kary_status, value_as_concept_id))

  return(df8)
}
