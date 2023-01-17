## age_gender_df function extracts a data frame containing age, gender, year_of_birth and age categories

#' Function to extract age and gender from the OMOP measurement table
#'
#' @param person The person data frame from the OMOP table
#' @param concept The concept data frame from the OMOP table
#' @param measurement The measurement data frame from the OMOP table
#'
#' @return A data frame containing age, gender, year of birth and age categories
#' @export
#'
age_gender_df <- function(person, concept, measurement){
  require(lubridate)
  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(concept, "data.frame") | !is(measurement, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }
  ## extract genderand year of birth
  gender <- person %>% select(person_id, gender_concept_id, year_of_birth) %>%
    merge(concept %>% select(concept_id, concept_name), by.x="gender_concept_id", by.y="concept_id", all.x=T) %>%
    select(-gender_concept_id) %>%
    dplyr::rename(gender="concept_name") %>%
    arrange(person_id)

  ## extract age, merge gender and add age_cat
  age_gender <- measurement %>% filter(grepl("3007016", measurement_concept_id)) %>%
    arrange(paste(person_id, measurement_concept_id, measurement_date)) %>%
    filter(!duplicated(paste(person_id, measurement_concept_id)) )  %>%
    mutate(age = value_as_number) %>%
    merge(gender, by = "person_id", all.x = T) %>%
    select(person_id, age, gender, year_of_birth) %>%
    mutate(age_cat = cut(age, breaks=c(-Inf, 5, 10, 15, 25, 40, Inf), labels=c("0-4", "5-9", "10-14", "15-24", "25-39", "40+"), right=FALSE))

  return(age_gender)
}
