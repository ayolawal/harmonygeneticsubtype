## extract_karyotype function extracts all the karyotype at diagnosis into a data frame

#' Function to extract karyotype from the OMOP measurement table
#'
#' @param measurement The measurement data frame from the OMOP table
#'
#' @return The extracted karyotype in a data frame format
#' @export
#'
extract_karyotype <- function(measurement) {
  ## check variable is of class data.frame
  if(!is(measurement, "data.frame")){
    stop("'measurement' must be of class 'data.frame'")
  }
  # Extract all individuals with karyotype filtered by (person_id, karyotype)
  kar <- measurement %>%
    filter( grepl("40765097", measurement_concept_id) & !grepl("Karyotype_after_transplant", measurement_source_value)) %>%
    select(person_id, value_source_value) %>%
    dplyr::rename(karyotype="value_source_value") %>%
    filter(!duplicated(paste(person_id, karyotype)))

  # Extract karyotype filtered by person_id - no duplicated person_id
  kar_1 <- kar %>% filter(!duplicated(person_id))

  # Extract karyotype for duplicated person_id
  kar_2 <- kar %>% filter(duplicated(person_id))

  # Extract karyotype in kary_1 without all person_id with double entry
  # Remove all double entry (FISH data)
  kar_red <- kar_1 %>% filter(person_id %ni% kar_2$person_id)

  # Ensuring only appropriate karyotype are extracted
  kar_refined <- rbind(kar_red, kar_2) %>% arrange(person_id) %>% as.data.frame

  return(kar_refined)
}
