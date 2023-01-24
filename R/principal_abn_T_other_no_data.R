
#' Function to extract T-other and No data principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#' @param condition_occurrence The condition_occurrence data frame from the OMOP table
#'
#' @return A data frame containing the T-other and No data principal genetic abnormality
#' @export
#'

principal_abn_T_other_no_data <- function(person, measurement, condition_occurrence) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame") | !is(condition_occurrence, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## Load the principal genetic abnormalities, B_other and B_other_plus
  df13_14 <- principal_abn_B_other_B_other_plus(person, measurement, condition_occurrence)

  df15_16 <- condition_occurrence %>% filter(grepl("4082462", condition_concept_id)) %>%
    filter(!duplicated(person_id)) %>%
    select(person_id) %>%
    mutate(value = person_id) %>%
    merge(df13_14, by = "person_id", all.y = T) %>%
    mutate(Immuno = ifelse(!is.na(value), "T_cell", Immuno)) %>% ## ADD IMMUNOPHENOTYPE
    mutate("T_other" = ifelse(is.na(value), "Absent", "Present")) %>%
    mutate(T_other = ifelse(is.na(gsubtype) & (T_other == "Present"), "Present", "Absent")) %>%
    mutate(gsubtype = ifelse(is.na(gsubtype) & (T_other == "Present"), "T_other", gsubtype)) %>%
    mutate("No_data" = ifelse(is.na(gsubtype), "Present", "Absent")) %>%
    mutate(gsubtype = ifelse(is.na(gsubtype), "No_data", gsubtype)) %>%
    select(person_id, ETV6_RUNX1, BCR_ABL1, KMT2A_AFF1, KMT2A_MLLT1, KMT2A_r, TCF3_PBX1, TCF3_HLF, iAMP21, heh, hap, hotr, complex_karyotype, B_other, B_other_plus, T_other, No_data, Immuno, gsubtype)

  return(df15_16)
}
