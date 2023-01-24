## extract TCF3::PBX1 principal genetic abnormality => t(1;19)(q23;p13)
#'
#' Function to extract TCF3::PBX1 principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#'
#' @return A data frame containing the TCF3::PBX1 principal genetic abnormality
#' @export
#'

tcf3_pbx1 <- function(person, measurement) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## extract karyotype
  kar_refined <- extract_karyotype(measurement)

  ## extract karyotype status
  karyotype_status <- kary_status(person, measurement)

  tcf3_pbx1_kar <- kar_refined %>%
    mutate(karyotype = gsub("((\\.).{0,5}[Ii][Ss][Hh]|\\.([Aa][Rr][Rr])|\\.([Mm][Oo][Ll])|(\\]\\.)).*$", "", karyotype)) %>%
    filter(grepl("t\\(1.19\\)\\(q23.p13|t\\(1.19\\)", karyotype) & !grepl("11.19|t\\(1.19\\)\\([pq](21|43|25|13|21|1|36)|t\\(1.19\\)\\(q23.p11", karyotype))

  df6 <- measurement %>% filter(grepl("3000296", measurement_concept_id) | grepl("2000108068", measurement_source_concept_id)) %>%
    filter(!duplicated(person_id)) %>%
    select( person_id, value_as_concept_id) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135, FALSE, TRUE)) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    mutate(TCF3_PBX1 = ifelse(is.na(value_as_concept_id), "Not_done", ifelse(value_as_concept_id == TRUE, "Present", "Absent"))) %>%
    mutate(TCF3_PBX1 = ifelse(person_id  %in% tcf3_pbx1_kar$person_id, "Present", TCF3_PBX1)) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(TCF3_PBX1 = ifelse(TCF3_PBX1 == "Not_done" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", TCF3_PBX1)) %>%
    select(-c(kary_status, value_as_concept_id))

  return(df6)
}
