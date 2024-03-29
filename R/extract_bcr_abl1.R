## extract BCR::ABL1 principal genetic abnormality

#' Function to extract BCR::ABL1 principal genetic abnormality
#'
#' @param person The person data frame from the OMOP table
#' @param measurement The measurement data frame from the OMOP table
#' @param visit_occurrence The visit_occurrence data frame from the OMOP table
#'
#' @return A data frame containing the BCR::ABL1 principal genetic abnormality
#' @export
#'

extract_bcr_abl1 <- function(person, measurement, visit_occurrence) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame") | !is(visit_occurrence, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## extract karyotype
  kar_refined <- extract_karyotype(measurement)

  ## Obttain karyotype status
  karyotype_status <- gen_kary_status(person, measurement)

  Visit_Date <- visit_occurrence %>% rename(VSD = "visit_start_date", VOI = "visit_occurrence_id", VSV = "visit_source_value") %>%
    filter(grepl("2000100001", visit_source_concept_id)) %>% select(person_id, VOI, VSD) %>%  arrange(person_id) %>%
    mutate(VSD = as.Date(ymd_hms(VSD)))

  bcr_abl1_kar <- kar_refined %>% mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    filter(grepl("9.22", karyotype) & !grepl("(dic|der)\\(9.22\\)|ABL\\-ve", karyotype))

  df2 <- measurement %>% filter(grepl("3028956|3011913|3016815|3013321|46235651|35977026", measurement_concept_id) | grepl("35977026", measurement_source_concept_id)) %>%
    select(person_id, value_as_concept_id, visit_occurrence_id) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    arrange(person_id, desc(value_as_concept_id)) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135, FALSE, TRUE)) %>%
    mutate(BCR_ABL1 = ifelse(is.na(value_as_concept_id), "Unknown", ifelse(value_as_concept_id == TRUE, "Present", "Absent"))) %>%
    mutate(BCR_ABL1 = ifelse(person_id  %in% bcr_abl1_kar$person_id, "Present", BCR_ABL1)) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(BCR_ABL1 = ifelse(BCR_ABL1 == "Unknown" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", BCR_ABL1)) %>%
    merge(Visit_Date, by = "person_id", all.x = T) %>%
    filter(!(BCR_ABL1 == "Present" & visit_occurrence_id != VOI)) %>%
    filter(!duplicated(person_id)) %>%
    select(person_id, BCR_ABL1)

  return(df2)
}
