## extract KMT2A_r principal genetic abnormality
#'
#' Function to extract KMT2A_r principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#'
#' @return  A data frame containing the KMT2A_r principal genetic abnormality
#' @export
#'

kmt2a_r <- function(person, measurement) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## extract karyotype
  kar_refined <- extract_karyotype(measurement)

  ## extract karyotype status
  karyotype_status <- kary_status(person, measurement)

  kmt2a_r_kar <- kar_refined %>%
    mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    mutate(karyotype = gsub("((\\.).{0,5}[Ii][Ss][Hh]|\\.([Aa][Rr][Rr])|\\.([Mm][Oo][Ll])|(\\]\\.)).*$", "", karyotype)) %>%
    filter(grepl("\\(9.11\\)\\(p([1-2][1-3])?.q23|t\\(11.1[0-9]\\)\\(q23|([0-9]|1[0-9]|2[0-2]).11\\)\\(.{0,8}q23|(inv|ins)\\(.{0,3}11\\)\\(.+q23", karyotype)) %>%
    filter(!(grepl("del.11.{0,6}q2.", karyotype) | grepl("4.11|1.11\\)\\(q11", karyotype) | grepl("11.19|9.22\\)\\(q34", karyotype) | grepl("del.11q23", karyotype)))

  ## REVISED AFTER REMOVING CONCEPT ID "3012425" FOR KMT2A::AFF1
  df5 <- measurement %>% filter(grepl("3022443|36017933|3037009|40764964", measurement_concept_id) | grepl("2000000269|2000000263", measurement_source_concept_id)) %>%
    select(person_id, value_as_concept_id, measurement_date) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    mutate(MD = as.numeric(measurement_date)/86400) %>%
    arrange(person_id, MD, desc(value_as_concept_id)) %>%             ## USED MEASUREMENT DATE TO SELECT THE EARLIEST MEASUREMENTS FOR INDIVIDUALS WITH MULTIPLE ENTRIES
    filter(!duplicated(person_id)) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135, FALSE, TRUE)) %>%
    mutate(KMT2A_r = ifelse(is.na(value_as_concept_id), "Not_done", ifelse(value_as_concept_id == TRUE, "Present", "Absent"))) %>%
    mutate(KMT2A_r = ifelse(person_id  %in% kmt2a_r_kar$person_id, "Present", KMT2A_r)) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(KMT2A_r = ifelse(KMT2A_r == "Not_done" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", KMT2A_r)) %>%
    select(-c(kary_status, value_as_concept_id, measurement_date, MD))

  return(df5)
}
