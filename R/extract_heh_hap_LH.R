## extract ANEUOPLOIDY genetic abnormality
#'
#' Function to extract ANEUOPLOIDY principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#'
#' @return A data frame containing the ANEUOPLOIDY principal genetic abnormality
#' @export
#'

extract_heh_hap_LH <- function(person, measurement) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## extract karyotype status
  karyotype_status <- gen_kary_status(person, measurement)

  ## extract karyotype
  kar_refined <- extract_karyotype(measurement)

  ## HeH - High hyperdiploidy extracted using concept_id
  heh_cid <- measurement %>% filter(grepl("36660734", measurement_concept_id) | grepl("2000000462", measurement_source_concept_id)) %>%
    select( person_id, value_as_concept_id) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135, "Absent", "Present")) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    arrange( person_id, desc(value_as_concept_id)) %>%
    filter(!duplicated(person_id)) %>%
    mutate(value_as_concept_id = ifelse(is.na(value_as_concept_id), "Unknown", value_as_concept_id)) %>%
    mutate( heh = value_as_concept_id) %>%
    select(-value_as_concept_id)

  heh_cid_red <- heh_cid %>% filter(heh == "Present")

  ## LH - Low hypodiploidy extracted using concept_id
  LH_cid <- measurement %>% filter((grepl("36660734", measurement_concept_id) | grepl("2000000463", measurement_source_concept_id)) & !grepl("high_hyperdiploidy",      measurement_source_value)) %>%
    select( person_id, value_as_concept_id) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135, "Absent", "Present")) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    arrange( person_id, desc(value_as_concept_id)) %>%
    filter(!duplicated(person_id)) %>%
    mutate(value_as_concept_id = ifelse(is.na(value_as_concept_id), "Unknown", value_as_concept_id)) %>%
    mutate( LH = value_as_concept_id) %>%
    select(-value_as_concept_id)

  LH_cid_red <- LH_cid %>% filter(LH == "Present")

  ## Hap - Haploidy/Near Haploidy extracted using concept_id
  hap_cid <- measurement %>% filter(grepl("2000000464", measurement_source_concept_id)) %>%
    select( person_id, value_as_concept_id) %>%
    mutate(value_as_concept_id = ifelse(value_as_concept_id == 4132135,"Absent", "Present")) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    arrange( person_id, desc(value_as_concept_id)) %>%
    filter(!duplicated(person_id)) %>%
    mutate(value_as_concept_id = ifelse(is.na(value_as_concept_id), "Unknown", value_as_concept_id)) %>%
    mutate(hap = value_as_concept_id) %>%
    select(-value_as_concept_id)

  hap_cid_red <- hap_cid %>% filter(hap == "Present")

  heh_kar1 <- kar_refined %>% mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    mutate(karyotype = gsub("\\|", ",", karyotype)) %>%                      # replace pipe (|) with coma
    mutate(kary1 = gsub(",.*$", "", karyotype)) %>%                          # Create another karyotype with everything after the 1st coma removed
    mutate(karyotype = gsub("((\\.).{0,5}[Ii][Ss][Hh]|\\.([Aa][Rr][Rr])|\\.([Mm][Oo][Ll])|(\\]\\.)).*$", "", karyotype)) %>%
    filter((grepl("^(5[1-9])|\\/(5[1-9])", karyotype) & !grepl("\\/(5[1-9])\\]", karyotype)) | grepl("[~-](5[1-9])", kary1)) %>%
    filter(person_id %ni% c(heh_cid_red$person_id, hap_cid_red$person_id, LH_cid_red$person_id)) %>%
    select(-kary1)

  heh_kar2 <- kar_refined %>% mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    mutate(karyotype = gsub("\\|", ",", karyotype)) %>%                      # replace pipe (|) with coma
    mutate(kary1 = gsub(",.*$", "", karyotype)) %>%                          # Create another karyotype with everything after the 1st coma removed
    mutate(karyotype = gsub("((\\.).{0,5}[Ii][Ss][Hh]|\\.([Aa][Rr][Rr])|\\.([Mm][Oo][Ll])|(\\]\\.)).*$", "", karyotype)) %>%
    filter((grepl("^(6[0-7])|\\/(6[0-7])", karyotype) & !grepl("\\/(6[0-7])\\]", karyotype)) | grepl("[~-](6[0-7])", kary1)) %>%
    filter(person_id %ni% c(heh_cid_red$person_id, hap_cid_red$person_id, LH_cid_red$person_id, heh_kar1$person_id)) %>%
    filter(!(grepl("\\+1", karyotype) & !grepl("\\+1[0-9]|\\+1[~-]2mar", karyotype))) %>%
    filter(!(grepl("XX[XY]|3[Nn]", karyotype) & grepl("\\+11|\\+18|\\+19", karyotype))) %>%
    filter(!grepl("(\\+11.{0,5}){2,}|(\\+18.{0,5}){2,}|(\\+19.{0,5}){2,}", karyotype)) %>%
    filter(!grepl("\\)x2", karyotype)) %>%
    select(-kary1)

  df9 <- heh_cid %>%
    mutate(heh = ifelse(person_id %in% c(heh_kar1$person_id, heh_kar2$person_id), "Present", heh)) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(heh = ifelse(heh == "Unknown" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", heh)) %>%
    select(-kary_status)

  hap_kar <- kar_refined %>% mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    mutate(karyotype = gsub("\\|", ",", karyotype)) %>%                      # replace pipe (|) with coma
    mutate(kary1 = gsub(",.*$", "", karyotype)) %>%                          # Create another karyotype with everything after the 1st coma removed
    mutate(karyotype = gsub("((\\.).{0,5}[Ii][Ss][Hh]|\\.([Aa][Rr][Rr])|\\.([Mm][Oo][Ll])|(\\]\\.)).*$", "", karyotype)) %>%
    filter((grepl("^(2[0-9])|\\/(2[0-9])", karyotype) & !grepl("\\/(2[0-9]).{0,2}\\]", karyotype)) | grepl("[~-](2[0-9])", kary1)) %>%
    filter(person_id %ni% c(heh_cid_red$person_id, hap_cid_red$person_id, LH_cid_red$person_id, heh_kar1$person_id, heh_kar2$person_id))  %>%
    select(-kary1)

  df10 <- hap_cid %>%
    mutate(hap = ifelse(person_id %in% hap_kar$person_id, "Present", hap)) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(hap = ifelse(hap == "Unknown" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", hap)) %>%
    select(-kary_status)

  p1 <- "\\+1"
  p2 <- "\\+1[0-9]|\\+1[~-]2mar"
  p3 <- "XX[XY]|3[Nn]"
  p4 <- "\\+11|\\+18|\\+19"
  p5 <- "(\\+11.{0,5}){2,}|(\\+18.{0,5}){2,}|(\\+19.{0,5}){2,}"
  p6 <- "\\)x2"

  LH_kar1 <- kar_refined %>% mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    mutate(karyotype = gsub("\\|", ",", karyotype)) %>%                      # replace pipe (|) with coma
    mutate(kary1 = gsub(",.*$", "", karyotype)) %>%                          # Create another karyotype with everything after the 1st coma removed
    mutate(karyotype = gsub("((\\.).{0,5}[Ii][Ss][Hh]|\\.([Aa][Rr][Rr])|\\.([Mm][Oo][Ll])|(\\]\\.)).*$", "", karyotype)) %>%
    filter((grepl("^(3[0-9])|\\/(3[0-9])", karyotype) & !grepl("\\/(3[0-9])\\]", karyotype)) | grepl("[~-](3[0-9])", kary1)) %>% ##  9
    filter(person_id %ni% c(heh_cid_red$person_id, hap_cid_red$person_id, LH_cid_red$person_id, heh_kar1$person_id, heh_kar2$person_id)) %>%
    select(-kary1)

  LH_kar2 <- kar_refined %>% mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    mutate(karyotype = gsub("\\|", ",", karyotype)) %>%                      # replace pipe (|) with coma
    mutate(kary1 = gsub(",.*$", "", karyotype)) %>%                          # Create another karyotype with everything after the 1st coma removed
    mutate(karyotype = gsub("((\\.).{0,5}[Ii][Ss][Hh]|\\.([Aa][Rr][Rr])|\\.([Mm][Oo][Ll])|(\\]\\.)).*$", "", karyotype)) %>%
    filter((grepl("^(6[8-9]|7[0-8])|\\/(6[8-9]|7[0-8])", karyotype) & !grepl("\\/(6[8-9]|7[0-8])\\]|4[Nn]|[~-]8[0-9]", karyotype)) | grepl("[~-](6[8-9]|7[0-8])", kary1)) %>% ##  18
    filter(person_id %ni% c(heh_cid_red$person_id, hap_cid_red$person_id, LH_cid_red$person_id, heh_kar1$person_id, heh_kar2$person_id, LH_kar1$person_id)) %>%
    select(-kary1)

  LH_kar3 <- kar_refined %>% mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    mutate(karyotype = gsub("\\|", ",", karyotype)) %>%                      # replace pipe (|) with coma
    mutate(kary1 = gsub(",.*$", "", karyotype)) %>%                          # Create another karyotype with everything after the 1st coma removed
    mutate(karyotype = gsub("((\\.).{0,5}[Ii][Ss][Hh]|\\.([Aa][Rr][Rr])|\\.([Mm][Oo][Ll])|(\\]\\.)).*$", "", karyotype)) %>%
    filter((grepl("^(6[0-7])|\\/(6[0-7])", karyotype) & !grepl("\\/(6[0-7])\\]", karyotype)) | grepl("[~-](6[0-7])", kary1)) %>%
    filter((grepl(p1, karyotype) & !grepl(p2, karyotype)) | (grepl(p3, karyotype) & grepl(p4, karyotype)) | grepl(p5, karyotype) | grepl(p6, karyotype)) %>%
    filter(person_id %ni% c(heh_cid_red$person_id, hap_cid_red$person_id, LH_cid_red$person_id, heh_kar1$person_id, heh_kar2$person_id, LH_kar1$person_id, LH_kar2$person_id)) %>%
    select(-kary1)

  df11 <- LH_cid %>%
    mutate(LH = ifelse(person_id %in% c(LH_kar1$person_id, LH_kar2$person_id, LH_kar3$person_id), "Present", LH)) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(LH = ifelse(LH == "Unknown" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", LH)) %>%
    select(-kary_status) %>%
    list(df9, df10) %>%
    reduce(inner_join, by = "person_id") %>%
    select(person_id, heh, hap, LH)

  return(df11)
}
