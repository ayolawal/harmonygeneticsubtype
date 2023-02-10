## extract kary_status as normal, failed, abnormal, and not_done

#' Function to extract karyotype status from the OMOP measurement table
#'
#' @param person The person data frame from the OMOP table
#' @param measurement The measurement data frame from the OMOP table
#'
#' @return A data frame containing the karyotype status
#' @export
#'
#'

gen_kary_status <- function(person, measurement) {
  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## extract karyotype
  kar_refined <- extract_karyotype(measurement)

  ## Normal karyotype
  kary_normal <- kar_refined %>% mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    mutate(karyotype = gsub("\\|", ",", karyotype)) %>%
    mutate(count_comma = str_count(karyotype, ",")) %>%
    filter(grepl("[Ii][Ss][Hh]", karyotype) | count_comma <= 1 | grepl("^(4[6-7].X[XY].?\\+?(21|X|Y)[Cc])", karyotype)) %>%
    filter(grepl("^(4[5-6].{0,2}[Xx][XxYy])", karyotype) | grepl("^(4[6-7].X[XY].?\\+?(21|X|Y)[Cc])", karyotype) & !grepl("[Ff][Aa][Ii][Ll]", karyotype) ) %>%
    filter(!grepl("^(4[5-6].{0,2}[Xx][XxYy].{0,2}\\(.+\\))", karyotype)) %>%
    filter(!grepl("^(4[5-6].{0,2}[Xx][XxYy]\\[([0-9]|1[0-9])\\])", karyotype)) %>%
    filter(!grepl("^(4[5-6].{0,2}[Xx][XxYy][XxYy]?.{0,2}(de[lr]|dup|add|dup|inv|ins|dic|idic|inc|idem|abn|(\\+|\\-)))", karyotype) | grepl("^(4[6-7].X[XY].?\\+?(21|X|Y)[Cc])", karyotype)) %>%
    filter(!grepl("^(4[5-6][Xx][XxYy][XY]?.{0,4}(\\+|\\-))", karyotype) ) %>%
    filter(!grepl("^(4[6-7].[Xx][XxYy].?\\+?(21|X|Y)[Cc](.inc)?\\[([0-9]|1[0-9])\\])|46.X[XY][XY][Cc].add", karyotype) ) %>%  ### DOWN SYNDROME 50 Individuals
    mutate(kary_status = "Normal") %>%
    select(person_id, kary_status)

  ## Failed karyotype
  kary_fail <- kar_refined %>% mutate(karyotype = gsub("\\s|\\?", "", karyotype)) %>%
    mutate(karyotype = gsub("\\|", ",", karyotype)) %>%
    mutate(count_comma = str_count(karyotype, ",")) %>%
    filter(grepl("[Ii][Ss][Hh]|inc", karyotype) | count_comma <= 1 | grepl("^(4[6-7].X[XY].?\\+([0-9]|1[0-9]|2[0-2]|X|Y)c)", karyotype)) %>%
    filter(grepl("^(46.{0,2}[Xx][XxYy]\\[([0-9]|1[0-9])\\])|^(46.{0,2}[Xx][XxYy].inc)|^(4[6-7].[Xx][XxYy].?\\+([0-9]|1[0-9]|2[0-2]|X|Y)c(.inc)?\\[([0-9]|1[0-9])\\])", karyotype) | grepl("[Ff][Aa][Ii][Ll]", karyotype)) %>%
    filter(!grepl("^(46.{0,2}X[XY].inc.{0,5}(4[0-9]|5[0-9]))", karyotype)) %>%
    filter(!grepl("^(46.{0,2}X[XY].{2,10}inc)", karyotype)) %>%
    filter(!(grepl("del", karyotype) & !grepl("[Ii][Ss][Hh]", karyotype))) %>%
    filter(!(grepl("^(4[6-7].X[XY].?\\+(21|X|Y)c)", karyotype) & grepl("4[6-7].X[XY].?\\-21", karyotype)) ) %>%
    merge((kar_refined %>% filter(person_id %in% c(110563, 246122))), by = "person_id", all = T) %>%     ### HARD-CODED "110563, 246122" TO THE PERSON_ID
    mutate(kary_status = "Failed") %>%
    select(person_id, kary_status)

  ## Abnormal karyotype
  kary_abnormal <- kar_refined %>% filter(person_id %ni% c(kary_normal$person_id, kary_fail$person_id)) %>%
    filter(!grepl("Not|Kphi", karyotype)) %>%
    mutate(kary_status = "Abnormal") %>%
    select(person_id, kary_status)

  ## Karyotype not done
  kary_none <- person %>% filter(person_id %ni% c(kary_normal$person_id, kary_fail$person_id, kary_abnormal$person_id)) %>%
    select(person_id) %>%
    mutate(kary_status = "Unknown")

  kary_status <- rbind(kary_normal, kary_fail, kary_abnormal, kary_none) %>% arrange(person_id)

  return(kary_status)
}
