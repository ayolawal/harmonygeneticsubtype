
#' Function to extract B-others and B-other-plus principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#' @param condition_occurrence The condition_occurrence data frame from the OMOP table
#'
#' @return A data frame containing the B-others and B-other-plus principal genetic abnormality
#' @export
#'

principal_abn_B_other_B_other_plus <- function(person, measurement, condition_occurrence) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame") | !is(condition_occurrence, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## Load the 12 principal genetic abnormalities
  df1 <- etv6_runx1(person, measurement)
  df2 <- bcr_abl1(person, measurement)
  df3 <- kmt2a_aff1(person, measurement)
  df4 <- kmt2a_mllt1(person, measurement)
  df5 <- kmt2a_r(person, measurement)
  df6 <- tcf3_pbx1(person, measurement)
  df7 <- tcf3_hlf(person, measurement)
  df8 <- iamp21(person, measurement)
  df9_11 <- heh_hap_hotr(person, measurement)
  df12 <- complex_karyotype(person, measurement)

  ## Construct dataset with the principal genetic abnormalities
  df_gsubtype <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9_11, df12) %>%
    reduce(inner_join, by = "person_id")

  #df_gsubtype$gsubtype <- df_gsubtype$heh

  for (i in seq(nrow(df_gsubtype))) {
    df_gsubtype$gsubtype[i] = colnames(df_gsubtype)[df_gsubtype[i, ]=="Present"][1]
  }

  ## extract karyotype status
  karyotype_status <- kary_status(person, measurement)

  # BCP-ALL # 4082461  56 Precursor B-cell lymphoblastic leukemia
  # 4173963 B-cell acute lymphoblastic leukemia
  df13_14 <- condition_occurrence %>% filter(grepl("4082461|4173963", condition_concept_id)) %>%
    filter(!duplicated(person_id)) %>%
    select(person_id) %>%
    mutate(value = person_id) %>%
    merge(df_gsubtype, by = "person_id", all.y = T) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate("Immuno" = ifelse(is.na(value), "No_data", "B_cell")) %>% ## ADD IMMUNOPHENOTYPE
    mutate("B_other" = ifelse(is.na(value), "Absent", "Present")) %>%
    mutate("B_other_plus" = B_other) %>%
    mutate(B_other = ifelse(is.na(gsubtype) & (B_other == "Present") & (kary_status=="Normal" | kary_status=="Abnormal"), "Present", "Absent")) %>%
    mutate(gsubtype = ifelse(is.na(gsubtype) & (B_other == "Present") & (kary_status=="Normal" | kary_status=="Abnormal"), "B_other", gsubtype)) %>%
    mutate(B_other_plus = ifelse(is.na(gsubtype) & (B_other_plus == "Present") & (kary_status=="Failed" | kary_status=="Not_done"), "Present", "Absent")) %>%
    mutate(gsubtype = ifelse(is.na(gsubtype) & (B_other_plus == "Present") & (kary_status=="Failed" | kary_status=="Not_done"), "B_other_plus", gsubtype)) %>%
    mutate(gsubtype = ifelse(gsubtype == "complex_karyotype", "complex", gsubtype)) %>%
    select(-c(value, kary_status))

  return(df13_14)
}
