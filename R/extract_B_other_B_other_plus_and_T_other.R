
#' Function to extract B-others and B-other-plus principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#' @param condition_occurrence The condition_occurrence data frame from the OMOP table
#'
#' @return A data frame containing the B-others and B-other-plus principal genetic abnormality
#' @export
#'

extract_B_other_B_other_plus_and_T_other <- function(person, measurement, condition_occurrence) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame") | !is(condition_occurrence, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## Load the 12 principal genetic abnormalities
  df1 <- extract_etv6_runx1(person, measurement)
  df2 <- extract_bcr_abl1(person, measurement)
  df3 <- extract_kmt2a_aff1(person, measurement)
  df4 <- extract_kmt2a_mllt1(person, measurement)
  df5 <- extract_kmt2a_r(person, measurement)
  df6 <- extract_tcf3_pbx1(person, measurement)
  df7 <- extract_tcf3_hlf(person, measurement)
  df8 <- extract_iamp21(person, measurement)
  df9_11 <- extract_heh_hap_hotr(person, measurement)
  df12 <- extract_complex_karyotype(person, measurement)

  ## Construct dataset with the principal genetic abnormalities
  df_gsubtype <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9_11, df12) %>%
    reduce(inner_join, by = "person_id")

  #df_gsubtype$gsubtype <- df_gsubtype$heh

  for (i in seq(nrow(df_gsubtype))) {
    df_gsubtype$gsubtype[i] = colnames(df_gsubtype)[df_gsubtype[i, ]=="Present"][1]
  }

  ## extract karyotype status
  karyotype_status <- gen_kary_status(person, measurement)

  # BCP-ALL # 4082461  56 Precursor B-cell lymphoblastic leukemia
  # 4173963 B-cell acute lymphoblastic leukemia
  df13_16 <- condition_occurrence %>% filter(grepl("4082461|4173963|4082462", condition_concept_id)) %>%
    filter(!duplicated(person_id)) %>%
    select(person_id, condition_concept_id) %>%
    mutate(Immuno = ifelse(condition_concept_id == 4082462, "T_cell", "B_cell")) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    mutate(Immuno = ifelse(is.na(Immuno), "B_cell", Immuno)) %>%
    select(-condition_concept_id) %>%
    merge(df_gsubtype, by = "person_id", all.y = T) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(B_other = ifelse((is.na(gsubtype) & (Immuno == "B_cell") & (kary_status=="Normal" | kary_status=="Abnormal")), "Present",
                            ifelse(!is.na(gsubtype) & (Immuno == "B_cell" | Immuno == "T_cell" ), "Absent", "Unknown"))) %>%
    mutate(gsubtype = ifelse(is.na(gsubtype) & B_other == "Present", "B_other", gsubtype)) %>%
    mutate(B_other_plus = ifelse((Immuno == "B_cell" & B_other == "Unknown" & is.na(gsubtype)), "Present",
                                 ifelse(!is.na(gsubtype) & (Immuno == "B_cell" | Immuno == "T_cell" ), "Absent", "Unknown"))) %>%
    mutate(gsubtype = ifelse(is.na(gsubtype) & B_other_plus == "Present", "B_other_plus", gsubtype)) %>%
    mutate(gsubtype = ifelse(gsubtype == "complex_karyotype", "complex", gsubtype)) %>%
    mutate(T_other = ifelse((is.na(gsubtype) & (Immuno == "T_cell") & (kary_status=="Normal" | kary_status=="Abnormal")), "Present",
                            ifelse(!is.na(gsubtype) & (Immuno == "T_cell" | Immuno == "B_cell"), "Absent", "Unknown"))) %>%
    mutate(gsubtype = ifelse(is.na(gsubtype) & T_other == "Present", "T_other", gsubtype)) %>%
    mutate("Unknown" = ifelse(is.na(gsubtype), "Present", "Absent")) %>%
    mutate(gsubtype = ifelse(is.na(gsubtype), "Unknown", gsubtype)) %>%
    select(person_id, ETV6_RUNX1, BCR_ABL1, KMT2A_AFF1, KMT2A_MLLT1, KMT2A_r, TCF3_PBX1, TCF3_HLF, iAMP21, heh, hap, hotr, complex_karyotype, B_other, B_other_plus, T_other, Unknown, Immuno, gsubtype)

  return(df13_16)
}
