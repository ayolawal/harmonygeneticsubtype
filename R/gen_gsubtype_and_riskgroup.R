
#' Function to extract gsubtype and riskgroup principal genetic abnormality
#'
#' @param person The person data frame from the OMOP table
#' @param measurement The measurement data frame from the OMOP table
#' @param condition_occurrence The condition_occurrence data frame from the OMOP table
#'
#' @return A data frame containing the gsubtype and riskgroup principal genetic abnormality
#' @export
#'

gen_gsubtype_and_riskgroup <- function(person, measurement, condition_occurrence) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame") | !is(condition_occurrence, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  riskgroup_gsubtype <- extract_B_other_B_other_plus_and_T_other(person, measurement, condition_occurrence) %>%
    select(person_id, gsubtype) %>%
    mutate(risk_group = ifelse((gsubtype=="ETV6_RUNX1" | gsubtype=="heh"), "good_risk",
                               ifelse(gsubtype=="TCF3_PBX1", "intermediate_risk",
                                      ifelse((gsubtype=="complex" | gsubtype=="B_other" | gsubtype=="B_other_plus" | gsubtype=="T_other" | gsubtype=="No_data"), gsubtype ,"poor_risk"))))
  return(riskgroup_gsubtype)
}
