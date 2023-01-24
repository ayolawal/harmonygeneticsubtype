## extract complex karyotype principal genetic abnormality
#'
#' Function to extract complex karyotype principal genetic abnormality
#'
#' @param person  The person data frame from the OMOP table
#' @param measurement  The measurement data frame from the OMOP table
#'
#' @return A data frame containing the complex karyotype principal genetic abnormality
#' @export
#'

complex_karyotype <- function(person, measurement) {

  ## check inputs are of class data.frame
  if(!is(person, "data.frame") | !is(measurement, "data.frame")){
    stop("All inputs must be of class 'data.frame'")
  }

  ## extract karyotype
  kar_refined <- extract_karyotype(measurement)

  ## extract karyotype status
  karyotype_status <- kary_status(person, measurement)

  complex_karyotype_list_1 <- c(342107, 193254, 113120, 594017, 358686, 515277, 225499, 108346, 158527, 371283, 514216, 406459, 491390, 520589, 395737, 587652, 470537, 341437, 240082, 197172, 239399, 283993, 535500, 187851, 549909, 434140, 192052, 300959, 499182, 137104, 337319, 130572, 535524, 173485, 388025, 470142, 538393, 520203, 495417, 316276, 191591, 328193, 486553, 582497, 185475, 153594, 235505, 377246, 198479, 340206, 371156, 272787, 117338, 458844, 243722, 568631, 272973, 481629, 404804, 492409, 305020, 358057, 199464, 357585, 107801, 328679, 372896, 576501, 188334, 183303, 187786, 389312, 392490, 547168, 288465, 314605, 522111, 583663, 584925, 395676, 222213, 246624, 350746, 235840, 462030, 242890, 483589, 353122, 241392, 273106, 216530, 567409, 584378, 396034, 183811, 210359, 400910, 150258, 169169, 308865, 477031, 163179, 374160, 258513, 141355, 533815, 447278, 346984, 564575, 339763, 104205, 181288, 109481, 369494, 382662, 415326, 590612, 412841, 335859, 488031, 293591, 558198, 144094, 339296, 250852, 576107, 160847, 150813, 218345, 202304, 543881, 413487, 319863, 355003, 362359, 166921, 252972, 560919, 336583, 106126, 187730, 538330, 473540, 418833, 353417, 195785, 529389, 112743, 115578, 376325, 334928, 583502, 429771, 238375, 410295, 550044, 385437, 481771, 523138, 307155, 348218, 295754)

  complex_karyotype_list_2 <- c("PRV_2_3dd571e0e17f38260328f3e4b99903384eab1a08", "PRV_2_0dee18dd6caa47c9f208c2dc4c37638eecf172f6", "PRV_2_8220204f96c7a737924cbe8704a162a72c6b11c9", "PRV_2_ca3f39d3f2dcd41fa6e03f4bc40eeaa0c897b1b6", "PRV_2_499c8f0fa768f2a336da88d3a3e898354241c126", "PRV_2_1658a84ba26c89b17bcaee115819d99ee559d1fa", "PRV_2_98857c096642aa7bdbd1a0765757a3370d38ddd6", "PRV_2_b422d88a92da3cf6eea75a757d6c7810d6cf2ee1", "PRV_2_9c43b78f60e735a4fb3b3dded1831ae0925a0059", "PRV_2_89cb73d32e627f917e50abb168f25f33816ad43e", "PRV_2_03964a268b79b1e4311afb49320cd8c1d5fe7cbf", "PRV_2_fcb21965b873e85cb07c0469a4239e31f97502aa", "PRV_2_29b1a9f1a602b209a63470d66d510c5772143119", "PRV_2_ef5b918d9de5eeccac14a7a01cf0027993cc2048", "PRV_2_8e2f09ab85ef3614026b3a3c0386c5ab88f54f3e", "PRV_2_f3246e5b902c0674f34cc3ba99988539b59caa9e", "PRV_2_a743781c9764ecf69e1bb695511cde022177af2f", "PRV_2_c8c02ed0df50b29a568a48d52e1ef9a0f6df2cf0", "PRV_2_192c674be6e1688319ba09fbc4da313609ba43d0", "PRV_2_2bf868f0cb4cb46b62da817ff2bd5cd2209d48e5", "PRV_2_8ca248c0ad4a7126a5c521e16eb2599c33a0803a", "PRV_2_c4eeb8e85eace480a2be6a5dfa4f903698a63cbe")

  df12 <- measurement %>% filter((grepl("2000000153", measurement_source_concept_id) & grepl("45877994", value_as_concept_id))) %>%
    select( person_id, value_as_concept_id) %>%
    merge(kar_refined, by = "person_id", all.y=T) %>%
    mutate(value_as_concept_id = ifelse(is.na(value_as_concept_id), FALSE, TRUE)) %>%
    merge((person %>% select(person_id)), by = "person_id", all.y = T) %>%
    mutate(value_as_concept_id = ifelse(is.na(value_as_concept_id), "Not_done", ifelse(value_as_concept_id == TRUE, "Present", "Absent"))) %>%
    mutate(value_as_concept_id = ifelse(person_id %in% c(complex_karyotype_list_1, complex_karyotype_list_2), "Present", value_as_concept_id)) %>%
    mutate("complex_karyotype" = value_as_concept_id) %>%
    merge(karyotype_status, by = "person_id", all=T) %>%
    mutate(complex_karyotype = ifelse(complex_karyotype == "Not_done" & (kary_status=="Normal" | kary_status=="Abnormal"), "Absent", complex_karyotype)) %>%
    select(-c(kary_status, value_as_concept_id, karyotype))

  return(df12)
}
