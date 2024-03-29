% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_age_gender_df.R
\name{gen_age_gender_df}
\alias{gen_age_gender_df}
\title{Function to extract age and gender from the OMOP measurement table}
\usage{
gen_age_gender_df(person, concept, measurement)
}
\arguments{
\item{person}{The person data frame from the OMOP table}

\item{concept}{The concept data frame from the OMOP table}

\item{measurement}{The measurement data frame from the OMOP table}
}
\value{
A data frame containing age, gender, diagnosis date, year of birth and age categories
}
\description{
Function to extract age and gender from the OMOP measurement table
}
