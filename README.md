# Harmony Genetic Subtype Package

Healthcare Alliance for Resourceful Medicines Offensive against Neoplasms in HematologY (HARMONY) is a European network of excellence that captures, integrates, analyzes, and harmonizes big data from high-quality multidisciplinary sources with the main purpose of unlocking valuable knowledge on various important hematologic malignancies (HMs).

HARMONY aims at assembling, connecting and analyzing big data from HM patients to define standard sets of outcome indicators to be measured and used for decision-making by key healthcare systems stakeholders.

HarmonyGeneticSubtype is an R package developed by the Leukaemia Research Cytogenetics Group (LRCG), Faculty of Medical Sciences at Newcastle University, specifically to extract some pre-defined chromosomal abnormalities from data collected from thousands of patients for downstream analyses. The package contains different tools that could be used to distil for instance individuals with certain genetic conditions like down syndrome, Intrachromosomal amplification of chromosome 21 (iAMP21), chromosomal translocations - ETV6-RUNX1, BCR-ABL1 etc from tens of thousands of high-quality anonymised harmonised hematologic malignancies datasets.

The package leverage on the representation of the datasets in Observational Medical Outcomes Partnership (OMOP) common data models (CDM) standardised tables format in the ekberfall database. More information on the OMOP common data models can be found on the [Observational Health Data Sciences and Informatics website](https://www.ohdsi.org/data-standardization/).

The package is continually being updated and the current tools that are available for use on the HARMONY platform are listed below. The input arguments to the tools are the pre-processed person, concept, measurement, and condition_occurrence  data frames obtained after executing the `preprocess_datasets.R` script available [here](https://github.com/ayolawal/harmonygeneticsubtype/blob/main/data-raw/preprocess_datasets.R).


## Tools Description

* `extract_karyotype`: **This function takes the pre-processed "measurement" data frame as input and outputs a data frame containing all the karyotypes at diagnosis for every patient with recorded karyotype on the HARMONY platform. The produced data frame has two columns - person_id and karyotype - and the number of rows equals the number of patients with recorded karyotype at diagnosis.**.
```
extract_karyotype(measurement = measurement)
```
* `gen_age_gender_df`: **The gen_age_gender_df function generates a data frame comprising 4 columns - age, gender, year_of_birth and age categories - for all the patients on the HARMONY platform. It takes the pre-processed person, concept, and measurement data frames as input arguments.**. 

~~~
gen_age_gender_df(person = person, concept = concept, measurement = measurement)
~~~

* `gen_kary_status`: **The gen_kary_status function generates a data frame with just two columns - person_id and kary_status. It essentially separates into 4 groups the karyotype status of all the patients on the HARMONY platform into normal karyotype, failed karyotype, abnormal karyotype, and unknown. The unknown refers to individuals with no recorded karyotype at diagnosis. It takes the pre-processed person and measurement data frames as inputs.**.
~~~
gen_kary_status(person = person, measurement = measurement)
~~~

* `extract_etv6_runx1` : **The extract_etv6_runx1 function extracts patients with ETV6::RUNX1 - t(12;21)(p13;q22.3) - translocation from the ekberfall database on the HARMONY platform, precisely from the measurement table using the concept_ids "3036903" and "35977015" from the "measurement_concept_id" and "measurement_source_concept_id" columns respectively to filter out the chromosomal abnormality. It takes the pre-processed person and measurement data frames as inputs. The output data frame has 2 columns - person_id and ETV6::RUNX1. The ETV6::RUNX1 column has 3 subgroups - Present indicating the presence of ETV6::RUNX1, Absent for absence and Unknown for no information on ETV6::RUNX1 for that individual.**.
~~~
extract_etv6_runx1(person = person, measurement = measurement)
~~~

* `extract_bcr_abl1` : **This function takes the pre-processed person and measurement data frames as input and outputs a two-column data frame comprising of the person_id and BCR::ABL1. The BCR::ABL1 column has Present, Absent, and Unknown to represent presence, absence and no information on the BCR::ABL1 genetic abnormality respectively. A combination of concept_ids for BCR::ABL1 and the presence of t(9;22)(q34;q11) pattern in the karyotype was used to extract the required information. The following concept_ids were used in the function to extract the required information from the measurement table - 3028956, 3011913, 3016815, 3013321, 46235651, and 35977026 from the measurement_concept_id column and 35977026 from the measurement_source_concept_id column of the measurement table respectively.**.
~~~
extract_bcr_abl1(person = person, measurement = measurement)
~~~

* `extract_kmt2a_aff1` : **This function extracts a two-column data frame - person_id and KMT2A::AFF1 columns for every patient on the HARMONY platform. It takes the pre-processed person and measurement data frames as inputs and uses 3012425 as measurement_concept_id and 37030733 and 2000109034 as measurement_source_concept_ids in the measurement table to produce the output. It also used the t(4;11)(q21;q23) pattern in the karyotype in addition to the concept_ids to ensure that every patient with KMT2A::AFF1 translocation is captured.**.
~~~
extract_kmt2a_aff1(person = person, measurement = measurement)
~~~

* `extract_kmt2a_mllt1` : **Like the extract_kmt2a_aff1, this function outputs a two-column data frame with the person_id and KMT2A::MLLT1 as its columns. the KMT2A::MLLT1 column has Present, Absent, and Unknown as its values representing the presence of KMT2A::MLLT1, absence of KMT2A::MLLT1 and no information on the genetic abnormality. It takes person and measurement data frames as input arguments, using 3002279 as measurement_concept_id, 37021664 as measurement_source_concept_id and the presence of t(11;19)(q23;p13.3) pattern in the karyotype to appropriately describe each patient concerning KMT2A::MLLT1 genetic abnormality.**.
~~~
extract_kmt2a_mllt1(person = person, measurement = measurement)
~~~

* `extract_kmt2a_r` : **This function extracts all other genetic abnormalities involving KMT2A other than KMT2A::AFF1 and KMT2A::MLLT1 into a two-column data frame. It takes the pre-processed person and measurement data frames as input arguments and used 3022443, 36017933, 3037009, and 40764964 as measurement_concept_id; and 2000000269 and 2000000263 as measurement_source_concept_id to produce the output data frame. This is in addition to a complex set of patterns used to extract the abnormalities from the karyotype.**.
~~~
extract_kmt2a_r(person = person, measurement = measurement)
~~~

* `extract_tcf3_pbx1` : **The extract_tcf3_pbx1 extracts TCF3::PBX1 chromosomal abnormality into a two-column data frame comprising of the persond_id and   TCF3::PBX1 columns. The combination of t(1;19)(q23;p13) pattern from the karyotype and 3000296 as measurement_concept_id and 2000108068 as measurement_source_concept_id were used to extract this particular genetic abnormality. The function takes the pre-processed person and measurement data frames as input arguments.**.
~~~
extract_tcf3_pbx1(person = person, measurement = measurement)
~~~

* `extract_tcf3_hlf` : **The extract_tcf3_hlf function extracts the rare TCF3::HLF genetic abnormality using only 42868760 as measurement_concept_id in the measurement table to produce a two-column data frame as output. Using t(17;19)(q22;p13) pattern in the karyotype did not produce additional information on the presence or otherwise of the genetic condition on all the patients on the HARMONY platform. The function takes the pre-processed person and measurement data frames as input arguments.**.
~~~
extract_tcf3_hlf(person = person, measurement = measurement)
~~~

* `extract_iamp21` : **The extract_iamp21 function extracts the presence or otherwise of the iAMP21 genetic abnormality for all the patients on the HARMONY platform into a two-column data frame taking the pre-processed person and measurement data frames as input arguments. Specifically, 35977099 as measurement_concept_id and 2000000466 as measurement_source_concept_id were used to extract the relevant information from the measurement table.**.
~~~
extract_iamp21(person = person, measurement = measurement)
~~~

* `extract_heh_hap_hotr` : **The extract_heh_hap_hotr function like most of the functions in the package takes the pre-processed person and measurement data frames as input arguments and outputs a four-column data frame comprising of the persond_id, HeH (High-hyperdiploidy), Hap (Haploidy and Near-haploidy), and HoTr (Low-hypodiploidy and Near-triploidy) columns. The function helps characterise the different aneuploidy subgroups in the HARMONY datasets. The HeH refers to the High-Hyperdiploidy subgroup having between 51-59 and 60-67 modal chromosome numbers and satisfying other pre-defined criteria. The Hap refers to the Haploidy or Near-Haploidy subgroup with between 20-29 modal chromosome numbers. The HoTr represents the Low-Hypodiploidy and Near-Triploidy subgroup with 30-39, 68-78 and 60-67 modal chromosome numbers with some pre-defined conditions. A combination of concept_ids and well-defined patterns in the karyotype was used to characterise each patient in the HARMONY datasets into the HeH, Hap and HoTr subgroups. A table containing all the concept_ids used to extract each of the genetic abnormalities considered in this package is in the data-raw folder and can be accessed by clicking [here](https://github.com/ayolawal/harmonygeneticsubtype/blob/main/data-raw/Concept_id_List.xlsx).**.
~~~
extract_heh_hap_hotr(person = person, measurement = measurement)
~~~

* `extract_complex_karyotype` : **The extract_complex_karyotype function extracts all individuals with five or more chromosomal abnormalities into a two-column data frame, taking the pre-processed person and measurement data frames as inputs. A combination of concept_ids - 2000000153 and 45877994 from the measurement_source_concept_id and the value_as_concept_id columns of the measurement table - were used to extract patients into 3 groups - Present, Absent and Unknown indicating the presence of complex karyotype, absence of complex karyotype and no information on complex karyotype for each patient in the HARMONY datasets. The unknown also refers to individuals with no recorded karyotype at diagnosis. Karyotypes of all the patients with recorded karyotype at diagnosis were examined very closely to identify patients with the presence of 5 or more chromosomal abnormalities and were added to those already captured by the concept_ids to produce a data frame that describes the presence or otherwise of complex karyotype in every patient in the HARMONY datasets.**.
~~~
extract_complex_karyotype(person = person, measurement = measurement)
~~~

* `extract_B_other_B_other_plus_and_T_other` : **The extract_B_other_B_other_plus_and_T_other function incorporates all the 12 previously described chromosomal abnormalities in addition to B_other, B_other_plus, T_other, and Unknown columns to fully describe all the patients in the HARMONY datasets concerning the chromosomal abnormalities of interest. B_other classifies patients with B-cells that have not been classified by the previous 12 chromosomal abnormalities in their order of preference into Present, Absent and Unknown. Present means the presence of B-cells with normal or abnormal karyotype status. Absent represents patients with B-cells that have already been described by at least one of the 12 previously defined genetic abnormalities, while Unknown represents patients with B-cells with no information to adequately classify them into any of the previously defined genetic subtypes. B_other_plus captures these individuals and further distils them into Present, Absent and Unknown. Present in the B_other_plus subgroup means the "Unknown" subgroup in B-other with B-cells that have not been classified into any of the previously defined genetic subtypes. T_other captures into Present subgroup patients with T-cells that have normal or abnormal karyotype status and are not already classified into any of the previously defined genetic subtypes. Absent in the T_other column represents patients with T-cells that have already been classified into at least one of the previously defined genetic subtypes, while the unknown simply represents all the T-cell patients with no adequate information to classify them into any of the genetic subtypes of interest.**.
~~~
extract_B_other_B_other_plus_and_T_other(person = person, measurement = measurement, condition_occurrence = condition_occurrence)
~~~

* `gen_gsubtype_and_riskgroup`: **The gen_gsubtype_and_riskgroup function produces a three-column data frame consisting of 3 columns - person_id, gsubtype and riskgroups. Essentially, the gsubtype classifies every patient into a genetic subtype in a particular order of preference. Though, there are individuals with multiple genetic abnormalities, the abnormality higher in the order of preference is allocated to such patients. The riskgroup groups the patients into mainly good_risk, intermediate_risk, and poor_risk using the relevant genetic subtypes that have been associated with each risk group. The function takes the pre-processed person, measurement, and condition_occurrence data frames as input.**
~~~
gen_gsubtype_and_riskgroup(person = person, measurement = measurement, condition_occurrence = condition_occurrence)
~~~

**The list of the concept ids used to extract all the genetic abnormalities are availabe [here](https://github.com/ayolawal/harmonygeneticsubtype/blob/main/data-raw/Concept_id_List.xlsx).
