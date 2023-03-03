- [Harmony Genetic Subtype Package](#harmony-genetic-subtype-package)
- [Tools Description](#tools-description)
  * [Extract Karyotype](#extract-karyotype)
  * [Generate Age and Gender Data Frame](#generate-age-and-gender-data-frame)
  * [Generate Karyotype Status](#generate-karyotype-status)
  * [Extract `ETV6::RUNX1`](#extract--etv6--runx1-)
  * [Extract `BCR::ABL1`](#extract--bcr--abl1-)
  * [Extract KMT2A::AFF1](#extract-kmt2a--aff1)
  * [Extract KMT2A::MLLT1](#extract-kmt2a--mllt1)
  * [Extract KMT2A Others](#extract-kmt2a-others)
  * [Extract TCF3::PBX1](#extract-tcf3--pbx1)
  * [Extract TCF3::HLF](#extract-tcf3--hlf)
  * [Extract iAMP21](#extract-iamp21)
  * [Extract High hyperdiploidy, Near Haploid/Haploidy, and Low hypodiploidy](#extract-high-hyperdiploidy--near-haploid-haploidy--and-low-hypodiploidy)
  * [Extract Complex Karyotype](#extract-complex-karyotype)
  * [Extract B_other, B_other_plus, and T_other Subgroups](#extract-b-other--b-other-plus--and-t-other-subgroups)
  * [List of concept IDs](#list-of-concept-ids)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>



# Harmony Genetic Subtype Package

Healthcare Alliance for Resourceful Medicines Offensive against Neoplasms in HematologY (HARMONY) is a European network of excellence that captures, integrates, analyzes, and harmonizes big data from high-quality multidisciplinary sources with the main purpose of unlocking valuable knowledge on various important hematologic malignancies (HMs).

HARMONY aims at assembling, connecting and analyzing big data from HM patients to define standard sets of outcome indicators to be measured and used for decision-making by key healthcare systems stakeholders.

HarmonyGeneticSubtype is an R package developed by the Leukaemia Research Cytogenetics Group (LRCG), Faculty of Medical Sciences at Newcastle University, specifically to extract some pre-defined chromosomal abnormalities from data collected from thousands of HM patients for downstream analyses. The package contains different tools that could be used to distil for instance individuals with certain genetic conditions like down syndrome, Intrachromosomal amplification of chromosome 21 (iAMP21), chromosomal translocations - ETV6::RUNX1, BCR::ABL1 etc from tens of thousands of high-quality anonymised harmonised hematologic malignancies datasets.

The package leverages the representation of the datasets in Observational Medical Outcomes Partnership (OMOP) common data models (CDM) standardised tables format in the European knowledge bank ... database. More information on the OMOP common data models can be found on the [Observational Health Data Sciences and Informatics website](https://www.ohdsi.org/data-standardization/).

The package is continually being updated and the tools that are available for use on the HARMONY platform are listed below. The input arguments to the tools are the pre-processed person, concept, measurement, and condition_occurrence  data frames obtained after executing the `preprocess_datasets.R` script available [here](https://github.com/ayolawal/harmonygeneticsubtype/blob/main/data-raw/preprocess_datasets.R).


# Tools Description

## Extract Karyotype
* `extract_karyotype`: **This function takes the pre-processed "measurement" data frame as input and outputs a data frame containing all the karyotypes at diagnosis for every patient with recorded karyotype on the HARMONY platform. The produced data frame has two columns - person_id and karyotype - and the number of rows equals the number of patients with recorded karyotype at diagnosis.**
```
extract_karyotype(measurement = measurement)
```
## Generate Age and Gender Data Frame
* `gen_age_gender_df`: **The gen_age_gender_df function generates a data frame comprising 4 columns - age, gender, year_of_birth and age categories - for all the patients on the HARMONY platform. It takes the pre-processed person, concept, and measurement data frames as input arguments.**

~~~
gen_age_gender_df(person = person, concept = concept, measurement = measurement)
~~~

## Generate Karyotype Status
* `gen_kary_status`: **The gen_kary_status function generates a data frame with just two columns - person_id and kary_status. It essentially separates into 4 groups the karyotype status of all the patients on the HARMONY platform into normal karyotype, failed karyotype, abnormal karyotype, and unknown. The unknown refers to individuals with no recorded karyotype at diagnosis. It takes the pre-processed person and measurement data frames as inputs.**
~~~
gen_kary_status(person = person, measurement = measurement)
~~~

## Extract `ETV6::RUNX1`
* `extract_etv6_runx1` : **The extract_etv6_runx1 function extracts patients with ETV6::RUNX1 - t(12;21)(p13;q22.3) - translocation from the database on the HARMONY platform using the relevant [concept_ids](#list-of-concept-ids) to filter out the chromosomal abnormality. It takes the pre-processed person and measurement data frames as inputs. The output data frame has 2 columns - person_id and ETV6::RUNX1. The ETV6::RUNX1 column has 3 subgroups - Present indicating the presence of ETV6::RUNX1, Absent for absence and Unknown for no information on ETV6::RUNX1 for that individual.**
~~~
extract_etv6_runx1(person = person, measurement = measurement)
~~~

## Extract `BCR::ABL1`
* `extract_bcr_abl1` : **This function takes the pre-processed person and measurement data frames as input and outputs a two-column data frame comprising of the person_id and BCR::ABL1. The BCR::ABL1 column has Present, Absent, and Unknown to represent presence, absence and no information on the BCR::ABL1 genetic abnormality respectively. A combination of relevant [concept_ids](#list-of-concept-ids) for BCR::ABL1 and the presence of t(9;22)(q34;q11) pattern in the karyotype was used to extract the required information.**
~~~
extract_bcr_abl1(person = person, measurement = measurement)
~~~

## Extract KMT2A::AFF1
* `extract_kmt2a_aff1` : **This function extracts a two-column data frame - person_id and KMT2A::AFF1 columns for every patient on the HARMONY platform. It takes the pre-processed person and measurement data frames as inputs and uses relevant [concept_ids](#list-of-concept-ids) to produce the output. It also used the t(4;11)(q21;q23) pattern in the karyotype in addition to the concept_ids to ensure that every patient with KMT2A::AFF1 translocation is captured.**
~~~
extract_kmt2a_aff1(person = person, measurement = measurement)
~~~

## Extract KMT2A::MLLT1
* `extract_kmt2a_mllt1` : **Like the extract_kmt2a_aff1, this function outputs a two-column data frame with the person_id and KMT2A::MLLT1 as its columns. the KMT2A::MLLT1 column has Present, Absent, and Unknown as its values representing the presence of KMT2A::MLLT1, absence of KMT2A::MLLT1 and no information on the genetic abnormality. It takes person and measurement data frames as input arguments, relevant [concept_ids](#list-of-concept-ids) and the presence of t(11;19)(q23;p13.3) pattern in the karyotype to appropriately describe each patient concerning KMT2A::MLLT1 genetic abnormality.**
~~~
extract_kmt2a_mllt1(person = person, measurement = measurement)
~~~

## Extract KMT2A Others
* `extract_kmt2a_r` : **This function extracts all other genetic abnormalities involving KMT2A other than KMT2A::AFF1 and KMT2A::MLLT1 into a two-column data frame. It takes the pre-processed person and measurement data frames as input arguments and used relevant [concept_ids](#list-of-concept-ids) to produce the output data frame. This is in addition to a complex set of patterns used to extract the abnormalities from the karyotype.**
~~~
extract_kmt2a_r(person = person, measurement = measurement)
~~~

## Extract TCF3::PBX1
* `extract_tcf3_pbx1` : **The extract_tcf3_pbx1 extracts TCF3::PBX1 chromosomal abnormality into a two-column data frame comprising of the persond_id and TCF3::PBX1 columns. The combination of t(1;19)(q23;p13) pattern from the karyotype and relevant [concept_ids](#list-of-concept-ids) were used to extract this particular genetic abnormality. The function takes the pre-processed person and measurement data frames as input arguments.**
~~~
extract_tcf3_pbx1(person = person, measurement = measurement)
~~~

## Extract TCF3::HLF
* `extract_tcf3_hlf` : **The extract_tcf3_hlf function extracts the rare TCF3::HLF genetic abnormality using the appropriate [concept_ids](#list-of-concept-ids) to produce a two-column data frame as output. Using t(17;19)(q22;p13) pattern in the karyotype did not produce additional information on the presence or otherwise of the genetic condition on all the patients on the HARMONY platform. The function takes the pre-processed person and measurement data frames as input arguments.**
~~~
extract_tcf3_hlf(person = person, measurement = measurement)
~~~

## Extract iAMP21
* `extract_iamp21` : **The extract_iamp21 function extracts the presence or otherwise of the iAMP21 genetic abnormality for all the patients on the HARMONY platform into a two-column data frame taking the pre-processed person and measurement data frames as input arguments, and using the appropriate [concept_ids](#list-of-concept-ids).**
~~~
extract_iamp21(person = person, measurement = measurement)
~~~

## Extract High hyperdiploidy, Near Haploid/Haploidy, and Low hypodiploidy
* `extract_heh_hap_LH` : **The extract_heh_hap_LH function like most of the functions in the package takes the pre-processed person and measurement data frames as input arguments and outputs a four-column data frame comprising of the persond_id, HeH (High-hyperdiploidy), Hap (Haploidy and Near-haploidy), and LH (Low-hypodiploidy and Near-triploidy) columns. The function helps characterise the different aneuploidy subgroups in the HARMONY datasets. The HeH refers to the High-Hyperdiploidy subgroup having between 51-59 and 60-67 modal chromosome numbers and satisfying other pre-defined criteria. The Hap refers to the Haploidy or Near-Haploidy subgroup with between 20-29 modal chromosome numbers. The LH represents the Low-Hypodiploidy and Near-Triploidy subgroup with 30-39, 68-78 and 60-67 modal chromosome numbers with some pre-defined conditions. A combination of concept_ids and well-defined patterns in the karyotype was used to characterise each patient in the HARMONY datasets into the HeH, Hap and LH subgroups. A table containing all the concept_ids used to extract each of the genetic abnormalities considered in this package is available [here](#list-of-concept-ids)**
~~~
extract_heh_hap_hotr(person = person, measurement = measurement)
~~~

## Extract Complex Karyotype
* `extract_complex_karyotype` : **The extract_complex_karyotype function extracts all individuals with five or more chromosomal abnormalities into a two-column data frame, taking the pre-processed person and measurement data frames as inputs. Relevant [concept_ids](#list-of-concept-ids) were used to extract patients into 3 groups - Present, Absent and Unknown indicating the presence of complex karyotype, absence of complex karyotype and no information on complex karyotype for each patient in the HARMONY datasets. The unknown also refers to individuals with no recorded karyotype at diagnosis. Karyotypes of all the patients with recorded karyotype at diagnosis were examined very closely to identify patients with the presence of 5 or more chromosomal abnormalities and were added to those already captured by the concept_ids to produce a data frame that describes the presence or otherwise of complex karyotype in every patient in the HARMONY datasets.**
~~~
extract_complex_karyotype(person = person, measurement = measurement)
~~~

## Extract B_other, B_other_plus, and T_other Subgroups
* `extract_B_other_B_other_plus_and_T_other` : **The extract_B_other_B_other_plus_and_T_other function combines all the 12 previously described chromosomal abnormalities with B_other, B_other_plus, T_other, and Unknown columns to fully describe all the patients in the HARMONY datasets with respect to the chromosomal abnormalities of interest. B_other classifies patients with B-cells that have not been classified by the previous 12 chromosomal abnormalities in their order of preference into Present, Absent and Unknown. Present means the presence of B-cells with normal or abnormal karyotype status. Absent represents patients with B-cells that have already been described by at least one of the 12 previously defined genetic abnormalities, while Unknown represents patients with B-cells with no information to adequately classify them into any of the previously defined genetic subtypes. B_other_plus captures these individuals and further distils them into Present, Absent and Unknown. Present in the B_other_plus subgroup means the "Unknown" subgroup in B-other with B-cells that have not been classified into any of the previously defined genetic subtypes. T_other captures into Present subgroup patients with T-cells that have normal or abnormal karyotype status and are not already classified into any of the previously defined genetic subtypes. Absent in the T_other column represents patients with T-cells that have already been classified into at least one of the previously defined genetic subtypes, while the unknown simply represents all the T-cell patients with no adequate information to classify them into any of the genetic subtypes of interest.**
~~~
extract_B_other_B_other_plus_and_T_other(person = person, measurement = measurement, condition_occurrence = condition_occurrence)
~~~


## List of concept IDs

**The concept ids used to extract all the genetic abnormalities are shown in the concept id table below.**
<!--- [here](https://github.com/ayolawal/harmonygeneticsubtype/blob/main/data-raw/Concept_id_List.xlsx) --->

<details><summary>Concept ID Table</summary>
<p>

| Abnormality | Concept_id | Column | Table | Concept_Name |
|:---------|:---------|:--------|:---------|:----------|
|karyotype| 40765097 | measurement_concept_id | Measurement | Chromosome analysis result in ISCN expression |
|Age | 3007016 | measurement_concept_id | Measurement | Age at cancer diagnosis |
| ETV6::RUNX1 | 3036903 <br> 35977015 | measurement_concept_id | Measurement | t(12;21)(p13;q22.3)(ETV6,RUNX1) fusion transcript in Blood or Tissue by Molecular genetics method <br> Karyotype t(12;21)(p13;q22) |
| BCR::ABL1 | 3011913 <br> 3013321 <br> 3016815 <br>3028956 <br> 46235651 <br>35977026  | measurement_concept_id <br> measurement_concept_id <br> measurement_concept_id <br> measurement_concept_id <br> measurement_concept_id <br> measurement_source_concept_id | Measurement | t(9;22)(q34.1;q11)(ABL1,BCR) fusion transcript, major, and minor break points in Blood or Tissue by Molecular genetics method <br> Karyotype t(9;22)(q34;q11.2) | 
| KMT2A::AFF1 | 2000109034 <br>  3012425 <br> 37030733 | measurement_source_concept_id  <br>   measurement_concept_id <br> measurement_source_concept_id  | Measurement  | Karyotype t(4;11)  <br> t(4;11)(q21.3;q23)(AFF1,MLL) fusion transcript in Blood or Tissue by Molecular genetics method |
| KMT2A::MLLT1 | 3002279  <br> 37021664 | measurement_concept_id <br>  measurement_source_concept_id | Measurement |  t(11;19)(q23;p13.3)(MLL,MLLT1) fusion transcript in Blood or Tissue  |
| KMT2A_r | 2000000263 <br> 2000000269 <br> 3022443 <br> 3037009  <br> 36017933 | measurement_source_concept_id <br> measurement_source_concept_id <br> measurement_concept_id <br> measurement_concept_id <br> measurement_concept_id | Measurement | MLL partial tandem duplication analysis <br> 11q23 rearrangement by cytogenetics <br> t(9;11)(p22;q23)(MLLT3,MLL) fusion transcript  in Blood <br> MLL gene rearrangements in Blood <br> t(9:11)(p22;q23);MLLT3-MLL |
| TCF3::PBX1 | 2000108068 <br> 3000296 | measurement_source_concept_id <br> measurement_concept_id | Measurement | Karyotype t(1;19)(q23;p13) <br> t(1;19)(q23.3;p13.3)(PBX1,TCF3) fusion transcript in Blood or Tissue |
| TCF3::HLF	| 42868760	| measurement_concept_id | Measurement | t(17;19)(q22;p13.3)(HLF,TCF3) fusion transcript in Blood or Tissue |
| iAMP21 | 2000000466 | measurement_source_concept_id	| Measurement | Intrachromosomal amplification of chromosome 21 measurement |
| HeH | 2000000462 <br> 36660734 | measurement_source_concept_id <br> measurement_concept_id | Measurement | Detection of high hyperdiploidy <br> Chromosome aneuploidy details in Blood or Tissue by Molecular genetics method Narrative |
| Hap	| 2000000464 | measurement_source_concept_id | Measurement | Detection of near haploidy measurement | 
| LH | 2000000463 <br> 36660734 | measurement_source_concept_id <br> measurement_concept_id | Measurement | Detection of low hypodiploidy <br> Chromosome aneuploidy details in Blood or Tissue by Molecular genetics method Narrative |
| Complex Karyotype | 2000000153 <br> 45877994 | measurement_source_concept_id <br> value_as_concept_id | Measurement | Complex karyotype <br> Yes |
| B-Cell | 4082461 <br> 4173963 | condition_concept_id | Condition_Occurrence | Precursor B-cell acute lymphoblastic leukaemia <br> B-cell acute lymphoblastic leukaemia |
| T-Cell | 4082464 | condition_concept_id	| Condition_Occurrence | T-cell acute lymphoblastic leukaemia |  

</p>
</details>
