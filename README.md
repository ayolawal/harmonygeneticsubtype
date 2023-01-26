# Harmony Genetic Subtype Package

Healthcare Alliance for Resourceful Medicines Offensive against Neoplasms in HematologY (HARMONY) is a European network of excellence that captures, integrates, analyzes, and harmonizes big data from high-quality multidisciplinary sources with the main purpose of unlocking valuable knowledge on various important hematologic malignancies (HMs).  

HARMONY aims at assembling, connecting and analyzing big data from HM patients to define standard sets of outcome indicators to be measured and used for decision-making by key healthcare systems stakeholders.

HarmonyGeneticSubtype is an R package developed by the Leukaemia Research Cytogenetics Group (LRCG), Faculty of Medical Sciences at Newcastle University, specifically to extract some pre-defined chromosomal abnormalities from data collected from thousands of patients for downstream analyses. The package contains different tools that could be used to distil for instance individuals with certain genetic conditions like down syndrome, Intrachromosomal amplification of chromosome 21 (iAMP21), chromosomal translocations - ETV6-RUNX1, BCR-ABL1 etc from tens of thousands of high-quality anonymised harmonised hematologic malignancies datasets.

The package leverage on the representation of the datasets in Observational Medical Outcomes Partnership (OMOP) common data models (CDM) standardised tables format in the ekberfall database. The package is continually being updated and the current tools that are available for use on the HARMONY platform are listed below. 


## Tools Description

* `extract_karyotype`: **Extracts karyotype from the OMOP measurement table**.
```
extract_karyotype(measurement = measurement)
```
* `age_gender_df`: **Extracts age and gender into a dataframe**. 

~~~
age_gender_df(person = person, concept = concept, measurement = measurement)
~~~

* `kary_status`: **Extracts karyotype status into normal, failed, abnormal, and not done**.
~~~
kary_status(person = person, measurement = measurement)
~~~

* `etv6_runx1` : **Extracts ETV6::RUNX1 principal genetic abnormality**.
~~~
etv6_runx1(person = person, measurement = measurement)
~~~

* `bcr_abl1` : **Extracts BCR::ABL1 principal genetic abnormality**.
~~~
bcr_abl1(person = person, measurement = measurement)
~~~

* `kmt2a_aff1` : **Extracts KMT2A::AFF1 principal genetic abnormality**.
~~~
kmt2a_aff1(person = person, measurement = measurement)
~~~

* `kmt2a_mllt1` : **Extracts KMT2A::MLLT1 principal genetic abnormality**.
~~~
kmt2a_mllt1(person = person, measurement = measurement)
~~~

* `kmt2a_r` : **Extracts KMT2A_r principal genetic abnormality**.
~~~
kmt2a_r(person = person, measurement = measurement)
~~~

* `tcf3_pbx1` : **Extracts TCF3::PBX1 principal genetic abnormality**.
~~~
tcf3_pbx1(person = person, measurement = measurement)
~~~

* `tcf3_hlf` : **Extracts TCF3::HLF principal genetic abnormality**.
~~~
tcf3_hlf(person = person, measurement = measurement)
~~~

* `iamp21` : **Extracts iAMP21 principal genetic abnormality**.
~~~
iamp21(person = person, measurement = measurement)
~~~

* `heh_hap_hotr` : **Extracts heh principal genetic abnormality**.
~~~
heh_hap_hotr(person = person, measurement = measurement)
~~~

* `complex_karyotype` : **Extracts complex karyotype**.
~~~
complex_karyotype(person = person, measurement = measurement)
~~~

* `principal_abn_B_other_B_other_plus` : **Extracts  B_other and B_other_plus with all the principal genetic abnormalities**.
~~~
principal_abn_B_other_B_other_plus(person = person, measurement = measurement, condition_occurrence = condition_occurrence)
~~~

* `principal_abn_T_other_no_data` : **Extracts everything in `principal_abn_B_other_B_other_plus` plus T_other and No data**.
~~~
principal_abn_T_other_no_data(person = person, measurement = measurement, condition_occurrence = condition_occurrence)
~~~

* `gsubtype_and_riskgroup`: **Extracts genetic subtype and risk groups for the identified chromosomal abnormalities**
~~~
gsubtype_and_riskgroup(person = person, measurement = measurement, condition_occurrence = condition_occurrence)
~~~
