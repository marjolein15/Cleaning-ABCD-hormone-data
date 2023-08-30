# Cleaning-ABCD-hormone-data
These scripts clean the salivary hormone data of the ABCD study. 

The script hormone_cleaning_baseline_2year_g.Rcontains several steps to clean the hormone data and confounds and correct the hormone levels for confounds.
It first follows the decision tree in Herting et al. 2021 (doi: 10.3389/fendo.2020.549928). Then it finds nonsensical values in time confounds and set these to NA. It adds in medication use variables (see below for explanation on how these are created). It calculates mentrual cycle regularity and phase. The hormone levels are log-transformed and regressed on confounds (hours since wake, collection duration, recent caffeine use, recent physical activity,  medication use).

The hormone cleaning script relies on the output of the medication script (categorize_medications_g.R).
Medication use in ABCD was categorized using the RXNORM database maintained by the NLM within the NIH. The medication codes in RXNORM are called RXCUI. 
To make the data useful as control variable(s), all RXCUIs were converted to MESHPA, another database maintained by NIH that has broader classifications based on physiological actions. This was done online using NIH's tool RxMix: https://mor.nlm.nih.gov/RxMix/. The output of this is called rxmix_output_meshpa_date.txt
This file is then used to categorize the ABCD data using the script categorize_medications_g.R
