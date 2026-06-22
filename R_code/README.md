

############################################################################################
# Assume you have already computed the TG alone,  ApoB_alone, and LDL-TG GWAS
#  
# For diabetes gwas outcomes: https://diagram-consortium.org/downloads.html
#
#
# Track A- Total effect (Univariable MR: LDL-TG)
#      Exposure: LDL-TG  
#      Outcome: CAD
#
# Hugh 12.1.2025 Updated
###########################################################################################



Step 1:   Reclump UKB GWAS data and output the summary statistics: clump for univariable LDL-TG total effect
Step 2:   Drop palindromic variants
Step 3:  Find the F-score for single exposure 
Step 4:   Harmonization with the CAD outcomes, SNP position : GRCh37
Step 5:   Univariable MR (Two-sample MR)
