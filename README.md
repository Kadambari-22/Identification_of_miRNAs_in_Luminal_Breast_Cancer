# Identification_of_miRNAs_in_Luminal_Breast_Cancer
Conducted a comprehensive bioinformatics analysis to identify expression of microRNAs in luminal breast cancer, providing insights into their potential  regulatory roles and therapeutic implications.

Based on clinical parameters and molecular profiling, Breast cancer is one of the most common malignancies, with multiple subtypes. Clinically, breast cancers are defined as ER+/PR+, HER+, or triple-negative. The ER+/PR+ breast cancers are predominately luminal A/B, which make up up to 70% of all BC patients. miRNAs, the non-coding RNAs that regulate gene expression by degrading mRNAs or inhibiting translation. miRNAs can play a dual role in breast cancer, with some acting as oncogenes and others acting as tumor suppressors. By assimilating and analyzing vast datasets from diverse sources, Bioinformatics tools can efficiently extract meaningful features and help identify miRNA signatures, predict targets, and predict tumor-specific biomarkers. 

Microarray data of microRNA expression was obtained from the Gene Expression Omnibus (GEO) repository using the GEOquery package in RStudio. Techniques including log transformation, normalization, clustering, principal component analysis and DEA facilitated the identification of 59 dysregulated microRNAs, with 45 upregulated microRNAs and 14 downregulated microRNAs (LogFC>1.5 & adj.p-value<0.05), offering insights into underlying data patterns and variability. 
DEA_code.R

MIENTURNET, an online tool was used for Over-Representation Analysis (ORA) & Functional enrichment analysis to prioritize miRNA-target interactions & gain insights into the biological processes underlying target gene activity.  The tool revealed specific microRNAs' potential roles in critical pathways and their interactions with key target genes implicated in tumorigenesis. 
http://userver.bio.uniroma1.it/apps/mienturnet/

The study employs a pioneering approach by integrating machine learning and bioinformatics in breast cancer research.  This study helped me understand the intricate regulatory roles of miRNAs in luminal BC, shedding light on potential therapeutic targets and pathways for further investigation for the development of targeted therapies aimed at improving patient outcomes. 


