# tax_QA

This is data and code that contributed to the following dissertation chapter, which is in preparation for publication:

**A robust template for increasing taxonomic quality assurance in an era of declining taxonomic capacity** 

*Erica T. Jarvis Mason<sup>1</sup>, Andrew R. Thompson<sup>2</sup>, Brice X. Semmens<sup>1</sup>*

<sup>1</sup>Scripps Institution of Oceanography, University of California San Diego, California, USA  
<sup>2</sup>NOAA Southwest Fisheries Science Center, California, USA  

**R-script files:**

1. JarvisMasonetal_AccPrec.R (Model probabilities of accurate and precise species classifications, and add an effect of taxonomist)
   - Compare and contrast taxonomist skill (Bayesian binomial model, with options for adding a fixed effect or random effect of taxonomist)

3. JarvisMasonetal_TaxMorph.R (Explore the utility of a suite of characters for species discrimination, and identify the most important characters)
   - Identify potential areas of taxonomist bias, subjectivity in interpreting characters (Bayesian multinomial logistic regression)
   - Identify which characters are most important (random forest)
