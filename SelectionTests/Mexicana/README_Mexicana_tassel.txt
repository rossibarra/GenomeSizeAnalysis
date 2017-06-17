Tassel Notebook

Raw file with the randomly imputed genotypes:

	Mexicana_Highelev_RndImp.txt

Read into tassel, selecting out only the individuals we want in the high elevation mexicana trial.  This file has single bp genotypes that we will need to convert to -1 00 01.

	(All mexicana populations except AM M TZ)

When we were dealing with the imputed data, we would also remove (da6 tx12) for missing data, but they are already absent in the raw randomly imputed file.
With the populations selected, export out of tassel to the file GBS_forSNPconversion.txt
Remove the header line, and run our perl script to convert snps to -1 0 01

	sed 1d GBS_forSNPconversion.txt > GBS_forSNPconversion2.txt
	perl ConvertSNP2.pl GBS_forSNPconversion2.txt rrblup_highelmex.csv
	
This file can be read into R and converted to the kinship matrix via rrBLUP.
For further analyses, see Mexicana_selectiontests.R
