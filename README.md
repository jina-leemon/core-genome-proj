# Core Genome Identification code

This code can be used to identify a set of core genes conserved in high percentages across a specified list of phages (or other organisms).
For analysis, you need: <br />
1. The result of running the proteins of the phage in psiblast (can be done with a local installation of PsiBlast, output format must be 9). I habe not tried non-PsiBlast results but theoretically it should work just as long as it's the right format
2. Text file with all the phages of interest.
The code extracts the homologs from phages of interest and outputs a .csv file with all the proteins in the query for the PsiBlast and the homologs above a certain threshold. There will be no overlaps (accession # unless the same protein has two accession numbers)
It then outputs a list of proteins that are conserved above a certain percentage across the phages of interest.
