# YFinder
Call Yfull (v11.02.0) Y haplogroups of multiple samples from a vcf file of chr Y

**TESTING**<br />
Tested on 1200+ modern males from 1000genomes project with 100% success rate (verified with Yfull (https://www.yfull.com/full-genomes/)<br />
Yfull has more resolution as it used whole genomes, whereas the 1000GP VCF file for chrY has lesser data.<br />
Also tested on 50 ancient male samples from V54.1 Harvard Allen Database. (converted Geno to Bed to VCF for select samples only for chr Y).<br />
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz<br />


**INSTALLATION** <br />
Works only on R <br />
install.packages("devtools") <br />
devtools::install_github("SahilArvind/YFinder") <br />

**LOAD PACKAGE**<br />
library(YFinder)<br />

**REQUIREMENTS**<br />
1. Download the tree_11.02.0.json and YFullpositions.txt file from the /data folder. They will be needed as input<br />
2. VCF data with Y Chr only (only hg19 reference genome supported in v1.0). Can use plink --chr 24 --recode vcf to convert bed/bim/fam to vcf only for Y <br />

**USAGE**<br />
output <- YFinder(vcffile,yfullpositionsfile,indivlistfile,
                          yfulljsonfile,transversionsonly=F,depth=1,
                          callquality=0.9,numsnptotal=2,viewrealtime=FALSE)<br />

**WHERE**<br />
1. vcffile, yfullpositionsfile,yfulljsonfile = /path/to/file/name of those respective files<br />
2. indivlistfile = txt file with list of samples to test. One per line and last line empty. Names should match those in vcffile exactly or they will be ignored. Plink conversion to VCF concatenates label name and sampleid so take care of that while making indivlist.txt <br />
3. transversionsonly: Default FALSE. Ignores C/T and G/A positions. Keep TRUE for ancient samples as deamination causes false C/T and G/A calls.<br />
4. numsnptotal: Default 2. Checks number of SNPs available for a subclade and assigns subclade +ve only if num snps >= numsnptotal and derived > callquality. Do not keep 1 as there will be many false derived calls. Higher the number of numsnptotal, more conservative is the estimate of haplogroup <br />
5. callquality: Default 0.9. Assigns subclade +ve only if % of derived positions > callquality. eg 4/5 positions derived for a subclade, but subclade will be marked Ancestral because 80%<90%.
6. depth: Used in case a particular SNP is derived, but only 1 SNP available for that subclade for the sample. This could be a false call. So ancestor status is checked. Depth tells us how many ancestors to check to give a derived status. If depth=1, one ancestor will be checked. If that subclade is Ancestral, then subclade in question is marked ancestral. If ancestor is Derived, then our subclade in question is also marked Derived. If ancestral subclade SNPs are missing, ancestor 1 step above is checked. Default depth 1, higher the number, more conservative is the haplogroup estimate.<br />
7. viewrealtime: Default FALSE. Keep TRUE if you want realtime update of output. Useful in case of large number of samples with a long runtime of function.<br />

**OUTPUT**<br />
List object <br />
view(output$Final_Assignment) will show output table<br />
view(output$CallsTable) will show all ancestral/derived status for all available SNPs in a table for manual cross checking<br />
Outputs can be written to file using write.table() function in R.
