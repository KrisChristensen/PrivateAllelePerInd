# PrivateAllelePerInd
A script to identify the number of private alleles between populations/groups and count the number of private alleles each individual has for that group (useful to identify if paricular individuals or sites contribute more private alleles between populations/groups, or to identify gradients).  Heterozygous and homozygous genotypes are counted as one occurence of the private allele instead of two for homozygous genotypes.

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#requirements">Requirements</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>

<!-- requirements -->
## Requirements

This script has been tested with python 3.
The script requires the following files:

&nbsp;&nbsp;&nbsp;VCF file that you wish to identify private alleles from<br />
&nbsp;&nbsp;&nbsp;A population file with each line having the format: individual population site<br />
&nbsp;&nbsp;&nbsp;The site in the population file can be made up.  It isn't used in the script except for reporting.<br />

<!-- usage -->
## Usage

1) Output the count of private alleles per individual (format: individual population site PrivateAlleleCountPerIndividual):<br /><br />
&nbsp;&nbsp;&nbsp;python PrivateAlleleIdentifierDistributionPerIndividual.v1.0.py -vcf file.vcf.gz -pop population.txt > PrivateAlleles.txt<br /><br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python PrivateAlleleIdentifierDistributionPerIndividual.v1.0.py -h

<!-- license -->
## License 

Distributed under the MIT License.
