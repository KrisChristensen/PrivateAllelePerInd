##########################################################
### Import Necessary Modules #############################

import argparse                       #provides options at the command line
import sys                            #take command line arguments and uses it in the script
import gzip                           #allows gzipped files to be read
import re                             #allows regular expressions to be used

##########################################################
### Command-line Arguments ###############################

parser = argparse.ArgumentParser(description="A script to identify the distribution of private alleles within groups/populations per individual.  Private alleles are identified between groups and they are counted per individual in that group.  This is output so the user can find the average number of private alleles per individual per site or population.")
parser.add_argument("-vcf", help = "The location of the vcf file", default=sys.stdin, required=True)
parser.add_argument("-pop", help = "The location of the population file (Individual<tab>Population<tab>Site)", default=sys.stdin, required=True)
parser.add_argument("-min", help = "The minimum number of individuals with the private allele in a population before it is considered in this analysis, default=1", default=1)
args = parser.parse_args()

#########################################################
###Variables ############################################

class Variables():
   population = {}  ### Individual => Population
   site = {}        ### Individual => Sampling site
   indCounts = {}   ### Individual => Count of private alleles from the population found for this individual 
   numIndividuals = 0 # The count of individuals in the population file (used to check the count in the vcf file)

#########################################################
### Body of script ######################################

class OpenFile():
    def __init__ (self, f, typ, occ):
        """Opens a file (gzipped) accepted"""
        if re.search(".gz$", f):
            self.filename = gzip.open(f, 'rb')
        else:
            self.filename = open(f, 'r') 
        if typ == "vcf":
            sys.stderr.write("\nOpened vcf file: {}\n".format(occ))
            OpenVcf(self.filename,occ)
        elif typ == "pop":
            sys.stderr.write("\nOpened pop file: {}\n".format(occ))
            OpenPop(self.filename,occ)

class OpenVcf():
    def __init__ (self,f,o):
        """Reads a vcf file to identify private alleles and count per group and site"""
        self.individuals = []   
        for self.line in f:
            ### Allows gzipped files to be read ###
            try:
                self.line = self.line.decode('utf-8')
            except:
                pass
            self.line = self.line.rstrip('\n')             
            if not re.search("^#", self.line):
                self.chr, self.pos, self.id, self.ref, self.alt, self.qual, self.filt, self.info, self.fmt = self.line.split()[0:9]
                self.individualGenotypes = self.line.split()[9:]
                self.privateTest = {}
                self.privateTest["0"] = {}
                self.privateTest["1"] = {}
                self.siteRef = {}
                self.siteAlt = {}
                for self.position, self.indGeno in enumerate(self.individualGenotypes):
                    self.indName = self.individuals[self.position]
                    try:
                        self.indPop = Variables.population[self.indName]
                    except:
                        continue ###VCF file can have individuals not assigned in population file.  These are ignored.
                    try:
                        self.indSite = Variables.site[self.indName]
                    except:
                        continue
                    self.indGenotype = self.indGeno.split(":")[0].split("/")
                    if self.indGenotype[0] == "." or self.indGenotype[1] == ".":
                        continue  ### Missing genotypes are ignored.
                    if int(self.indGenotype[0]) == 0 or int(self.indGenotype[1]) == 0:
                        if self.indPop in self.privateTest["0"]:
                            self.privateTest["0"][self.indPop] += 1
                        else:
                            self.privateTest["0"][self.indPop] = 1
                        if self.indName in self.siteRef:
                            self.siteRef[self.indName] += 1
                        else:
                            self.siteRef[self.indName] = 1
                    if int(self.indGenotype[0]) == 1 or int(self.indGenotype[1]) == 1:
                        if self.indPop in self.privateTest["1"]:
                            self.privateTest["1"][self.indPop] += 1
                        else:
                            self.privateTest["1"][self.indPop] = 1
                        if self.indName in self.siteAlt:
                            self.siteAlt[self.indName] += 1
                        else:
                            self.siteAlt[self.indName] = 1                        
                if len(self.privateTest["0"].keys()) == 1:
                    for self.pop in self.privateTest["0"]:
                        if int(self.privateTest["0"][self.pop]) >= int(args.min):
                            for self.indName in self.siteRef:
                                if self.indName in Variables.indCounts:
                                    Variables.indCounts[self.indName] += 1
                                else:
                                    Variables.indCounts[self.indName] = 1
                elif len(self.privateTest["1"].keys()) == 1:
                    for self.pop in self.privateTest["1"]:
                        if int(self.privateTest["1"][self.pop]) >= int(args.min):
                            for self.indName in self.siteAlt:
                                if self.indName in Variables.indCounts:
                                    Variables.indCounts[self.indName] += 1
                                else:
                                    Variables.indCounts[self.indName] = 1                                                                   
            elif re.search("^#CHROM", self.line):
                self.individuals = self.line.split()[9:]
                self.numInds = len(self.individuals)
                if int(self.numInds)  != int(Variables.numIndividuals):
                    sys.stderr.write("Warning, population and vcf files have different number of individuals. VCF: {}, Pop: {}\n".format(int(self.numInds), int(Variables.numIndividuals)))
        for self.ind in Variables.indCounts:
            self.indPop = Variables.population[self.ind]
            self.indSite = Variables.site[self.ind]
            print("{}\t{}\t{}\t{}".format(self.ind, self.indPop, self.indSite, Variables.indCounts[self.ind]))
        f.close()
        
class OpenPop():
    def __init__ (self,f,o):
        """Reads a population file to identify which individual goes to which population/site"""
        self.pops = {}
        self.sites = {}
        for self.line in f:
            ### Allows gzipped files to be read ###
            try:
                self.line = self.line.decode('utf-8')
            except:
                pass
            self.line = self.line.rstrip('\n')            
            if not re.search("^#", self.line):
                self.individual, self.pop, self.site = self.line.split()
                if self.individual in Variables.population:
                    sys.stderr.write("\tWarning: {} already defined for population {}, replacing with population {}\n\n".format(self.individual, Variables.population[self.individual], self.pop))
                Variables.population[self.individual] = self.pop
                if self.individual in Variables.site:
                    sys.stderr.write("\tWarning: {} already defined for site {}, replacing with site {}\n\n".format(self.individual, Variables.site[self.individual], self.site))
                Variables.site[self.individual] = self.site
                if self.pop in self.pops:
                    self.pops[self.pop] += 1
                else:
                    self.pops[self.pop] = 1
                if self.site in self.sites:
                    self.sites[self.site] += 1
                else:
                    self.sites[self.site] = 1
                Variables.numIndividuals += 1
        for self.pop in sorted(self.pops):
            sys.stderr.write("\tIdentified population {} with {} samples\n".format(self.pop, self.pops[self.pop]))
        for self.site in sorted(self.sites):
            sys.stderr.write("\tIdentified site {} with {} samples\n".format(self.site, self.sites[self.site]))
        f.close()

### Order of script ####
if __name__ == '__main__':
    Variables()
    open_aln = OpenFile(args.pop, "pop", args.pop)
    open_aln = OpenFile(args.vcf, "vcf", args.vcf)
