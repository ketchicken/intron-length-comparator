# Classes used within the system

import matplotlib.pyplot as plt
import re
import statistics
from scipy.stats import linregress
from scipy.stats import zscore
import math
from ortholog_dictionary.intron_dict_creator import IntronTranscriptIDDictionary as intronDict

class BEDReader:
    """ 
    Docstring
    Reads the BED files, annotated with UCSC refFormat
    """
    # gene symbol matching
    # ex. 7_intron_21_flank0_chr1_
    # can use this format to regex something
    
    REG_INTRON_ID_0 = r"\d+_intron_"
    REG_INTRON_ID_1 = r"_flank\d+_chr\d+_"
    REG_INTRON_ID = r"\d+_intron_\d+_flank\d+_chr\d+_"
    REG_ENSEMBL_ID = r"ENS([A-Z]+)?T[0-9]{11}"
    

    def __init__(self, bedFile, species_name=None):
        self.intronlengths = {}
        with open(bedFile, 'r') as bed:

            # assume no header first
            header = bed.readline()
            if "track name" not in header:
                # Then there was no file header included in the BED file
                # and the first line is BED line
                intronID, start, stop = self.parse_line(header)
                if intronID:
                    self.add_intron(intronID, start, stop)  
            
            for line in bed:
                intronID, start, stop = self.parse_line(line)
                if intronID:
                    self.add_intron(intronID, start, stop)

        if species_name is None:
            # Grab the last parsed line, and get the ensembl species code
            intron_info = line.split("\t")
            # [chromosome (0), start (1), stop (2), name (3), score (4), strand direction (5)]
            species_name = re.search(self.REG_ENSEMBL_ID, intron_info[3])

            if species_name and species_name.group(1):
                species_name = species_name.group(1)
            else:
                # nothing was found, so assume it is a human 
                species_name = "HUM"
        
        self.species_name = species_name
            

    def add_intron(self, intronID, start, stop):
        if intronID not in self.intronlengths.keys():
            self.intronlengths[intronID] = [math.log10(abs(stop - start))]
        else:
            self.intronlengths[intronID].append(math.log10(abs(stop - start)))
        

    def parse_line(self, line):
        intron_info = line.split("\t")
        # This results in a list with the following entries:
        #   [chromosome (0), start (1), stop (2), name (3), score (4), strand direction (5)]
        # BED files downloaded from the UCSC Genome Browser Table Browser will have region names
        # in the Ensembl Format. The code follows as suit.

        # the 2nd number found in intronID is actually the part that the intron has been split to
        
        intronID = re.search(self.REG_ENSEMBL_ID, intron_info[3])
        if intronID:
            return intronID.group(), int(intron_info[1]), int(intron_info[2])
        else:
            return None, None, None
    
    def get_intron_lengths(self):
        return self.intronlengths
    
    def get_species_name(self):
        return self.species_name
    
    def get_average_intron_lengths(self):
        avg_length = {}
        for ensid in self.intronlengths.keys():
            avg_length[ensid] = statistics.mean(self.intronlengths[ensid])
        return avg_length


class IntronComp:
    """ Docstring """

    def __init__(self, bedFile1, bedFile2, id_dict, std):
        self.bed1 = BEDReader(bedFile1)
        self.bed2 = BEDReader(bedFile2)

        self.std = float(std)

        self.main_species_name = self.bed1.get_species_name()
        self.other_species_name = self.bed2.get_species_name()
        self.intron_dict = id_dict[self.main_species_name]
    
    def pair_avg_ortholog_lengths(self, s1_intron_lengths, s2_intron_lengths):
        x_val = []
        y_val = []

        for enst_id in s1_intron_lengths.keys():
            orthologs_for_intron = self.intron_dict.get(enst_id, None)
            if orthologs_for_intron:
                ortholog_id = orthologs_for_intron.get(self.other_species_name, None)
                if s2_intron_lengths.get(ortholog_id, None):
                    x_val.append(s1_intron_lengths[enst_id])
                    y_val.append(s2_intron_lengths[ortholog_id])

        return x_val, y_val

    def pair_ortholog_lengths(self, s1_introns, s2_introns):
        x_val = []
        y_val = []
        
        for enst_id in s1_introns.keys():
            s1_orths = self.intron_dict.get(enst_id, None)
            if s1_orths:
                s2_orth_id = s1_orths.get(self.other_species_name, None)
                if s2_introns.get(s2_orth_id, None):
                    # pair up the ones 
                    for i in range(min(len(s1_introns[enst_id]), len(s2_introns[s2_orth_id]))):
                        x_val.append(s1_introns[enst_id][i])
                        y_val.append(s2_introns[s2_orth_id][i])

                    while len(x_val) < len(y_val):
                        x_val.append(0)     # fill with 0s
                    
                    while len(y_val) < len(x_val):
                        y_val.append(0)     # fill with 0s

        return x_val, y_val
    
    def remove_outliers(self, x, y):
        print("Removing Outliers")
        # get the mean, get the std
        # go through values and pop out any pairs that are not part of the thing
        if self.std == -1:
            # no need to process
            return x, y
        
        x_zscores = zscore(x).tolist()
        y_zscores = zscore(y).tolist()

        i = 0
        while i < len(x_zscores):
            if abs(x_zscores[i]) > self.std:
                x.pop(i)
                y.pop(i)
                x_zscores.pop(i)
                y_zscores.pop(i)
            elif abs(y_zscores[i]) > self.std:
                x.pop(i)
                y.pop(i)
                x_zscores.pop(i)
                y_zscores.pop(i)
            else:
                i+=1
        
        return x, y
    
    def plot_graph(self, outfile_path=None, average=False, **kwargs):
        '''
        Saves an image comparing intron lengths to specified path

        Parameters:
        :param: outfile_path: the path to save the image. if unspecified, it saves in the current directory.
        :param: average: whether or not to plot based on the average length of all grouped introns of a gene
        :param: kwargs: formatting arguments for graph
        '''
        print("Starting Plot")

        if outfile_path is None:
            outfile_path = "."
        
        if ".png" not in outfile_path:
            outfile_path += f"/intron_comparator_{self.bed1.get_species_name()}_{self.bed2.get_species_name()}.png"

        fig, ax = plt.subplots(figsize=kwargs.get("figsize", (7, 7)))

        if average:
            x_val, y_val = self.pair_avg_ortholog_lengths(self.bed1.get_average_intron_lengths(), self.bed2.get_average_intron_lengths())
        else:
            x_val, y_val = self.pair_ortholog_lengths(self.bed1.get_intron_lengths(), self.bed2.get_intron_lengths())

        if not x_val or not y_val:
            print(" No orthologous introns found ")
            return 
        
        x_clean, y_clean = self.remove_outliers(x_val, y_val)

        slope, intercept, r_value, _, _ = linregress(x_clean, y_clean) # ignore outliers
            
        ax.scatter(x_val, y_val, s=kwargs.get("size", 2.5), c=kwargs.get("color", "blue"), marker='o')
        
        ax.plot(x_val, [intercept + slope * x for x in x_val], linestyle='-', color='red')

        ax.set_xlabel(f"{self.bed1.get_species_name()} Total Intron length (log10bp)")
        ax.set_ylabel(f"{self.bed2.get_species_name()} Total Intron length (log10bp)")                                     
        ax.set_title(f"{self.bed1.get_species_name()} v. {self.bed2.get_species_name()} Intron Length Correlation\nr-value = {r_value} (n = {len(x_val)})")

        fig.savefig(outfile_path)
        print("Completed successfully")
        
# Command line flags
import argparse as arg
class ArgParser:
    '''
        Implement a parser to interpret the command line argv string using argparse.
    '''
    def __init__(self, inOpts=None):
        
        ap = arg.ArgumentParser(description = ' Compares the length of paralog genes given two BED files. \n BED files are taken from UCSC genome browser format', 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s <bed1> <bed2> [options] <tsv>')
        ap.add_argument('bed1', help='Path to bed files')
        ap.add_argument('bed2', help='Path to bed files')
        ap.add_argument('-o', '-outdir', action='store', nargs='?', default=None, help='path to output file')
        ap.add_argument('-q', '-quiet', action='store', nargs='?', default=False, help='Toggle for output log printing')
        ap.add_argument('-std', '-standardDeviations', action='store', nargs="?", default=-1, help='Which outliers to consider - how many standard deviations away. default -1')
        ap.add_argument('-s1', '-species1', action='store', nargs="?", default=None, help="Name of first species (If not specified, will list Ensembl ID from BED file)")
        ap.add_argument('tsv', help="path to tsv file of orthologs of species of bed1")

        if inOpts is None :
            self.args = ap.parse_args()
        else :
            self.args = ap.parse_args(inOpts)

def main():
    argp = ArgParser()

    # create the ortholog dictionary
    id_dict = intronDict()
    id_dict.add_species(argp.args.tsv)

    ic = IntronComp(argp.args.bed1, argp.args.bed2, id_dict.get_entire_dict(), argp.args.std)

    ic.plot_graph(outfile_path=argp.args.o)

    return


if __name__ == '__main__':
    main()