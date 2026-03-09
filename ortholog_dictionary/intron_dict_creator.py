import re

class IntronTranscriptIDDictionary():

    REG_TRANSCRIPT_ID = r'ENS([A-Z]+)?T[0-9]{11}'
    # Group 1 represents the Species ID
    # Group 0 is the Ensembl Transcript ID

    def __init__(self):

        self.transcript_dict = {}
        # {MAIN SPECIES ID : {MAIN SPECIES ENSEMBL TRANSCRIPT ID : {OTHER SPECIES ID : ORTHOLOGOUS TRANSCRIPT ID}}}

    def add_species(self, tsv_file):

        with open(tsv_file, 'r') as tsv:
            # TSV file is formatted as:
            #   <GENOMENAME>-<GENEID>@<TRANSCRIPTID>-intron<NUMBER>(???) OtherSpecies,OtherSpecies,OtherSpecies,OtherSpecies
            for line in tsv:
                # first, isolate the main species intron id
                head = line.split("\t")[0]
                other_species = line.split("\t")[1].split(",")

                match = re.search(self.REG_TRANSCRIPT_ID, head)
                if match:
                    # a match was found
                    selfTID = str(match.group(0))    # Transcript ID of main species
                    if match.group(1):
                        selfSPID = str(match.group(1))   # Species Symbol of Main Species
                    else:
                        # Humans don't get symbols because they're so special 
                        selfSPID="HUM"

                    if selfSPID not in self.transcript_dict.keys():
                        self.transcript_dict[selfSPID] = {}

                    self.transcript_dict[selfSPID][selfTID] = {}

                    # add self to dictionary
                    self.transcript_dict[selfSPID][selfTID][selfSPID] = selfTID

                    # parse orthologs
                    for species in other_species:
                        match = re.search(self.REG_TRANSCRIPT_ID, species)
                        if match:
                            spTID = str(match.group(0))
                            spSPID = str(match.group(1))

                            self.transcript_dict[selfSPID][selfTID][spSPID] = spTID

    def get_ortholog_dict_of(self, speciesID):
        return self.transcript_dict[speciesID]

    def get_entire_dict(self):
        return self.transcript_dict
    

def main():
    tid = IntronTranscriptIDDictionary()
    tid.add_species('ortholog_dictionary/GRCh38_U12.tsv')
    #tid.add_species('FinalProject/ortholog_dictionary/GRCm38_U2.tsv')
    print(tid.get_ortholog_dict_of("HUM"))

if __name__ == '__main__':
    main()