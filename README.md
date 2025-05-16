import re
import doctest
kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}

aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

class Seq:

    def __init__(self,sequence,gene,species):
        self.sequence=sequence.strip().replace(" ","").upper() # This part cleans the sequence
        self.gene=gene
        self.species=species
        self.kmers = []



    def make_kmers(self,k=3):
        
        """
        >>> d = DNA("GCGC", "gene", "species", "id")
        >>> d.analysis()
        4
        """
    
        for i in range(0,len(self.sequence)-k + 1):
            self.kmers.append(self.sequence[i:i+k])
        return self.kmers
        

    def print_record(self):
        print(self.species + " " + self.gene + ": " + "\n" +self.sequence + "\nkmers:" + str(self.kmers))
    
    def fasta(self):
        print(">" + self.species + " " + self.gene+ "\n" +self.sequence )
        
    

class DNA(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence=re.sub('[^ATGCU]','N',self.sequence)
        self.geneid=geneid
        self.reverse_comp = None
        self.all_6_frames= None

 
    def analysis(self):

        """
        >>> d = DNA("GCGC", "gene", "species", "id")
        >>> d.analysis()
        4
        """
        
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return gc

       


    def print_info(self):
        print(" " + self.species + "\n" + self.gene + ": " + "\n" + self.geneid)

    
    def reverse_complement(self):

        """
        >>> d = DNA("ATGC", "gene", "species", "id")
        >>> d.reverse_complement()
        'GCAT'
        """

        complement="" # Empty string that will have the compliment added to
        for base in self.sequence:
            if base =="A":
                complement+="T" # If the base is A, it will add a T to the compliment empty string. The "+=" part adds T to the compliment string.
            elif base =="T":
                complement+="A"
            elif base =="G":
                complement+="C"
            elif base =="C":
                complement+="G"
                
        self.reverse_comp=complement[::-1]
        return self.reverse_comp




  
    def six_frames(self):

        """
        >>> d = DNA("ATGC", "gene", "species", "id")
        >>> d.reverse_complement()
        'GCAT'
        >>> d.six_frames()
        ['ATGC', 'TGC', 'GC', 'GCAT', 'CAT', 'AT']
        """


        frames=[]
        frames.append(self.sequence)
        frame2=self.sequence[1:]
        frames.append(frame2)
        frame3=self.sequence[2:]
        frames.append(frame3)
        frames.append(self.reverse_comp)
        reverse2=self.reverse_comp[1:]
        frames.append(reverse2)
        reverse3=self.reverse_comp[2:]
        frames.append(reverse3)
        return frames


class RNA(DNA):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species,geneid)
        self.sequence=sequence.strip().replace(" ","").upper()
        self.sequence=self.sequence.replace("T","U")
        self.sequence=re.sub('[^ATGCU]','N',self.sequence)
        self.geneid=geneid
        self.codons= []
     

    def print_info(self):
        print(" " + self.species + "\n" +self.sequence + "\n" + self.gene + ": " + "\n" + self.geneid + "\n" + str(self.codons))



    def make_codons(self):
        """
        >>> r = RNA("AUGGCC", "gene", "species", "id")
        >>> r.make_codons()
        >>> r.codons
        ['AUG', 'GCC']
        """
        for i in range(0,len(self.sequence),3):
            self.codons.append(self.sequence[i:i+3])



    def translate(self):

        """
        >>> r = RNA("AUGGCU", "gene", "species", "id")
        >>> r.make_codons()
        >>> r.translate()
        'MA'
        """
        protien=""
        for codon in self.codons:# checks if the codon in standard code is in self.codons.
            if codon in standard_code:# If the codon from standard code and self.codon match, add the AA (the value) to the empty string "protein"
                protien+=standard_code[codon]
        return protien



class Protein(Seq):

    def __init__(self,sequence,gene,species,protid,**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence=sequence.strip().replace(" ","").upper()
        self.sequence=re.sub('[^A-Z]','X',self.sequence)
        self.protid=protid

    def print_info(self):
        print(" " + self.species + "\n" +self.sequence + "\n" + self.gene + ": " + "\n" + self.protid)

    def total_hydro(self):
        """
        >>> p = Protein("ACDE", "gene", "species", "id")
        >>> round(p.total_hydro(), 1)
        -3.7
        """
        hydro=0 # makes the hydrophobicity equal to 0
        for seq in self.sequence: # loops through the sequence aa by aa
            if seq in kyte_doolittle: # checks to see if the aa is in the kyle.dolittle dictionary
                hydro+=kyte_doolittle[seq] # if the aa is infact in the dictionary, it grabs the associated value with kyle_dolittle[seq] which is then added to hydro. Thereon, the repitition of this process with the loop will calculate the total hydro.
        return hydro

    

    def mol_weight(self):

        """
        >>> p = Protein("ACDE", "gene", "species", "id")
        >>> round(p.mol_weight(), 2)
        490.47
        """
        molecular_weight=0
        for seq in self.sequence:
            if seq in aa_mol_weights:
                molecular_weight+=aa_mol_weights[seq]
        return molecular_weight





if __name__ == "__main__":
    import doctest
    doctest.testmod()

    myseq = Seq("  gATATAGGACctttaGGACCAC ", "my_gene", "H.sapiens")
    myseq.print_record()
    myseq.make_kmers(3)
    print(myseq.print_record())
    print(myseq.fasta())

    d = DNA(" -tcaaaGCGGCGGATCTCCCaaatga\n", "my_dna", "D.terebrans", "AX5667")
    d.print_info()
    rc = d.reverse_complement()
    print(rc)
    all_6_frames = d.six_frames()
    print(all_6_frames)

    r = RNA(" g?ATATAGGACctttaGGACCAC  ", "my_rna", "G.gallus", "R5990999")
    r.print_info()
    print(r)
    r.make_codons()
    print(r.codons)
    trans = r.translate()
    print(trans)

    p = Protein(" WCVALKKKCCYhhhhh-yyrsQ\t", "my_prot", "D.melanogaster", "56008009")
    p.print_info()
    hydro_value = p.total_hydro()
    print(hydro_value)
    prot_mol_weight = p.mol_weight()
    print(prot_mol_weight)

