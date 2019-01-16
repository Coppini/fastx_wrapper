#!/usr/bin/python3

import sys
import math
from hashlib import sha1
# from datetime import datetime

class Sequence(object):
    ''' This class is used for a sequence, determining it's header and sequence.'''
    def __init__(self, header, seq):
        self.header = header # Header of the sequence (it's "name")
        self.seq = seq # The sequence itself (nucleotides/aminoacids)
    
    def __repr__(self):
        return f"{self.header}\n{self.seq}" # Will return the header, a new line, and then the sequence. Notice that the > sign must be part of the header itself.


''' # This function may be used in order to load the sequences onto the memory in the proper format before starting the job
def load_fasta_file(path):
    #'#''Loads a fasta file onto the memory (usually a multi-fasta). Will then load each sequence as a Sequence object inside an array.
    Notice that newlines inside the sequence will be discarded, so that each sequence has only two lines:
        HEADER [the greater sign (">") is part of the Header]
        SEQUENCE [Everything that is not in the header's line, and before the next header]
    As everything is loaded onto the memory, computers with low memory may have trouble with big files.'#''
    contents = open(path).read().split("\n") # Opens the file and loads it onto the memory
    sequences = [] 
    current_sequence = Sequence("", "")
    for line in contents: ### Parses through each line of the file loaded onto the memory, determining if they are Headers or Sequences, and loading them as such
        if line.startswith(">"): # If the line starts with ">", consider it a header
            sequences.append(current_sequence)
            #'#''A header is the start of a new sequence, so before starting a new one, we first must append the previous sequence to our array.
            Notice that this will generate an empty sequence ("","") at the first position in the array [0]. We will remove that later.'#''
            current_sequence = Sequence("","") # Clearing the variable to load the next sequence in
            current_sequence.header = line # Assigns the line as a header
        else: # If the line doesn't start with ">", considers it as part of the sequence
            current_sequence.seq += line # Adding the sequence to any previous lines that weren't headers (in the case of newlines inside the sequence itself)
    sequences.append(current_sequence)
    sequences[:1]="" # Removes the first sequence of the array, the empty sequence ("","") mentioned before.
    contents = "" # Clears the variable from the memory. No need to keep our memory full.
    return sequences '''


def norepeat_splitter(bigseq, min_len):
    '''This function splits fasta sequences (bigseq) into smaller sequences with a minimum defined length (min_len),
    without repeating sequences (with no overlaps between each generated sequence).
    Notice that the original sequences are split into smaller sequences with AT LEAST min_len, so if the original 
    sequence (bigseq) is not a multiple of the min_len for the splitted sequences, the function will automatically
    increase the length to make similar sized sequences. If the length of the bigseq is close to min_len, this can
    yield splitted sequences with a maximum of 2x the given min_len.'''
    norepeat_split_sequences = []
    nseqs = int(len(bigseq.seq) / min_len) # Defines the number of smaller sequences generated from the given sequence (to have at least 'length' bases)
    if nseqs == 0:
        nseqs = 1
        min_len = len(bigseq.seq)
    remainder = len(bigseq.seq) % min_len # Calculates how many bases would remain after splitting into 'nseqs' with 'length' bases
    while remainder / nseqs >= 1: ### While the remainder is bigger than the number of split sequences to generate, add bases to the length of every sequence
        remainder = len(bigseq.seq) % min_len
        min_len = int(min_len + ( remainder / nseqs ))  ### Adds the remainders to every sequence, to even them out. Notice that there still may be some remainder
    end = 0
    for i in range(nseqs - remainder): ### Creates sequences with the min_len defined by the previous loop
        start = end # At which point the newly generated sequence will start
        end = start + min_len # And where it will end
        currentseq = Sequence(
            f"{bigseq.header}_{i}_[{start}:{end}]", # Adds to the header of the sequence a number, as well as the start and end base in the original sequence
            bigseq.seq[start:end]
        )
        print(currentseq)
        ''' # To use instead of print if you want to hold it in the memory instead of printing it as soon as the sequence is processed.
        norepeat_split_sequences.append(currentseq) '''
    for i in range(nseqs - remainder , nseqs): # Adds an extra base to the last sequences generated, in order to include every sequence of the original bigseq
        start = end
        end = start + min_len + 1 # Adds 1 to the end, to include the remainders
        currentseq = Sequence(
            f"{bigseq.header}_{i}_[{start}:{end}]",
            bigseq.seq[start:end]
        )
        print(currentseq)
        ''' # To use instead of print if you want to hold it in the memory instead of printing it as soon as the sequence is processed.
        norepeat_split_sequences.append(currentseq) '''
    total = [s.seq for s in norepeat_split_sequences] # This will contain the entire sequence, in order to compare with the bigseq
    assert "".join(total) == bigseq.seq , f"Error. Split sequences are not equal to original sequence for {bigseq.header}." # The comparison asserts every base is accounted for
    return norepeat_split_sequences

def overlap_splitter(bigseq, fixed_len, overlap=0, coverage=0):
    '''This function splits fasta sequences (bigseq) into smaller sequences of a fixed length (fixed_len),
    making overlaps between each generated sequence (overlap). The user can give a specific amount of bases that
    will overlap between one generated sequence and the next; or decide a coverage for the sequences, so that
    the overlaps will be automatically calculated in order to generate the desired coverage.'''
    index=set()
    totalbpcount = len(bigseq.seq) # Calculates the total size of the sequence
    if coverage > 1:
        if coverage >= fixed_len: # In case the coverage is too big, fixes the 
            coverage = fixed_len - 0.0001
            gap = 1
            overlap = fixed_len - gap
        else:
            overlap = fixed_len - int(fixed_len / coverage)
    elif overlap > 0:
        if overlap >= fixed_len:
            overlap = fixed_len - 1
            coverage = fixed_len - 0.001
        else:
            coverage = fixed_len / (fixed_len - overlap)
    elif overlap <= 0 and coverage <= 1:
        overlap = 0
        coverage = 1
    overlap_split_sequences = []
    gap = fixed_len - overlap
    if fixed_len >= totalbpcount:
        currentseq = Sequence(
            f"{bigseq.header}_{1}_[0:{totalbpcount}]",
            bigseq.seq
        )
    else:
        start = 0
        end = fixed_len
        i=1
        Nseqs = math.ceil((totalbpcount - fixed_len) / gap)+1
        restofrest = ((totalbpcount - fixed_len) % (Nseqs-1))
        fullgap=(totalbpcount - fixed_len)/(Nseqs - 1)
        while end < totalbpcount:
            while (i <= (restofrest + 1)) and (end < totalbpcount):
                currentseq = Sequence(
                    f"{bigseq.header}_{i}_[{start}:{end}]", # Adds to the header of the sequence a number, as well as the start and end base in the original sequence
                    bigseq.seq[start:end]
                )
                if sha1(currentseq.seq.encode('utf-8')).hexdigest() not in index: ### Avoids redundancy by testing if an identical sequence has already been generated
                    print(currentseq)
                    ''' # To use instead of print if you want to hold it in the memory instead of printing it as soon as the sequence is processed.
                    overlap_split_sequences.append(currentseq) '''
                    index.add(sha1(currentseq.seq.encode('utf-8')).hexdigest())
                start = start + math.ceil((totalbpcount - fixed_len)/(Nseqs - 1))
                end = start + fixed_len
                i+=1
            currentseq = Sequence(
                f"{bigseq.header}_{i}_[{start}:{end}]", # Adds to the header of the sequence a number, as well as the start and end base in the original sequence
                bigseq.seq[start:end]
            )
            if sha1(currentseq.seq.encode('utf-8')).hexdigest() not in index: ### Avoids redundancy by testing if an identical sequence has already been generated
                print(currentseq)
                ''' # To use instead of print if you want to hold it in the memory instead of printing it as soon as the sequence is processed.
                overlap_split_sequences.append(currentseq) '''
                index.add(sha1(currentseq.seq.encode('utf-8')).hexdigest())
            start = start + int((totalbpcount - fixed_len)/(Nseqs - 1))
            end = start + fixed_len
            i+=1
        currentseq = Sequence(
            f"{bigseq.header}_{i}_[{totalbpcount - fixed_len}:{totalbpcount}]", # Adds to the header of the sequence a number, as well as the start and end base in the original sequence
            bigseq.seq[start:end]
        )
    if sha1(currentseq.seq.encode('utf-8')).hexdigest() not in index: ### Avoids redundancy by testing if an identical sequence has already been generated
        print(currentseq)
        ''' # To use instead of print if you want to hold it in the memory instead of printing it as soon as the sequence is processed.
        overlap_split_sequences.append(currentseq) '''
        index.add(sha1(currentseq.seq.encode('utf-8')).hexdigest())

    ''' # To use instead of print if you want to hold it in the memory instead of printing it as soon as the sequence is processed.
    return overlap_split_sequences '''

if __name__ == "__main__":
#    start=datetime.now()
    coverage = 0
    overlap = 0
    length = 0
    sequences = []
    if len(sys.argv) > 2:
        length = int(sys.argv[2]) # The desired length for your split sequences (min_len)

        if len(sys.argv) > 3:       # Tests if an overlap or coverage value was given
            if sys.argv[3].startswith("-o="):
                overlap = sys.argv[3][3:]
                coverage = 0
            elif sys.argv[3].startswith("--o="):
                overlap = sys.argv[3][4:]
                coverage = 0
            elif sys.argv[3].startswith("-overlap="):
                overlap = sys.argv[3][9:]
                coverage = 0
            elif sys.argv[3].startswith("--overlap="):
                overlap = sys.argv[3][10:]
                coverage = 0
            elif sys.argv[3].startswith("-c="):
                overlap = 0
                coverage = sys.argv[3][3:]
            elif sys.argv[3].startswith("--c="):
                overlap = 0
                coverage = sys.argv[3][4:]
            elif sys.argv[3].startswith("--cov="):
                overlap = 0
                coverage = sys.argv[3][6:]
            elif sys.argv[3].startswith("-cov="):
                overlap = 0
                coverage = sys.argv[3][5:]
            elif sys.argv[3].startswith("-coverage="):
                overlap = 0
                coverage = sys.argv[3][10:]
            elif sys.argv[3].startswith("--coverage="):
                overlap = 0
                coverage = sys.argv[3][11:]
            else:
                overlap = sys.argv[3]
            if len(sys.argv) > 4: 
                if sys.argv[4].startswith("-o="):
                    overlap = sys.argv[4][3:]
                elif sys.argv[4].startswith("--o="):
                    overlap = sys.argv[4][4:]
                elif sys.argv[4].startswith("-overlap="):
                    overlap = sys.argv[4][9:]
                elif sys.argv[4].startswith("--overlap="):
                    overlap = sys.argv[4][10:]
                elif sys.argv[4].startswith("-c="):
                    coverage = sys.argv[4][3:]
                elif sys.argv[4].startswith("--c="):
                    coverage = sys.argv[4][4:]
                elif sys.argv[4].startswith("--cov="):
                    coverage = sys.argv[4][6:]
                elif sys.argv[4].startswith("-cov="):
                    coverage = sys.argv[4][5:]
                elif sys.argv[4].startswith("-coverage="):
                    coverage = sys.argv[4][10:]
                elif sys.argv[4].startswith("--coverage="):
                    coverage = sys.argv[4][11:]
                else:
                    coverage = sys.argv[4]
        else:
            overlap = 0
            coverage = 0
    else:
        print("You need to provide a `Fasta File` and a `Length` to work with.\nUsage:\n\tpython3 Splitter.py <Fasta File> <Length> <Overlap> <Coverage>")
        sys.exit()

    overlap=int(overlap)
    coverage=float(coverage)

    # print(f"Length={length},\nOverlap={overlap},\nCoverage={coverage})")

#'#'#'# If you want to withhold sequences from output until all of them are processed, replace the section of code directly below with the next one.####
    contents = open(sys.argv[1]).read().split("\n") # Opens the file and loads it onto the memory

    sequences = [] 
    current_sequence = Sequence("", "")
    for line in contents: ### Parses through each line of the file loaded onto the memory, determining if they are Headers or Sequences, and loading them as such
        if line.startswith(">"): # If the line starts with ">", consider it a header
            if current_sequence.header and current_sequence.seq: # If there was a sequence already loaded
                overlap_splitter(current_sequence, length, overlap, coverage)
            current_sequence = Sequence("","") # Clearing the variable to load the next sequence in
            current_sequence.header = line # Assigns the line as a header
        else: # If the line doesn't start with ">", considers it as part of the sequence
            current_sequence.seq += line # Adding the sequence to any previous lines that weren't headers (in the case of newlines inside the sequence itself)
    overlap_splitter(current_sequence, length, overlap, coverage)   
    contents = ""
#### If you want to withhold sequences from output until all of them are processed, replace the section of code directly above with the next (below) one.#'#'#'#

    ''' # Uses extra memory. Not recommended. Use this if you're using the load_fasta_file function.
    This will withhold the sequences until all of them are processed, before sending them to the standard output. You may use this in a program if you want
    all the sequences to be written to output at the same time, instead of one by one as they are processed.
    # sequences = load_fasta_file(sys.argv[1]) # Loads your fasta file (see function load_fasta_file), with the original sequences (bigseq's), into an array of Sequences. # 
    for curseq in sequences:
        # splitted = norepeat_splitter(curseq,length) # Splits one of the original sequences, making an array with multiple smaller sequences
        splitted = overlap_splitter(curseq,length,overlap,coverage)
        for s in splitted:
            print(s) # Each split sequence in the generated array is printed
        splitted=[] '''
#    print(datetime.now() - start)
