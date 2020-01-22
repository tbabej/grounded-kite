# main Twister script, which scores a list of input sequences

"""
Usage:
    TwistSequences.py <input_sequences_filepath> <lookup_table_filepath> <output_prefix>

Example:
    python3 TwistSequences.py Punam_30.txt Library/EnergyTable_IncGen0.1Helix.8Termini.polyAla output
"""

import json
import time
from docopt import docopt

# take a sequence and mask it to a fragment (e.g. KRTVEWYKR -> KR--V--KR)
def MaskSequenceToFragment(sequenceWindow, location, length):

    if location == "CENTRAL":
        fragment = sequenceWindow[:2] + "--" + sequenceWindow[4] + "--" + sequenceWindow[7:9]

    elif location == "NTERM":
        if length == 8:
            fragment = sequenceWindow[0] + "--" + sequenceWindow[3] + "--" + sequenceWindow[6:8]
        elif length == 7:
            fragment = "--" + sequenceWindow[2] + "--" + sequenceWindow[5:7]
        elif length == 6:
            fragment = "-" + sequenceWindow[1] + "--" + sequenceWindow[4:6]
        else:
            fragment = sequenceWindow[0] + "--" + sequenceWindow[3:5]

    elif location == "CTERM":
        if length == 8:
            fragment = sequenceWindow[:2] + "--" + sequenceWindow[4] + "--" + sequenceWindow[7]
        elif length == 7:
            fragment = sequenceWindow[:2] + "--" + sequenceWindow[4] + "--"
        elif length == 6:
            fragment = sequenceWindow[:2] + "--" + sequenceWindow[4] + "-"
        else:
            fragment = sequenceWindow[:2] + "--" + sequenceWindow[4]

    return fragment


# main function that scores a sequence
def ScoreSequence(sequence):
    
    rosettaEnergy = 0.0
    rosettaExtraEnergy = 0.0

    # slide over the N-terminal portion of the sequence
    for i in range(5,9):
        [fragmentRosettaEnergy, fragmentRosettaExtraEnergy] = EnergyTable[ MaskSequenceToFragment(sequence[:i], "NTERM", i) + "_NTERM_" + i ]
        rosettaEnergy += fragmentRosettaEnergy
        rosettaExtraEnergy += fragmentRosettaExtraEnergy

    # slide over the central part of the sequence
    for i in range(0,20):
        [fragmentRosettaEnergy, fragmentRosettaExtraEnergy] = EnergyTable[ MaskSequenceToFragment(sequence[i: i+9], "CENTRAL", 9) + "_CENTRAL_9"]
        rosettaEnergy += fragmentRosettaEnergy
        rosettaExtraEnergy += fragmentRosettaExtraEnergy

    # slide over the C-terminal portion of the sequence
    for i in range(20,24):
        [fragmentRosettaEnergy, fragmentRosettaExtraEnergy] = EnergyTable[ MaskSequenceToFragment(sequence[i:], "CTERM", 28-i) + "_CTERM_" + str(28-i)]
        rosettaEnergy += fragmentRosettaEnergy
        rosettaExtraEnergy += fragmentRosettaExtraEnergy

    # convert the energy to helicity once we have fit values for the tool on Gen0 at least
    rosettaHelicity = rosettaEnergy * LinearFitValues["Rosetta_Slope"] + LinearFitValues["Rosetta_Intercept"]
    rosettaExtraHelicity = rosettaExtraEnergy * LinearFitValues["Rosetta_Extra_Slope"] + LinearFitValues["Rosetta_Extra_Intercept"]

    return [rosettaEnergy, rosettaExtraEnergy, rosettaHelicity, rosettaExtraHelicity]


# define the linear fit values for getting helicity
LinearFitValues = {"Rosetta_Slope":-0.419199368287855, "Rosetta_Intercept":-2.66325889012837, "Rosetta_Extra_Slope":-0.710699066520611, "Rosetta_Extra_Intercept":-18.7554614780813}


################ MAIN ###############################    

if __name__ == "__main__":
    arguments = docopt(__doc__)

    # load the sequences
    IDs = []
    Sequences = []
    sequenceFile = open(arguments["<input_sequences_filepath>"])
    for i, line in enumerate(sequenceFile.readlines()):
        Entry = line.rstrip().split()
        if len(Entry) == 1:
            IDs.append(str(i))
            Sequences.append(Entry[0])
        else:
            IDs.append(Entry[0])
            Sequences.append(Entry[1])
    sequenceFile.close()


    # load the energy table
    EnergyTable = json.load( open(arguments["<lookup_table_filepath>"]) )

    # score the sequences
    Scores = []
    for sequence in Sequences:
        Scores.append(ScoreSequence(sequence))


    # output the scores
    #output = "ID\tSequence\tRosetta_Energy\tRosetta_Helicity\tRosetta_Extra_Energy\tRosetta_Extra_Helicity\n"
    Output = []
    for i in range(len(Scores)):
        #output += f"{i}\t{Sequences[i]}\t{Scores[i][0]:.3f}\t{Scores[i][1]:.2f}\t{Scores[i][2]:.3f}\t{Scores[i][3]:.2f}\n"
        Output.append( {"sequence":Sequences[i], "rosetta_energy":Scores[i][0], "rosetta_helicity":Scores[i][2], "rosetta_energy_extra":Scores[i][1], "rosetta_extra_helicity":Scores[i][3]} )

    with open(arguments["<output_prefix>"] + ".json") as fileHandle:
        json.dump(Output, fileHandle)
