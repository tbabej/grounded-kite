# generate all sequences for a particular design rule

# TODO pull design rule from the command
# ----------XXX-X--XX--XX-----
X = ["A", "V", "L", "F", "Y", "M", "S"]
SequenceTemplate = [ ["S"], ["A"], ["E"], ["E"], ["E"], ["K"], ["R"], ["K"], ["A"], ["E"], X, X, X, ["R"], X, ["A"], ["E"], X, X, ["K"], ["R"], X, X, ["E"], ["E"], ["E"], ["K"], ["W"] ]

# nested looping function for permuting the design rules
def IterateOverDesignSpace(currentSequenceFragment):

	# determine which position we're at
	positionIndex = len(currentSequenceFragment)

	# loop over the options at this level
	for aminoAcid in SequenceTemplate[positionIndex]:

		# add the next amino acid
		#print(currentSequenceFragment, positionIndex)
		thisSequenceFragment = currentSequenceFragment + aminoAcid

		# if we've reached end of the sequence, append it
		if len(thisSequenceFragment) == len(SequenceTemplate):
			Sequences.append(thisSequenceFragment)
			continue

		# otherwise, move into the next level of the loop
		else:
			IterateOverDesignSpace(thisSequenceFragment)

	# if we get here we're done with this level, so collapse it so we can move back up
	return True


# expand the design rule by nested looping
Sequences = []
IterateOverDesignSpace("")	


# output all sequences
output = ""
for sequence in Sequences:
	output += sequence + "\n"

outputFile = open("test_sequences.txt", "w")
outputFile.write(output)
outputFile.close()

# sanity-check
print(Sequences[0])
print(Sequences[-1])
print(len(Sequences))
