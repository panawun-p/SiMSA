#!/usr/bin/env python
import sys
import os.path
from os import system, remove
from itertools import combinations
from time import strftime, ctime
from collections import Counter

#####################################################################
#####                                                           #####
#####                 GLOBALS and CONFIGURATION                 #####
#####                                                           #####
#####################################################################
# @LOG_PATH[string]: full path to folder containing the log file
# @LOG_FILE[file]: file object containing log file information
# @FASTA_EXTENSION[list]: list of input sequence extension that the program accepts
# @SIMILARITY_MATRIX[list]: list of similarity matrix available in the program
# @PROTEIN_GROUPS[dictionary]: dictionary of protein groups derived from similarity matrix

LOG_PATH = 'log/'
LOG_FILE = None
FASTA_EXTENSION = ['.fasta', '.fa', '.fas', '.fsa', '.seq']
SIMILARITY_MATRIX = ['NONE','PAM60','PAM250','BLOSUM40','BLOSUM62','BLOSUM80']
PROTEIN_GROUPS = {'PAM250': \
	            [('?'), ('M','I','V','L'), ('D','N','H','Q','E'), ('F','I','L'), ('S','P','A'), \
	            ('S','A','G'), ('Q','K','N'), ('R','H','Q'), ('S','T'), ('R','K','Q'), ('S','N'), \
	            ('F','W'), ('R','W'), ('G','D')], \
	            'PAM60': \
	            [('?'), ('S','A','T'), ('H','N'), ('S','N'), ('N','D'), ('D','E'), ('H','Q'), \
	            ('Q','E'), ('Y','F'), ('R','K'), ('M','I'), ('L','M'), ('I','V')], \
	            'BLOSUM40': \
	            [('?'), ('S','T'), ('S','A'), ('S','Q','N'), ('H','Y','M'), ('D','E'), ('N','D'), \
	            ('H','N'), ('W','Y','F'), ('E','Q','K'), ('K','R'), ('R','Q'), ('L','M','I','V'), \
	            ('E','K'), ('A','G'), ('L','F','I'), ('V','T')], \
	            'BLOSUM62': \
	            [('?'), ('S','T'), ('S','A'), ('S','N'), ('H','Y'), ('D','E'), ('N','D'), ('H','N'), \
	            ('W','Y','F'), ('E','Q','K'), ('K','R','Q'),('L','M','I','V')], \
	            'BLOSUM80': \
	            [('?'), ('Q','R','K'), ('Q','E','K'), ('E','D'), ('D','N'), ('Q','H'), ('Y','H'), \
	            ('Y','W'), ('Y','F'), ('S','T'), ('S','A'), ('M','V','L','I')]}

# Create log folder if not exists
if not os.path.exists(LOG_PATH):
	system("mkdir -p " + LOG_PATH)
#####################################################################



######################## FUNCTION log_write #########################
# write a message to the log file with timestamp
# ============================== INPUT ==============================
# @message[string]: message to write to log file
# ============================= OUTPUT ==============================
# function write a line of message to the log file
#####################################################################
def log_write(message):
	LOG_FILE.write('[' + strftime("%Y-%m-%d %H:%M:%S") + '] ' + message + '\n')

######################## FUNCTION terminate #########################
# terminate SiMSA software
#####################################################################
def terminate():
	log_write('Program terminated.')
	sys.exit()

####################### FUNCTION check_input ########################
# check input file format
# ============================= OUTPUT ==============================
# @input_name[string]: the input name as specified in the command arguments
#####################################################################
def check_input():
	# ======================= LOCAL VARIABLES =======================
	# @extension[string]: temporary to hold each file extension in @FASTA_EXTENSION
	# ===============================================================
	if '-i' in sys.argv:
		# If input argument is specified
		# Retrieve input name from the argument
		input_name = sys.argv[sys.argv.index('-i')+1]

		# Check if input file extension is acceptible
		if not any(extension in input_name for extension in FASTA_EXTENSION):
			# If the extension is not specified in the input name
			print 'Input format is not correct.'
			# Write log and terminate program
			log_write('Input ERROR: Input format is not correct.')
			terminate()
	else:
		# input name is not specified, and the program cannot be run
		print 'Please specify input file.'
		print 'For help plase call SiMSA with "-h" after filename.'
		# Write log and terminate the program
		log_write('Argument ERROR: Input file name not specified.')
		terminate()
	# If everything is okay, then return the specified input name
	return input_name

####################### FUNCTION check_output #######################
# check if the output file with user-specified name already exists
# ============================= OUTPUT ==============================
# @output_name[string]: the output name as specified in the command arguments
#####################################################################
# check and assign output name
def check_output():
	if '-o' in sys.argv:
		# If output name argument is specified
		# Get output name from the argument
		output_name = sys.argv[sys.argv.index('-o')+1]

		# Check if the output file with that name already exist or not
		if os.path.isfile(output_name+'_rpt.txt') or os.path.isfile(output_name+'_mat.txt'):
			# If file already exist
			# =================== LOCAL VARIABLES ===================
			# @overwrite[string]: contain answer of the user whether to
			#					  overwrite the output file or not
			# =======================================================
			# Prompt to confirm whether to overwrite the file
			overwrite = raw_input('This file already exists. Overwrite? (y/n) : ')
			# If the user does not agree to overwrite
			if overwrite.lower() not in ['y','yes']:
				# Write log and terminate
				log_write('Argument ERROR: Output file of that name already exists.')
				terminate()
	else:
		# If the output name is not specified,
		# Auto-generate using current timestamp
		output_name = 'report[' + strftime("%Y-%m-%d %H-%M") + ']'

	# Write log and return output name
	log_write('Output files: ' + output_name)
	return output_name

######################## FUNCTION read_fasta ########################
# read sequences from the input fasta file (taken from BioPython)
# ============================== INPUT ==============================
# @input_file_obj[file]: contain information of the input fasta file
# ============================= OUTPUT ==============================
# function return the name of the sequence and the sequence itself
#####################################################################
def read_fasta(input_file_obj):
	# ======================= LOCAL VARIABLES =======================
	# @name[string]: contain sequence name as specify in the input file
	# @sequence[string]: contain amino acid or nucleotide sequence
	# ===============================================================

	# Initialize @name and @sequence
    name, sequence = None, []

    # Loop through each line of the input file and read data
    for line in input_file_obj:
		# ======================= LOCAL VARIABLES =======================
		# @line[string]: a temporary contain each line in the file
		# ===============================================================    	

		# Strip out all whitespace characters
        line = line.rstrip()
        # Check whether the line is header line or not
        if line.startswith(">"):
        	# If this line is header line
        	# Return the sequence in case that name has already been declared
            if name: yield (name, ''.join(sequence))
            # If name is not yet declared, that use this line as the name
            # And initialize the list to contain sequence
            name, sequence = line, []
        else:
        	# If the line is not header line - it contain sequence information
        	# Append the line to sequence list
            sequence.append(line)
    # Return the last sequence in the file
    if name: yield (name, ''.join(sequence))

##################### FUNCTION fetch_sequences ######################
# get sequence from the given input file name
# ============================== INPUT ==============================
# @input_name[string]: name of the input file
# ============================= OUTPUT ==============================
# @input_content[list]: a list of dictionary with name and sequence information
#						the dictionary contain following key:
#						'name': name of the sequence
#						'seq': protein or DNA sequence
#####################################################################
def fetch_sequences(input_name):
	# ======================= LOCAL VARIABLES =======================
	# @seq_count[int]: use to count number of sequence in the file
	# ===============================================================

	# Initiate important local variables
	input_content = [] # list of dictionary | dictionary contain 'name' and 'seq'
	seq_count = 0

	# Print console message and write log
	print 'Reading file: ' + input_name
	log_write('Reading file: ' + input_name)

	# Open input fasta file and parse its content
	with open(input_name) as input_file_obj:
		# ===================== LOCAL VARIABLES =====================
		# @input_file_obj[file]: contain information of the input fasta file
		# ===========================================================
	    for name, seq in read_fasta(fp):
	    	# =================== LOCAL VARIABLES ===================
			# @name[string]: a temp for name of each sequence in input file
			# @seq[string]: a sequence for each sequence in input file
			# @input_data[dictionary]: a temporary dict to be added to @input_content
			# =======================================================

			# Construct tmp dictionary and eliminate '>' from sequence name
	    	input_data = {'name':name.translate(None, '>'), 'seq':seq}
	    	# Append sequence information to @input_content
	        input_content.append(input_data)
	        # Update number of sequences
	        seq_count += 1

	# If in total there's 0 or 1 sequence, then terminate program
	# Since the similarity and identity cannot be calculated
	if seq_count < 2:
		# Write message to console and log file, then exit program
		print '2 or more sequences are required.'
		log_write('Input ERROR: 2 or more sequences are required.')
		terminate()

	# If the read is complete, write to console and log
	print 'Read completed.'
	log_write('Read completed.')

	# Then return all the input sequences
	return input_content

###################### FUNCTION check_sequence ######################
# check the input sequences for its type (protein or DNA) and also
# for its length
# ============================== INPUT ==============================
# @input_content[list]: a list of dictionary with name and sequence information
#						the dictionary contain following key:
#						'name': name of the sequence
#						'seq': protein or DNA sequence
# ============================= OUTPUT ==============================
# function return the type of sequence. 'P' for protein and 'N' for DNA
#####################################################################
def check_sequence(input_content):
	# ======================= LOCAL VARIABLES =======================
	# @types[list]: contain sequence types of all input sequences
	#				the value in each can be 'P' or 'N'
	# ===============================================================

	# Initiate @types to contain the result of sequence type checking
	types = []

	# Print progress message to console and log file
	print 'Checking sequences.'
	log_write('Checking sequences.')

	# Check that every input sequences have the same length
	# Then check the sequence type
	if all( len(i['seq']) == len(input_content[0]['seq']) for i in input_content ):
		# If all sequences has the same length
		# Check sequence type of all sequence
		for each_input in input_content:
			# =================== LOCAL VARIABLES ===================
			# @each_input[dictionary]: contain sequence informations
			#						   'name': name of the sequence
			#						   'seq': protein or DNA sequence
			# =======================================================

			# Check the sequence types of sequence in current iteration
			if any(c.upper() in 'EFILPQZ' for c in each_input['seq']):
				# sequence is protein if contain any of 'EFILPQZ'
				types.append('P')
			else:
				# recheck the sequence to ensure it really is DNA
				# ================= LOCAL VARIABLES =================
				# @base_count[list]: contain 2 most common bases in the sequence
				# ===================================================

				# Get 2 most common base from the sequence
				base_count = Counter(each_input['seq'].upper()).most_common(2)

				# Recheck the sequence type using the most common base in @base_count
				if all(x[0] in 'ATCGN-.' for x in base_count):
					# base ATCG is majority of the sequence
					types.append('N')
				else:
					types.append('P')

		# Check if every input sequences is of the same type
		if all(each == types[0] for each in types):
			# If every sequences is the same type
			# Write message to console and log file
			print 'Sequence check completed.'
			log_write('Sequence check completed.')
			# Then return the sequence type as single character
			return types[0]
		else:
			# If not all sequence is the same type
			# Write error message to console and log file
			print 'All inputs must consist of same kind of sequences.'
			log_write('Input ERROR: All inputs must consist of same kind of sequences.')
			# Then exit the program
			terminate()
	else:
		# If not all sequences have the same length
		# Write error message to console and log file
		print 'All inputs must have the same length.'
		log_write('Input ERROR: All inputs must have the same length.')
		# Then exit the program
		terminate()

######################### FUNCTION pairing ##########################
# generate a list of sequence pair combination
# ============================== INPUT ==============================
# @input_content[list]: a list of dictionary with name and sequence information
# ============================= OUTPUT ==============================
# @pair[list]: list of dictionary contain 'pair' and 'length' where:
#			   'pair': list of 2 sequence index
#			   'length': length of the sequence after check the 
#####################################################################
def pairing(input_content):
	# ======================= LOCAL VARIABLES =======================
	# @all_combination[list]: list of all combinations with no repeat
	# ===============================================================

	# Initiate @pair list
	pair = []

	# generate combination with no repeat
	all_combination = list(combinations(range(len(input_content)), 2))

	# Loop through each combination
	for each_pair in all_combination:
		# ===================== LOCAL VARIABLES =====================
		# @each_pair[tuple]: a temporary for a combination in @all_combination
		# @tmp_dict[dictionary]: a temporary dictionary to hold data for each combination
		# ===========================================================

		# Initialize pair dictionary with the pair and it's total length
		# Initialize 'length' as 0, the value will be assigned after base pruning
		tmp_dict = {'pair':each_pair, 'length':0} 

		# Append new pair to @pair
		pair.append(tmp_dict)

	return pair

#################### FUNCTION get_matrix_choice #####################
# get similarity matrix choice from the user through console
# ============================= OUTPUT ==============================
# @matrix_choice[string]: matrix choice that the user selected
#####################################################################
def get_matrix_choice():
	# Print prompt message and options to console
	print '\nChoose protein grouping choice.\n'

	# Print a line for each similarity matrix with number label
	for i in range(len(SIMILARITY_MATRIX)):
		# ===================== LOCAL VARIABLES =====================
		# @i[int]: temporary to hold index for each similarity matrix option
		# ===========================================================
		if i == 0:
			# Print full text for the first choice
			print 'Enter ' + str(i+1) + ' for ' + SIMILARITY_MATRIX[i] + '\n'
		else:
			# Left out some text for other choices
			print '      ' + str(i+1) + '     ' + SIMILARITY_MATRIX[i] + '\n'

	# Ask for the input from user
	while True: 
		# While keyboard input is not valid, continue to ask again
		# ===================== LOCAL VARIABLES =====================
		# @choice_index[int]: the index of the choice user choose
		# ===========================================================
		
		# Print the prompt on the console
		choice_index = int(raw_input('Your choice: '))

		# Check if user's input is valid or not
		if choice_index in range(1, len(SIMILARITY_MATRIX) + 1):
			# If user's input is valid
			# Then get the matrix choice from @SIMILARITY_MATRIX list
			matrix_choice = SIMILARITY_MATRIX[choice_index - 1]

			# Then print message on the console and log file
			print 'Using ' + matrix_choice + ' for protein grouping.'
			log_write('Using ' + matrix_choice + ' for protein grouping.')
			# Then break the loop to return the value
			break
		else:
			# If user's input is invalid, then print the error message
			# And continue to ask for the valid input
			print 'Choice not valid. Please enter your choice again.'

	return matrix_choice

###################### FUNCTION fetch_grouping ######################
# fetch protein groups from the specified matrix choice
# ============================== INPUT ==============================
# @matrix_choice[string]: matrix choice that the user selected
# ============================= OUTPUT ==============================
# @p_groups[list]: each member contain a tuple of proteins in the same group
#####################################################################
def fetch_grouping(matrix_choice):
	# If user does not wish to calculate the similarity
	if matrix_choice == 'NONE':
		return 'NONE'

	# Else return the list contain protein groups accordingly
	p_groups = PROTEIN_GROUPS[matrix_choice]

	return p_groups

###################### FUNCTION prune_sequence ######################
# prune some ambiguous base positions from the sequence, since we
# cannot calculate similarity and identity for this base positions
# ============================== INPUT ==============================
# @seq1[string]: one of the sequence in the pair
# @seq2[string]: another sequence in the pair to be calculate with @seq1
# @seq_type[string]: sequence type, the value is either 'P' or 'N'
# ============================= OUTPUT ==============================
# @pruned_index[list]: list of sequence index to be disregard from the calculation
#####################################################################
def prune_sequence(seq1, seq2, seq_type):
	# ======================= LOCAL VARIABLES =======================
	# @seq_length[int]: length of the sequence before pruning
	# ===============================================================

	# Get the length of the sequence
	seq_length = len(seq1)
	# Initialize a list to contain positions to prune out
	pruned_index = []
	
	for i in range(len(seq1)):
		# ===================== LOCAL VARIABLES =====================
		# @i[int]: a temporary to keep track of sequence position
		# ===========================================================
		if (seq1[i] in '-.?' and seq2[i] in '-.?') or (seq_type == 'P' and 'X?' in [seq1[i], seq2[i]]) or (seq_type == 'N' and 'N?' in [seq1[i], seq2[i]]):
			# prune if position is (gap with gap, ? with ?) 
			# or (one protein sequence contain X or ?) 
			# or (one DNA sequence contain N or ?)
			pruned_index.append(i)

	return pruned_index

#################### FUNCTION parewise_identity #####################
# calculate pairwise identity between 2 sequences
# ============================== INPUT ==============================
# @seq1[string]: one of the sequence in the pair
# @seq2[string]: another sequence in the pair to be calculate with @seq1
# @prune_index[list]: list of sequence index to ignore in the calculation
# ============================= OUTPUT ==============================
# @percent_identity[float]: identity value of the sequence pair in percentage
#####################################################################
def parewise_identity(seq1, seq2, prune_index):
	# ======================= LOCAL VARIABLES =======================
	# @seq_length[int]: length of the sequence before pruning
	# @matched_position[int]: count of the position with identical base
	# ===============================================================

	# if seq1 and seq2 is exactly the same, then identity equals 100%
	if seq1 == seq2:
		return 100

	# Get the length of sequence and initiate @matched_position
	seq_length = len(seq1)
	matched_position = 0

	# Loop through each position in the sequence and see if it is match
	for i in range(seq_length):
		# ===================== LOCAL VARIABLES =====================
		# @i[int]: a temp to keep track of position in the sequence
		# ===========================================================
		if i not in prune_index and seq1[i] == seq2[i]:
			# If the position is not to be pruned and the base is identical
			# Then add it to the matched position
			matched_position += 1

	# Calculate the identity
	percent_identity = float(matched_position)/float(len(seq1)-len(prune_index))*100.00

	return percent_identity

################### FUNCTION pairwise_similarity ####################
# calculate pairwise identity between 2 sequences
# ============================== INPUT ==============================
# @seq1[string]: one of the sequence in the pair
# @seq2[string]: another sequence in the pair to be calculate with @seq1
# @p_groups[string]: list of protein groups from the selected similarity matrix
# @prune_index[list]: list of sequence index to ignore in the calculation
# ============================= OUTPUT ==============================
# @percent_similarity[float]: similarity value of the sequence pair in percentage
#####################################################################
def pairwise_similarity(seq1, seq2, p_groups, prune_index):
	# ======================= LOCAL VARIABLES =======================
	# @seq_length[int]: length of the sequence before pruning
	# @matched_position[int]: count of the position with identical base
	# ===============================================================

	# if seq1 and seq2 is exactly the same, then identity equals 100%
	if seq1 == seq2:
		return 100

	# Get the length of sequence and initiate @matched_position
	seq_length = len(seq1)
	matched_position = 0

	# Loop through each position in the sequence and see if it is match
	for i in range(seq_length):
		# ===================== LOCAL VARIABLES =====================
		# @i[int]: a temp to keep track of position in the sequence
		# ===========================================================
		if i not in prune_index and (any(set([seq1[i].upper(), seq2[i].upper()]).issubset(group) for group in p_groups) or seq1[i] == seq2[i]):
			# If the position is not to be pruned and (both are in the same protein group, or is exactly the same)
			# Then add it to the matched position
			matched_position += 1

	# Calculate the similarity
	percent_similarity = float(matched_position)/float(len(seq1)-len(prune_index))*100.00

	return percent_similarity

######################## FUNCTION statistic #########################
# calculate statistical values: max, min, median, mode, mean for any given @value_list
# ============================== INPUT ==============================
# @value_list[list]: list of value to get statistics from
# @middle_value[boolean]: whether to calculate median, mode, mean, or not
# @non_percentage[boolean]: whether the value in @value_list is percentage value
# ============================= OUTPUT ==============================
# @return_stats[dictionary]: stored the calculated statistic
#####################################################################
def statistic(value_list, middle_value = True, non_percentage = False):
	# ======================= LOCAL VARIABLES =======================
	# @tmp[list]: temporary variable to hold the list of values
	# @most_common[list]: list of tuples, contain 1 most common value in @value_list
	#					  the tuple contains the value at [0] and number of occurance at [1]
	# ===============================================================

	# Return the function with empty value if there's nothing to calculate from
    if value_list == []:
        return None

    # Initiate @return_stats dictionary
    return_stats = {}

    # Sort input value list for mode calculation
    tmp = [] + value_list
    tmp.sort()

    # Calculate median, mode, and mean if @middle_value is True
    if middle_value:
        # Calculate MEAN
        return_stats['Mean'] = sum(tmp)/len(tmp)
        
        # Calculate MEDIAN
        # Check if the value to calculate have odd or even number of members
        if len(tmp) % 2:
            # If the list of value have odd number of members
            if len(tmp) == 1:
            	# If there's only 1 member in the list
                return_stats['Median'] = tmp[0]
            else:
            	# If there's 3, 5, 7, ... members in the list
                return_stats['Median'] = tmp[(len(tmp)/2) + 1]
        else:
            # If the list of value have even number of members
            return_stats['Median'] = (tmp[len(tmp)/2] + tmp[len(tmp)/2 + 1])/2.0

        # Get MODE
        # Find the most common member in list
        most_common = Counter(tmp).most_common(1)
        # Check if there is most common member or not
        if most_common[0][1] == 1:
        	# If there is no most common member, there's no mode
            return_stats['Mode'] = 'NONE'
        else:
        	# If there is most common member, return the value with number of occurances
            return_stats['Mode'] = str(round(most_common[0][0], 2)) + '%  with ' + str(most_common[0][1]) + ' occurances.'

    # Get MAX and MIN
    return_stats['Max'] = max(tmp)
    return_stats['Min'] = min(tmp)

    # If the input value is suppose to be percentage
    if not non_percentage:
        # Iterate through the calculated statistics and format
        for key, value in return_stats.iteritems():
        	# =================== LOCAL VARIABLES ===================
			# @key[string]: a temp for the key to each @return_stats
			# @value[string/float]: a temp for value i each @return_stats
			# =======================================================
            if key != 'Mode':
            	# Format each value to string in [xx.xx%] format
                return_stats[key] = str(round(value, 2)) + '%'

    return return_stats

####################### FUNCTION init_matrix ########################
# initiate the empty matrix to be written into output file
# ============================== INPUT ==============================
# @header[list]: list of headers to each row and column of the matrix
# ============================= OUTPUT ==============================
# @matrix[list]: a 2 dimensional list, each member contain a row of values
#####################################################################
def init_matrix(header):
	# ======================= LOCAL VARIABLES =======================
	# @total_line[int]: total number of lines this matrix should have
	# ===============================================================

	# Initiate the matrix
    matrix = []
    # Calculate total line include header line, seperator line, and value lines
    total_line = len(header) + 2

    # Initialize the matrix with header and filled value cells with '.'
    for line_count in range(total_line):
    	# ===================== LOCAL VARIABLES =====================
		# @line_count[int]: a temp to locate the current row in matrix
		# ===========================================================

        # Construct each line of the matrix
        if line_count == 0:
            # If this is the first line of matrix
            # Then it is the header line
            tmp = ['|'] + header

        elif line_count == 1:
            # If this is the second line of matrix
            # Then it is the seperator line
            tmp = ['='*7 for i in range(total_line-1)]

        else:
            # The rest of lines are the value lines
            # Add header at the left and the rest is values
            tmp = [str(header[line_count-2]) + ' |'] + ['.' for i in range(total_line - 2)] # initialize value as '.'

        # Append the line to @matrix
        matrix.append(tmp)

    return matrix

##################### FUNCTION construct_matrix #####################
# construct the string matrix to be printed as an output
# ============================== INPUT ==============================
# @header[list]: list of headers to each row and column of the matrix
# @seq_pair[list]: list of sequence pairs
# @value[list]: list of values to add to the matrix
# @non_percentage[boolean]: whether the value in @value is percentage value
# ============================= OUTPUT ==============================
# @matrix[list]: a 2 dimensional list, each member contain a row of values
#####################################################################
def construct_matrix(header, seq_pair, value, non_percentage = False):
	# If the value list is empty, then there's nothing to construct
    if value == []:
        return None

    # Initiate the matrix
    matrix = init_matrix(header)

    # Loop through each sequence pairing and add values to matrix
    for each_pair in seq_pair:
        # ===================== LOCAL VARIABLES =====================
        # @each_pair[list]: a temp for each sequence pairing
		# @i1[int]: the greater sequence number label
		# @i2[int]: the lesser sequence number label
		# ===========================================================

		# assign each index value so matrix forms into triangle
        i1 = max(each_pair['pair'][0], each_pair['pair'][1])
        i2 = min(each_pair['pair'][0], each_pair['pair'][1])

        # Add the value to matrix position
        if not non_percentage:
        	# If the value should be the percentage, add '%' and format
            matrix[i1+2][i2+1] = '%.2f' % round(value[seq_pair.index(each_pair)], 2) + '%'
        else:
        	# If the value should not be the percentage, add the value as it is
            matrix[i1+2][i2+1] = '%d' % round(value[seq_pair.index(each_pair)], 2)

    return matrix

########################## FUNCTION report ##########################
# print all the reports of the calculated values
# ============================== INPUT ==============================
# @input_name[string]: input file name
# @output_name[string]: output file name
# @input_content[list]: a list of dictionary with name and sequence information
# @seq_type[string]: sequence type, the value is either 'P' or 'N'
# @pair[list]: list of dictionary contain 'pair' and 'length'
# @identity[list]: contain identity value for all pair of sequences
# @similarity[list]: contain similarity value for all pair of sequences
# @gap_count[list]: contain gap count value for all pair of sequences
# ============================= OUTPUT ==============================
# the function print all types of report files into file system
#####################################################################
def report(input_name, output_name, input_content, seq_type, pair, identity, similarity, gap_count):
	# ======================= LOCAL VARIABLES =======================
    # @report_output[file]: report file to be written
	# @matrix_output[file]: matrix file to be written
	# @gap_output[file]: gap file to be written
	# ===============================================================

	#################################################################
	####################### WRITE REPORT FILE #######################
	#################################################################
	# Print message to console and log file
	print 'Writing report file.'
	log_write('Writing report file.')

	# Open the file object to be written
	report_output = open(output_name + '_rpt.txt', 'w')
	# Write the file header
	report_write_header(report_output, input_content, seq_type, input_name)

	# Write identity statistic information
	report_write_stats('Identity', report_output, identity)
	# Write similarity statistic information if the sequence type is protein and the value is calculated
	if seq_type == 'P' and similarity != []:
		report_write_stats('Similarity', report_output, similarity)

	# Close the written file
	report_output.close()

	#################################################################
	####################### WRITE MATRIX FILE #######################
	#################################################################
	# Print message to console and log file
	print 'Writing matrix file.'
	log_write('Writing matrix file.')

	# Open the file object to be written
	matrix_output = open(output_name + '_mat.txt', 'w')
	# Write the file header
	report_write_header(matrix_output, input_content, seq_type, input_name)

	# Write full identity information matrix
	report_write_matrix('Identity', matrix_output, input_content, pair, identity)
	# Write full similarity information matrix if the sequence type is protein and the value is calculated
	if seq_type == 'P' and similarity != []:
		report_write_matrix('Similarity', matrix_output, input_content, pair, similarity)

	# Close the written file
	matrix_output.close()

	#################################################################
	##################### WRITE GAP REPORT FILE #####################
	#################################################################
	# Print message to console and log file
	print 'Writing Gap Information file.'
	log_write('Writing Gap Information file.')

	# Open the file object to be written
	gap_output = open(output_name + '_gap.txt', 'w')
	# Write the file header
	report_write_header(gap_output, input_content, seq_type, input_name)
	# Write gap statistic information and the full gap information matrix
	report_write_stats('Gap Information', gap_output, gap_count)
	report_write_matrix('Gap Information', gap_output, input_content, pair, gap_count)

	# Close the written file
	gap_output.close()

################### FUNCTION report_write_header ####################
# write the header section to the specified output file
# ============================== INPUT ==============================
# @file_data[file]: a file object of the file to write this header to
# @input_content[list]: a list of dictionary with name and sequence information
# @seq_type[string]: sequence type, the value is either 'P' or 'N'
# @input_name[string]: input file name
# ============================= OUTPUT ==============================
# this function write header content to the specified file object
#####################################################################
def report_write_header(file_data, input_content, seq_type, input_name):
	# Write the input file and its path
	file_data.write('Input File: ' + str(input_name) + '\n')

	# Write the type of this sequence
	file_data.write('Sequence Type: ')
	if seq_type == 'P':
		file_data.write('Protein\n')
	else:
		file_data.write('Nucleotide\n')

	# Write names of the input sequences
	file_data.write('\nSequences:\n')
	for each_seq in input_content:
		# ===================== LOCAL VARIABLES =====================
	    # @each_seq[dictionary]: a temp contain informations of each sequence
		# ===========================================================

		# Write number labels and sequence name to the file
		file_data.write('    ' + str(input_content.index(each_seq)+1) + '. ' + each_seq['name'] + '\n')

#################### FUNCTION report_write_stats ####################
# write the statistics section to the specified output file
# ============================== INPUT ==============================
# @identifier[string]: indicate which data this statistical value is from
# @file_data[file]: a file object of the file to write statistics to
# @value_list[string]: list of value to write statistics from
# ============================= OUTPUT ==============================
# this function write statistics to the specified file object
#####################################################################
def report_write_stats(identifier, file_data, value_list):
	# ======================= LOCAL VARIABLES =======================
	# @stats[dictionary]: contain all stats value calculated
	# ===============================================================

	# Write the header for the section
	file_data.write('\n' + identifier +' Statistics:\n')
	# Compute the statistics
	if not identifier == 'Gap Information':
		# If the value is not gap information
		# Then calculate with default
		stats = statistic(value_list)
	else:
		# If the value is gap information
		# Then calculate with @middle_value = False and @non_percentage = False
		stats = statistic(value_list, False, True)

	# print label and value for each
	for key, value in stats.iteritems():
		# ===================== LOCAL VARIABLES =====================
	    # @key[string]: a temp for each key in @stats
	    # @value[string/float]: a temp for the value of each @stats
		# ===========================================================
		file_data.write('    ' + key + ': ' + str(value) + '\n')

################### FUNCTION report_write_matrix ####################
# write the full value matrix to the specified output file
# ============================== INPUT ==============================
# @identifier[string]: indicate which data this matrix is from
# @file_data[file]: a file object of the file to write matrix to
# @input_content[list]: a list of dictionary with name and sequence information
# @pair[list]: list of dictionary contain 'pair' and 'length'
# @value_list[string]: list of value to write matrix from
# ============================= OUTPUT ==============================
# this function write full value matrix to the specified file object
#####################################################################
def report_write_matrix(identifier, file_data, input_content, pair, value_list):
	# ======================= LOCAL VARIABLES =======================
	# @matrix[list]: a 2 dimensional list to be written
	# ===============================================================

	# Write the header for the section
	file_data.write('\n' + identifier +' Matrix:\n\n')
	# Build string matrix from the given @value_list
	if not identifier == 'Gap Information':
		# If the value is not gap information
		# Then calculate with default
		matrix = construct_matrix(range(1, len(input_content)+1), pair, value_list)
	else:
		# If the value is gap information
		# Write value to @matrix as non-percentage values
		matrix = construct_matrix(range(1, len(input_content)+1), pair, value_list, True)

	# Write the @matrix to file
	for each_row in matrix:
		# ===================== LOCAL VARIABLES =====================
		# @each_row[list]: a temp for each row in the @matrix
		# ===========================================================
		for each_column in each_row:
			# =================== LOCAL VARIABLES ===================
			# @each_column[string]: a temp for each column in @each_row in @matrix
			# =======================================================
			file_data.write('{:>7s}'.format(str(each_column)) + '\t')

		file_data.write('\n')








#####################################################################
#####                                                           #####
#####                           MAIN                            #####
#####                                                           #####
#####################################################################
def main(modulecall = False, input_name = 'NONE', output_name = 'NONE', matrix_choice = 'NONE'):
	# ======================= LOCAL VARIABLES =======================
	# @modulecall[boolean]: indicate whether this is the direct program run
	#						or the script has been called as a module from other program
	# @input_name[string]: input file name
	# @output_name[string]: output file name
	# @matrix_choice[string]: matrix choice that the user selected
	# @input_content[list]: a list of dictionary with name and sequence information
	# @seq_type[string]: sequence type, the value is either 'P' or 'N'
	# @identity[list]: contain identity value for all pair of sequences
	# @similarity[list]: contain similarity value for all pair of sequences
	# @gap_count[list]: contain gap count value for all pair of sequences
	# @pair[list]: list of dictionary contain 'pair' and 'length'
	# @p_groups[list]: each member contain a tuple of proteins in the same group
	# ===============================================================

	# Initiate log writing
	if not os.path.exists(LOG_PATH):
		# If @LOG_PATH does not exists, then create the folder
		os.makedirs(LOG_PATH)
	# Open @LOG_FILE in append mode
	LOG_FILE = open(LOG_PATH + 'log[' + strftime("%Y-%m-%d") + '].txt','a')
	LOG_FILE.write('\n')

	# Print starting message to console and log file
	print '='*20
	print '\nSiMSA start.'
	log_write('Program start.')

	#################################################################
	########################### CHECKING ############################
	#################################################################
	# Checking arguments, input, and output
	if not modulecall:
		# If this is not module calling
		# Then check and get input and output form system arguments
		input_name = check_input()
		output_name = check_output()
	# Then get all sequence information from the input file
	input_content = fetch_sequences(input_name)
	# And check for the sequence type
	seq_type = check_sequence(input_content)

	#################################################################
	########################## CALCULATION ##########################
	#################################################################
	# Initiate variables to contain data to be calculated	
	similarity = []
	identity = []
	gap_count = []
	# Generate sequence pairing
	pair = pairing(input_content)
	# Get protein groups in case of protein sequence input
	if seq_type == 'P':
		# Prompt user to give matrix choice if this is not module calling
		if not modulecall:
			matrix_choice = get_matrix_choice()
		# Get protein groups from user-specified matrix choice
		p_groups = fetch_grouping(matrix_choice)

	# Print message to console and log file
	print 'Calculating similarity and identity...'
	log_write('Calculating similarity and identity...')

	# Loop through each pair of the input sequence
	# And calculate identity, similarity, and gap for the pair
	for each_pair in pair:
		# ===================== LOCAL VARIABLES =====================
		# @each_pair[dictionary]: a temp contain information of a pair of sequence
		#						  'pair': list of 2 sequence index
		#						  'length': length of the sequence after check the 
		# @seq1[string]: one of the sequence in the pair
		# @seq2[string]: another sequence in the pair to be calculate with @seq1
		# @pruned_index[list]: list of sequence index to be disregard from the calculation
		# ===========================================================

		# Get the sequence of each pair from the input sequences
		seq1 = input_content[each_pair['pair'][0]]['seq']
		seq2 = input_content[each_pair['pair'][1]]['seq']

		# Prune ambiguous positions out from the sequence
		prune_index = prune_sequence(seq1, seq2, seq_type)
		# Update the sequence length after pruning
		each_pair['length'] = len(seq1) - len(prune_index)
		# Append data of pruned position as gap position
		gap_count.append(len(prune_index))

		# Calculate similarity for protein sequence
		if seq_type == 'P' and p_groups != 'NONE':
			similarity.append(pairwise_similarity(seq1, seq2, p_groups, prune_index))
		# Calculate identity
		identity.append(parewise_identity(seq1, seq2, prune_index))

	# Print message to console and log file
	print 'Calculation completed.'
	log_write('Calculation completed.')

	#################################################################
	########################### REPORTING ###########################
	#################################################################
	# Print out all report to file system
	report(input_name, output_name, input_content, seq_type, pair, identity, similarity, gap_count)

	# Print ending message to console and log file
	print 'SiMSA ended.\n'
	print '='*20
	log_write('Program ended.')
	# Close the log file
	LOG_FILE.close()



#####################################################################
#####                           HELP                            #####
#####################################################################
def help():
	helpText = """
=======================================================================
||                             SiMSA.py                              ||
||                         SiMSA Help Menu                           ||
=======================================================================
||                                                                   ||
|| Calling script:                                                   ||
||     python SiMSA.py (options)                                     ||
||                                                                   ||
|| Options:                                                          ||
||     -i      Specify input files. (FASTA)                          ||
||     -o      Specify output files name, exclude file extension.    ||
||                                                                   ||
|| Supported FASTA format:                                           ||
||     ".fasta", ".fa", ".fas", ".fsa", and ".seq"                   ||
||                                                                   ||
|| Supported protein grouping:                                       ||
||     Protein in same bracket will be considered as similar         ||
||     ------------------------------------------------------------  ||
||     PAM60    (S, A, T) (H, N) (S, N) (N, D) (D, E) (H, Q) (Q, E)  ||
||              (Y, F) (R, K) (M, I) (L, M) (I, V)                   ||
||     ------------------------------------------------------------  ||
||     PAM250   (M, I, V, L) (D, N, H, Q, E) (F, I, L) (S, P, A)     ||
||              (S, A, G) (Q, K, N) (R, H, Q) (S, T) (R, K, Q)       ||
||              (S, N) (F, W) (R, W) (G, D)                          ||
||     ------------------------------------------------------------  ||
||     BLOSUM40 (S, T) (S, A) (S, Q, N) (H, Y, M) (D, E) (N, D)      ||
||              (H, N) (W, Y, F) (E, Q, K) (K, R) (R, Q)             ||
||              (L, M, I, V) (E, K) (A, G) (L, F, I) (V, T)          ||
||     ------------------------------------------------------------  ||
||     BLOSUM62 (S, T) (S, A) (S, N) (H, Y) (D, E) (N, D) (H, N)     ||
||              (W, Y, F) (E, Q, K) (K, R, Q) (L, M, I, V)           ||
||     ------------------------------------------------------------  ||
||     BLOSUM80 (Q, R, K) (Q, E, K) (E, D) (D, N) (Q, H) (Y, H)      ||
||              (Y, W) (Y, F) (S, T) (S, A) (M, V, L, I)             ||
||     ------------------------------------------------------------  ||
||                                                                   ||
=======================================================================
||          Laboratory of Microbial Diversity and Evolution          ||
||             Faculty of Sciences, Mahidol University               ||
=======================================================================
||             Last modified: """+ctime(os.path.getmtime('SiMSA.py'))+"""               ||
======================================================================="""
	print helpText



#####################################################################
#####                  EXTERNAL CALL FUNCTION                   #####
#####################################################################
def call_SiMSA(input_name, output_name, matrix_choice = 'NONE'):
	main(True, input_name, output_name, matrix_choice)



#####################################################################
#####                          GLOBAL                           #####
#####################################################################
if __name__ == "__main__":
	if '-h' in sys.argv:
		help() # check if help menu is requested
	else: # proceed to the main program
		main()