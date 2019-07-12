# validateFasta.py

"""
Set of functions used to validate FASTA file of taxonomic reference set
analysis object before submission to ENA.
"""

import unittest
import re
import gzip


def validate_txmb_fasta(file_name, table_identifiers):
	"""Main function for validating the FASTA component of a taxonomic
	reference set

	Keyword arguments:
	file_name -- str, name of file to be checked
	table_identifiers -- list of identifiers taken from TSV

	Returns:
	fasta_errors -- list of errors found with file
	"""

	fasta_errors = []
	linecount = 0
	line_flag = ''
	remaining_ids = table_identifiers.copy()

	try:
		fasta_file = gzip.open(file_name, 'rt')
	except FileNotFoundError as e:
		fasta_errors.append(e)
		return fasta_errors
	except OSError as e:
		fasta_errors.append(e)
		return fasta_errors

	for line in fasta_file:
		linecount += 1

		if line[0] == '>':
			if line_flag == 'se' or line_flag =='':
				line_flag = 'id'
				id_line_errors, id_index = check_identifier(line, remaining_ids, linecount)
				remaining_ids.pop(id_index)
				fasta_errors.extend(id_line_errors)
			elif line_flag == 'id':
				message = ("Two consecutive ID lines in FASTA at line {0}"
						   .format(linecount))
				fasta_errors.append(message)
		else:
			if line_flag == 'id':
				line_flag = 'se'
				seq_line_errors = check_sequence(line, linecount)
				fasta_errors.extend(seq_line_errors)
			elif line_flag == 'se':
				message = ("Two consecutive sequence lines in FASTA at line {0}"
						   .format(linecount))
				fasta_errors.append(message)

	if remaining_ids:
		message = ("The following identifiers found in the metadata table "
				   "were not found in the FASTA file:\n")
		for identifier in remaining_ids:
			message = message + identifier + "\n"
			fasta_errors.append(message)

	if (linecount % 2 != 0):
		message = ("Input FASTA does not contain an even number of lines")
		fasta_errors.append(message)

	fasta_file.close()

	return fasta_errors


def check_identifier(id_line, remaining_ids, linecount):
	"""Takes a line recognised as an ID line and checks for an identifier
	which matches one found in the metadata TSV

	Keyword arguments:
	id_line -- str, full ID line found in FASTA file
	remaining_ids -- list, identifiers not yet found in file
	linecount -- current line in file

	Returns:
	id_line_errors
	identifier_index
	"""

	id_line_errors = []

	line_content = id_line.split('|')
	identifier_section = line_content[0]
	identifier = identifier_section[1:]
	id_index = 0

	try:
		id_index = remaining_ids.index(identifier)
	except ValueError:
		message = ("FASTA identifier '{0}' at line {1} does not match anything"
				   " in metadata table".format(identifier, linecount))
		id_line_errors.append(message)

	return id_line_errors, id_index


def check_sequence(sequence_line, linecount):
	"""Takes a line recognised as sequence line and confirms it matches
	the nucleotide regex

	Keyword arguments:
	sequence_line -- str, a sequence line from FASTA
	linecount -- current line in file

	Returns:
	sequence_errors -- list of any errors found to be added to main error log
	"""

	sequence_errors = []
	nucleotide_regex = re.compile("^[ATCGactgRYSWKMBDHVryswkmbdhvNn]*$")
	non_nucleotide_regex = re.compile("[efijlopquxzEFIJLOPQUXZ]")

	if re.match(nucleotide_regex, sequence_line):
		return sequence_errors

	if re.match(r"\s", sequence_line):
		message = ("Sequence at line {0} contains whitespace".format(linecount))
		sequence_errors.append(message)

	if re.match(non_nucleotide_regex, sequence_line):
		message = ("Sequence at line {0} contains non-nucleotide characters"
			       .format(linecount))
		sequence_errors.append(message)

	if not sequence_line:
		message = ("Sequence at line {0} is null"
			       .format(linecount))
		sequence_errors.append(message)

	if not sequence_errors:
		message = ("Sequence errors at line {0}, could not be diagnosed"
				   .format(linecount))
		sequence_errors.append(message)

	return sequence_errors


class Test(unittest.TestCase):
	"""Test suite for module validateFasta.py"""

	local_identifiers = ['ITS1DB00887249',
						 'ITS1DB00944432',
						 'ITS1DB01095025',
						 'ITS1DB01019240',
						 'ITS1DB00588026',
						 'ITS1DB00588027',
						 'ITS1DB00588024',
						 'ITS1DB00588025',
						 'ITS1DB00588022',
						 'ITS1DB00588023']

	local_identifiers_missing = ['ITS1DB00887249',
								 'ITS1DB00944432',
								 'ITS1DB01095025',
								 'ITS1DB01019240',
								 'ITS1DB00588026',
								 'ITS1DB00588027',
								 'ITS1DB00588024',
								 'ITS1DB00588025',
								 'ITS1DB00588022']

	# validate_txmb_fasta tests
	def test_fasta_validation_valid(self):
		fasta_val_case_1 = validate_txmb_fasta('Test_FASTAs/valid.fasta.gz', self.local_identifiers)
		assert(not fasta_val_case_1)

	def test_fasta_validation_valid_missing_list_entry(self):
		fasta_val_case_2 = validate_txmb_fasta('Test_FASTAs/valid.fasta.gz', self.local_identifiers_missing)
		assert(fasta_val_case_2[0])

	def test_fasta_validation_2_seq_lines(self):
		fasta_val_case_3 = validate_txmb_fasta('Test_FASTAs/two_seqs.fasta.gz', self.local_identifiers)
		assert(fasta_val_case_3[0])

	def test_fasta_validation_2_id_lines(self):
		fasta_val_case_4 = validate_txmb_fasta('Test_FASTAs/two_ids.fasta.gz', self.local_identifiers)
		assert(fasta_val_case_4[0])

	def test_fasta_validation_missing_entry(self):
		fasta_val_case_5 = validate_txmb_fasta('Test_FASTAs/missing_entry.fasta.gz', self.local_identifiers)
		assert(fasta_val_case_5[0])

	def test_fasta_validation_empty_seq_line(self):
		fasta_val_case_6 = validate_txmb_fasta('Test_FASTAs/empty_line_seq.fasta.gz', self.local_identifiers)
		assert(fasta_val_case_6[0])

	def test_fasta_validation_empty_id_line(self):
		fasta_val_case_7 = validate_txmb_fasta('Test_FASTAs/empty_line_id.fasta.gz', self.local_identifiers)
		assert(fasta_val_case_7[0])

	def test_fasta_validation_empty_id(self):
		fasta_val_case_8 = validate_txmb_fasta('Test_FASTAs/empty_id.fasta.gz', self.local_identifiers)
		assert(fasta_val_case_8[0])


	# check_identifier tests
	def test_valid_identifier_line(self):
		test_id_line_case_1 = check_identifier('>ITS1DB00588023|175245|uncultured fungus|ITS1 located by HMM annotation, 142bp', self.local_identifiers, 1)
		assert(not test_id_line_case_1[0])
		assert(test_id_line_case_1[1])

	def test_missing_pipe_full_line(self):
		test_id_line_case_2 = check_identifier('>ITS1DB00588023175245uncultured fungusITS1 located by HMM annotation, 142bp', self.local_identifiers, 1)
		assert(test_id_line_case_2[0])
		assert(not test_id_line_case_2[1])

	def test_missing_pipe_only_id(self):
		test_id_line_case_3 = check_identifier('>ITS1DB00588023', self.local_identifiers, 1)
		assert(not test_id_line_case_3[0])
		assert(test_id_line_case_3[1])

	def test_only_id(self):
		test_id_line_case_4 = check_identifier('>ITS1DB00588023|', self.local_identifiers, 1)
		assert(not test_id_line_case_4[0])
		assert(test_id_line_case_4[1])

	def test_id_not_in_list(self):
		test_id_line_case_5 = check_identifier('>ITS1DB99999999|175245|uncultured fungus|ITS1 located by HMM annotation, 142bp', self.local_identifiers, 1)
		assert(test_id_line_case_5[0])
		assert(not test_id_line_case_5[1])


	# check_sequence tests
	def test_seq_check_valid(self):
		test_seq_line_case_1 = check_sequence('gtcgaaccctgcgatagcagacgacccggtaacatgtaaacacatcgggagagcatcaaggggcatcgatcccttgaccgcgaaccctaggttggtgtccgcttagattctaaggggacgccgattgacataaccatccgggcgcggcatgcgccaaggaacttaaaatcgaattgtatgttcgcttcccgttatcgggcagcagcgtcattccaaaaaacaacg', 1)
		assert(not test_seq_line_case_1)

	def test_seq_check_bad_chars(self):
		test_seq_line_case_2 = check_sequence('abcdefghijklmnopqrstuvwxyz!£$%^&*()', 1)
		assert(test_seq_line_case_2[0])

	def test_seq_check_newline(self):
		test_seq_line_case_3 = check_sequence('actg\nactg', 1)
		assert(test_seq_line_case_3[0])

	def test_seq_check_spaces(self):
		test_seq_line_case_4 = check_sequence('actg actg', 1)
		assert(test_seq_line_case_4[0])

	def test_seq_check_null(self):
		test_seq_line_case_5 = check_sequence('', 1)
		assert(test_seq_line_case_5[0])


if (__name__ == "__main__"):
	unittest.main()