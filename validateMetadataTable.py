# validateMetadataRecord.py

"""
Set of functions used to validate metadata table of taxonomic reference set
analysis object before submission to ENA.
"""

import unittest
import requests
import re

character_regex = re.compile("^[A-Za-z0-9_]+$")
field_length_limit = 50

def validate_number_of_columns(mandatory_headers, input_headers, record_custom_columns):
	"""Checks metadata table contains correct number of columns

	Keyword arguments:
	mandatory_headers -- list of internally defined required headers
	input_headers -- full list of headers in given table
	record_custom_columns -- dictionary of custom column names given in metadata record

	Returns:
	column_number_errors -- list of any errors found with number of columns
	"""

	column_number_errors = []
	message = ''

	num_mandatory_headers = len(mandatory_headers)
	num_input_headers = len(input_headers)
	num_custom_columns = len(record_custom_columns)

	num_columns_expected = num_mandatory_headers + num_custom_columns

	try:
		assert(num_columns_expected == num_input_headers)
		message = ''
	except AssertionError:
		message = ("Metadata table contains {0} columns, {1} were expected. Please "
				 "ensure all custom columns are defined in the metadata record.".format(
				 num_input_headers, num_columns_expected))

	column_number_errors.append(message)

	return column_number_errors



def validate_mandatory_headers(input_headers, mandatory_headers):
	"""Checks all mandatory headers are present in the input table

	Keyword arguments:
	input_headers -- full list of headers in given table
	mandatory_headers -- list of mandatory header names

	Returns:
	custom_headers -- list of any non-mandatory headers identified
	mandatory_header_errors -- list of any errors found with header names
	"""

	mandatory_header_errors = []
	custom_headers = []

	for mand_index in range(0, len(mandatory_headers)):
		try:
			input_index = input_headers.index(mandatory_headers[mand_index])
		except ValueError:
			message = ("A mandatory header, '{0}', could not be found in the input"
					   " metadata table".format(mandatory_headers[mand_index]))
			mandatory_header_errors.append(message)
			continue
		mandatory_headers = mandatory_headers.pop(mand_index)
		input_headers = input_headers.pop(input_index)

	if (mandatory_headers):
		message = ("One or more mandatory headers were not found in the input "
				   "table: {0}".format(mandatory_headers))
		mandatory_header_errors.append(message)

	custom_headers = input_headers

	return mandatory_header_errors, custom_headers


def validate_custom_columns(table_custom_columns, record_custom_columns):
	"""Checks custom columns from metadata record all match with custom columns
	from metadata table

	Keyword arguments:
	table_custom_columns -- list, custom columns found in table
	record_custom_columns -- dict, custom columns specified in metadata record

	Returns:
	custom_column_errors -- list of any errors found with custom columns
	"""

	custom_column_errors = []

	record_headers = list(record_custom_columns.keys())

	sorted_record_headers = record_headers.sort()
	sorted_table_headers = record_headers.sort()

	if (len(sorted_record_headers) != len(sorted_table_headers)):
		message = ("Number of custom headers is mismatched between metadata "
				   "record and metadata table.")
		custom_column_errors.append(message)

	for record_index in range(0, len(sorted_record_headers)):
		try:
			header_to_check = sorted_record_headers[record_index]
			table_index = sorted_table_headers.index(header_to_check)
		except ValueError:
			message = ("A custom header specified in the metadata record, does"
					   " not exist in the metadata table: {0}".format(header_to_check))
			custom_column_errors.append(message)
			continue

		sorted_record_headers = sorted_record_headers.pop(record_index)
		sorted_table_headers = sorted_table_headers.pop(table_index)

	if sorted_table_headers:
		message = ("One or more headers used in the metadata table were "
				   "not defined in the metadata record: {0}".format(sorted_table_headers))
		custom_column_errors.append(message)

	return custom_column_errors


def validate_identifier(sequence_identifier):
	"""Checks whether an internal sequence identifier is valid. Not null.

	Keyword arguments:
	sequence_identifier -- string

	Returns:
	identifier_errors -- list of errors found with identifier
	"""

	identifier_errors = []

	if (not sequence_identifier):
		message = ("Sequence identifier is null")
		identifier_errors.append(message)

	if (len(sequence_identifier) > field_length_limit):
		message = ("Sequence identifier is too long (>50): {0}".format(sequence_identifier))
		identifier_errors.append(message)

	if (not re.match(character_regex, sequence_identifier)):
		message = ("Sequence identifier {0} does not match regular expression {1}"
				   .format(sequence_identifier, character_regex.pattern))

	return identifier_errors


def validate_insdc_sequence_accession(insdc_sequence_accession):
	"""Minimal checks on INSDC sequence accessions. Null input allowed.

	Keyword arguments:
	insdc_sequence_accession -- string accession

	Returns:
	accession_errors -- list of errors found with accession
	accession_present -- bool, whether accession is present or not
						 invalid accession still counts as accession
	"""

	accession_errors = []
	accession_present = False
	insdc_acc_regex = re.compile(r'^[A-Z]{1,6}[0-9]{5,8}(\.[0-9])?$')

	if (not insdc_sequence_accession):
		return accession_errors, accession_present
	else:
		accession_present = True

	valid_accession = re.match(insdc_acc_regex, insdc_sequence_accession)

	if not valid_accession:
		message = ("{0} does not appear to be a valid INSDC sequence accession"
				   .format(insdc_sequence_accession))
		accession_errors.append(message)

	return accession_errors, accession_present


def validate_insdc_sequence_range(insdc_sequence_range, accession_present):
	"""Checks that an INSDC sequence range is in appropriate format. Null input allowed.

	Keyword arguments:
	insdc_sequence_range -- string sequence range
	accession_present -- bool, whether an accession was given

	Returns:
	range_errors -- list of errors found with range
	"""

	range_errors = []
	range_regex = re.compile(r'^\d+\.\.\d+$')

	if not insdc_sequence_range:
		return range_errors
	elif (insdc_sequence_range and not accession_present):
		message = ("An accession range {0} was given for a record with no "
				   "accession specified".format(insdc_sequence_range))
		range_errors.append(message)

	if (not str(insdc_sequence_range)):
		message = ("A given accession range could not be read as a string")
		range_errors.append(message)

	valid_range = re.match(range_regex, insdc_sequence_range)

	if not valid_range:
		message = ("{0} is not a valid sequence range".format(insdc_sequence_range))
		range_errors.append(message)

	return range_errors


def validate_local_organism_name(local_organism_name, ncbi_tax):
	"""Checks that given organism name is valid. Not null.
	Could be modified to check NCBI tax browser

	Keyword arguments:
	local_organism_name -- string, organism name
	ncbi_tax -- bool, whether we are using NCBI tax or not

	Returns:
	organism_name_errors -- list of errors found with name
	expected_ncbi_tax_id -- int list, tax IDs determined if using NCBI tax
	"""

	organism_name_errors = []
	expected_ncbi_tax_id = 0
	tax_suggest_url = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/suggest-for-submission/'

	if not local_organism_name:
		message = ("No organism name has been given for this record")
		organism_name_errors.append(message)
		return organism_name_errors, expected_ncbi_tax_id

	if ncbi_tax:
		final_url = (tax_suggest_url + local_organism_name)


def validate_local_lineage(local_lineage):
	"""Minimal checks of organism lineage format. Null input allowed.
	Doesn't check that lineage is appropriate to organism name.

	Keyword arguments:
	local_lineage -- string, input linage of organism

	Returns:
	local_lineage_errors -- list of errors found
	"""


def validate_ncbi_tax_id(input_ncbi_tax_id, expected_ncbi_tax_id, ncbi_tax):
	"""Checks given Tax ID matches expected Tax ID, if NCBI tax in use

	Keyword arguments:
	input_ncbi_tax_id -- int, tax id from table
	expected_ncbi_tax_id -- int, tax id from function validate_local_organism_name()
	ncbi_tax -- bool, whether NCBI tax is being used for this submission

	Returns:
	ncbi_tax_id_errors -- list of errors found
	"""


class Test(unittest.TestCase):
	"""Test suite for module validateMetadataTable.py"""

	empty_list = ()
	boolean = True
	string = 'string'
	integer = 6
	too_long_name = 'reallylongnameofatleast50charactersisthisenoughofthemyet'
	unacceptable_characters_name = '!"£$%^&*()'
	mandatory_headers = ('Local Identifier',
						 'INSDC Sequence Accession',
						 'INSDC Sequence Range',
						 'Local Organism Name',
						 'Local Lineage',
						 'NCBI Taxon ID')
	correct_headers = ('Local Identifier',
					   'INSDC Sequence Accession',
					   'INSDC Sequence Range',
					   'Local Organism Name',
					   'Local Lineage',
					   'NCBI Taxon ID')
	missing_headers = ('Local Identifier',
					   'INSDC Sequence Accession',
					   'Local Lineage',
					   'NCBI Taxon ID')
	correct_headers_extras = ('Local Identifier',
					 		  'INSDC Sequence Accession',
					 		  'INSDC Sequence Range',
					 		  'Extra Header One',
					 	      'Local Organism Name',
					 		  'Local Lineage',
					 		  'Extra Header Two',
					 		  'NCBI Taxon ID')
	custom_columns = ('custom_column_1',
					  'custom_column_desc_1',
					  'custom_column_2',
					  'custom_column_desc_2')


	# validate_number_of_columns tests
	def test_num_columns_correct(self):
		num_columns_case_1 = validate_number_of_columns(self.mandatory_headers, self.correct_headers, self.custom_columns)
		assert(not num_columns_case_1[0])

	def test_num_columns_no_custom(self):
		num_columns_case_2 = validate_number_of_columns(self.mandatory_headers, self.correct_headers, self.empty_list)
		assert(not num_columns_case_2[0])

	def test_num_columns_missing_w_custom(self):
		num_columns_case_3 = validate_number_of_columns(self.mandatory_headers, self.missing_headers, self.custom_columns)
		assert(num_columns_case_3[0])

	def test_num_columns_missing_no_custom(self):
		num_columns_case_4 = validate_number_of_columns(self.mandatory_headers, self.missing_headers, self.empty_list)
		assert(num_columns_case_4[0])

	def test_num_columns_wrong_kw_1(self):
		num_columns_case_5 = validate_number_of_columns(self.mandatory_headers, self.integer, self.custom_columns)
		assert(num_columns_case_5[0])

	def test_num_columns_wrong_kw_2(self):
		num_columns_case_6 = validate_number_of_columns(self.mandatory_headers, self.correct_headers, self.integer)
		assert(num_columns_case_6[0])

	# validate_mandatory_headers tests
	def test_validate_mandatory_correct(self):
		mandatory_cols_case_1 = validate_mandatory_headers(self.correct_headers, self.mandatory_headers)
		assert(not mandatory_cols_case_1[0])
		assert(not mandatory_cols_case_1[1])

	def test_validate_mandatory_missing(self):
		mandatory_cols_case_2 = validate_mandatory_headers(self.missing_headers, self.mandatory_headers)
		assert(mandatory_cols_case_1[0])
		assert(not mandatory_cols_case_2[1])

	def test_validate_mandatory_reverse_mismatch(self):
		mandatory_cols_case_3 = validate_mandatory_headers(self.mandatory_headers, self.missing_headers)
		assert(not mandatory_cols_case_3[0])
		assert(mandatory_cols_case_3[1])

	def test_validate_mandatory_wrong_type_1(self):
		mandatory_cols_case_4 = validate_mandatory_headers(self.integer, self.mandatory_headers)
		assert(mandatory_cols_case_4[0])
		assert(not mandatory_cols_case_4[1])

	def test_validate_mandatory_wrong_type_2(self):
		mandatory_cols_case_5 = validate_mandatory_headers(self.correct_headers, integer)
		assert(mandatory_cols_case_5[0])
		assert(not mandatory_cols_case_5[1])

	# validate_identifier tests
	def test_validate_identifier_valid(self):
		validate_identifer_case_1 = validate_identifier(self.string)
		assert(not validate_identifer_case_1[0])

	def test_validate_identifer_empty(self):
		validate_identifer_case_2 = validate_identifier('')
		assert(validate_identifer_case_2[0])

	def test_validate_identifier_too_long(self):
		validate_identifer_case_3 = validate_identifier(self.too_long_name)
		assert(validate_identifer_case_3[0])

	def test_validate_identier_bad_chars(self):
		validate_identifer_case_4 = validate_identifier(self.unacceptable_characters_name)
		assert(validate_identifer_case_4[0])

	# validate_insdc_sequence_accession tests
	def test_validate_accession_valid(self):
		validate_accession_case_1 = validate_identifier('LR590077')
		assert(not validate_accession_case_1[0])
		assert(validate_accession_case_1[1])

	def test_validate_accession_invalid_string(self):
		validate_accession_case_2 = validate_identifier(self.string)
		assert(validate_accession_case_2[0])
		assert(validate_accession_case_2[1])

	def test_validate_accession_bad_insdc(self):
		validate_accession_case_3 = validate_identifier('58LR0077')
		assert(validate_accession_case_3[0])
		assert(validate_accession_case_3[1])

	def test_validate_accession_integer(self):
		validate_accession_case_4 = validate_identifier(self.integer)
		assert(validate_accession_case_4[0])
		assert(validate_accession_case_4[1])

	def test_validate_accession_bool(self):
		validate_accession_case_5 = validate_identifier(self.boolean)
		assert(validate_accession_case_5[0])
		assert(validate_accession_case_5[1])

	def test_validate_accession_null(self):
		validate_accession_case_6 = validate_identifier('')
		assert(not validate_accession_case_6[0])
		assert(not validate_accession_case_6[1])

	# validate_insdc_sequence_range tests
	def test_validate_sequence_range_valid_wanted(self):
		validate_seq_range_case_1 = validate_insdc_sequence_range('21..3523', True)
		assert(not validate_seq_range_case_1[0])

	def test_validate_sequence_range_valid_unwanted(self):
		validate_seq_range_case_2 = validate_insdc_sequence_range('21..3523', False)
		assert(validate_seq_range_case_2[0])

	def test_validate_sequence_range_invalid_string_wanted(self):
		validate_seq_range_case_3 = validate_insdc_sequence_range(self.string, True)
		assert(validate_seq_range_case_3[0])

	def test_validate_sequence_range_invalid_string_unwanted(self):
		validate_seq_range_case_4 = validate_insdc_sequence_range(self.string, False)
		assert(validate_seq_range_case_4[0])

	def test_validate_sequence_range_off_format(self):
		validate_seq_range_case_5 = validate_insdc_sequence_range('21,.3523', True)
		assert(validate_seq_range_case_5[0])

	def test_validate_sequence_range_integer(self):
		validate_seq_range_case_6 = validate_insdc_sequence_range(self.integer, True)
		assert(validate_seq_range_case_6[0])

	def test_validate_sequence_range_bool(self):
		validate_seq_range_case_7 = validate_insdc_sequence_range(self.boolean, True)
		assert(validate_seq_range_case_7[0])

	def test_validate_sequence_range_null_wanted(self):
		validate_seq_range_case_8 = validate_insdc_sequence_range('', True)
		assert(validate_seq_range_case_8[0])

	def test_validate_sequence_range_null_unwanted(self):
		validate_seq_range_case_9 = validate_insdc_sequence_range('', False)
		assert(not validate_seq_range_case_9[0])

	# validate_local_organism_name tests
	def test_validate_organism_name_human_ncbi(self):
		validate_organism_name_case_1 = validate_local_organism_name('Homo sapiens', True)
		assert(not validate_organism_name_case_1[0])
		assert(validate_organism_name_case_1[1] == '9606')

	def test_validate_organism_name_human_not_ncbi(self):
		validate_organism_name_case_2 = validate_local_organism_name('Homo sapiens', False)
		assert(not validate_organism_name_case_2[0])
		assert(not validate_organism_name_case_2[1])

	def test_validate_organism_name_stram_ncbi(self):
		validate_organism_name_case_3 = validate_local_organism_name('Stramenopiles sp. MAST-7 TOSAG23-7', True)
		assert(not validate_organism_name_case_3[0])
		assert(validate_organism_name_case_3[1] == '2590674')

	def test_validate_organism_name_stram_not_ncbi(self):
		validate_organism_name_case_4 = validate_local_organism_name('Stramenopiles sp. MAST-7 TOSAG23-7', False)
		assert(not validate_organism_name_case_4[0])
		assert(not validate_organism_name_case_4[1])

	def test_validate_organism_name_fake_ncbi(self): # fake name with ncbi validation
		validate_organism_name_case_5 = validate_local_organism_name('Schmlub gablub', True)
		assert(validate_organism_name_case_5[0])
		assert(not validate_organism_name_case_5[1])

	def test_validate_organism_name_fake_non_ncbi(self): # fake name without ncbi validation
		validate_organism_name_case_6 = validate_local_organism_name('Schmlub gablub', False)
		assert(not validate_organism_name_case_6[0])
		assert(not validate_organism_name_case_6[1])

	def test_validate_organism_name_integer(self):
		validate_organism_name_case_7 = validate_local_organism_name(self.integer, True)
		assert(validate_organism_name_case_7[0])
		assert(not validate_organism_name_case_7[1])

	def test_validate_organism_name_bool(self):
		validate_organism_name_case_8 = validate_local_organism_name(self.boolean, True)
		assert(validate_organism_name_case_8[0])
		assert(not validate_organism_name_case_8[1])

	def test_validate_organism_name_null(self):
		validate_organism_name_case_9 = validate_local_organism_name('', True)
		assert(validate_organism_name_case_9[0])
		assert(not validate_organism_name_case_9[1])

	# validate_local_lineage tests
	def test_validate_lineage_human(self):
		validate_lineage_case_1 = validate_local_lineage('cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; Euarchontoglires; Primates; Haplorrhini; Simiiformes; Catarrhini; Hominoidea; Hominidae; Homininae; Homo')
		assert(not validate_lineage_case_1[0])

	def test_validate_lineage_ecoli(self):
		validate_lineage_case_2 = validate_local_lineage('cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia')
		assert(not validate_lineage_case_2[0])

	def test_validate_lineage_string(self):
		validate_lineage_case_3 = validate_local_lineage(self.string)
		assert(not validate_lineage_case_3[0])

	def test_validate_lineage_integer(self):
		validate_lineage_case_4 = validate_local_lineage(self.integer)
		assert(validate_lineage_case_4[0])

	def test_validate_lineage_bool(self):
		validate_lineage_case_5 = validate_local_lineage(self.boolean)
		assert(not validate_lineage_case_5[0])

	def test_validate_lineage_null(self):
		validate_lineage_case_6 = validate_local_lineage('')
		assert(not validate_lineage_case_6[0])

	# validate_ncbi_tax_id tests
	def test_validate_tax_id_valid_w_ncbi(self):
		validate_tax_id_case_1 = validate_ncbi_tax_id(self.integer, self.integer, True)
		assert(validate_tax_id_case_1[0])

	def test_validate_tax_id_valid_wo_ncbi(self):
		validate_tax_id_case_2 = validate_ncbi_tax_id(self.integer, self.integer, False)
		assert(not validate_tax_id_case_2[0])

	def test_validate_tax_id_invalid_w_ncbi(self):
		validate_tax_id_case_3 = validate_ncbi_tax_id(self.integer, (self.integer + 1, True))
		assert(validate_tax_id_case_3[0])

	def test_validate_tax_id_invalid_wo_ncbi(self):
		validate_tax_id_case_4 = validate_ncbi_tax_id(self.integer, (self.integer + 1, False))
		assert(validate_tax_id_case_4[0])

	def test_validate_tax_id_null_1(self):
		validate_tax_id_case_5 = validate_ncbi_tax_id(0, self.integer, True)
		assert(validate_tax_id_case_5[0])

	def test_validate_tax_id_null_2(self):
		validate_tax_id_case_ = validate_ncbi_tax_id(self.integer, 0, True)
		assert(validate_tax_id_case_[0])

	def test_validate_tax_id_null_3(self):
		validate_tax_id_case_ = validate_ncbi_tax_id(self.integer, self.integer, 0)
		assert(validate_tax_id_case_[0])

	def test_validate_tax_id_null_all(self):
		validate_tax_id_case_ = validate_ncbi_tax_id(0, 0, 0)
		assert(not validate_tax_id_case_[0])