# validateMetadataRecord.py

"""
Set of functions used to validate metadata record of taxonomic reference set
analysis object before submission to ENA.
"""

import unittest
import re

character_regex = re.compile("^[A-Za-z0-9_]+$")
character_regex_with_space = re.compile("^[A-Za-z0-9_ ]+$")
field_length_limit = 50

def validate_file_content(metadata_record_dict, mandatory_record_content):
	"""Checks for list of required fields in manifest file

	Keyword arguments:
	metadata_record_dict -- dict, content of manifest file
	mandatory_record_content -- list, required tags

	Returns:
	metadata_record_file_errors -- list of errors found
	"""

	metadata_record_file_errors = []
	input_fields = list(metadata_record_dict.keys())

	for required_field in mandatory_record_content:
		if required_field in input_fields:
			pass
		else:
			message = ("Required field '{0}' not found in ")
			metadata_record_file_errors.append(message)

	return metadata_record_file_errors


def validate_local_taxonomy(local_taxonomy):
	"""Validate definition of taxonomic system

	Keyword arguments:
	local_taxonomy -- string giving name of taxonomic system

	Returns:
	taxonomy_validation_errors -- text of any errors found
	ncbi_tax -- True or False
	"""

	field_name = 'local_taxonomy'
	ncbi_synonyms = ('NCBI', 'ncbi')
	ncbi_tax = False
	taxonomy_validation_errors = []

	try:
		assert(local_taxonomy)
	except AssertionError:
		taxonomy_validation_errors.append(
									get_empty_mandatory_value_error(field_name))
		return taxonomy_validation_errors

	try:
		assert(str(local_taxonomy))
	except AssertionError:
		taxonomy_validation_errors.append(
					get_data_type_error_text(field_name, local_taxonomy, 'str'))

	try:
		assert(re.match(character_regex, str(local_taxonomy)))
	except AssertionError:
		taxonomy_validation_errors.append(
			get_regex_mismatch_error_text(field_name, character_regex))

	try:
		assert(len(str(local_taxonomy)) <= field_length_limit)
	except AssertionError:
		taxonomy_validation_errors.append(
			get_field_length_error_text(field_name))

	if (local_taxonomy in ncbi_synonyms):
		ncbi_tax = True

	return (taxonomy_validation_errors, ncbi_tax)


def validate_local_taxonomy_version(local_taxonomy_version):
	"""Validate version information of taxonomic system

	Keyword arguments:
	local_taxonomy_version -- user-provided version of taxonomic system

	Returns:
	local_taxonomy_version_errors -- text of any errors found
	"""

	field_name = "local_taxonomy_version"
	tax_version_regex = re.compile("^[A-Za-z0-9_.]+$")
	local_taxonomy_version_errors = []

	try:
		assert(local_taxonomy_version == '')
		return local_taxonomy_version_errors
	except AssertionError:
		pass

	try:
		assert(str(local_taxonomy_version))
	except AssertionError:
		local_taxonomy_version_errors.append(
					get_data_type_error_text(field_name, local_taxonomy_version,
											 'str'))
		return local_taxonomy_version_errors

	try:
		assert(re.match(tax_version_regex, str(local_taxonomy_version)))
	except AssertionError:
		local_taxonomy_version_errors.append(
			get_regex_mismatch_error_text(field_name,
										  tax_version_regex))

	try:
		assert(len(str(local_taxonomy_version)) <= field_length_limit)
	except AssertionError:
		local_taxonomy_version_errors.append(get_field_length_error_text(
													  				field_name))

	return local_taxonomy_version_errors


def validate_dataset_name(reference_dataset_name):
	"""Validate user-define name of analysis object

	Keyword arguments:
	reference_dataset_name -- user-defined name for this submission

	Returns:
	dataset_name_errors -- text of any errors found
	"""

	field_name = "reference_dataset_name"
	reference_dataset_name_errors = []

	try:
		assert(re.match(character_regex, reference_dataset_name))
	except AssertionError:
		reference_dataset_name_errors.append(get_regex_mismatch_error_text(
												field_name, character_regex))

	try:
		assert(len(reference_dataset_name) <= field_length_limit)
	except AssertionError:
		reference_dataset_name_errors.append(get_field_length_error_text(
																	field_name))

	return reference_dataset_name_errors


def generate_custom_col_dict(raw_custom_columns):
	"""Takes custom column definitions from metadata record, parses into dict
	with structure {'col name' : 'col description'}

	Keyword arguments:
	raw_custom_columns -- dict of user-provided custom fields in manifest

	Returns:
	custom_col_generation_errors -- list of errors found in parsing input
	record_custom_columns -- dict of custom columns and their descriptions
	"""

	custom_col_generation_errors = []
	record_custom_columns = {}

	num_extra_cols = len(raw_custom_columns) / 2

	if num_extra_cols.is_integer():
		pass
	else:
		message = ('Custom columns do not all have definitions')
		custom_col_generation_errors.append(message)

	for i in range(0, int(num_extra_cols)):
		key_key = 'CUSTOMCOLUMNHEADER' + str(i+1)
		val_key = 'CUSTOMCOLUMNHEADER' + str(i+1) + 'DESCRIPTION'

		try:
			record_custom_columns[raw_custom_columns[key_key]] = raw_custom_columns[val_key]
		except KeyError as e:
			message = ('Expected to find line \'{0}\' in manifest file'.format(e))
			custom_col_generation_errors.append(message)

	return custom_col_generation_errors, record_custom_columns


def validate_custom_columns(record_custom_columns):
	"""Confirms custom columns meet some minimal standards including character
	restriction, name length, and pairing

	Keyword arguments:
	record_custom_columns -- dict of user-provided column names and descriptions

	Returns:
	custom_column_name_errors -- text of any errors found
	"""

	custom_column_name_errors = []

	for k, v in record_custom_columns.items():
		field_name = "Custom Column: " + k
		field_desc = "Column Description: " + v

		if not re.match(character_regex_with_space, k):
			message = get_regex_mismatch_error_text(field_name, character_regex)
			custom_column_name_errors.append(message)

		if not re.match(character_regex_with_space, v):
			message = get_regex_mismatch_error_text(field_desc, character_regex_with_space)
			custom_column_name_errors.append(message)

		if (len(k) > field_length_limit):
			message = get_field_length_error_text(field_name)
			custom_column_name_errors.append(message)

		if (len(v) > field_length_limit):
			message = get_field_length_error_text(field_desc)
			custom_column_name_errors.append(message)

	return custom_column_name_errors


def get_regex_mismatch_error_text(field_name, source_regex):
	"""Provides error message when entered data does not match regex defined
	above

	Keyword arguments:
	field_name -- name of field where error was found

	Returns
	string containing error message
	"""

	return("Value entered for '{0}' does not match regex '{1}'"
		   .format(field_name, source_regex.pattern))


def get_field_length_error_text(field_name):
	"""Provides error message when entered data is too long

	Keyword arguments:
	field_name -- name of field where error was found

	Returns
	string containing error message
	"""

	return("Value entered for '{0}' exceeds character length limit of {1}"
		   .format(field_name, str(field_length_limit)))


def get_data_type_error_text(field_name, field_value, type_name):
	"""Provides error message when an entered value cannot be converted to
	required type

	Keyword arguments:
	field_name -- name of field where error was found
	field_value -- value which was found to be in error
	type_name -- type of value which is required

	Returns
	string containing error message
	"""

	message = ''

	try:
		message = ("Value '{0}' entered for '{1}' could not be parsed as a valid {2}"
				   .format(str(field_value),field_name,type_name))
	except TypeError:
		message = ("A value entered for '{0}' could not be read".format(field_name))

	return message


def get_empty_mandatory_value_error(field_name):
	"""Provides error message for when a mandatory field is provided an empty
	value

	Keyword arguments;
	field_name -- str name of field which was not provided

	Returns
	string containing error message
	"""

	message = ("No value was given for mandatory field '{0}'".format(field_name))

	return message


class Test(unittest.TestCase):
	"""Test suite for module validateMetadataRecord.py"""
	too_long_name = 'reallylongnameofatleast50charactersisthisenoughofthemyet'
	unacceptable_characters_name = '!"£$%^&*()'

	mandatory_record_content = ["LOCALTAXONOMY", "REFERENCEDATASETNAME",
								"FASTA", "TABLE"]
	mandatory_record_content_no_table = ["LOCALTAXONOMY", "LOCALTAXONOMYVERSION",
										 "REFERENCEDATASETNAME", "FASTA"]
	mdata_record_dict = {"LOCALTAXONOMY" : "NCBI", "LOCALTAXONOMYVERSION" : "1",
						 "REFERENCEDATASETNAME" : "test_name", "FASTA" : "file.fasta.gz",
						 "TABLE" : "file.tsv.gz"}
	mdata_record_dict_no_table = {"LOCALTAXONOMY" : "NCBI", "LOCALTAXONOMYVERSION" : "1",
								  "REFERENCEDATASETNAME" : "test_name",
								  "FASTA" : "file.fasta.gz"}

	# Metadata record input tests
	def test_mdata_record_valid(self):
		record_validation_case_1 = validate_file_content(self.mdata_record_dict, self.mandatory_record_content)
		assert(not record_validation_case_1)

	def test_mdata_record_dict_missing_one(self):
		record_validation_case_1 = validate_file_content(self.mdata_record_dict_no_table, self.mandatory_record_content)
		assert(record_validation_case_1[0])

	def test_mdata_record_missing_mandatory(self):
		record_validation_case_1 = validate_file_content(self.mdata_record_dict, self.mandatory_record_content_no_table)
		assert(not record_validation_case_1)


	# Taxonomy name tests
	def test_local_taxonomy_length(self):
		local_taxonomy_case_1 = validate_local_taxonomy(self.too_long_name)
		assert(local_taxonomy_case_1[0])
		assert(not local_taxonomy_case_1[1])

	def test_local_taxonomy_regex(self):
		local_taxonomy_case_2 = validate_local_taxonomy(self.unacceptable_characters_name)
		assert(local_taxonomy_case_2[0])
		assert(not local_taxonomy_case_2[1])

	def test_local_taxonomy_good_input(self):
		local_taxonomy_case_3 = validate_local_taxonomy('A_Taxonomy_Database')
		assert(not local_taxonomy_case_3[0])
		assert(not local_taxonomy_case_3[1])

	def test_local_taxonomy_upper_ncbi_input(self):
		local_taxonomy_case_4 = validate_local_taxonomy('NCBI')
		assert(not local_taxonomy_case_4[0])
		assert(local_taxonomy_case_4[1])

	def test_local_taxonomy_lower_ncbi_input(self):
		local_taxonomy_case_5 = validate_local_taxonomy('ncbi')
		assert(not local_taxonomy_case_5[0])
		assert(local_taxonomy_case_5[1])

	def test_local_taxonomy_integer(self):
		local_taxonomy_case_6 = validate_local_taxonomy(6)
		assert(not local_taxonomy_case_6[0])
		assert(not local_taxonomy_case_6[1])

	def test_local_taxonomy_empty(self):
		local_taxonomy_case_7 = validate_local_taxonomy('')
		assert(local_taxonomy_case_7[0])
		assert(len(local_taxonomy_case_7) == 1)

	# Taxonomy version tests
	def test_local_taxonomy_version_length(self):
		local_taxonomy_version_case_1 = validate_local_taxonomy_version(self.too_long_name)
		assert(local_taxonomy_version_case_1[0])

	def test_local_taxonomy_version_regex(self):
		local_taxonomy_version_case_2 = validate_local_taxonomy_version(self.unacceptable_characters_name)
		assert(local_taxonomy_version_case_2[0])

	def test_local_taxonomy_version_integer(self):
		local_taxonomy_version_case_3 = validate_local_taxonomy_version(3)
		assert(not local_taxonomy_version_case_3)

	def test_local_taxonomy_version_empty(self):
		local_taxonomy_version_case_4 = validate_local_taxonomy_version('')
		assert(not local_taxonomy_version_case_4)

	# Dataset name tests
	def test_dataset_name_length(self):
		dataset_name_case_1 = validate_dataset_name(self.too_long_name)
		assert(dataset_name_case_1[0])

	def test_dataset_name_regex(self):
		dataset_name_case_2 = validate_dataset_name('$%^&*()')
		assert(dataset_name_case_2[0])

	# Custom column dict creation tests
	input_custom_cols = {'CUSTOMCOLUMNHEADER1' : 'Annotation',
						 'CUSTOMCOLUMNHEADER1DESCRIPTION' : 'Source of annotation',
						 'CUSTOMCOLUMNHEADER2' : 'ITSoneDB URL',
						 'CUSTOMCOLUMNHEADER2DESCRIPTION' : 'URL within ITSoneDB'}
	processed_custom_cols = {'Annotation' : 'Source of annotation',
							 'ITSoneDB URL' : 'URL within ITSoneDB'}
	input_custom_cols_misnumbered = {'CUSTOMCOLUMNHEADER1' : 'Annotation',
						 			 'CUSTOMCOLUMNHEADER2DESCRIPTION' : 'Source of annotation',
						 			 'CUSTOMCOLUMNHEADER3' : 'ITSoneDB URL',
						 			 'CUSTOMCOLUMNHEADER4DESCRIPTION' : 'URL within ITSoneDB'}
	input_custom_cols_misnamed  = {'CUSTOMCOLUMNHEADER1' : 'Annotation',
						 		   'CUSTOMCOLHEADER2DESCRIPTION' : 'Source of annotation',
						 		   'CUSTOMCOLUMNDER3' : 'ITSoneDB URL',
						 		   'CUSTOMCOLUMNHEADER4DESCRIPTION' : 'URL within ITSoneDB'}
	input_custom_cols_long = {'Annotation lotsofextracharacterstomakethistoolongandthereforeinvalid' : 'Source of annotation',
							  'ITSoneDB URL' : 'URL within ITSoneDB'}
	input_custom_cols_bad_chars = {'Annotation!!!' : '$$Source of annotation',
								   'ITSoneDB URL@@' : '##URL within ITSoneDB'}

	def test_custom_col_gen_valid(self):
		col_gen_case_1 = generate_custom_col_dict(self.input_custom_cols)
		self.assertFalse(col_gen_case_1[0])
		self.assertTrue(col_gen_case_1[1] == self.processed_custom_cols)

	def test_custom_col_gen_mismatched(self):
		col_gen_case_2 = generate_custom_col_dict(self.input_custom_cols_misnumbered)
		self.assertTrue(col_gen_case_2[0])
		self.assertFalse(col_gen_case_2[1])

	def test_custom_col_gen_mispelled(self):
		col_gen_case_3 = generate_custom_col_dict(self.input_custom_cols_misnamed)
		self.assertTrue(col_gen_case_3[0])
		self.assertFalse(col_gen_case_3[1])

	# Custom column validation tests
	long_columns = [too_long_name, too_long_name]
	bad_name_columns = [unacceptable_characters_name, unacceptable_characters_name]
	odd_no_of_columns = ['col1', 'desc1', 'col2']

	def test_columns_valid(self):
		columns_name_case_3 = validate_custom_columns(self.processed_custom_cols)
		self.assertFalse(columns_name_case_3)

	def test_columns_long(self):
		columns_name_case_1 = validate_custom_columns(self.input_custom_cols_long)
		self.assertTrue(columns_name_case_1[0])

	def test_columns_bad_name(self):
		columns_name_case_2 = validate_custom_columns(self.input_custom_cols_bad_chars)
		self.assertTrue(columns_name_case_2[0])



if __name__ == '__main__':
	unittest.main()