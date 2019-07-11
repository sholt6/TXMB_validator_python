# validateMetadataRecord.py

"""
Set of functions used to validate metadata record of taxonomic reference set
analysis object before submission to ENA.
"""

import unittest
import re

character_regex = re.compile("^[A-Za-z0-9_]+$")
field_length_limit = 50

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
		assert(local_taxonomy_version)
	except AssertionError:
		local_taxonomy_version_errors.append(
									get_empty_mandatory_value_error(field_name))
		return local_taxonomy_version_errors

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


def validate_custom_columns(input_custom_columns):
	"""Confirms custom columns meet some minimal standards including character
	restriction, name length, and pairing

	Keyword arguments:
	input_custom_columns -- list of user-provided column names and descriptions

	Returns:
	custom_column_name_errors -- text of any errors found
	"""

	custom_column_name_errors = []

	try:
		assert(len(input_custom_columns) % 2 == 0)
	except AssertionError:
		column_no_message = ("Number of custom columns and descriptions is not "
		"even, indicating that some are not matched.")
		custom_column_name_errors.append(column_no_message)

	for name in input_custom_columns:
		field_name = "custom_column " + name

		try:
			assert(re.match(character_regex, name))
		except AssertionError:
			custom_column_name_errors.append(get_regex_mismatch_error_text(
												   field_name, character_regex))

		try:
			assert(len(name) <= field_length_limit)
		except AssertionError:
			custom_column_name_errors.append(get_field_length_error_text(
																	field_name))

	return custom_column_name_errors


def get_regex_mismatch_error_text(field_name, source_regex):
	"""Provides error message when entered data does not match regex defined
	above

	Keyword arguments:
	field_name -- name of field where error was found

	Returns
	string containing error message
	"""

	return("Value entered for " + field_name +
			" does not match regex " + source_regex.pattern)


def get_field_length_error_text(field_name):
	"""Provides error message when entered data is too long

	Keyword arguments:
	field_name -- name of field where error was found

	Returns
	string containing error message
	"""

	return("Value entered for " + field_name +
			" exceeds character length limit of " + str(field_length_limit))


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
		message = ("Value " + str(field_value) + " entered for " + field_name +
		   " could not be parsed as a valid " + type_name)
	except TypeError:
		message = "A value entered for " + field_name + " could not be read"

	return message


def get_empty_mandatory_value_error(field_name):
	"""Provides error message for when a mandatory field is provided an empty
	value

	Keyword arguments;
	field_name -- str name of field which was not provided

	Returns
	string containing error message
	"""

	message = "No value was given for mandatory field " + field_name

	return message


class Test(unittest.TestCase):
	"""Test suite for module validateMetadataRecord.py"""
	too_long_name = 'reallylongnameofatleast50charactersisthisenoughofthemyet'
	unacceptable_characters_name = '!"£$%^&*()'

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
		assert(local_taxonomy_version_case_4[0])

	# Dataset name tests
	def test_dataset_name_length(self):
		dataset_name_case_1 = validate_dataset_name(self.too_long_name)
		assert(dataset_name_case_1[0])

	def test_dataset_name_regex(self):
		dataset_name_case_2 = validate_dataset_name('$%^&*()')
		assert(dataset_name_case_2[0])

	# Custom column tests
	long_columns = [too_long_name, too_long_name]
	bad_name_columns = [unacceptable_characters_name, unacceptable_characters_name]
	odd_no_of_columns = ['col1', 'desc1', 'col2']

	def test_columns_long(self):
		columns_name_case_1 = validate_custom_columns(self.long_columns)
		assert(columns_name_case_1[0])

	def test_columns_bad_name(self):
		columns_name_case_2 = validate_custom_columns(self.bad_name_columns)
		assert(columns_name_case_2[0])

	def test_columns_odd_no(self):
		columns_name_case_3 = validate_custom_columns(self.odd_no_of_columns)
		assert(columns_name_case_3[0])

if __name__ == '__main__':
	unittest.main()