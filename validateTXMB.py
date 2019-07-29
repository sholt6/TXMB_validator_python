import unittest
import sys
import gzip
import validateMetadataRecord
import validateMetadataTable
import validateFasta

def validate_metadata_record(metadata_record_filename):
	"""Validate metadata record from manifest file

	Keyword arguments:
	metadata_record_filename -- str, name of manifest file

	Returns:
	metadata_record_errors -- list of errors found with metadata record
	metadata_record -- dict, content of record
	custom_columns -- dict, custom columns with values
	"""

def generate_custom_col_dict(notyetsure):
	"""Generate dictionary of custom column names and descriptions
	### from what exactly?

	Keyword arguments:
	not yet sure

	Returns:
	custom_columns -- dictionary of custom column names and descriptions
	"""

def validate_metadata_table(table_file_name, record_custom_columns, ncbi_tax):
	"""Validate table of sequence metadata

	Keyword arguments:
	table_file_name -- str, file to be validated
	record_custom_columns -- dict, custom columns derived from metadata record
	ncbi_tax -- bool, whether NCBI taxonomy is in use

	Returns:
	metadata_table_errors -- list of errors found in metadata table
	local_identifiers -- list of sequence identifiers found in metadata table
	"""

def validate_fasta(fasta_filename, table_identifiers):
	"""Validate sequence fasta file

	Keyword arguments:
	fasta_filename -- str, name of fasta to be validated
	table_identifiers -- list of identifiers found in

	Returns:
	fasta_errors -- list of errors found in fasta_file
	"""


def report_errors(error_messages, dataset_name):
	"""Report on errors found in dataset

	Keyword arguments:
	error_messages -- list of errors found in all files
	dataset_name -- str, name of dataset as found in manifest file
	"""

def validate_txmb(manifest_filename):
	"""Main function to validate input TXMB dataset

	Keyword arguments:
	manifest_filename - str, name of manifest file

	Returns:
	validation_result - bool, True if no errors found
	"""

class Test(unittest.TestCase):

	# validate_metadata_record tests

	# generate_custom_col_dict tests

	# validate_metadata_table

	# validate_fasta tests

	# report_errors tests

	# validate_txmb tests
	def test_valid(self):
		test_valid_result = validate_txmb('valid.txt')
		assertTrue(test_valid_result)

	def test_valid_w_customs(self):
		valid_w_customs_result = validate_txmb('valid_w_customs.txt')
		assertTrue(valid_w_customs_result)

	def test_two_ids_f_valid_t(self):
		two_ids_f_valid_t_result = validate_txmb('two_ids_f_valid_t.txt')
		assertFalse(two_ids_f_valid_t_result)

	def test_two_seqs_f_valid_t(self):
		two_seqs_f_valid_t_result = validate_txmb('two_seqs_f_valid_t.txt')
		assertFalse(two_seqs_f_valid_t_result)

	def test_empty_id_f_valid_t(self):
		empty_id_f_valid_t_result = validate_txmb('empty_id_f_valid_t_result.txt')
		assertFalse(empty_id_f_valid_t_result)

	def test_empty_line_id_f_valid_t(self):
		empty_line_id_f_valid_t_result = validate_txmb('empty_line_id_f_valid_t.txt')
		assertFalse(empty_line_id_f_valid_t_result)

	def test_empty_line_seq_f_valid_t(self):
		empty_line_seq_f_valid_t_result = validate_txmb('empty_line_seq_f_valid_t.txt')
		assertFalse(empty_line_seq_f_valid_t_result)

	def test_missing_entry_f_valid_t(self):
		missing_entry_f_valid_t_result = validate_txmb('missing_entry_f_valid_t.txt')
		assertFalse(missing_entry_f_valid_t_result)

	def test_valid_f_empty_id_col_t(self):
		valid_f_empty_id_col_t_result = validate_txmb('valid_f_empty_id_col_t.txt')
		assertFalse(valid_f_empty_id_col_t_result)

	def test_valid_f_empty_row_t(self):
		valid_f_empty_row_t_result = validate_txmb('valid_f_empty_row_t.txt')
		assertFalse(valid_f_empty_row_t_result)

	def test_valid_f_missing_id_col_t(self):
		valid_f_missing_id_col_t_result = validate_txmb('valid_f_missing_id_col_t.txt')
		assertFalse(valid_f_missing_id_col_t_result)

	def test_valid_f_missing_identifier_t(self):
		valid_f_missing_identifier_t_result = validate_txmb('valid_f_missing_identifier_t.txt')
		assertFalse(valid_f_missing_identifier_t_result)

	def test_valid_f_missing_insdc_acc_and_range_t(self):
		valid_f_missing_insdc_acc_and_range_t_result = validate_txmb('valid_f_missing_insdc_acc_and_range_t.txt')
		assertTrue(valid_f_missing_insdc_acc_and_range_t_result)

	def test_valid_f_missing_insdc_acc_t(self):
		valid_f_missing_insdc_acc_t_result = validate_txmb('valid_f_missing_insdc_acc_t.txt')
		assertFalse(valid_f_missing_insdc_acc_t_result)

	def test_valid_f_missing_insdc_range_tself):
		valid_f_missing_insdc_range_t_result = validate_txmb('valid_f_missing_insdc_range_t.txt')
		assertTrue(valid_f_missing_insdc_range_t_result)

	def test_valid_f_missing_lineage_t(self):
		valid_f_missing_lineage_t_result = validate_txmb('valid_f_missing_lineage_t.txt')
		assertTrue(valid_f_missing_lineage_t_result)

	def test_valid_f_missing_row_t(self):
		valid_f_missing_row_t_result = validate_txmb('valid_f_missing_row_t.txt')
		assertFalse(valid_f_missing_row_t_result)

	def test_valid_f_missing_tax_id_t(self):
		valid_f_missing_tax_id_t_result = validate_txmb('valid_f_missing_tax_id_t.txt')
		assertFalse(valid_f_missing_tax_id_t_result)

	def test_valid_f_missing_taxon_name_t(self):
		valid_f_missing_taxon_name_t_result = validate_txmb('valid_f_missing_taxon_name_t.txt')
		assertFalse(valid_f_missing_taxon_name_t_result)

	def test_non_ncbi_w_tax_ids(self):
		non_ncbi_w_tax_ids_result = validate_txmb('non_ncbi_w_tax_ids.txt')
		assertFalse(non_ncbi_w_tax_ids_result)

	def test_non_ncbi_no_tax_ids(self):
		non_ncbi_no_tax_ids_result = validate_txmb('non_ncbi_no_tax_ids.txt')
		assert(non_ncbi_no_tax_ids_result)

	def test_invalid_dataset_name(self):
		invalid_dataset_name_result = validate_txmb('invalid_dataset_name.txt')
		assert(invalid_dataset_name_result)

	def test_different_tax_system_w_taxids(self):
		different_tax_system_w_taxids_result = validate_txmb('different_tax_system_w_taxids.txt')
		assert(different_tax_system_w_taxids_result)

	def test_different_tax_system_no_taxids(self):
		different_tax_system_no_taxids_result = validate_txmb('different_tax_system_no_taxids.txt')
		assert(different_tax_system_no_taxids_result)