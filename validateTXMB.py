#!/usr/bin/env python3

import unittest
import sys
import gzip
import validateMetadataRecord as vmr
import validateMetadataTable as vmt
import validateFasta as vFa
import pandas as pd


def validate_metadata_record(metadata_record_filename):
	"""Validate metadata record from manifest file

	Keyword arguments:
	metadata_record_filename -- str, name of manifest file

	Returns:
	metadata_record_errors -- list of errors found with metadata record
	metadata_record -- dict, content of record
	record_custom_columns -- dict, custom columns with values
	ncbi_tax -- bool, whether NCBI taxonomy has been specified for this submission
	"""

	metadata_record_errors = []
	metadata_record = {}
	raw_custom_columns = {}
	record_custom_columns = {}
	ncbi_tax = False

	mandatory_record_content = ['LOCALTAXONOMY', 'REFERENCEDATASETNAME',
								'FASTA', 'TABLE']

	try:
		with open(metadata_record_filename) as metadata_record_file:
			record_content = metadata_record_file.readlines()
	except FileNotFoundError:
		message = ('Could not find ' + metadata_record_filename)
		metadata_record_errors.append(message)
		return (metadata_record_errors, metadata_record,
				record_custom_columns, ncbi_tax)

	for line in record_content:
		line_content = line.split(None, 1)
		try:
			assert(len(line_content) == 2)
		except AssertionError:
			if line_content[0] not in mandatory_record_content:
				continue
			else:
				message = ('The following metadata record line could not be '
						   'processed. This may be malformed, or a mandatory'
						   ' value was not supplied:\n' + line)
				metadata_record_errors.append(message)
				return (metadata_record_errors, metadata_record,
						record_custom_columns, ncbi_tax)

		if line_content[0] in mandatory_record_content:
			metadata_record[line_content[0]] = line_content[1]
		else:
			raw_custom_columns[line_content[0]] = line_content[1]

	metadata_record_errors.extend(
		vmr.validate_file_content(metadata_record, mandatory_record_content))

	metadata_record_errors.extend(
		vmr.validate_dataset_name(metadata_record['REFERENCEDATASETNAME']))

	tax_validation = vmr.validate_local_taxonomy(metadata_record['LOCALTAXONOMY'])
	metadata_record_errors.extend(tax_validation[0])
	ncbi_tax = tax_validation[1]

	if 'LOCALTAXONOMYVERSION' in metadata_record:
		metadata_record_errors.extend(
			vmr.validate_local_taxonomy_version(metadata_record['LOCALTAXONOMYVERSION']))

	if raw_custom_columns:
		custom_col_gen_result = vmr.generate_custom_col_dict(raw_custom_columns)
		metadata_record_errors.extend(custom_col_gen_result[0])
		record_custom_columns = custom_col_gen_result[1]
		metadata_record_errors.extend(
			vmr.validate_custom_columns(record_custom_columns))

	return (metadata_record_errors, metadata_record,
			record_custom_columns, ncbi_tax)


def validate_metadata_table(table_filename, record_custom_columns, ncbi_tax):
	"""Validate table of sequence metadata

	Keyword arguments:
	table_file_name -- str, file to be validated
	record_custom_columns -- dict, custom columns derived from metadata record
	ncbi_tax -- bool, whether NCBI taxonomy is in use

	Returns:
	metadata_table_errors -- list of errors found in metadata table
	local_identifiers -- list of sequence identifiers found in metadata table
	"""

	mandatory_headers = ('local_identifier',
						 'insdc_sequence_accession',
						 'insdc_sequence_range',
						 'local_organism_name',
						 'local_lineage',
						 'ncbi_tax_id')

	metadata_table_errors = []
	local_identifiers = []

	try:
		sequence_metadata_table = pd.read_csv(table_filename,
										compression='gzip', header=0, sep='\t')
	except (FileNotFoundError, UnicodeDecodeError, OSError) as e:
		return e, local_identifiers

	input_headers = list(sequence_metadata_table.columns.values)

	metadata_table_errors.extend(
		vmt.validate_number_of_columns(mandatory_headers, input_headers,
									   record_custom_columns))

	mandatory_header_validation = vmt.validate_mandatory_headers(input_headers,
															 mandatory_headers)
	metadata_table_errors.extend(mandatory_header_validation[0])
	table_custom_headers = mandatory_header_validation[1]

	metadata_table_errors.extend(vmt.validate_custom_columns(
		table_custom_headers, record_custom_columns))

	total_records = sequence_metadata_table.shape[0]
	current_record = 1

	for row_index in range(0, total_records):
		row_errors = []

		sequence_identifier = sequence_metadata_table.loc[row_index, mandatory_headers[0]]
		insdc_sequence_accession = sequence_metadata_table.loc[row_index, mandatory_headers[1]]
		insdc_sequence_range = sequence_metadata_table.loc[row_index, mandatory_headers[2]]
		local_organism_name = sequence_metadata_table.loc[row_index, mandatory_headers[3]]
		local_lineage = sequence_metadata_table.loc[row_index, mandatory_headers[4]]
		input_ncbi_tax_id = sequence_metadata_table.loc[row_index, mandatory_headers[5]]

		row_errors.extend(vmt.validate_identifier(sequence_identifier))
		local_identifiers.append(sequence_identifier)

		accession_errors, accession_present = vmt.validate_insdc_sequence_accession(
			insdc_sequence_accession)
		row_errors.extend(accession_errors)

		row_errors.extend(vmt.validate_insdc_sequence_range(
									insdc_sequence_range, accession_present))

		organism_errors, expected_tax_ids = vmt.validate_local_organism_name(
												local_organism_name, ncbi_tax)
		row_errors.extend(organism_errors)

		row_errors.extend(vmt.validate_local_lineage(local_lineage))

		row_errors.extend(vmt.validate_ncbi_tax_id(input_ncbi_tax_id,
							expected_tax_ids, ncbi_tax, local_organism_name))

		for i in range(0, len(row_errors)):
			row_errors[i] = ("Metadata Table Line {0}: {1}".
							 format(current_record, row_errors[i]))

		metadata_table_errors.extend(row_errors)

		current_record += 1

	if len(local_identifiers) != len(set(local_identifiers)):
		message = ("Metadata table contains duplicate local identifiers")
		metadata_table_errors.append(message)

	return metadata_table_errors, local_identifiers


def validate_fasta(fasta_filename, table_identifiers):
	"""Validate sequence fasta file

	Keyword arguments:
	fasta_filename -- str, name of fasta to be validated
	table_identifiers -- list of identifiers found in

	Returns:
	fasta_errors -- list of errors found in fasta_file
	"""

	fasta_errors = vFa.validate_txmb_fasta(fasta_filename, table_identifiers)

	return fasta_errors


def report_errors(report_filename, error_messages):
	"""Report on errors found in dataset

	Keyword arguments:
	report_file -- str, name for error report
	error_messages -- list of errors found in all files
	"""

	report_file = open(report_filename, 'w')
	message_start = 'ERROR: '
	message_end = '\n'

	for error in error_messages:
		report_file.write(message_start + error + message_end)

	report_file.close()

	return True


def validate_txmb(manifest_filename):
	"""Main function to validate input TXMB dataset

	Keyword arguments:
	manifest_filename - str, name of manifest file

	Returns:
	validation_result - bool, True if no errors found
	"""

	metadata_record_errors = []
	metadata_table_errors = []
	fasta_errors = []

	record_custom_columns = {}
	local_identifiers = []
	ncbi_tax = False

	metadata_record_errors, metadata_record, record_custom_columns, ncbi_tax = \
									 validate_metadata_record(manifest_filename)

	report_filename = (metadata_record['REFERENCEDATASETNAME'] + ".report")

	if metadata_record_errors:
		report_errors(report_filename, metadata_record_errors)
		sys.exit('Could not validate metadata record, please view error ' +
				 'messages in ' + report_filename)

	table_filename = metadata_record('TABLE')

	metadata_table_errors, local_identifiers = validate_metadata_table(
							table_filename, record_custom_columns, ncbi_tax)

	if metadata_table_errors:
		report_errors(report_filename, metadata_table_errors)
		sys.exit('Could not validate ' + table_filename + ', please view ' +
				 'error messages in ' + report_filename)

	fasta_filename = metadata_record('FASTA')

	fasta_errors = validate_fasta(fasta_filename, local_identifiers)

	if fasta_errors:
		report_errors(report_filename, fasta_errors)
		sys.exit('Could not validate ' + fasta_filename + ', please view' +
				 'error messages in ' + report_filename)

	return metadata_record['REFERENCEDATASETNAME']


def main():

	manifest_filename = sys.argv[1]

	submission_valid = validate_txmb(manifest_filename)

	if submission_valid:
		print('Taxonomic reference set \'{0}\' validated successfully'
			  .format(submission_valid))

if __name__ == '__main__':
	main()


class Test(unittest.TestCase):
	valid_record = {'REFERENCEDATASETNAME' : 'valid_submission',
					'LOCALTAXONOMY' : 'NCBI',
					'LOCALTAXONOMYVERSION' : '',
					'FASTA' : 'valid.fasta.gz',
					'TABLE' : 'valid.tsv.gz'}

	valid_record_w_customs = {'REFERENCEDATASETNAME' : 'valid_submission',
						      'LOCALTAXONOMY' : 'NCBI',
							  'LOCALTAXONOMYVERSION' : '',
							  'FASTA' : 'valid.fasta.gz',
							  'TABLE' : 'valid_w_customs.tsv.gz'}

	custom_columns = {'Annotation' : 'Source of annotation',
					  'ITSoneDB URL' : 'URL within ITSoneDB'}

	table_identifiers = ['ITS1DB00887249', 'ITS1DB00944432', 'ITS1DB01095025',
					     'ITS1DB01019240', 'ITS1DB00588026', 'ITS1DB00588027',
				   	 	 'ITS1DB00588024', 'ITS1DB00588025', 'ITS1DB00588022',
				   		 'ITS1DB00588023']

	# validate_metadata_record tests
	def test_mdata_record_val_valid(self):
		mdata_record_val_result_valid = validate_metadata_record('Test_Files/valid.txt')
		self.assertFalse(mdata_record_val_result_valid[0])
		self.assertTrue(mdata_record_val_result_valid[1] == self.valid_record)
		self.assertFalse(mdata_record_val_result_valid[2])

	def test_mdata_record_val_valid_w_customs(self):
		mdata_record_val_result_valid_w_customs = validate_metadata_record('Test_Files/valid_w_customs.txt')
		self.assertFalse(mdata_record_val_result_valid_w_customs[0])
		self.assertTrue(mdata_record_val_result_valid_w_customs[1] == self.valid_record_w_customs)
		self.assertTrue(mdata_record_val_result_valid_w_customs[2] == self.custom_columns)

	# generate_custom_col_dict tests


	# validate_metadata_table
	def test_mdata_tab_val_valid(self):
		mdata_tab_val_result_valid = validate_metadata_table('Test_Files/valid.tsv.gz', {}, True)
		self.assertFalse(mdata_tab_val_result_valid[0])
		self.assertTrue(mdata_tab_val_result_valid[1] == self.table_identifiers)

	def test_mdata_tab_val_valid_w_customs(self):
		mdata_tab_val_result_valid = validate_metadata_table('Test_Files/valid_w_customs.tsv.gz', self.custom_columns, True)
		self.assertFalse(mdata_tab_val_result_valid_w_customs[0])
		self.assertTrue(mdata_tab_val_result_valid_w_customs[1] == self.table_identifiers)

	def test_mdata_tab_val_missing_taxon_name(self):
		mdata_tab_val_result_missing_taxon_name = validate_metadata_table('Test_Files/missing_taxon_name.tsv.gz', {}, True)
		self.assertTrue(mdata_tab_val_result_missing_taxon_name[0])
		self.assertTrue(mdata_tab_val_result_missing_taxon_name[1] == self.table_identifiers)

	def test_mdata_tab_val_missing_tax_id(self):
		mdata_tab_val_result_missing_tax_id = validate_metadata_table('Test_Files/missing_tax_id.tsv.gz', {}, True)
		self.assertTrue(mdata_tab_val_result_missing_tax_id[0])
		self.assertTrue(mdata_tab_val_result_missing_tax_id[1] == self.table_identifiers)

	def test_mdata_tab_val_missing_row(self):
		mdata_tab_val_result_missing_row = validate_metadata_table('Test_Files/missing_row.tsv.gz', {}, True)
		self.assertTrue(mdata_tab_val_result_missing_row[0])
		self.assertFalse(mdata_tab_val_result_missing_row[1] == self.table_identifiers)

	def test_mdata_tab_val_missing_lineage(self):
		mdata_tab_val_result_missing_lineage = validate_metadata_table('Test_Files/missing_lineage.tsv.gz', {}, True)
		self.assertFalse(mdata_tab_val_result_[0])
		self.assertTrue(mdata_tab_val_result_[1] == self.table_identifiers)

	def test_mdata_tab_val_missing_insdc_range(self):
		mdata_tab_val_result_missing_insdc_range = validate_metadata_table('Test_Files/missing_insdc_range.tsv.gz', {}, True)
		self.assertFalse(mdata_tab_val_result_missing_insdc_range[0])
		self.assertTrue(mdata_tab_val_result_missing_insdc_range[1] == self.table_identifiers)

	def test_mdata_tab_val_missing_insdc_acc(self):
		mdata_tab_val_result_missing_insdc_acc = validate_metadata_table('Test_Files/missing_insdc_acc.tsv.gz', {}, True)
		self.assertTrue(mdata_tab_val_result_missing_insdc_acc[0])
		self.assertTrue(mdata_tab_val_result_missing_insdc_acc[1] == self.table_identifiers)

	def test_mdata_tab_val_missing_insdc_acc_and_range(self):
		mdata_tab_val_result_missing_insdc_acc_and_range = validate_metadata_table('Test_Files/missing_insdc_acc_and_range.tsv.gz', {}, True)
		self.assertFalse(mdata_tab_val_result_missing_insdc_acc_and_range[0])
		self.assertTrue(mdata_tab_val_result_missing_insdc_acc_and_range[1] == self.table_identifiers)

	def test_mdata_tab_val_missing_identifier(self):
		mdata_tab_val_result_missing_identifier = validate_metadata_table('Test_Files/missing_identifier.tsv.gz', {}, True)
		self.assertTrue(mdata_tab_val_result_missing_identifier[0])
		self.assertFalse(mdata_tab_val_result_missing_identifier[1] == self.table_identifiers)

	def test_mdata_tab_val_missing_id_column(self):
		mdata_tab_val_result_missing_id_column = validate_metadata_table('Test_Files/missing_id_column.tsv.gz', {}, True)
		self.assertTrue(mdata_tab_val_result_missing_id_column[0])
		self.assertFalse(mdata_tab_val_result_missing_id_column[1] == self.table_identifiers)

	def test_mdata_tab_val_empty_row(self):
		mdata_tab_val_result_empty_row = validate_metadata_table('Test_Files/empty_row.tsv.gz', {}, True)
		self.assertTrue(mdata_tab_val_result_empty_row[0])
		self.assertFalse(mdata_tab_val_result_empty_row[1] == self.table_identifiers)

	def test_mdata_tab_val_empty_id_column(self):
		mdata_tab_val_result_empty_id_column = validate_metadata_table('Test_Files/empty_id_column.tsv.gz', {}, True)
		self.assertTrue(mdata_tab_val_result_empty_id_column[0])
		self.assertFalse(mdata_tab_val_result_empty_id_column[1] == self.table_identifiers)

	# validate_fasta tests
	def test_fasta_val_valid(self):
		fasta_val_result_valid = validate_fasta('Test_Files/valid.fasta.gz', self.table_identifiers)
		self.assertFalse(fasta_val_result_valid)

	def test_fasta_val_(self):
		fasta_val_result_two_seqs = validate_fasta('Test_Files/two_seqs.fasta.gz', self.table_identifiers)
		self.assertTrue(fasta_val_result_two_seqs)

	def test_fasta_val_two_ids(self):
		fasta_val_result_two_ids = validate_fasta('Test_Files/two_ids.fasta.gz', self.table_identifiers)
		self.assertTrue(fasta_val_result_two_ids)

	def test_fasta_val_missing_entry(self):
		fasta_val_result_missing_entry = validate_fasta('Test_Files/missing_entry.fasta.gz', self.table_identifiers)
		self.assertTrue(fasta_val_result_missing_entry)

	def test_fasta_val_empty_line_seq(self):
		fasta_val_result_empty_line_seq = validate_fasta('Test_Files/empty_line_seq.fasta.gz', self.table_identifiers)
		self.assertTrue(fasta_val_result_empty_line_seq)

	def test_fasta_val_empty_line_id(self):
		fasta_val_result_empty_line_id = validate_fasta('Test_Files/empty_line_id.fasta.gz', self.table_identifiers)
		self.assertTrue(fasta_val_result_empty_line_id)

	def test_fasta_val_empty_id(self):
		fasta_val_result_empty_id = validate_fasta('Test_Files/empty_id.fasta.gz', self.table_identifiers)
		self.assertTrue(fasta_val_result_empty_id)

	# validate_txmb tests
	def test_txmb_val_valid(self):
		test_txmb_val_valid_result = validate_txmb('valid.txt')
		self.assertTrue(test_txmb_val_valid_result)

	def test_txmb_val_valid_w_customs(self):
		valid_w_customs_result = validate_txmb('valid_w_customs.txt')
		self.assertTrue(valid_w_customs_result)

	def test_txmb_val_two_ids_f_valid_t(self):
		two_ids_f_valid_t_result = validate_txmb('two_ids_f_valid_t.txt')
		self.assertFalse(two_ids_f_valid_t_result)

	def test_txmb_val_two_seqs_f_valid_t(self):
		two_seqs_f_valid_t_result = validate_txmb('two_seqs_f_valid_t.txt')
		self.assertFalse(two_seqs_f_valid_t_result)

	def test_txmb_val_empty_id_f_valid_t(self):
		empty_id_f_valid_t_result = validate_txmb('empty_id_f_valid_t_result.txt')
		self.assertFalse(empty_id_f_valid_t_result)

	def test_txmb_val_empty_line_id_f_valid_t(self):
		empty_line_id_f_valid_t_result = validate_txmb('empty_line_id_f_valid_t.txt')
		self.assertFalse(empty_line_id_f_valid_t_result)

	def test_txmb_val_empty_line_seq_f_valid_t(self):
		empty_line_seq_f_valid_t_result = validate_txmb('empty_line_seq_f_valid_t.txt')
		self.assertFalse(empty_line_seq_f_valid_t_result)

	def test_txmb_val_missing_entry_f_valid_t(self):
		missing_entry_f_valid_t_result = validate_txmb('missing_entry_f_valid_t.txt')
		self.assertFalse(missing_entry_f_valid_t_result)

	def test_txmb_val_valid_f_empty_id_col_t(self):
		valid_f_empty_id_col_t_result = validate_txmb('valid_f_empty_id_col_t.txt')
		self.assertFalse(valid_f_empty_id_col_t_result)

	def test_txmb_val_valid_f_empty_row_t(self):
		valid_f_empty_row_t_result = validate_txmb('valid_f_empty_row_t.txt')
		self.assertFalse(valid_f_empty_row_t_result)

	def test_txmb_val_valid_f_missing_id_col_t(self):
		valid_f_missing_id_col_t_result = validate_txmb('valid_f_missing_id_col_t.txt')
		self.assertFalse(valid_f_missing_id_col_t_result)

	def test_txmb_val_valid_f_missing_identifier_t(self):
		valid_f_missing_identifier_t_result = validate_txmb('valid_f_missing_identifier_t.txt')
		self.assertFalse(valid_f_missing_identifier_t_result)

	def test_txmb_val_valid_f_missing_insdc_acc_and_range_t(self):
		valid_f_missing_insdc_acc_and_range_t_result = validate_txmb('valid_f_missing_insdc_acc_and_range_t.txt')
		self.assertTrue(valid_f_missing_insdc_acc_and_range_t_result)

	def test_txmb_val_valid_f_missing_insdc_acc_t(self):
		valid_f_missing_insdc_acc_t_result = validate_txmb('valid_f_missing_insdc_acc_t.txt')
		self.assertFalse(valid_f_missing_insdc_acc_t_result)

	def test_txmb_val_valid_f_missing_insdc_range_t(self):
		valid_f_missing_insdc_range_t_result = validate_txmb('valid_f_missing_insdc_range_t.txt')
		self.assertTrue(valid_f_missing_insdc_range_t_result)

	def test_txmb_val_valid_f_missing_lineage_t(self):
		valid_f_missing_lineage_t_result = validate_txmb('valid_f_missing_lineage_t.txt')
		self.assertTrue(valid_f_missing_lineage_t_result)

	def test_txmb_val_valid_f_missing_row_t(self):
		valid_f_missing_row_t_result = validate_txmb('valid_f_missing_row_t.txt')
		self.assertFalse(valid_f_missing_row_t_result)

	def test_txmb_val_valid_f_missing_tax_id_t(self):
		valid_f_missing_tax_id_t_result = validate_txmb('valid_f_missing_tax_id_t.txt')
		self.assertFalse(valid_f_missing_tax_id_t_result)

	def test_txmb_val_valid_f_missing_taxon_name_t(self):
		valid_f_missing_taxon_name_t_result = validate_txmb('valid_f_missing_taxon_name_t.txt')
		self.assertFalse(valid_f_missing_taxon_name_t_result)

	def test_txmb_val_non_ncbi_w_tax_ids(self):
		non_ncbi_w_tax_ids_result = validate_txmb('non_ncbi_w_tax_ids.txt')
		self.assertFalse(non_ncbi_w_tax_ids_result)

	def test_txmb_val_non_ncbi_no_tax_ids(self):
		non_ncbi_no_tax_ids_result = validate_txmb('non_ncbi_no_tax_ids.txt')
		assert(non_ncbi_no_tax_ids_result)

	def test_txmb_val_invalid_dataset_name(self):
		invalid_dataset_name_result = validate_txmb('invalid_dataset_name.txt')
		assert(invalid_dataset_name_result)

	def test_txmb_val_different_tax_system_w_taxids(self):
		different_tax_system_w_taxids_result = validate_txmb('different_tax_system_w_taxids.txt')
		assert(different_tax_system_w_taxids_result)

	def test_txmb_val_different_tax_system_no_taxids(self):
		different_tax_system_no_taxids_result = validate_txmb('different_tax_system_no_taxids.txt')
		assert(different_tax_system_no_taxids_result)