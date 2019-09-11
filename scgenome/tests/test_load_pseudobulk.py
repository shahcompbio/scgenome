import sys
import click
import logging
import collections

import dbclients.tantalus
import datamanagement.transfer_files

from scgenome.loaders.snv import load_cached_snv_data
from scgenome.loaders.breakpoint import load_cached_breakpoint_data
from scgenome.loaders.allele import load_cached_haplotype_allele_data

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


dtypes_check = {
    'snv_data': {
        'chrom': 'category',
        'coord': 'int64',
        'ref': 'category',
        'alt': 'category',
        'alt_counts_sum': 'int64',
        'ref_counts_sum': 'int64',
        'mappability': 'float64',
        'is_cosmic': 'object',
        'gene_name': 'category',
        'effect': 'category',
        'effect_impact': 'category',
        'amino_acid_change': 'category',
        'tri_nucleotide_context': 'category',
        'strelka_score': 'float64',
        'museq_score': 'float64',
        'num_cells': 'int64',
    },
    'snv_count_data': {
        'chrom': 'category',
        'coord': 'int64',
        'ref': 'category',
        'alt': 'category',
        'ref_counts': 'int64',
        'alt_counts': 'int64',
        'cell_id': 'category',
        'total_counts': 'int64',
        'sample_id': 'category',
    },
    'breakpoint_data': {
        'prediction_id': 'int64',
        'chromosome_1': 'object',
        'strand_1': 'object',
        'position_1': 'int64',
        'chromosome_2': 'object',
        'strand_2': 'object',
        'position_2': 'int64',
        'library_id': 'object', 
        'sample_id': 'object',
    },
    'breakpoint_count_data': {
        'prediction_id': 'int64',
        'cell_id': 'object',
        'read_count': 'int64',
        'library_id': 'object',
        'sample_id': 'object',
    },
    'allele_counts': {
        'allele_id': 'int64',
        'cell_id': 'category',
        'chromosome': 'category',
        'end': 'int64',
        'hap_label': 'int64',
        'readcount': 'int64',
        'start': 'int64',
    },
}

@click.command()
@click.option('--local_cache_directory')
@click.option('--local_storage_name')
@click.option('--single_ticket_id')
def test_load_cached_pseudobulk_data(local_cache_directory, local_storage_name, single_ticket_id):
    if local_cache_directory is not None and local_storage_name is not None:
        raise ValueError('local_cache_directory and local_storage_name are mutually exclusive')

    if local_cache_directory is None and local_storage_name is None:
        raise ValueError('require one of local_cache_directory and local_storage_name')

    tantalus_api = dbclients.tantalus.TantalusApi()

    if single_ticket_id is None:
        pseudobulk_analyses = tantalus_api.list('analysis', analysis_type__name='pseudobulk')

        version_tickets = collections.defaultdict(list)
        for analysis in pseudobulk_analyses:
            version_tickets[analysis['version']].append(analysis['jira_ticket'])
        
        ticket_ids = []
        for version in version_tickets:
            ticket_ids.append(version_tickets[version][-1])

    else:
        ticket_ids = (single_ticket_id,)

    for ticket_id in ticket_ids:
        logging.info(ticket_id)

        if local_cache_directory is not None:
            ticket_results = tantalus_api.list('results', analysis__jira_ticket=ticket_id)

            for results in ticket_results:
                filepaths = datamanagement.transfer_files.cache_dataset(
                    tantalus_api,
                    results['id'],
                    'resultsdataset',
                    'singlecellresults',
                    local_cache_directory,
                )

            local_results_directory = local_cache_directory

        elif local_storage_name is not None:
            local_results_directory = tantalus_api.get('storage', name=local_storage_name)['storage_directory']

        snv_results_tables = load_cached_snv_data(
            ticket_id,
            local_results_directory,
        )

        breakpoint_results_tables = load_cached_breakpoint_data(
            ticket_id,
            local_results_directory,
        )

        haplotype_results_tables = load_cached_haplotype_allele_data(
            ticket_id,
            local_results_directory,
        )

        for results_tables in (snv_results_tables, breakpoint_results_tables, haplotype_results_tables):
            for table_name, table_data in results_tables.items():
                logging.info(f'table {table_name} has size {len(table_data)}')
                for column_name, dtype_name in dtypes_check[table_name].items():
                    column_dtype = str(results_tables[table_name][column_name].dtype)
                    if not column_dtype == dtype_name:
                        raise Exception(f'{column_name} has dtype {column_dtype} not {dtype_name}')


if __name__ == '__main__':
    logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)
    test_load_cached_pseudobulk_data()
