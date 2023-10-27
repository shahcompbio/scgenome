import os
import sys
import logging

import scgenome.loaders.snv
import dbclients.tantalus
import datamanagement.transfer_files

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(format=LOGGING_FORMAT, stream=sys.stderr, level=logging.INFO)

local_cache_directory = os.environ['TANTALUS_CACHE_DIR']

tantalus_api = dbclients.tantalus.TantalusApi()

genotyping_ticket_id = 'SC-3341'

genotyping_results = tantalus_api.get(
    'resultsdataset',
    analysis__jira_ticket=genotyping_ticket_id,
    results_type='snv_genotyping')

genotyping_analysis = tantalus_api.get(
    'analysis',
    id=genotyping_results['analysis'],
)

variant_calling_results = []
cache_results_ids = [genotyping_results['id']]

for results_id in genotyping_analysis['input_results']:
    results = tantalus_api.get('resultsdataset', id=results_id)
    if results['results_type'] != 'variant_calling':
        continue
    variant_calling_results.append(results)
    cache_results_ids.append(results_id)

for results_id in cache_results_ids:
    filepaths = datamanagement.transfer_files.cache_dataset(
        tantalus_api,
        results_id,
        'resultsdataset',
        'singlecellresults',
        local_cache_directory,
    )

snv_data = []

for results in variant_calling_results:
    analysis = tantalus_api.get('analysis', id=results['analysis'])
    ticket_id = analysis['jira_ticket']
    ticket_directory = os.path.join(local_cache_directory, ticket_id)
    snv_results = scgenome.loaders.snv.load_snv_data(ticket_directory)
    snv_data.append(snv_results['snv_data'])

snv_data = scgenome.utils.concat_with_categories(snv_data, ignore_index=True)

positions = snv_data[['chrom', 'coord', 'ref', 'alt']].drop_duplicates()

ticket_id = genotyping_analysis['jira_ticket']
ticket_directory = os.path.join(local_cache_directory, ticket_id)
snv_results = scgenome.loaders.snv.load_snv_data(ticket_directory, positions=positions)

snv_count_data =  snv_results['snv_count_data']

scgenome.utils.union_categories([
    snv_data,
    snv_count_data,
])

print(snv_data.head())
print(snv_data.dtypes)

print(snv_count_data.head())
print(snv_count_data.dtypes)

