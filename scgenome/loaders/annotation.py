import os
from collections import defaultdict

import pandas as pd
import scgenome.csvutils
import scgenome.loaders.utils
import scgenome.utils
import yaml

_categorical_cols = [
    'cell_id',
    'sample_id',
    'library_id',
]

_table_suffixes_v0_2_25 = (
    ('annotation_metrics', '_metrics.csv.gz'),
)

table_suffixes = defaultdict(lambda: _table_suffixes_v0_2_25, {
    'v0.2.25': _table_suffixes_v0_2_25,
    'v0.3.0': _table_suffixes_v0_2_25,
    'v0.3.1': _table_suffixes_v0_2_25,
})


def load_annotation_data_from_file(filepath, table_name="annotation_metrics"):
    results_tables = {}

    if table_name == 'annotation_metrics':
        data = process_annotation_file(filepath, is_anno_metrics_table=True)
    else:
        data = process_annotation_file(filepath, is_anno_metrics_table=False)

    results_tables[table_name] = data

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def load_annotation_data(
        results_dir,
):
    """ Load copy number tables
    
    Args:
        results_dir (str): results directory to load from.
    
    Returns:
        dict: pandas.DataFrame tables keyed by table name
    """

    analysis_dirs = scgenome.loaders.utils.find_results_directories(
        results_dir)

    if 'annotation' not in analysis_dirs:
        raise ValueError(f'no annotation found for directory {results_dir}')

    annotation_results_dir = analysis_dirs['annotation']
    assert len(annotation_results_dir) == 1
    annotation_results_dir = annotation_results_dir[0]

    manifest_filename = os.path.join(annotation_results_dir, 'metadata.yaml')
    manifest = yaml.safe_load(open(manifest_filename))

    # KLUDGE: 0.3.1 -> v0.3.1
    if not manifest['meta']['version'].startswith('v'):
        manifest['meta']['version'] = 'v' + manifest['meta']['version']

    version = manifest['meta']['version']

    results_tables = {}

    for table_name, suffix in table_suffixes[version]:
        filenames = scgenome.loaders.utils.find_filenames(manifest['filenames'], suffix)

        if len(filenames) != 1:
            raise ValueError(f'found filenames {filenames} for suffix {suffix}')

        filename = filenames[0]

        filepath = os.path.join(annotation_results_dir, filename)

        if table_name == 'annotation_metrics':
            data = process_annotation_file(filepath, is_anno_metrics_table=True)
        else:
            data = process_annotation_file(filepath, is_anno_metrics_table=False)

        results_tables[table_name] = data

    scgenome.utils.union_categories(results_tables.values())

    return results_tables


def process_annotation_file(filepath, is_anno_metrics_table=False):
    csv_input = scgenome.csvutils.CsvInput(filepath)

    dtypes_override = None
    if is_anno_metrics_table:
        dtypes_directory = os.path.join(os.path.dirname(__file__), 'dtypes')
        dtypes_filename = os.path.join(dtypes_directory, 'metrics_column_defs.yaml')
        dtypes_override = yaml.safe_load(open(dtypes_filename))
        dtypes_override = {a['name']: a['dtype'] for a in dtypes_override}

    data = csv_input.read_csv(dtypes_override=dtypes_override)

    data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
    data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

    for col in _categorical_cols:
        if col in data:
            data[col] = pd.Categorical(data[col])

    return data
