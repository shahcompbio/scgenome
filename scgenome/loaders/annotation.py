import logging
import yaml
import os
import pandas as pd

import scgenome.utils
import scgenome.loaders.utils


_categorical_cols = [
    'cell_id',
    'sample_id',
    'library_id',
]


_table_suffixes_v0_2_25 = [
    ('annotation_metrics', '_metrics.csv.gz'),
]


table_suffixes = {
    'v0.2.25': _table_suffixes_v0_2_25,
    'v0.3.0': _table_suffixes_v0_2_25,
    'v0.3.1': _table_suffixes_v0_2_25,
}


def _table_fixes_v0_2_25(results_tables):
    pass # TODO


_table_fixes = {
    'v0.2.25': _table_fixes_v0_2_25,
    'v0.3.0': _table_fixes_v0_2_25,
    'v0.3.1': _table_fixes_v0_2_25,
}


def load_annotation_data(
        annotation_results_dir,
    ):
    """ Load copy number tables
    
    Args:
        annotation_results_dir (str): path to hmmcopy results.
    
    Returns:
        dict: pandas.DataFrame tables keyed by table name
    """

    manifest_filename = os.path.join(annotation_results_dir, 'metadata.yaml')
    manifest = yaml.load(open(manifest_filename))

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

        data = pd.read_csv(filepath)

        data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
        data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

        for col in _categorical_cols:
            if col in data:
                data[col] = pd.Categorical(data[col])

        results_tables[table_name] = data

    scgenome.utils.union_categories(results_tables.values())

    _table_fixes[version](results_tables)

    return results_tables

