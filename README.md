# A repository of code for analyzing single cell genomes

## Installation

It is recommended that you install all prerequisites with pip in a virtual environment:

```
virtualenv venv
source venv/bin/activate
pip install numpy cython
pip install -r requirements.txt
python setup.py develop
```

Note that you will have to install numpy and cython prior to other requirements.

## Analyses

### Cell cycle

The `compute_cell_cycle_state.py` script can be used to run the cell cycle classifier and upload results into tantalus.

### Extract cellenone features

The `extract_cellenone_features.py` script searches for unprocessed cellenone data and extracts a table of features (Diameter, Elongation, Circularity) from the cellenone tables.

### Clonal inference

The `infer_clones.py` script can be used to run the full clonal analysis on a library or set of libraries.

#### Usage

The `infer_clones.py` script runs in 3 stages: `retrieve-cn`, `cluster-cn`, and `pseudobulk-analysis`, to be run in that order.

The CLI requires you to select the stage, and also specify the results prefix and the local storage directory for caching files.  Results of the analysis will be stored in files named with the results prefix and tantalus data will be cached in the local storage directory.

The `rretrieve-cn` stage requests metadata from tantalus and caches the data locally.  This stage requires one or more library ids and sample ids.  Copy number tables will be stored with the results prefix.

The `cluster-cn` stage runs the current copy number clustering, produces tables including the cluster labels, and distance from each cell to each cluster.  Additional plots will be output to files with the results prefix.

The `pseudobulk-analysis` stage runs the current pseudobulk analyses on a given raw pseudobulk run.  This includes:
- an analysis of the data as bulk with plots that include mutation signatures and genome wide SNV frequencies
- SNV phylogenetic analysis of clone specific SNVs
- Inference of allele specific copy number
Tables and plots are output with the results prefix.

## Scripts

Filter cells and bins:
```
filter_copynumber copynumber_matrix.tsv copynumber_matrix_filt.tsv \
    --cell-scores all_metrics_summary_classified.csv \
    --qdnaseq-blacklist QDNAseq_1.14.0_500kb_blacklist.tsv
```

Same, but also remove consecutive bins with the same copy number across cells:
```
filter_copynumber copynumber_matrix.tsv copynumber_matrix_filt.tsv \
    --cell-scores all_metrics_summary_classified.csv \
    --qdnaseq-blacklist QDNAseq_1.14.0_500kb_blacklist.tsv \
    --filter-contig-dup-bins
```

Cluster cells:
```
cluster_cells copynumber_matrix_filt.tsv \
    copynumber_cell_clusters.tsv \
    --plot copynumber_cell_clusters.pdf
```

![cell cluster scatterplot](https://user-images.githubusercontent.com/381464/45980923-56f2b300-c021-11e8-9b0e-9dcf4b53f9c7.png)

There is a sample Snakemake file included in the pipelines directory. You can run it like this:
```
snakemake --config \
    copynumber_matrix=cn_matrix.csv \
    classified=all_metrics_summary_classified.csv \
    qdnaseq_blacklist=QDNAseq_1.14.0_500kb_blacklist.tsv
```
where classified is a file with cell classifications, and qndaseq_blacklist is a tsv file generated from QDNAseq with the following columns:
* chromosome
* start
* end
* bases
* gc
* mappability
* blacklist
* residual
* use

## API

The API allows access to both HMMCopy and Pseudobulk data stored in blob and managed by tantalus.

### Prerequisites

#### Software

You should set up an environment with the requirements from `requirements.txt` including sisyphus.

#### Accounts

You must have accounts for colossus and tantalus, as well as access to azure blob storage. The list of credentials below are required in your environment:
```
TANTALUS_API_PASSWORD
TANTALUS_API_USERNAME
CLIENT_ID
TENANT_ID
SECRET_KEY
COLOSSUS_API_PASSWORD
COLOSSUS_API_USERNAME
```

### HMMCopy Data

The following example code snippet will provide access to HMMCopy data for the OV cell line data:

```
import dbclients
import scgenome.utils

from scgenome.loaders.qc import load_qc_data
from scgenome.db.qc import cache_qc_results


tantalus_api = dbclients.tantalus.TantalusApi()

hmmcopy_tickets = [
    'SC-1935',
    'SC-1936',
    'SC-1937',
]

sample_ids = [
    'SA1090',
    'SA921',
    'SA922',
]

local_cache_directory = '/Your/Local/Cache'

cn_data = []
segs_data = []
metrics_data = []
align_metrics_data = []

for jira_ticket in hmmcopy_tickets:
    analysis = tantalus_api.get(
        'analysis',
        analysis_type__name='hmmcopy',
        jira_ticket=jira_ticket)

    cache_qc_results(jira_ticket, local_cache_directory)
    ticket_directory = os.path.join(local_cache_directory, ticket_id)
    hmmcopy_data = load_qc_data(ticket_directory)

    cn_data.append(hmmcopy_data['hmmcopy_reads'])
    segs_data.append(hmmcopy_data['hmmcopy_segs'])
    metrics_data.append(hmmcopy_data['hmmcopy_metrics'])
    align_metrics_data.append(hmmcopy_data['align_metrics'])

cn_data = scgenome.utils.concat_with_categories(cn_data)
segs_data = scgenome.utils.concat_with_categories(segs_data)
metrics_data = scgenome.utils.concat_with_categories(metrics_data)
align_metrics_data = scgenome.utils.concat_with_categories(align_metrics_data)
```

### Pseudobulk Data

The following example code snippet will provide access to pseudobulk SNV, allele and breakpoint data for the OV cell line data:

```
import scgenome.snvdata
import scgenome.loaders.snv
import scgenome.loaders.allele
import scgenome.loaders.breakpoint

import dbclients.tantalus
import datamanagement.transfer_files

ticket_id = 'SC-1939'

results_prefix = './results'

local_cache_directory = '/Your/Local/Cache'

# Download the results

tantalus_api = dbclients.tantalus.TantalusApi()

ticket_results = tantalus_api.list('results', analysis__jira_ticket=ticket_id)

for results in ticket_results:
    filepaths = datamanagement.transfer_files.cache_dataset(
        tantalus_api,
        results['id'],
        'resultsdataset',
        'singlecellresults',
        local_cache_directory,
    )

# Load from cache

museq_score_threshold = None
strelka_score_threshold = None
snvs_num_cells_threshold = 2
snvs_sum_alt_threshold = 2

ticket_directory = os.path.join(local_cache_directory, ticket_id)

snv_results = scgenome.loaders.snv.load_snv_data(
    ticket_directory,
    museq_filter=museq_score_threshold,
    strelka_filter=strelka_score_threshold,
)

allele_results = scgenome.loaders.allele.load_haplotype_allele_data(
    ticket_directory,
)

breakpoint_results = scgenome.loaders.breakpoint.load_breakpoint_data(
    ticket_directory,
)

```

For additional filtering and annotation see `scgenome.analyses.infer_clones.retrieve_pseudobulk_data`.

### Testing on ceto/juno

```
(venv) -bash-4.2$ bsub -Is -R "rusage[mem=50]select[type==CentOS7]" python scgenome/tests/test_load_qc.py test-cached-single-ticket SC-2140 --local_storage_name juno
Job <23793165> is submitted to default queue <general>.
<<Waiting for dispatch ...>>
<<Starting on ja10>>
2019-11-01 14:27:36,999 - INFO - table align_metrics has size 782
2019-11-01 14:27:37,000 - INFO - table hmmcopy_reads has size 4853092
2019-11-01 14:27:37,002 - INFO - table hmmcopy_segs has size 413918
2019-11-01 14:27:37,003 - INFO - table hmmcopy_metrics has size 782
2019-11-01 14:27:37,004 - INFO - table annotation_metrics has size 775
2019-11-01 14:27:37,009 - WARNING - fastqscreen_grch37_multihit not in table annotation_metrics
2019-11-01 14:27:37,010 - WARNING - fastqscreen_grch37 not in table annotation_metrics
2019-11-01 14:27:37,011 - WARNING - fastqscreen_mm10_multihit not in table annotation_metrics
2019-11-01 14:27:37,011 - WARNING - fastqscreen_mm10 not in table annotation_metrics
2019-11-01 14:27:37,012 - WARNING - fastqscreen_nohit not in table annotation_metrics
2019-11-01 14:27:37,013 - WARNING - fastqscreen_salmon_multihit not in table annotation_metrics
2019-11-01 14:27:37,014 - WARNING - fastqscreen_salmon not in table annotation_metrics
2019-11-01 14:27:37,014 - WARNING - grch37_multihit not in table annotation_metrics
2019-11-01 14:27:37,015 - WARNING - grch37 not in table annotation_metrics
2019-11-01 14:27:37,018 - WARNING - is_contaminated not in table annotation_metrics
2019-11-01 14:27:37,020 - WARNING - mm10_multihit not in table annotation_metrics
2019-11-01 14:27:37,021 - WARNING - mm10 not in table annotation_metrics
2019-11-01 14:27:37,022 - WARNING - nohit not in table annotation_metrics
2019-11-01 14:27:37,023 - WARNING - order_corrupt_tree not in table annotation_metrics
2019-11-01 14:27:37,026 - WARNING - salmon_multihit not in table annotation_metrics
2019-11-01 14:27:37,026 - WARNING - salmon not in table annotation_metrics
2019-11-01 14:27:37,030 - INFO - successfully loaded results from /juno/work/shah/tantalus/SC-2140
```

```
(venv) -bash-4.2$ bsub -Is -R "rusage[mem=50]select[type==CentOS7]" python scgenome/tests/test_load_pseudobulk.py test-cached-single-ticket SC-2658 --local_storage_name juno
Job <23795221> is submitted to default queue <general>.
<<Waiting for dispatch ...>>
<<Starting on ja10>>
2019-11-01 15:19:29,456 - INFO - starting load
Traceback (most recent call last):
  File "scgenome/tests/test_load_pseudobulk.py", line 189, in <module>
    cli()
  File "/home/vatrtwaa/scgenome/venv/lib/python3.7/site-packages/click/core.py", line 764, in __call__
    return self.main(*args, **kwargs)
  File "/home/vatrtwaa/scgenome/venv/lib/python3.7/site-packages/click/core.py", line 717, in main
    rv = self.invoke(ctx)
  File "/home/vatrtwaa/scgenome/venv/lib/python3.7/site-packages/click/core.py", line 1137, in invoke
    return _process_result(sub_ctx.command.invoke(sub_ctx))
  File "/home/vatrtwaa/scgenome/venv/lib/python3.7/site-packages/click/core.py", line 956, in invoke
    return ctx.invoke(self.callback, **ctx.params)
  File "/home/vatrtwaa/scgenome/venv/lib/python3.7/site-packages/click/core.py", line 555, in invoke
    return callback(*args, **kwargs)
  File "scgenome/tests/test_load_pseudobulk.py", line 151, in test_cached_single_ticket
    local_storage_name=local_storage_name,
  File "scgenome/tests/test_load_pseudobulk.py", line 132, in test_load_stored_pseudobulk_data
    test_load_local_pseudobulk_data(ticket_directory)
  File "scgenome/tests/test_load_pseudobulk.py", line 88, in test_load_local_pseudobulk_data
    results_dir,
  File "/home/vatrtwaa/scgenome/scgenome/loaders/snv.py", line 274, in load_snv_data
    strelka_filter=strelka_filter)
  File "/home/vatrtwaa/scgenome/scgenome/loaders/snv.py", line 158, in load_snv_annotation_results
    mappability = load_snv_annotation_table(pseudobulk_dir, 'mappability')
  File "/home/vatrtwaa/scgenome/scgenome/loaders/snv.py", line 82, in load_snv_annotation_table
    for sample_id, library_id, filepath in scgenome.loaders.utils.get_pseudobulk_files(pseudobulk_dir, f'snv_{table_name}.csv.gz'):
  File "/home/vatrtwaa/scgenome/scgenome/loaders/utils.py", line 85, in get_pseudobulk_files
    raise ValueError(f'found {len(sample_lib_filenames)} {suffix} files for {sample_id}, {library_id}, {results_dir}')
ValueError: found 0 snv_mappability.csv.gz files for SA1090, A96213A, /juno/work/shah/tantalus/SC-2658/results/variants
```
