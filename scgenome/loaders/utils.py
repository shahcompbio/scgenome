import os
import yaml


def find_manifest_filenames(
        ticket_id,
        local_cache_directory,
    ):
    ticket_directory = os.path.join(local_cache_directory, ticket_id)

    for dirpath, dirnames, filenames in os.walk(ticket_directory):
        for filename in filenames:
            if not filename == 'metadata.yaml':
                continue

            manifest_filename = os.path.join(dirpath, filename)

            yield manifest_filename


def _find_manifest_filenames(
        ticket_id,
        local_cache_directory,
    ):
    ticket_directory = os.path.join(local_cache_directory, ticket_id)

    for dirpath, dirnames, filenames in os.walk(ticket_directory):
        for filename in filenames:
            if not filename == 'metadata.yaml':
                continue

            manifest_filename = os.path.join(dirpath, filename)

            yield manifest_filename


def find_results_directories(
        ticket_id,
        local_cache_directory,
    ):
    results_directories = {}

    for manifest_filename in _find_manifest_filenames(ticket_id, local_cache_directory):
        manifest = yaml.load(open(manifest_filename))

        results_type = manifest['meta']['type']

        if results_type in results_directories:
            raise Exception('found {} and {} with results type {}'.format(
                results_directories[results_type],
                os.path.dirname(manifest_filename),
                results_type))

        results_directories[results_type] = os.path.dirname(manifest_filename)

    return results_directories
