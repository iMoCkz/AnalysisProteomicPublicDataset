import csv
import shutil
import urllib.request as request
from collections import OrderedDict
from contextlib import closing
from pyteomics import fasta
import requests
from pathlib import Path
import pandas as pd
import gzip
import subprocess
import os
import re

from search_engine import SearchEngineInterface
from comet import Comet


class PeptideIdentificationPipeline:
    # option to disable console output from external .exe
    _FNULL = open(os.devnull, 'w')
    # current working directory
    _cwd = Path(__file__).parent.absolute()
    _fasta_url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/' \
                 'knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz'
    _fasta_name = 'UP000005640_9606.fasta'
    _fasta_decoy_name = None
    sdrf_files = 0

    def __init__(self, accession: str, search_engine: SearchEngineInterface,
                 thermorawfileparser_path: str, search_engine_specific_fdr=False,
                 separate_sdrf_entries=False):
        self._accession = accession
        self._separate_sdrf_entries = separate_sdrf_entries
        self._search_engine = search_engine
        # tested ThermoRawFileParser: v11.3.2 (https://github.com/compomics/ThermoRawFileParser)
        self._thermorawfileparser_path = thermorawfileparser_path
        self._use_search_engine_specific_fdr = search_engine_specific_fdr

        assert issubclass(search_engine, SearchEngineInterface), \
            'Search engine does does not implement search engine interface!'

        self._search_engine.separate_sdrf_entries = separate_sdrf_entries

    def fdr(self):
        # TODO: Implement FDR on .mzid oder .mztab (e.g. from pyOpenMS or pyteomics)
        pass

    def start(self):
        # check weather ThermoRawFileParser.exe is available
        assert Path(self._thermorawfileparser_path).is_file(), \
            'Executable to ThermoRawFileParser not found.'

        # access PRIDE repository by accession
        req = requests.get(f'https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={self._accession}',
                           headers={'Accept': 'application/json'})

        assert req.status_code == 200, \
            f'Unsuccessful PRIDE access via accession (HTTP response status code {req.status_code}).'

        # find .sdrf file(s)
        req = req.json()
        files = [file['value'] for accession in req
                 for file in accession['publicFileLocations']
                 if 'ftp' in file['value'] and 'sdrf' in file['value'].lower()]

        assert len(files) > 0, 'PRIDE accession does not contain SDRF file(s).'

        self.sdrf_files = len(files)

        # create directory for acession
        Path(self._accession).mkdir(parents=False, exist_ok=True)

        # download database (.fasta) for search engine
        if not Path(self._fasta_name).is_file():
            with request.urlopen(self._fasta_url) as response, open(self._fasta_name, 'wb') as out_file:
                with gzip.GzipFile(fileobj=response) as uncompressed:
                    shutil.copyfileobj(uncompressed, out_file)

        # add decoys to database
        self._fasta_decoy_name = self._fasta_name.replace('.fasta', '_decoy.fasta')
        fasta.write_decoy_db(source=self._fasta_name, output=self._fasta_decoy_name)

        # set current working directory for search engine
        self._search_engine.cwd = self._cwd

        # download and iterate every .sdrf file found for appropriate accession
        for idx, file in enumerate(files):
            print(f'Processing .sdrf file ({idx + 1}/{len(files)})')
            with closing(request.urlopen(file)) as r:
                # download .sdrf
                sdrf = '{}/{}'.format(self._accession, file.split('/')[-1])
                with open(sdrf, 'wb') as f:
                    shutil.copyfileobj(r, f)
                with open(sdrf, 'r') as csv_file:
                    entries = [entry for entry in csv.DictReader(csv_file, delimiter='\t')]
                    for entry_idx, entry in enumerate(entries):
                        sdrf_infos = self._read_config_sdrf(entry)
                        # if needed, create new directory
                        if self._separate_sdrf_entries:
                            experiment_name = '/{}'.format(sdrf_infos['name'])
                            Path(f'{self._accession}{experiment_name}').mkdir(parents=False, exist_ok=True)
                        else:
                            experiment_name = ''
                        # download .raw file
                        file_name = sdrf_infos['file name']
                        print('Processing {} ({}/{})'.format(sdrf_infos['name'],
                                                              entry_idx + 1, len(entries)))
                        with closing(request.urlopen(sdrf_infos['uri'])) as r:
                            with open(f'{self._accession}{experiment_name}/{file_name}', 'wb') as f:
                                shutil.copyfileobj(r, f)
                        # convert .raw to .mgf using ThermoRawFileParser.exe
                        arguments = f'{self._thermorawfileparser_path} ' \
                                    f'-i={self._accession}{experiment_name}/{file_name} ' \
                                    f'-o={self._accession}{experiment_name} ' \
                                    f'-f=0'
                        subprocess.call(arguments, stdout=self._FNULL, stderr=self._FNULL, shell=False)

                        mgf_file = '{}{}/{}'.format(self._accession,
                                                    experiment_name,
                                                    file_name.replace('raw', 'mgf'))
                        self._search_engine.search(database=self._fasta_decoy_name,
                                                   sdrf_entry=sdrf_infos,
                                                   mgf_file=mgf_file)

                        if self._use_search_engine_specific_fdr:
                            self._search_engine.fdr()
                        else:
                            self.fdr()

    @staticmethod
    def _read_config_sdrf(information: OrderedDict) -> dict:
        return {
            'experiment name': information['source name'],
            'organism': information['characteristics[organism]'],
            'name': information['assay name'],
            'file name': information['comment[data file]'],
            'uri': information['comment[file uri]'],
            'fraction id': information['comment[fraction identifier]'],
            'label': information['comment[label]'],
            'instrument': information['comment[instrument]'],
            'enzyme': re.search('NT=(.*?);', information['comment[cleavage agent details]']).group(1),
            'precursor mass tolerance': information['comment[precursor mass tolerance]'].replace(' ppm', ''),
            'fragment mass tolerance': information['comment[fragment mass tolerance]'].replace(' Da', ''),
            'missed cleavages': information['comment[number of missed cleavages]'],
        }


if __name__ == '__main__':
    comet = Comet('Executables/SearchEngines/comet.exe')
    pep_ident_pipeline = PeptideIdentificationPipeline(accession='PXD002171',
                                                       search_engine=comet,
                                                       thermorawfileparser_path=
                                                       'Executables/ThermoRawFileParser/ThermoRawFileParser.exe',
                                                       search_engine_specific_fdr=True,
                                                       separate_sdrf_entries=True)
    pep_ident_pipeline.start()