import csv
import shutil
import urllib.request as request
from collections import OrderedDict
from contextlib import closing
from pyteomics import fasta
import pandas as pd
import requests
from pathlib import Path
import gzip
import subprocess
import os
import re

from search_engine import SearchEngineInterface
from comet import Comet


class PeptideIdentificationPipeline:
    """
    This class represents the peptide identification pipeline beginning
    with locating the .sdrf file from PRIDE accession and finishing with
    determining fdr on search engine output.
    """
    # option to disable console output from external .exe
    # is used for thermorawfileparser.exe
    _FNULL = open(os.devnull, 'w')
    # path to .fasta overview file (README) from uniprot
    UNIPROT_FASTA_BASIS_PATH = \
        'https://ftp.uniprot.org/pub/databases/uniprot/' \
        'current_release/knowledgebase/reference_proteomes'
    # path to .fasta map file
    FASTA_MAP_FILE = 'fasta_map.csv'
    # using same .fasta database for every sample in .sdrf file
    use_same_fasta = False
    # current working directory
    _cwd = Path(__file__).parent.absolute()
    # .sdrf file count that is most of the times one
    sdrf_files = 0
    # .sdrf files necessarily need several information to execute and reproduce
    # search engine results
    prerequisite_sdrf_cols = ['comment[cleavage agent details]',
                              'comment[precursor mass tolerance]',
                              'comment[fragment mass tolerance]',
                              'comment[number of missed cleavages]',
                              'characteristics[organism]']
    fasta_map = None

    def __init__(self, accession: str, search_engine: SearchEngineInterface,
                 thermorawfileparser_path: str, search_engine_specific_fdr=False,
                 separate_sdrf_entries=False):
        self._accession = accession
        # possibility to split different 'source name's in separate directories
        self._separate_sdrf_entries = separate_sdrf_entries
        # used search engine for this pipeline
        self._search_engine = search_engine
        # tested ThermoRawFileParser: v11.3.2 (https://github.com/compomics/ThermoRawFileParser)
        self._thermorawfileparser_path = thermorawfileparser_path
        # Use FDR calculation given in this pipeline or given by the search engine.
        # A general approach to calculate the FDR based on .mztab or .mzid would be
        # recommended, but needs a current implementation on this file types. That
        # should point to implementations e.g. by pyOpenMS or pyteomics.
        self._use_search_engine_specific_fdr = search_engine_specific_fdr
        # The given search engine must implement methods given by SearchEngineInterface,
        # which provides the search() method from any search engine.
        assert issubclass(search_engine, SearchEngineInterface), \
            'Search engine does does not implement search engine interface!'

    def fdr(self):
        """
        Performs False Discovery Rate (FDR) on .mzid or .mztab file
        given by search engine.
        """
        # TODO: Implement FDR on .mzid or .mztab (i.e. from pyOpenMS or pyteomics)
        # https://pyteomics.readthedocs.io/en/latest/api/mzid.html
        # did not work for me: mzid.filter()
        pass

    def create_fasta_map(self, readme_url: str):
        """

        :param readme_url:
        """
        # create directory for databases
        Path('databases').mkdir(parents=False, exist_ok=True)

        with closing(request.urlopen(readme_url)) as r:
            with open('fasta_cache.txt', 'wb') as f:
                shutil.copyfileobj(r, f)

        with open('fasta_cache.txt', 'r') as f:
            content = f.readlines()

        os.remove('fasta_cache.txt')

        start_idx = -1
        for idx, line in enumerate(content):
            if 'Proteome_ID\t' in line:
                start_idx = idx
                break

        end_idx = -1
        for idx in range(start_idx + 1, len(content)):
            if content[idx] == '\n':
                end_idx = idx
                break

        data = pd.DataFrame(columns=[col.replace('\n', '')
                                     for col in content[start_idx].split('\t')],
                            data=[line.replace('\n', '').split('\t')
                                  for line in content[start_idx + 1:end_idx]])

        data.to_csv('databases/fasta_map.csv', sep=';', index=False)

    def _receive_fasta_database(self, organism: str) -> str:
        # find organism in .fasta map while ignoring case
        relevant_fasta_databases = \
            self.fasta_map[self.fasta_map['Species Name'].str.contains(f'(?i){organism}')]
        relevant_fasta_databases.reset_index(inplace=True)

        if len(relevant_fasta_databases) > 0:
            while True:
                print('\nChoose your FASTA database by given index:')
                for idx, entry in relevant_fasta_databases.iterrows():
                    print('\t({})\t Proteome ID: {}\t Taxonomy ID: {}\tsuperregnum: {}\tspecies name: {}'.format(
                        idx + 1, entry['Proteome_ID'], entry['Tax_ID'], entry['SUPERREGNUM'], entry['Species Name']))

                selection = input('Selection: ')
                if selection.isnumeric() and 0 < int(selection) <= len(relevant_fasta_databases):
                    selection = int(selection) - 1
                    break
                else:
                    print('Please only choose a valid index!')

            fasta_db_info = relevant_fasta_databases.iloc[selection]
            fasta_name = 'databases/{}'.format(fasta_db_info['Proteome_ID'])
            fasta_url  = '{}/{}/{}/{}_{}.fasta.gz'.format(
                self.UNIPROT_FASTA_BASIS_PATH, fasta_db_info['SUPERREGNUM'].capitalize(),
                fasta_db_info['Proteome_ID'], fasta_db_info['Proteome_ID'],
                fasta_db_info['Tax_ID'])

            if not Path(f'{fasta_name}_decoy.fasta').is_file():
                if not Path(f'{fasta_name}.fasta').is_file():
                    with request.urlopen(fasta_url) as response, \
                            open(f'{fasta_name}.fasta', 'wb') as out_file:
                        with gzip.GzipFile(fileobj=response) as uncompressed:
                            shutil.copyfileobj(uncompressed, out_file)
                # remove old decoy database (otherwise it appends data)
                # and add decoys to database
                if Path(f'{fasta_name}_decoy.fasta').is_file():
                    os.remove(f'{fasta_name}_decoy.fasta')
                fasta.write_decoy_db(source=f'{fasta_name}.fasta',
                                     output=f'{fasta_name}_decoy.fasta')

            while True:
                answer = \
                    input('Do you want to keep this FASTA database for every sample (y/n)? ')
                if answer in ['y', 'n']:
                    if answer == 'y':
                        self.use_same_fasta = True
                    break

            return f'{fasta_name}_decoy.fasta'
        else:
            return None

    def start(self):
        """
        Main methode that starts the pipeline to go from PRIDE accession to
        FDR calculation on search engine results.
        """
        # check whether ThermoRawFileParser.exe is available
        assert Path(self._thermorawfileparser_path).is_file(), \
            'Executable to ThermoRawFileParser not found.'

        # access PRIDE repository by accession
        req = requests.get(f'https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={self._accession}',
                           headers={'Accept': 'application/json'})

        # check whether API request contains success notification
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
        # load .fasta map to acquire according .fasta darabase(s)
        if not Path(f'databases/{self.FASTA_MAP_FILE}').is_file():
            self.create_fasta_map(f'{self.UNIPROT_FASTA_BASIS_PATH}/README')
        self.fasta_map = pd.read_csv(f'databases/{self.FASTA_MAP_FILE}', sep=';')

        # download and iterate every .sdrf file found for appropriate accession
        for idx, file in enumerate(files):
            print(f'Processing .sdrf file ({idx + 1}/{len(files)}) from {self._accession}')
            with closing(request.urlopen(file)) as r:
                # download .sdrf
                sdrf = '{}/{}'.format(self._accession, file.split('/')[-1])
                with open(sdrf, 'wb') as f:
                    shutil.copyfileobj(r, f)
                with open(sdrf, 'r') as csv_file:
                    entries = [entry for entry in csv.DictReader(csv_file, delimiter='\t')]

                    if not all([col in list(entries[0].keys())
                                for col in self.prerequisite_sdrf_cols]):
                        print('ERROR: SDRF file does not provide prerequisite information '
                              '(column names). SDRF file is skipped!')
                        break

                col_names = list(entries[0].keys())
                for entry_idx, entry in enumerate(entries):
                    sdrf_infos = self._read_config_sdrf(entry, col_names)
                    #
                    if not self.use_same_fasta or entry_idx == 0:
                        fasta_db = self._receive_fasta_database(sdrf_infos['organism'])
                    if fasta_db:
                        # if needed, create new directory
                        if self._separate_sdrf_entries:
                            sample_name = '/{}'.format(sdrf_infos['name'])
                            Path(f'{self._accession}{sample_name}').mkdir(parents=False, exist_ok=True)
                        else:
                            sample_name = ''
                        file_name = sdrf_infos['file name']
                        print('\nProcessing {} ({}/{})'.format(sdrf_infos['name'],
                                                             entry_idx + 1, len(entries)))
                        # download .raw file
                        with closing(request.urlopen(sdrf_infos['uri'])) as r:
                            with open(f'{self._accession}{sample_name}/{file_name}', 'wb') as f:
                                shutil.copyfileobj(r, f)
                        # convert .raw to .mgf using ThermoRawFileParser.exe
                        arguments = f'{self._thermorawfileparser_path} ' \
                                    f'-i={self._accession}{sample_name}/{file_name} ' \
                                    f'-o={self._accession}{sample_name} ' \
                                    f'-f=0'
                        subprocess.call(arguments, stdout=self._FNULL, stderr=self._FNULL, shell=False)
                        # determine path to created .mgf file
                        mgf_file = '{}{}/{}'.format(self._accession,
                                                    sample_name,
                                                    file_name.replace('raw', 'mgf'))
                        # start search engine search
                        self._search_engine.search(cwd=self._cwd,
                                                   database=fasta_db,
                                                   sdrf_entry=sdrf_infos,
                                                   mgf_file=mgf_file)
                        # perform FDR on results
                        if self._use_search_engine_specific_fdr:
                            self._search_engine.fdr()
                        else:
                            self.fdr()
                    else:
                        print('There is no associated FASTA database.')

    @staticmethod
    def _read_config_sdrf(information: OrderedDict, col_names: list) -> dict:
        """
        Reads and attracts configuration from .sdrf entry.
        :param information: Information from .sdrf entry.
        :param col_names: Actual column names from given .sdrf file.
        :return: Processed configuration from .sdrf entry.
        """
        return {
            # these four entries are mandatory (because they are necessary for search engine)
            'enzyme':
                re.search('NT=(.*?);', information['comment[cleavage agent details]']).group(1),
            'precursor mass tolerance':
                information['comment[precursor mass tolerance]'].replace(' ppm', ''),
            'fragment mass tolerance':
                information['comment[fragment mass tolerance]'].replace(' Da', ''),
            'missed cleavages':
                information['comment[number of missed cleavages]'],
            'organism':
                information['characteristics[organism]'],
            # these are additional information
            'experiment name':
                information['source name']
                if 'source name' in col_names else 'Undefined',
            'name':
                information['assay name']
                if 'assay name' in col_names else 'Undefined',
            'file name':
                information['comment[data file]']
                if 'comment[data file]' in col_names else 'Undefined',
            'uri':
                information['comment[file uri]']
                if 'comment[file uri]' in col_names else 'Undefined',
            'fraction id':
                information['comment[fraction identifier]']
                if 'comment[fraction identifier] ' in col_names else 'Undefined',
            'label':
                information['comment[label]']
                if 'comment[label]' in col_names else 'Undefined',
            'instrument':
                information['comment[instrument]']
                if 'comment[instrument]' in col_names else 'Undefined'
        }


if __name__ == '__main__':
    comet = Comet('executables/search_engines/comet.exe')
    pep_ident_pipeline = PeptideIdentificationPipeline(accession='PXD002171',
                                                       search_engine=comet,
                                                       thermorawfileparser_path=
                                                       'executables/ThermoRawFileParser/ThermoRawFileParser.exe',
                                                       search_engine_specific_fdr=True,
                                                       separate_sdrf_entries=True)
    pep_ident_pipeline.start()
