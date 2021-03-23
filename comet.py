from pathlib import Path
from search_engine import SearchEngineInterface
import subprocess
import os
import pandas as pd
import numpy as np


class Comet(SearchEngineInterface):
    _cwd = None
    _separate_sdrf_entries = False
    # option to disable console output from external .exe
    _FNULL = open(os.devnull, 'w')
    _param_file = None
    _mgf_file = None

    def __init__(self, path_to_executable):
        self._path_to_executable = path_to_executable

        assert Path(self._path_to_executable).is_file(), \
            'Executable to comet search engine not found. Make sure to name it \'comet.exe\'.'

    @property
    def cwd(self):
        return self._cwd

    @cwd.setter
    def cwd(self, value):
        self._cwd = value

    @property
    def separate_sdrf_entries(self):
        return self._separate_sdrf_entries

    @separate_sdrf_entries.setter
    def separate_sdrf_entries(self, value):
        self._separate_sdrf_entries = value

    @staticmethod
    def _adjust_param_file(params: list, search_params: dict) -> list:
        params_left_to_change = len(search_params)
        for idx in range(len(params)):
            for key in search_params.keys():
                if key in params[idx]:
                    params[idx] = '{} = {}\n'.format(key, search_params[key])
                    params_left_to_change -= 1
                    if params_left_to_change == 0:
                        break

        return params

    @staticmethod
    def _enzym_number(enzyme: str) -> str:
        return {
            'Trypsin': '1',
            'Trypsin/P': '2',
            'Lys_C': '3',
            'Lys_N': '4',
            'Arg_C': '5',
            'Asp_N': '6',
            'CNBr': '7',
            'Glu_C': '8',
            'PepsinA': '9',
            'Chymotrypsin': '10',
        }[enzyme]

    def _create_comet_param_file(self):
        # create comet parameter file and move/rename it (see README.txt from comet)
        # this parameter file will be created once and edited for every experiment
        arguments = f'{self._path_to_executable} -p'
        subprocess.call(arguments, stdout=self._FNULL, stderr=self._FNULL, shell=False)

        self._param_file = f'{self._cwd}/comet.params'
        if Path(self._param_file).is_file():
            os.remove(self._param_file)
        os.rename(f'{self._cwd}/comet.params.new', self._param_file)

    def search(self, database: str, sdrf_entry: dict, mgf_file: str):
        self._mgf_file = mgf_file
        self._create_comet_param_file()

        with open(self._param_file, 'r') as f:
            content = f.readlines()

        params = self._adjust_param_file(content,
                                         {'database_name': f'{self._cwd}/{database}',
                                          'peptide_mass_tolerance': sdrf_entry['precursor mass tolerance'],
                                          'search_enzyme_number': self._enzym_number(sdrf_entry['enzyme']),
                                          'allowed_missed_cleavage': sdrf_entry['missed cleavages'],
                                          'remove_precursor_tolerance': sdrf_entry['fragment mass tolerance'],
                                          'output_txtfile': '1',
                                          'output_pepxmlfile': '0',
                                          'output_mzidentmlfile': '1'})

        with open(self._param_file, 'w') as f:
            f.writelines(params)

        # e.g. comet.exe C:/Users/Max/PycharmProjects/ComputationalProteomics/PXD002171/OEI06439.mgf
        # e.g. comet.exe C:/Users/../ComputationalProteomics/PXD002171/OEI06439.mgf
        arguments = f'{self._cwd}/{self._path_to_executable} {self._mgf_file}'
        subprocess.call(arguments, shell=False)

    def _fdr_on_txt_output(self):
        output_file = self._mgf_file.replace('.mgf', '')

        data = pd.read_csv(f'{output_file}.txt', delimiter='\t', skiprows=[0])
        # see 'A deeper look into Comet – implementation and features'
        data.sort_values(['e-value'], inplace=True)
        #
        data.reset_index(inplace=True)
        # add column that indicates whether protein is decoy or not
        data['is_decoy'] = np.where(data['protein'].str.contains('DECOY'), 1, 0)

        is_decoy = data['is_decoy'].to_numpy()
        decoy_count = np.zeros(len(is_decoy))

        for entry_idx in range(1, len(data)):
            decoy_count[entry_idx] = decoy_count[entry_idx - 1] + is_decoy[entry_idx]

        data['decoy_count'] = decoy_count
        data['fdr'] = data['decoy_count'] / (data.index + 1 - data['decoy_count'])

        data.to_csv(f'{output_file}_fdr.csv', sep=';')

    def fdr(self):
        # another option that did not work for:
        # https://pyteomics.readthedocs.io/en/latest/api/mzid.html
        self._fdr_on_txt_output()
