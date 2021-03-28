from pathlib import Path
from search_engine import SearchEngineInterface
import subprocess
import os
import pandas as pd
import numpy as np


class Comet(SearchEngineInterface):
    # option to disable console output from external .exe
    _FNULL = open(os.devnull, 'w')
    _param_file = None
    _mgf_file = None

    def __init__(self, path_to_executable):
        self._path_to_executable = path_to_executable

        assert Path(self._path_to_executable).is_file(), \
            'Executable to comet search engine not found. Make sure to name it \'comet.exe\'.'

    @staticmethod
    def _adjust_param_file(params: list, search_params: dict) -> list:
        """
        Takes a newly created comet.params (to be more accurate: comet.params.new)
        file as list (one row == one entry) and adjusts the parameter file
        according to the given .sdrf entry.
        :param params: List of rows from parameter file
        :param search_params: Parameters to be adjusted in the original parameter file.
        :return: List of adjusted parameters (removed comments behind '#')
        """
        params_left_to_change = len(search_params)
        for idx in range(len(params)):
            for key in search_params.keys():
                if key in params[idx]:
                    params[idx] = f'{key} = {search_params[key]}\n'
                    params_left_to_change -= 1
                    if params_left_to_change == 0:
                        break

        return params

    @staticmethod
    def _enzym_number(enzyme: str) -> str:
        """
        Matches cleavage enzyme to number in comet.exe (see bottom of comet.params)
        :param enzyme: cleavage enzyme
        :return: number of cleavage enzyme for comet.params
        """
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

    def _create_comet_param_file(self, cwd):
        """
        Creates comet parameter file that is mandatory to execute search on comet.
        """
        # create comet parameter file (see README.txt from comet)
        arguments = f'{self._path_to_executable} -p'
        subprocess.call(arguments, stdout=self._FNULL, stderr=self._FNULL, shell=False)

        # comet.params.new is created in directory of PeptideIdentificationPipeline, so
        # the parameter file will be renamed 'comet.params' as expected by comet.exe
        self._param_file = f'{cwd}/comet.params'
        if Path(self._param_file).is_file():
            os.remove(self._param_file)
        os.rename(f'{cwd}/comet.params.new', self._param_file)

    def search(self, cwd: str, database: str, sdrf_entry: dict, mgf_file: str):
        """
        Start comet search with specific .fasta database, .sdrf entry and .mgf file.
        :param cwd: current working directory
        :param database: .fasta database suitable for given PRIDE accession.
        :param sdrf_entry: Entry from .sdrf file with meta information.
        :param mgf_file: .mgf file that contains data to be investigated.
        """
        self._mgf_file = mgf_file
        self._create_comet_param_file(cwd)

        # read comet.params
        with open(self._param_file, 'r') as f:
            content = f.readlines()

        # determine adjusted parameter and rewrite to content
        params = self._adjust_param_file(content,
                                         {'database_name': f'{cwd}/{database}',
                                          'peptide_mass_tolerance': sdrf_entry['precursor mass tolerance'],
                                          'search_enzyme_number': self._enzym_number(sdrf_entry['enzyme']),
                                          'allowed_missed_cleavage': sdrf_entry['missed cleavages'],
                                          'remove_precursor_tolerance': sdrf_entry['fragment mass tolerance'],
                                          'output_txtfile': '1',
                                          'output_pepxmlfile': '0',
                                          'output_mzidentmlfile': '1'})

        # overwrite parameter file
        with open(self._param_file, 'w') as f:
            f.writelines(params)

        # e.g. comet.exe C:/Users/../ComputationalProteomics/PXD002171/OEI06439.mgf
        arguments = f'{cwd}/{self._path_to_executable} {self._mgf_file}'
        subprocess.call(arguments, shell=False)

    def _fdr_on_txt_output(self):
        """
        Perform FDR on .txt output from comet.exe (must be set in parameter file)
        """
        output_file = self._mgf_file.replace('.mgf', '')

        if Path(f'{output_file}.txt').is_file():
            data = pd.read_csv(f'{output_file}.txt', delimiter='\t', skiprows=[0])
            # see 'A deeper look into Comet – implementation and features'
            data.sort_values(['e-value'], inplace=True)
            data.reset_index(inplace=True)
            # add column that indicates whether protein is decoy or not
            data['is_decoy'] = np.where(data['protein'].str.contains('DECOY'), 1, 0)

            # perform row-wise calculation on numpy array (for best performance)
            is_decoy = data['is_decoy'].to_numpy()
            decoy_count = np.zeros(len(is_decoy))

            # cumulate decoys for each PSM
            for entry_idx in range(1, len(data)):
                decoy_count[entry_idx] = decoy_count[entry_idx - 1] + is_decoy[entry_idx]

            data['decoy_count'] = decoy_count
            # calculate FDR
            data['fdr'] = data['decoy_count'] / (data.index + 1 - data['decoy_count'])

            # write dataframe after adding colums and calculating FDR to a new file
            data.to_csv(f'{output_file}_fdr.csv', sep=';', index=False)
            # write cleaned up data (remove decoys and PSM with FDR > 0.01)
            data = data[(~data['protein'].str.contains('DECOY')) & (data['fdr'] <= 0.01)]
            data.to_csv(f'{output_file}_fdr_cleaned.csv', sep=';', index=False)
        else:
            print(f'{output_file} not found. FDR is skipped.')

    def fdr(self):
        """
        Perform search engine specific FDR
        :return:
        """
        # an search engine specific fdr option performed on the outputted .txt
        self._fdr_on_txt_output()
