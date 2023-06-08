import numpy as np
import pandas as pd
import Data
from os import sep as sep
import os

root_path = r'Y:\Lior&Einav\Experiments'
exp_details = pd.read_excel(root_path + sep + 'Experiments_details.xlsx', engine='openpyxl')


def get_exp_path(exp_num, root_path):
    folderlist = os.listdir(root_path)
    exp_folder = [x for x in folderlist if x.startswith('experiment'+str(exp_num))]
    return root_path + sep + exp_folder[0]


def get_PC_ratios_dict(exp_details, row_idx):
    PC_ratios_dict = {}
    for color in ['yellow', 'red']:
        PC_ratios_dict[color] = [exp_details.loc[row_idx, 'P_'+color], exp_details.loc[row_idx, 'C_'+color]]
    return PC_ratios_dict


def get_bdata_filename(exp_path):
    bdata_filename = [filename for filename in os.listdir(exp_path+sep+r'with food\blob analysis normalized by white paper')
                      if filename.startswith(r'bdata')][0]
    return bdata_filename


def get_experiment_data(exp_num, exp_details, root_path, files=('bdata', 'fdata', 'tdata', 'crops',
                                                                'conversion_factors', 'transparency', 'ants', 'PC_ratios')):

    exp_path = get_exp_path(exp_num, root_path)
    exp_idx = exp_details.index[exp_details.Experiment==exp_num][0]

    data_dict = {}
    for file in files:
        match file:
            case 'bdata':
                bdata_filename = get_bdata_filename(exp_path)
                data = pd.read_csv(exp_path + sep + bdata_path + sep + bdata_filename)
            case 'fdata':
                data = pd.read_csv(exp_path + sep + 'forager_table_with_feeding_sizes_ul_transparency_corrected.csv')
            case 'tdata':
                data = pd.read_csv(exp_path + sep + r'clean_trophallaxis_table_transparency_corrected.csv')
            case 'crops':
                data = pd.read_csv(exp_path + sep + r'clean_crops_transparency_corrected.csv', header=[0, 1], index_col=[0]).swaplevel(axis=1)
            case 'conversion_factors':
                data = pd.read_csv(exp_path + sep + r'conversion_factors_by_weight_and_feeding_sum.csv')
            case 'transparency':
                data = pd.read_csv(exp_path + sep + r'transparency_table.csv')
            case 'ants':
                data = pd.read_csv(exp_path + sep + r'ants_list.csv')
                data_dict['foragers'] = [ant['ant_id'] for idx, ant in ants.iterrows() if ant['is_forager']]
            case 'PC_ratios':
                data = get_PC_ratios_dict(exp_details, exp_idx)

        data_dict[file] = data

    return exp_idx, exp_path, data_dict
