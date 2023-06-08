import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import sep as sep
import Data
import os
import seaborn as sns
import matplotlib.lines as mlines


def get_exp_path(exp_num, root_path):
    folderlist = os.listdir(root_path)
    exp_folder = [x for x in folderlist if x.startswith('experiment'+str(exp_num))]
    return root_path + sep + exp_folder[0]


def get_PC_ratios_dict(exp_details,row_idx):
    PC_ratios_dict = {}
    for color in ['yellow', 'red']:
        PC_ratios_dict[color] = [exp_details.loc[row_idx,'P_'+color], exp_details.loc[row_idx,'C_'+color]]
    return PC_ratios_dict


def get_exp_data(exp_num, exp_details, root_path):
    exp_path = get_exp_path(exp_num, root_path)
    exp_idx = exp_details.index[exp_details.Experiment==exp_num][0]
    PC_ratios_dict = get_PC_ratios_dict(exp_details,exp_idx)
    return exp_idx, PC_ratios_dict, exp_path

