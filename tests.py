# import Data
# import pandas as pd
# from os import sep as sep
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Data import AntData, ForagerData, ExperimentData
from os import sep as sep

exp1 = ExperimentData(24, bdata_path='blob analysis normalized by white paper')
# enriched_interactions_df = exp1.enrich_interactions_df(write_data=True, filename='trophallaxis_table_enriched_temp_with_conf')
# clean_interactions_df = exp1.make_clean_interaction_table(write_data=True, filename='clean_trophallaxis_table_temp')
transparency_table = pd.read_csv(exp1.exp_path+sep+'transparency_table.csv', index_col='ant')
exp1.correct_data_for_transparency(transparency_table, save_corrected_files=True)

# A = AntData(1004, bdata_df=exp20.bdata, interactions_df=exp20.interactions_df)
# A.plot_raw_timeline()


# exp11 = Data.ExperimentData(11, bdata_path='blob analysis normalized by old norm_mats')
# ant_id = 289
# ant = Data.ForagerData(ant_id, feedings_df=exp11.feedings_df, bdata_df=exp11.bdata, interactions_df=exp11.interactions_df)
# jj = ant.clean_crop()

# exp_path = r"Y:\Lior&Einav\Experiments\experiment16_250820"
# feeding_filename = r"forager_feeding_table.xlsx"
# bdata_path = r"with food\blob analysis normalized by white paper"
# trophallaxis_filename = r"trophallaxis_table.csv"
# bdata_filename = "bdata_16_250820.csv"
#
#
#
# bdata = pd.read_csv(exp_path+sep+bdata_path+sep+bdata_filename)
# bdata = bdata.drop(bdata.columns[0],axis=1)
# #bdata['frame'] = bdata.index
# start_frame = bdata.frame[0]
# #bdata.index = bdata.index.to_series().apply(lambda x: x-start_frame)
#
#
# #bdata_try = pd.read_csv(exp_path+sep+bdata_path+sep+bdata_filename)
# fdata = pd.read_excel(exp_path+sep+feeding_filename)
# tdata = pd.read_csv(exp_path+sep+trophallaxis_filename)
# tdata[['general_start_frame','general_end_frame']] = tdata[['general_start_frame','general_end_frame']].apply(lambda x: x-start_frame)
# tdata = tdata.query('general_start_frame >= 0')
# foragers_list = fdata['ant_id'].unique()
#
# # forager_id = 540
# for forager_id in foragers_list:
#     F = ForagerData(ant_id=forager_id, bdata_df=bdata, feedings_df=fdata, interactions_df=tdata)
#     # F.plot_raw_timeline()
#     F.plot_clean_timeline(times_s=bdata[['time']],conversion_factors=pd.DataFrame({'red':[1], 'yellow':[1]}))



a = 1