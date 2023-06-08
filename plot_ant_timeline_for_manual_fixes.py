import Data
import pandas as pd
from os import sep as sep
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

exp11 = Data.ExperimentData(20, bdata_path='blob analysis normalized by white paper')
# exp11.correct_data_for_transparency()

ant_id = 424

if ant_id in exp11.foragers_list:
    ant = Data.ForagerData(ant_id, feedings_df=exp11.feedings_df, bdata_df=exp11.bdata,
                           interactions_df=exp11.interactions_df)

else:
    ant = Data.AntData(ant_id, bdata_df=exp11.bdata, interactions_df=exp11.interactions_df)


fig = ant.plot_clean_timeline(times_s=exp11.bdata[['time']], conversion_factors=exp11.conversion_factors_by_weights_df, show=True, x_axis='frame')
plt.show()
