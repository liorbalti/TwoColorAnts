import Data
import pandas as pd
from os import sep as sep
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

exp11 = Data.ExperimentData(11, bdata_path='blob analysis normalized by old norm_mats')
ant_id = 289
ant = Data.ForagerData(ant_id, feedings_df=exp11.feedings_df, bdata_df=exp11.bdata, interactions_df=exp11.interactions_df)
jj = ant.clean_crop()
a=1