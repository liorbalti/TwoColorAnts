import Data
import pandas as pd
from os import sep as sep
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

exp11 = Data.ExperimentData(11, bdata_path='blob analysis normalized by white paper')

# ant_id = 535
pdf_foragers = PdfPages(exp11.exp_path + sep + 'forager_timelines.pdf')
pdf_workers = PdfPages(exp11.exp_path + sep + 'nonforager_timelines.pdf')
for ant_id in exp11.ants_list:
    if ant_id in exp11.foragers_list:
        ant = Data.ForagerData(ant_id, feedings_df=exp11.feedings_df, bdata_df=exp11.bdata,
                               interactions_df=exp11.interactions_df)
        pdf = pdf_foragers
    else:
        ant = Data.AntData(ant_id, bdata_df=exp11.bdata, interactions_df=exp11.interactions_df)
        pdf = pdf_workers

    fig = ant.plot_clean_timeline(times_s=exp11.bdata[['time']], conversion_factors=exp11.conversion_factors_by_weights_df, show=False)
    #plt.show()
    pdf.savefig(fig)
pdf_foragers.close()
pdf_workers.close()
#plt.savefig(f"Y:\\Lior&Einav\\Experiments\\experiment11_140720\\ant_{ant_id}.tif", format="tiff")