from Data import AntData, ForagerData, ExperimentData
from os import sep as sep

exp19 = ExperimentData(19, bdata_path='blob analysis normalized by white paper')

# for a given ant
# find frames between events
# search for windows where the smoothed derrivative has a constant sign