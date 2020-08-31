import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import sep as sep
from copy import copy


def food_volume_to_PC_amounts(vol_ul, PC_ratio, concentration_mg_per_ml=100):
    """
    :param vol_ul:
    :param PC_ratio: list [P,C], for example: for P:C 1:3 --> [1,3]
    :param concentration_mg_per_ml: default 100 mg/ml
    :return: P_mg, C_mg
    """
    vol_ml = vol_ul/1000
    total_nutrients_mg = concentration_mg_per_ml*vol_ml
    P_mg, C_mg = total_nutrients_mg*np.array(PC_ratio)/np.sum(PC_ratio)
    return P_mg, C_mg


class AntData:
    def __init__(self, ant_id, bdata_df, interactions_df=None):
        self.ant_id = ant_id
        self.start_frame = bdata_df.index[0]
        self.crop_dict_raw = self.get_raw_crop(bdata_df)
        self.crop_dict_clean = {'red': None, 'yellow': None}
        self.is_forager = False
        self.x_raw, self.raw_y = self.get_raw_xy(bdata_df)
        self.x_interp = None
        self.y_interp = None
        self.interactions_df = pd.DataFrame()
        self.transparency = None
        self.ymax = {'red': 800000, 'yellow': 800000}

    def get_raw_crop(self, bdata_df):
        red_crop = bdata_df['a'+str(self.ant_id)+'-crop_intensity'][bdata_df.acquisition == 'GLRF']
        yellow_crop = bdata_df['a'+str(self.ant_id)+'-crop_intensity'][bdata_df.acquisition == 'BLGF']
        red_crop[red_crop == -1] = np.nan
        yellow_crop[yellow_crop == -1] = np.nan
        return {'red': red_crop, 'yellow': yellow_crop}

    def get_raw_xy(self, bdata_df):
        x = bdata_df['a'+str(self.ant_id)+'-x']
        y = bdata_df['a'+str(self.ant_id)+'-y']
        return x, y

    def clean_crop(self):
        pass

    def plot_raw_timeline(self):
        plt.figure(figsize=[18,4])
        plt.plot(self.crop_dict_raw['yellow'], '.y', alpha=0.4)
        plt.plot(self.crop_dict_raw['red'], '.r', alpha=0.4)
        plt.ylim([0, max(self.ymax['yellow'], self.ymax['red'])])
        plt.xlabel('frame')
        plt.ylabel('Fluorescence')
        plt.title(self.ant_id)

    def plot_heatmap(self):
        pass


class ForagerData(AntData):
    def __init__(self, ant_id, feedings_df, bdata_df=None, interactions_df=None):
        super().__init__(ant_id,bdata_df,interactions_df)
        self.is_forager = True
        self.feedings_dict = self.get_feedings_dict(feedings_df)

    def get_feedings_dict(self,feedings_df):
        feedings = feedings_df[feedings_df.ant_id == self.ant_id]
        feedings[['feeding_start', 'feeding_end', 'last_interaction_before_end', 'first_interaction_after_start']] += \
            self.start_frame
        red_feedings = feedings[feedings.food_source == 'red']
        yellow_feedings = feedings[feedings.food_source == 'yellow']
        return {'red': red_feedings, 'yellow': yellow_feedings}

    def get_feeding_sizes_intensity(self):
        for food_source in ['red', 'yellow']:
            ymax = self.ymax[food_source]
            crop = self.crop_dict_raw[food_source]
            if len(self.feedings_dict[food_source]) > 0:
                self.feedings_dict[food_source]['crop_before'] = self.feedings_dict[food_source].apply(
                    lambda row: ForagerData.get_crop_before(row,crop=crop,ymax=ymax), axis=1)
                self.feedings_dict[food_source]['crop_after'] = self.feedings_dict[food_source].apply(
                    lambda row: ForagerData.get_crop_after(row, crop=crop, ymax=ymax), axis=1)
                self.feedings_dict[food_source]['feeding_size_intensity'] = \
                    self.feedings_dict[food_source]['crop_after'] - self.feedings_dict[food_source]['crop_before']
            else:
                self.feedings_dict[food_source] = pd.DataFrame(columns=
                    ['ant_id', 'feeding_start', 'feeding_end', 'food_source', 'last_interaction_before_end',
                     'last_interaction_before_partner', 'first_interaction_after_start',
                     'first_interaction_after_partner', 'crop_before', 'crop_after', 'feeding_size_intensity'])

    @staticmethod
    def get_crop_before(feed_row, crop, ymax):
        crop_before = crop.loc[feed_row['last_interaction_before_end']:(feed_row['feeding_start']-1)]
        crop_before95 = np.nanpercentile(crop_before[crop_before <= ymax], 95)
        return crop_before95

    @staticmethod
    def get_crop_after(feed_row, crop, ymax):
        crop_after = crop.loc[feed_row['feeding_end']:feed_row['first_interaction_after_start']]
        crop_after95 = np.nanpercentile(crop_after[crop_after <= ymax], 95)
        return crop_after95

    def get_feeding_sizes_ul(self, feeding_df, conversion_factors=None):
        pass

    def plot_raw_timeline(self):
        super().plot_raw_timeline()
        for start_frame, end_frame in self.feedings_dict['yellow'][['feeding_start', 'feeding_end']].itertuples(index=False):
            plt.axvspan(start_frame, end_frame, facecolor='y', alpha=0.3)
        for start_frame, end_frame in self.feedings_dict['red'][['feeding_start', 'feeding_end']].itertuples(index=False):
            plt.axvspan(start_frame, end_frame, facecolor='r', alpha=0.3)

    def plot_timeline_around_feedings(self):
        self.plot_raw_timeline()
        colors = {'red':'#8B0000', 'yellow':'#B8860B'}
        for food_source, feeding_table in zip(self.feedings_dict.keys(), self.feedings_dict.values()):
            x = feeding_table[['last_interaction_before_end', 'feeding_start', 'feeding_end',
                               'first_interaction_after_start']].to_numpy().flatten()
            y = feeding_table[['crop_before','crop_before','crop_after','crop_after']].to_numpy().flatten()
            plt.plot(x, y, colors[food_source])
            interactions = feeding_table[['last_interaction_before_end', 'first_interaction_after_start']].to_numpy().flatten()
            plt.vlines(x=interactions, ymin=0, ymax=max(self.ymax['red'], self.ymax['yellow']),
                       linestyles='dotted', colors='k')

    def manual_correction(self,experiment_path):
        corrections = pd.read_excel(experiment_path+sep+'manual_fixes_forager_timelines.xlsx')
        forager_corrections = corrections[corrections['ant_id']==self.ant_id]
        for food_source, crop_data in self.crop_dict_raw.items():
            corrected_crop = copy(crop_data)
            frames_to_replace = forager_corrections['frame_to_replace'][forager_corrections['color'] == food_source]
            frames_of_input_values = forager_corrections['frame_to_take_value_from'][forager_corrections['color'] == food_source]
            corrected_crop[frames_to_replace] = crop_data[frames_of_input_values]
            self.crop_dict_raw[food_source] = corrected_crop


class InteractionData:
    def __init__(self):
        pass
