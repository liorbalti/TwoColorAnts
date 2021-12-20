import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import sep as sep
from copy import copy
import os.path
import csv


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


def get_general_frame(frame_in_video, video_index, video_length=1001):
    general_frame = (video_index-1)*video_length + frame_in_video
    return general_frame


class AntData:
    def __init__(self, ant_id, bdata_df, interactions_df=None):
        self.ant_id = ant_id
        self.start_frame = bdata_df.loc[0,'frame']# bdata_df.index[0]
        self.crop_dict_raw = self.get_raw_crop(bdata_df)
        self.crop_dict_clean = {'red': None, 'yellow': None}
        self.is_forager = False
        self.x_raw, self.raw_y = self.get_raw_xy(bdata_df)
        self.x_interp = None
        self.y_interp = None
        self.interactions_df = self.get_interactions(interactions_df)
        self.transparency = None
        self.ymax = {'red': 800000, 'yellow': 800000}

    def get_interactions(self, all_interactions_df):
        if all_interactions_df is None:
            return None
        ant_interactions = all_interactions_df[
            (all_interactions_df.actual_ant1 == self.ant_id) | (all_interactions_df.actual_ant2 == self.ant_id)]
        return ant_interactions

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

    def crop_to_ul(self, conversion_factors):
        pass

    def unify_overlapping_interactions(self):
        end_frame = max(max(self.crop_dict_raw['red'].index), max(self.crop_dict_raw['yellow'].index))
        raw_interaction_starts = self.interactions_df.general_start_frame
        raw_interaction_ends = self.interactions_df.general_end_frame
        unified_interaction_starts = []
        unified_interaction_ends = [0]
        for f_st, f_end in zip(raw_interaction_starts, raw_interaction_ends):
            if f_st < unified_interaction_ends[-1]:
                unified_interaction_ends[-1] = max(unified_interaction_ends[-1], f_end)
            else:
                unified_interaction_starts.append(f_st)
                unified_interaction_ends.append(f_end)
        unified_interaction_starts.append(end_frame)
        return unified_interaction_starts, unified_interaction_ends

    @staticmethod
    def get_percentile_in_intervals(measurements, interval_starts, interval_ends, max_frame, percentile=90, set_first_chunk=True):
        out_dict = {}
        for ii, (i_st, i_end) in enumerate(zip(interval_starts, interval_ends)):
            frame_range = np.arange(i_st, i_end)
            if i_end == 0:
                i_end = 1
            if i_end > max_frame:
                print(f'get_percentile_in_intervals: frame {i_end} out of bounds \n'
                      f'cutting or skipping interval')
                i_end = max_frame
                if i_st <= max_frame:
                    frame_range = np.arange(i_st, i_end)
                else:
                    continue
            if i_end-i_st == 1 and i_st not in measurements.index:
                measurements_chunk = measurements[measurements.index.to_series().ge(i_st)].iloc[0]
            else:
                measurements_chunk = measurements.loc[frame_range]
            if set_first_chunk and (ii == 0):
                p = 0
            else:
                p = np.nanpercentile(measurements_chunk, percentile)
            out_dict.update({frame: p for frame in frame_range})
        return out_dict

    def clean_crop(self, percentile=90):
        end_frame = max(max(self.crop_dict_raw['red'].index), max(self.crop_dict_raw['red'].index))
        interaction_starts, interaction_ends = self.unify_overlapping_interactions()

        filtered_yellow = self.crop_dict_raw['yellow'].where(lambda x: x < self.ymax['yellow'])
        filtered_red = self.crop_dict_raw['red'].where(lambda x: x < self.ymax['red'])

        medfiltered_yellow = filtered_yellow.rolling(15,min_periods=1,center=True).median()
        medfiltered_red = filtered_red.rolling(15,min_periods=1,center=True).median()

        yellow_dict = AntData.get_percentile_in_intervals(medfiltered_yellow, interaction_ends,
                                                          interaction_starts, max_frame=end_frame, percentile=percentile)
        red_dict = AntData.get_percentile_in_intervals(medfiltered_red, interaction_ends,
                                                       interaction_starts, max_frame=end_frame, percentile=percentile)

        clean_crop_df = pd.DataFrame({'red': red_dict, 'yellow': yellow_dict})
        clean_crop_df = clean_crop_df.reindex(range(end_frame), fill_value=np.NaN)
        clean_crop_df = clean_crop_df.apply(lambda x: x.interpolate(method='polynomial', order=1))
        return clean_crop_df

    def plot_raw_timeline(self):
        fig = plt.figure(figsize=[18, 4])
        plt.plot(self.crop_dict_raw['yellow'].index,self.crop_dict_raw['yellow'], '.y', alpha=0.4)
        plt.plot(self.crop_dict_raw['red'].index,self.crop_dict_raw['red'], '.r', alpha=0.4)
        plt.ylim([0, max(self.ymax['yellow'], self.ymax['red'])])
        plt.xlabel('frame')
        plt.ylabel('Fluorescence')
        plt.title(f'Ant {self.ant_id}')
        for start_frame, end_frame in self.interactions_df[['general_start_frame', 'general_end_frame']].itertuples(index=False):
            plt.axvspan(start_frame, end_frame, facecolor='0.5', alpha=0.2, zorder=-100)
        return fig

    def plot_clean_timeline(self, times_s, conversion_factors, include_raw=True, show=False):

        clean_crop_df = self.clean_crop()

        raw_crop_ul = {}
        clean_crop_ul = pd.DataFrame()
        for c in ['yellow', 'red']:
            raw_crop_ul[c] = self.crop_dict_raw[c]/conversion_factors.loc[0, c]
            clean_crop_ul.loc[:, c] = clean_crop_df[c]/conversion_factors.loc[0, c]

        fig = plt.figure(figsize=[18, 4])
        plt.plot(raw_crop_ul['yellow'].index * 2.8 / 60, raw_crop_ul['yellow'], '.y', alpha=0.4)
        plt.plot(raw_crop_ul['red'].index * 2.8 / 60, raw_crop_ul['red'], '.r', alpha=0.4)
        plt.ylim([0, max(self.ymax['yellow'] / conversion_factors.loc[0, 'yellow'],
                         self.ymax['red'] / conversion_factors.loc[0, 'red'])])
        plt.xlabel('Time [min]')
        plt.ylabel(r'Crop load [$\mu$l]')
        plt.title(f'Ant {self.ant_id}')

        for start_frame, end_frame in self.interactions_df[['general_start_frame', 'general_end_frame']].itertuples(
                index=False):
            plt.axvspan(start_frame*2.8/60, end_frame*2.8/60, facecolor='0.5', alpha=0.2, zorder=-100)

        clean_crop_ul.set_index(clean_crop_ul.index * 2.8 / 60, inplace=True)
        clean_crop_ul.plot(ax=fig.axes[0], color=[(184 / 255, 134 / 255, 11 / 255), (139 / 255, 0, 0), ], legend=None)

        if show:
            plt.show()

        return fig

    def plot_heatmap(self):
        pass


class ForagerData(AntData):
    def __init__(self, ant_id, feedings_df, bdata_df=None, interactions_df=None):
        super().__init__(ant_id, bdata_df, interactions_df)
        self.is_forager = True
        self.feedings_dict = self.get_feedings_dict(feedings_df)

    def insert_feedings_for_crop_clean(self, interaction_starts, interaction_ends, crop_color):
        feeding_starts = self.feedings_dict[crop_color]['feeding_start'].tolist()
        feeding_ends = self.feedings_dict[crop_color]['feeding_end'].tolist()
        event_starts = interaction_starts + feeding_starts
        event_ends = interaction_ends + feeding_ends
        all_events = [event_starts, event_ends]
        sorted_events = np.sort(all_events)
        return sorted_events[0], sorted_events[1]

    def clean_crop(self, percentile=90):
        end_frame = max(max(self.crop_dict_raw['red'].index), max(self.crop_dict_raw['red'].index))
        interaction_starts, interaction_ends = self.unify_overlapping_interactions()

        filtered_yellow = self.crop_dict_raw['yellow'].where(lambda x: x < self.ymax['yellow'])
        filtered_red = self.crop_dict_raw['red'].where(lambda x: x < self.ymax['red'])

        # TODO: remove data from feeding events of other color
        red_feedings_start = self.feedings_dict['red'].feeding_start
        red_feedings_end = self.feedings_dict['red'].feeding_end
        yellow_feedings_start = self.feedings_dict['yellow'].feeding_start
        yellow_feedings_end = self.feedings_dict['yellow'].feeding_end

        for t_st, t_end in zip(red_feedings_start, red_feedings_end):
            filtered_yellow.loc[t_st:t_end] = np.nan

        for t_st, t_end in zip(yellow_feedings_start, yellow_feedings_end):
            filtered_red.loc[t_st:t_end] = np.nan


        medfiltered_yellow = filtered_yellow.rolling(15,min_periods=1,center=True).median()
        medfiltered_red = filtered_red.rolling(15,min_periods=1,center=True).median()

        clean_crop_dict = {}
        for color, crop_measurement in zip(['red', 'yellow'], [medfiltered_red, medfiltered_yellow]):
            event_starts, event_ends = self.insert_feedings_for_crop_clean(interaction_starts, interaction_ends, color)
            clean_crop_dict[color] = AntData.get_percentile_in_intervals(crop_measurement, event_ends, event_starts,
                                                                         max_frame=end_frame, percentile=percentile)

        clean_crop_df = pd.DataFrame(clean_crop_dict)
        clean_crop_df = clean_crop_df.reindex(range(end_frame), fill_value=np.NaN)
        clean_crop_df = clean_crop_df.apply(lambda x: x.interpolate(method='polynomial', order=1))
        return clean_crop_df

    def get_feedings_dict(self, feedings_df):
        #feedings_df[['feeding_start','feeding_end']] = feedings_df[['feeding_start','feeding_end']].apply(lambda x: x+self.start_frame)
        feedings = feedings_df[feedings_df.ant_id == self.ant_id]
        red_feedings = feedings[feedings.food_source == 'red']
        yellow_feedings = feedings[feedings.food_source == 'yellow']
        return {'red': red_feedings, 'yellow': yellow_feedings}

    def get_feeding_sizes_intensity_old(self):
        for food_source in ['red', 'yellow']:
            ymax = self.ymax[food_source]
            crop = self.crop_dict_raw[food_source]
            if len(self.feedings_dict[food_source]) > 0:
                self.feedings_dict[food_source]['crop_before'] = self.feedings_dict[food_source].apply(
                    lambda row: ForagerData.get_crop_before_row(row, crop=crop, ymax=ymax), axis=1)
                self.feedings_dict[food_source]['crop_after'] = self.feedings_dict[food_source].apply(
                    lambda row: ForagerData.get_crop_after_row(row, crop=crop, ymax=ymax), axis=1)
                self.feedings_dict[food_source]['feeding_size_intensity'] = \
                    self.feedings_dict[food_source]['crop_after'] - self.feedings_dict[food_source]['crop_before']
            else:
                self.feedings_dict[food_source] = pd.DataFrame(columns=
                    ['ant_id', 'feeding_start', 'feeding_end', 'food_source', 'last_interaction_before_end',
                     'last_interaction_before_partner', 'first_interaction_after_start',
                     'first_interaction_after_partner', 'crop_before', 'crop_after', 'feeding_size_intensity'])

    def get_feeding_sizes_intensity(self):
        for food_source in ['red', 'yellow']:
            clean_crop = self.clean_crop().loc[:, [food_source]]
            if len(self.feedings_dict[food_source]) > 0:
                self.feedings_dict[food_source]['crop_before'] = self.feedings_dict[food_source].apply(
                    lambda row: ForagerData.get_crop_before(row, clean_crop=clean_crop, color=food_source), axis=1)
                self.feedings_dict[food_source]['crop_after'] = self.feedings_dict[food_source].apply(
                    lambda row: ForagerData.get_crop_after(row, clean_crop=clean_crop), axis=1)
                self.feedings_dict[food_source]['feeding_size_intensity'] = \
                    self.feedings_dict[food_source]['crop_after'] - self.feedings_dict[food_source]['crop_before']
            else:
                self.feedings_dict[food_source] = pd.DataFrame(columns=
                    ['ant_id', 'feeding_start', 'feeding_end', 'food_source', 'last_interaction_before_end',
                     'last_interaction_before_partner', 'first_interaction_after_start',
                     'first_interaction_after_partner', 'crop_before', 'crop_after', 'feeding_size_intensity'])

    @staticmethod
    def get_crop_before(feed_row, clean_crop, color):
        if feed_row['feeding_start'] == 0:
            crop_before = pd.Series({color: 0})
        else:
            crop_before = clean_crop.iloc[feed_row['feeding_start']-1]
        if type(crop_before) is np.ndarray:
            crop_before = crop_before[0]
        return crop_before

    @staticmethod
    def get_crop_after(feed_row, clean_crop):
        if feed_row['feeding_end'] >= max(clean_crop.index):
            crop_after = clean_crop.iloc[-1]
        else:
            crop_after = clean_crop.iloc[feed_row['feeding_end']+1]
        return crop_after

    @staticmethod
    def get_crop_before_row(feed_row, crop, ymax):
        crop_before = crop.loc[feed_row['last_interaction_before_end']:(feed_row['feeding_start']-1)]
        crop_before95 = np.nanpercentile(crop_before[crop_before <= ymax], 95)
        return crop_before95

    @staticmethod
    def get_crop_after_row(feed_row, crop, ymax):
        crop_after = crop.loc[feed_row['feeding_end']:feed_row['first_interaction_after_start']]
        crop_after95 = np.nanpercentile(crop_after[crop_after <= ymax], 95)
        return crop_after95

    def get_feeding_sizes_ul(self, feeding_df, conversion_factors=None):
        pass

    def plot_raw_timeline(self):
        fig = super().plot_raw_timeline()
        for start_frame, end_frame in self.feedings_dict['yellow'][['feeding_start', 'feeding_end']].itertuples(index=False):
            plt.axvspan(start_frame, end_frame, facecolor='y', alpha=0.3)
        for start_frame, end_frame in self.feedings_dict['red'][['feeding_start', 'feeding_end']].itertuples(index=False):
            plt.axvspan(start_frame, end_frame, facecolor='r', alpha=0.3)
        plt.title(f'Forager {self.ant_id}')
        return fig

    def plot_clean_timeline(self, times_s, conversion_factors, include_raw=True, show=False):

        clean_crop_df = self.clean_crop()

        raw_crop_ul = {}
        clean_crop_ul = pd.DataFrame()
        for c in ['yellow', 'red']:
            raw_crop_ul[c] = self.crop_dict_raw[c] / conversion_factors.loc[0, c]
            clean_crop_ul.loc[:, c] = clean_crop_df[c] / conversion_factors.loc[0, c]

        fig = plt.figure(figsize=[18, 4])
        plt.plot(raw_crop_ul['yellow'].index*2.8/60, raw_crop_ul['yellow'], '.y', alpha=0.4)
        plt.plot(raw_crop_ul['red'].index*2.8/60, raw_crop_ul['red'], '.r', alpha=0.4)
        plt.ylim([0, max(self.ymax['yellow'] / conversion_factors.loc[0, 'yellow'], self.ymax['red'] / conversion_factors.loc[0, 'red'])])
        plt.xlabel('Time [min]')
        plt.ylabel(r'Crop load [$\mu$l]')
        plt.title(f'Ant {self.ant_id}')

        for start_frame, end_frame in self.interactions_df[['general_start_frame', 'general_end_frame']].itertuples(
                index=False):
            plt.axvspan(start_frame*2.8/60, end_frame*2.8/60, facecolor='0.5', alpha=0.2, zorder=-100)
        for start_frame, end_frame in self.feedings_dict['yellow'][['feeding_start', 'feeding_end']].itertuples(index=False):
            plt.axvspan(start_frame*2.8/60, end_frame*2.8/60, facecolor='y', alpha=0.3)
        for start_frame, end_frame in self.feedings_dict['red'][['feeding_start', 'feeding_end']].itertuples(index=False):
            plt.axvspan(start_frame*2.8/60, end_frame*2.8/60, facecolor='r', alpha=0.3)

        clean_crop_ul.set_index(clean_crop_ul.index*2.8/60, inplace=True)
        clean_crop_ul.plot(ax=fig.axes[0], color=[(184 / 255, 134 / 255, 11 / 255), (139 / 255, 0, 0), ], legend=None)

        if show:
            plt.show()
        return fig

    def plot_timeline_around_feedings(self):
        self.plot_raw_timeline()
        colors = {'red': '#8B0000', 'yellow': '#B8860B'}
        for food_source, feeding_table in zip(self.feedings_dict.keys(), self.feedings_dict.values()):
            x = feeding_table[['last_interaction_before_end', 'feeding_start', 'feeding_end',
                               'first_interaction_after_start']].to_numpy().flatten()
            y = feeding_table[['crop_before', 'crop_before', 'crop_after', 'crop_after']].to_numpy().flatten()
            plt.plot(x, y, colors[food_source])
            interactions = feeding_table[['last_interaction_before_end', 'first_interaction_after_start']].to_numpy().flatten()
            plt.vlines(x=interactions, ymin=0, ymax=max(self.ymax['red'], self.ymax['yellow']),
                       linestyles='dotted', colors='k')

    def manual_correction(self, experiment_path):
        corrections = pd.read_excel(experiment_path+sep+'manual_fixes_forager_timelines.xlsx')
        forager_corrections = corrections[corrections['ant_id'] == self.ant_id]
        for food_source, crop_data in self.crop_dict_raw.items():
            corrected_crop = copy(crop_data)
            frames_to_replace = forager_corrections['frame_to_replace'][forager_corrections['color'] == food_source]
            frames_of_input_values = forager_corrections['frame_to_take_value_from'][forager_corrections['color'] == food_source]
            corrected_crop[frames_to_replace] = crop_data[frames_of_input_values]
            self.crop_dict_raw[food_source] = corrected_crop


class InteractionData:
    def __init__(self, interactions_df_row, bdata, clean_crops, final_frame):
        self.ants = [str(int(interactions_df_row.actual_ant1)), str(int(interactions_df_row.actual_ant2))]
        self.start_frame, self.end_frame = interactions_df_row.general_start_frame, \
                                           min(interactions_df_row.general_end_frame, final_frame-1)
        # TODO: self.group_id = self.get_group_id(interactions_df_row)
        self.ant1_x, self.ant1_y = self.get_ant_location(self.ants[0], bdata)
        self.ant2_x, self.ant2_y = self.get_ant_location(self.ants[1], bdata)
        self.ant1_crop_before, self.ant1_got = self.get_ant_measurement(self.ants[0], clean_crops)
        self.ant2_crop_before, self.ant2_got = self.get_ant_measurement(self.ants[1], clean_crops)
        self.x = np.mean([self.ant1_x, self.ant2_x])
        self.y = np.mean([self.ant1_y, self.ant2_y])
        # TODO: self.ant1_confidence = self.rate_ant_confidence(self.ants[0], clean_crops)
        # TODO: self.ant2_confidence = self.rate_ant_confidence(self.ants[1], clean_crops)
        # TODO: self.size, self.giver, self.receiver = self.get_interaction_volume_and_direction()

    def get_ant_measurement(self, ant, clean_crops):
        if ant == '-1':
            return pd.Series({'red': np.nan, 'yellow': np.nan}), pd.Series({'red': np.nan, 'yellow': np.nan})
        ant_crop_before = clean_crops.loc[self.start_frame-1, ant]
        ant_crop_after = clean_crops.loc[self.end_frame+1, ant]
        ant_got = ant_crop_after - ant_crop_before
        return ant_crop_before, ant_got

    # TODO rate confidence
    def rate_ant_confidence(self, ant_id, bdata, clean_crops):
        ant = AntData(ant_id, bdata)
        clean_crop = clean_crops[ant]

        # find intervals before and after based on clean crop
        raw_crops = ant.crop_dict_raw
        pass

    def get_ant_location(self, ant, bdata):
        if ant == '-1':  # if ant_id is unknown
            return np.nan, np.nan
        x_data = bdata.loc[(self.start_frame - 1):(self.end_frame + 1), 'a' + ant + '-x']
        y_data = bdata.loc[(self.start_frame - 1):(self.end_frame + 1), 'a' + ant + '-x']
        x_detections = [x for x in x_data if x != -1]  # take values where ant was actually detected
        y_detections = [y for y in y_data if y != -1]
        if not x_detections:  # if ant was not detected
            return np.nan, np.nan
        ant_x = np.mean(x_detections)
        ant_y = np.mean(y_detections)
        return ant_x, ant_y

    # TODO: get unique group id
    def get_group_id(self, interactions_df_row):
        pass

    # TODO: get interaction volume
    def get_interaction_volume_and_direction(self):
        pass


class ExperimentData:
    def __init__(self, exp_num, root_path=r'Y:\Lior&Einav\Experiments', condition='with food', bdata_path='blob analysis'):
        # general experiment information
        self.exp_num = exp_num
        self.colony = None
        self.larvae = None
        self.food_ratios_dict = None
        self.food_concentration = None

        # paths
        self.exp_path = self.get_exp_path(root_path)
        self.condition = condition
        self.bdata_path = bdata_path

        # bdata - crops, locations and angles
        self.bdata_filename = self.get_bdata_filename()
        self.bdata = pd.read_csv(self.exp_path+sep+self.condition+sep+self.bdata_path + sep + self.bdata_filename)  # , index_col='frame')

        self.start_frame = self.bdata.frame[0]  # index[0]

        # interactions data
        self.interactions_df = self.load_or_create(r'trophallaxis_table.csv', self.get_interactions_df, write_data=True)

        # ants data
        self.num_ants, self.ants_list, self.foragers_list = self.get_ants_info(write_data=True)

        # weights data
        self.weights_dict = self.get_consumed_weights()  # how much was consumed from each food type by weight (grams)

        self.final_intake_by_weights = None

        # foragers feeding data
        self.feedings_df = self.load_or_create('forager_table_with_feeding_sizes_ul.csv',
                                               self.convert_feeding_sizes_to_ul,
                                               write_data=True)
        self.conversion_factors_by_weights_df = self.load_or_create('conversion_factors_by_weight_and_feeding_sum.csv',
                                                                    self.get_conversion_factors_by_weights,
                                                                    write_data=True)

        # clean crop data
        self.clean_crops = self.load_or_create('clean_crops_initial.csv', self.clean_bdata_initial, write_data=True,
                                               header=[0, 1])

    def get_ants_info(self, write_data=False):
        ants_info_df = self.load_or_create('ants_list.csv', self.create_ant_info_file,write_data=write_data)
        num_ants = len(ants_info_df)
        ants_list = ants_info_df.ant_id.tolist()
        foragers_list = ants_info_df.ant_id.loc[ants_info_df['is_forager']].tolist()
        return num_ants, ants_list, foragers_list

    def create_ant_info_file(self, min_detections=100, certainty_thres=0, write_data=False):
        # get recognized ants
        errors = self.bdata.filter(regex='error$', axis=1)
        all_ants = [int(x[1:-6]) for x in errors.columns]
        high_certainty = errors.where(lambda x: x <= certainty_thres)
        recognized = high_certainty.count() >= min_detections
        recognized_ants = np.array(all_ants)[recognized.to_numpy()]

        # find the foragers
        raw_feedings_file = pd.read_excel(self.exp_path+sep+"forager_feeding_table.xlsx")
        foragers_list = raw_feedings_file['ant_id'].unique()
        is_forager = np.isin(recognized_ants, foragers_list)

        # create the file
        ants_info = pd.DataFrame({'ant_id': recognized_ants, 'is_forager': is_forager})
        if write_data:
            ants_info.to_csv(self.exp_path + sep + r'ants_list.csv', index=False)

        return ants_info

    def visualize_ant_recognitions(self):
        errors = self.bdata.filter(regex='error$', axis=1)
        mean_error = errors.mean()
        n_detections = errors.count()
        plt.scatter(n_detections, mean_error)
        plt.xlabel('number of detections')
        plt.ylabel('mean detection error')
        plt.show()

    def get_foragers_list(self, write_data=False):
        foragers_list = self.feedings_df['ant_id'].unique().tolist()
        if write_data:
            np.savetxt(self.exp_path + sep + r'foragers_list.csv', foragers_list,delimiter=',')
        return foragers_list

    def get_exp_path(self, root_path):
        folderlist = os.listdir(root_path)
        exp_folder = [x for x in folderlist if x.startswith('experiment'+str(self.exp_num))]
        return root_path + sep + exp_folder[0]

    def get_experiment_info(self):
        pass

    def get_bdata_filename(self):
        filenames = os.listdir(self.exp_path+sep+self.condition+sep+self.bdata_path)
        bdata_filename = [x for x in filenames if x.startswith('bdata')]
        return bdata_filename[0]

    def load_or_create(self, filename, creator_method, write_data=False, load_as='df', header=0):
        # check if file exists
        # if exists load it
        if os.path.isfile(self.exp_path + sep + filename):
            if filename.endswith('csv'):
                if load_as == 'df':
                    df = pd.read_csv(self.exp_path + sep + filename, header=header)
                elif load_as == 'list':
                    with open(self.exp_path + sep + filename) as f:
                        reader = csv.reader(f)
                        df = list(reader)
                elif load_as == 'flat_list':
                    with open(self.exp_path + sep + filename) as f:
                        reader = csv.reader(f)
                        df = [int(row[0]) for row in reader]
            elif filename.endswith('xls') or filename.endswith('xlsx'):
                df = pd.read_excel(self.exp_path + sep + filename)
            else:
                df = None
                print('unrecognized file extension. Currently supports csv, xls, xlsx')
        # else create it using method
        else:
            df = creator_method(write_data=write_data)

        # convert trophallaxis table frame indexes
        if filename == 'trophallaxis_table.csv':
            df[['general_start_frame','general_end_frame']] = \
                df[['general_start_frame','general_end_frame']].apply(lambda x: x-self.start_frame)
            df = df.query('general_start_frame >= 0')
        return df

    def get_interactions_df(self, write_data=False):
        event_tables_path = self.exp_path+sep+self.condition+sep+'trophallaxis detection'+sep+'analyzed event tables'
        filelist = [f for f in os.listdir(event_tables_path) if f.endswith('xlsx')]
        df_list = []
        for filename in filelist:
            if not filename.startswith('~$'):
                p = pd.read_excel(event_tables_path+sep+filename)
                p = p.dropna(subset=['actual_ant1'])
                p = p.query("actual_ant1 != 'X' & actual_ant1 != 'x'")
                split1 = filename.split('_')
                split2 = split1[-1].split('.')
                vidnum = float(split2[0])
                p = p.assign(vidnum=vidnum)
                df_list.append(p[['vidnum', 'id', 'actual_ant1', 'actual_ant2', 'actual_start', 'actual_end']])
        all_interactions_df = pd.concat(df_list, ignore_index=True)
        all_interactions_df = all_interactions_df.assign(
            general_start_frame=lambda x: get_general_frame(x.actual_start, x.vidnum),
            general_end_frame=lambda x: get_general_frame(x.actual_end, x.vidnum))
        all_interactions_df.sort_values(by='general_start_frame', inplace=True)
        all_interactions_df.reset_index(inplace=True, drop=True)

        if write_data:
            all_interactions_df.to_csv(self.exp_path+sep+r'trophallaxis_table.csv', index=False)

        # all_interactions_df[['general_start_frame','general_end_frame']] = \
        #     all_interactions_df[['general_start_frame','general_end_frame']].apply(lambda x: x-self.start_frame)
        # all_interactions_df = all_interactions_df.query('general_start_frame >= 0')

        return all_interactions_df

    def enrich_interactions_df(self, write_data=False, filename='trophallaxis_table_enriched'):
        ant1_got_dict = {}
        ant2_got_dict = {}
        ant1_x_dict = {}
        ant1_y_dict = {}
        ant2_x_dict = {}
        ant2_y_dict = {}
        ant1_crop_before_dict = {}
        ant2_crop_before_dict = {}
        final_frame = max(self.clean_crops.index)
        for idx, trop_row in self.interactions_df.iterrows():
            trop = InteractionData(trop_row, self.bdata, self.clean_crops, final_frame)
            ant1_got_dict[idx] = trop.ant1_got
            ant2_got_dict[idx] = trop.ant2_got
            ant1_crop_before_dict[idx] = trop.ant1_crop_before
            ant2_crop_before_dict[idx] = trop.ant2_crop_before

            ant1_x_dict[idx] = trop.ant1_x
            ant2_x_dict[idx] = trop.ant2_x
            ant1_y_dict[idx] = trop.ant1_y
            ant2_y_dict[idx] = trop.ant2_y

        ant1_got_df = pd.concat(ant1_got_dict, axis=1).T
        ant2_got_df = pd.concat(ant2_got_dict, axis=1).T
        ant1_got_df.rename(columns={'red': 'ant1_got_red', 'yellow': 'ant1_got_yellow'}, inplace=True)
        ant2_got_df.rename(columns={'red': 'ant2_got_red', 'yellow': 'ant2_got_yellow'}, inplace=True)

        ant1_crop_before_df = pd.concat(ant1_crop_before_dict, axis=1).T
        ant2_crop_before_df = pd.concat(ant2_crop_before_dict, axis=1).T
        ant1_crop_before_df.rename(columns={'red': 'ant1_crop_before_red', 'yellow': 'ant1_crop_before_yellow'}, inplace=True)
        ant2_crop_before_df.rename(columns={'red': 'ant2_crop_before_red', 'yellow': 'ant2_crop_before_yellow'}, inplace=True)

        loc_df = pd.DataFrame({'ant1_x': ant1_x_dict, 'ant1_y': ant1_y_dict, 'ant2_x': ant2_x_dict, ' ant2_y': ant2_y_dict})

        enriched_interactions_df = self.interactions_df.join([ant1_got_df, ant2_got_df, ant1_crop_before_df, ant2_crop_before_df, loc_df])

        if write_data:
            enriched_interactions_df.to_csv(self.exp_path + sep + filename + '.csv', index=False)

        return enriched_interactions_df

    def enrich_interactions_df_size_and_loc_old(self, write_data=False):
        ant1_got_dict = {}
        ant2_got_dict = {}
        ant1_x_dict = {}
        ant1_y_dict = {}
        ant2_x_dict = {}
        ant2_y_dict = {}
        ant1_crop_before_dict = {}
        ant2_crop_before_dict = {}
        final_frame = max(self.clean_crops.index)
        for idx, trop in self.interactions_df.iterrows():
            ant1, ant2 = str(int(trop.actual_ant1)), str(int(trop.actual_ant2))
            start_frame, end_frame = trop.general_start_frame, trop.general_end_frame
            if end_frame > final_frame:
                end_frame = final_frame-1
            ant1_crop_before = self.clean_crops.loc[start_frame-1, ant1]
            ant1_crop_after = self.clean_crops.loc[end_frame+1, ant1]
            ant1_got = ant1_crop_after - ant1_crop_before
            ant1_x = self.bdata.loc[(start_frame - 1):(end_frame + 1), 'a' + ant1 + '-x'].mean()
            ant1_y = self.bdata.loc[(start_frame - 1):(end_frame + 1), 'a' + ant1 + '-y'].mean()
            if ant2 == '-1':
                ant2_got = pd.Series({'red': np.nan, 'yellow': np.nan})
                ant2_x = np.nan
                ant2_y = np.nan
                ant2_crop_before = pd.Series({'red': np.nan, 'yellow': np.nan})
            else:
                ant2_crop_before = self.clean_crops.loc[start_frame - 1, ant2]
                ant2_crop_after = self.clean_crops.loc[end_frame + 1, ant2]
                ant2_got = ant2_crop_after - ant2_crop_before
                ant2_x = self.bdata.loc[(start_frame - 1):(end_frame + 1), 'a' + ant2 + '-x'].mean()
                ant2_y = self.bdata.loc[(start_frame - 1):(end_frame + 1), 'a' + ant2 + '-y'].mean()
            ant1_got_dict[idx] = ant1_got
            ant2_got_dict[idx] = ant2_got
            ant1_crop_before_dict[idx] = ant1_crop_before
            ant2_crop_before_dict[idx] = ant2_crop_before

            ant1_x_dict[idx] = ant1_x
            ant2_x_dict[idx] = ant2_x
            ant1_y_dict[idx] = ant1_y
            ant2_y_dict[idx] = ant2_y

        ant1_got_df = pd.concat(ant1_got_dict, axis=1).T
        ant2_got_df = pd.concat(ant2_got_dict, axis=1).T
        ant1_got_df.rename(columns={'red': 'ant1_got_red', 'yellow': 'ant1_got_yellow'}, inplace=True)
        ant2_got_df.rename(columns={'red': 'ant2_got_red', 'yellow': 'ant2_got_yellow'}, inplace=True)

        ant1_crop_before_df = pd.concat(ant1_crop_before_dict, axis=1).T
        ant2_crop_before_df = pd.concat(ant2_crop_before_dict, axis=1).T
        ant1_crop_before_df.rename(columns={'red': 'ant1_crop_before_red', 'yellow': 'ant1_crop_before_yellow'}, inplace=True)
        ant2_crop_before_df.rename(columns={'red': 'ant2_crop_before_red', 'yellow': 'ant2_crop_before_yellow'}, inplace=True)

        loc_df = pd.DataFrame({'ant1_x': ant1_x_dict, 'ant1_y': ant1_y_dict, 'ant2_x': ant2_x_dict, ' ant2_y': ant2_y_dict})

        enriched_interactions_df = self.interactions_df.join([ant1_got_df, ant2_got_df, ant1_crop_before_df, ant2_crop_before_df, loc_df])
        return enriched_interactions_df

    def clean_bdata_initial(self, percentile=90, write_data=False):
        clean_crops_dict = {}
        for ant_id in self.ants_list:
            if ant_id in self.foragers_list:
                ant = ForagerData(ant_id, feedings_df=self.feedings_df, bdata_df=self.bdata,
                                  interactions_df=self.interactions_df)
            else:
                ant = AntData(ant_id, bdata_df=self.bdata, interactions_df=self.interactions_df)
            clean_crop_ant = ant.clean_crop()
            clean_crops_dict[ant_id] = clean_crop_ant
        clean_crops_df = pd.concat(clean_crops_dict, axis=1)

        if write_data:
            clean_crops_df.to_csv(self.exp_path + sep + r'clean_crops_initial.csv', index=False)
        return clean_crops_df

    def get_consumed_weights(self):
        weights_df = pd.read_excel(self.exp_path+sep+'weights.xlsx')
        weights_df['difference'] = weights_df['before_g'] - weights_df['after_g']
        control_rows = weights_df['type'] == 'control'
        treatment_rows = ~control_rows
        evap = np.mean(weights_df['difference'][control_rows])
        consumed_g = {}
        for c in ['yellow', 'red']:
            consumed_g[c] = weights_df.difference[treatment_rows & (weights_df['color'] == c)].to_numpy()[0] - evap
        return consumed_g

    def get_feeding_sizes_intensity(self, plot_timelines=False, write_data=False):
        # try:
        #     fdata = pd.read_excel(self.exp_path+sep+"forager_feeding_table_with_interaction_data.xlsx", usecols="A:D,F:I")
        # except FileNotFoundError:
        #     print('feeding_sizes_intensity is None (no file forager_feeding_table_with_interaction_data.xlsx)')
        #     return None
        bdata = self.bdata
        bdata = bdata.drop(bdata.columns[0], axis=1)

        tdata = self.interactions_df
        # tdata[['general_start_frame','general_end_frame']] = tdata[['general_start_frame','general_end_frame']].apply(lambda x: x-self.start_frame)
        # tdata = tdata.query('general_start_frame >= 0')

        fdata = pd.read_excel(self.exp_path + sep + r"forager_feeding_table.xlsx")

        framelist = []
        for ant_id in fdata.ant_id.unique():
            if ant_id == 535:
                print('535')
            ant = ForagerData(ant_id=ant_id, bdata_df=bdata, interactions_df=tdata, feedings_df=fdata)
            ant.manual_correction(self.exp_path)
            ant.get_feeding_sizes_intensity()
            framelist.extend([ant.feedings_dict['red'], ant.feedings_dict['yellow']])
            if plot_timelines:
                ant.plot_timeline_around_feedings()
        fdata_with_feed_sizes_intensity = pd.concat(framelist)
        fdata_with_feed_sizes_intensity.sort_index(inplace=True)
        fdata_with_feed_sizes_intensity['feeding_size_intensity'][
            fdata_with_feed_sizes_intensity['feeding_size_intensity'] < 0] = 0

        # fdata_with_feed_sizes_intensity.loc[:, ['feeding_start', 'feeding_end', 'last_interaction_before_end', 'first_interaction_after_start']] += self.start_frame

        if write_data:
            fdata_with_feed_sizes_intensity.to_csv(
                self.exp_path+sep+r'forager_table_with_feeding_sizes.csv', index=False)

        return fdata_with_feed_sizes_intensity

    def get_conversion_factors_by_weights(self, write_data=False):
        fdata = self.load_or_create('forager_table_with_feeding_sizes.csv', self.get_feeding_sizes_intensity, write_data=True)
        if fdata is None:
            print('conversion_factors_by_weights is None (no file forager_table_with_feeding_sizes.csv)')
            return None
        consumed_intensities = fdata['feeding_size_intensity'].groupby(fdata['food_source']).sum()
        conversion_factors = {}
        for c in ['red', 'yellow']:
            ul = self.weights_dict[c] * 1000
            conversion_factors[c] = consumed_intensities[c] / ul

        conversion_factors = pd.DataFrame(conversion_factors, index=[0])
        if write_data:
            conversion_factors.to_csv(
                self.exp_path+sep+r'conversion_factors_by_weight_and_feeding_sum.csv', index=False)

        return conversion_factors

    def convert_feeding_sizes_to_ul(self, write_data=False):
        fdata_with_feed_sizes_ul = self.load_or_create('forager_table_with_feeding_sizes.csv', self.get_feeding_sizes_intensity, write_data=True)
        if fdata_with_feed_sizes_ul is None:
            print('fdata_with_feed_sizes_ul is None (no file forager_table_with_feeding_sizes.csv)')
            return None
        self.conversion_factors_by_weights_df = self.load_or_create(filename='conversion_factors_by_weight_and_feeding_sum.csv',
                                                                    creator_method=self.get_conversion_factors_by_weights,
                                                                    write_data=True)
        fdata_with_feed_sizes_ul['feeding_size_ul'] = fdata_with_feed_sizes_ul.apply(
            lambda x: x['feeding_size_intensity'] / self.conversion_factors_by_weights_df[x['food_source']], axis=1)

        if write_data:
            fdata_with_feed_sizes_ul.to_csv(self.exp_path+sep+r'forager_table_with_feeding_sizes_ul.csv',index=False)

        return fdata_with_feed_sizes_ul

    def get_final_intake_by_weights(self):
        pass

    def plot_intake_trajectory(self):
        pass

    def plot_intake_trajectory_by_foragers(self):
        pass

    def crops_animation(self):
        pass
