"""
File:         filter_and_recode.py
Created:      06/05/2021
Last Changed: 07/06/2021
Author(s):    H.H.Wiersma

Code to create distribution plot and a mask image used in the research plan
Copyright (C) 2019 H.H.Wiersma

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

import sys
import os
import glob
import pandas as pd
import numpy as np
from datetime import datetime
import json


def main():
    input_data_path = sys.argv[1] # Qestuinnaire data input file path
    question_id_path = sys.argv[2] # question id file path
    recode_model_path = sys.argv[3] # recode information file path
    participant_selection = sys.argv[4] # participant selection filter file path
    participant__longitudinal_selection = sys.argv[5] # participant selection longitidional selection
    output_path = sys.argv[6] # output file path
    question_dir_path = sys.argv[7] # directory path with the question information
    ennumeration_dir_path = sys.argv[8] # directory path with the enummeration information files?
    question_file_paths = glob.glob(os.path.join(question_dir_path, "*.csv"))
    total_dfs = []

    # read questionnaire information files
    for quest_file_path in question_file_paths:
        df = pd.read_csv(quest_file_path, sep=",")
        total_dfs.append(df)
    df_variables = pd.concat(total_dfs, axis=0)
    df_variables = df_variables.reset_index(drop=True)
    df_variables = df_variables.loc[df_variables["VARIABLE_NAME"].str.startswith("covt"), :]

    # read enumeration information files
    enum_file_paths = glob.glob(os.path.join(ennumeration_dir_path, "*.csv"))
    total_enum_dfs = []
    for enum_file_path in enum_file_paths:
        df = pd.read_csv(enum_file_path, sep=",")
        total_enum_dfs.append(df)
    df_enumerations = pd.concat(total_enum_dfs, axis=0)
    df_enumerations = df_enumerations.reset_index(drop=True)

    # read other input files
    df_data = pd.read_pickle(input_data_path)
    df_question_ids = pd.read_csv(question_id_path, sep="\t", index_col=0)
    df_recode_info = pd.read_csv(recode_model_path, sep="\t", index_col=0)
    df_participant_selection = pd.read_csv(participant_selection, sep="\t", index_col=0)
    df_participant_longitudinal_selection = pd.read_csv(participant__longitudinal_selection, sep="\t")

    # create output dir
    create_dir(output_path)
    ###
    ### perform filtering
    ###
    df_data_filtered = df_data.loc[df_data.index.intersection(df_participant_selection.index), :]
    col_trans = {
        "X1.0" :"covt01",
        "X2.0" :"covt02",
        "X3.0" :"covt03",
        "X4.0" :"covt04",
        "X5.0" :"covt05",
        "X6.0" :"covt06",
        "X7.0" :"covt07",
        "X8.0" :"covt08",
        "X9.0" :"covt09",
        "X10.0" :"covt10",
        "X11.0" :"covt11",
        "X12.0" :"covt12",
        "X13.0" :"covt13",
        "X14.0" :"covt14",
        "X15.0" :"covt15",
        "X15.5" :"covt15b",
        "X16.0" :"covt16",
        "X16.5" :"covt16b",
        "X17.0" :"covt17",
    }
    df_participant_selection = df_participant_selection.rename(columns=col_trans)
    for col_index, column in df_participant_selection.iteritems():
        col_selection = df_data.columns[df_data.columns.str.startswith(col_index)]
        df_data_filtered.loc[column.index, col_selection] = df_data_filtered.loc[column.index, col_selection].where(column, np.nan)
    export_df(df_data_filtered, output_path, "participants_filtered")
    print("Filtering ready", df_data_filtered.shape)
    ###
    ### // perform filtering
    ###


    ###
    ### recode answers
    ###
    question_selection = df_recode_info.index.intersection(df_question_ids.index)
    selected_cols = ['1.0', '2.0', '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0', '10.0', '11.0', '12.0', '13.0', '14.0', '15.0', '15.5', '16.0', '16.5', '17.0']
    df_data_filtered_recoded = df_data_filtered.copy()
    total_df_info_per_question = []
    #loop through all the questions
    for question_id in question_selection:

        #get the recoding information
        replace_data = df_recode_info.loc[question_id, "Recode answers"]
        model_type = df_recode_info.loc[question_id, "Type"]

        quest_ids_for_replace = df_question_ids.loc[question_id, selected_cols]
        quest_ids_for_replace = quest_ids_for_replace.dropna()

        question_first_time_asked = True
        # loop through de single ids of a particular question
        for single_quest_id in quest_ids_for_replace:
            if single_quest_id in df_data_filtered.columns:
                label_en = ""
                if single_quest_id in df_variables.columns:
                    label_en = df_variables.loc[single_quest_id, "DEFINITION_EN"]

                # Get the the information for the recoding overview sheet
                question_enum_nl = get_enum_of_question(df_enumerations, single_quest_id, en=False)
                question_enum_nl_str = create_json_str(question_enum_nl)

                question_enum_en = get_enum_of_question(df_enumerations, single_quest_id, en=True)
                question_enum_en_str = create_json_str(question_enum_en)

                df_single_question = df_data_filtered.loc[:, single_quest_id]
                df_single_question = df_single_question.astype(float)
                value_counts_before = df_single_question.value_counts(dropna=False)
                if len(value_counts_before) > 20:
                    value_counts_before = pd.Series()

                value_counts_before_str = create_json_str(value_counts_before.to_dict())
                fin_recode_info = create_recode_dict(df_single_question, model_type, replace_data)
                if len(fin_recode_info) > 0:
                    #recode the data
                    recoded_values = df_single_question.replace(fin_recode_info)
                    df_data_filtered_recoded.loc[:, single_quest_id] = recoded_values

                    # calculate value counts after recode
                    value_counts_after = df_data_filtered_recoded.loc[:, single_quest_id].value_counts(dropna=False)
                    if len(value_counts_after) > 20:
                        value_counts_after = pd.Series()

                    # get enum info
                    question_enum_nl_recoded = recode_enum(question_enum_nl, fin_recode_info)
                    question_enum_nl_recoded_str = create_json_str(question_enum_nl_recoded)

                    question_enum_en_recoded = recode_enum(question_enum_en, fin_recode_info)
                    question_enum_en_recoded_str = create_json_str(question_enum_en_recoded)
                    fin_recode_info_str = create_json_str(fin_recode_info)
                    value_counts_after_str = create_json_str(value_counts_after.to_dict())

                    # save all the information for the recoding file
                    dict_data = {
                        "question_id": single_quest_id,
                        "is_first_time_asked": question_first_time_asked,
                        "question_label_nl": question_id,
                        "question_label_en": label_en,
                        "model_type": model_type,
                        "given_recode_json": replace_data,
                        "is_recoded": True,
                        "given_recode_human": human_redable_json(replace_data),
                        "total_recode_json": fin_recode_info_str,
                        "total_recode_human": human_redable_json(fin_recode_info_str),
                        "value_counts": value_counts_before_str,
                        "value_counts_human": human_redable_json(value_counts_before_str),
                        "value_counts_after_recode": value_counts_after_str,
                        "value_counts_after_recode_human": human_redable_json(value_counts_after_str),
                        "answer_options_nl_json": question_enum_nl_str,
                        "answer_options_nl_human": human_redable_json(question_enum_nl_str),
                        "answer_options_nl_after_recode_json": question_enum_nl_recoded_str,
                        "answer_options_nl_after_recode_human": human_redable_json(question_enum_nl_recoded_str),
                        "answer_options_en_json": question_enum_en_str,
                        "answer_options_en_human": human_redable_json(question_enum_en_str),
                        "answer_options_en_after_recode_json": question_enum_en_recoded_str,
                        "answer_options_en_after_recode_human": human_redable_json(question_enum_en_recoded_str),
                    }
                    total_df_info_per_question.append(dict_data)
                else:
                    dict_data = {
                        "question_id": single_quest_id,
                        "is_first_time_asked": question_first_time_asked,
                        "question_label_nl": question_id,
                        "question_label_en": label_en,
                        "model_type": model_type,
                        "given_recode_json": replace_data,
                        "is_recoded": False,
                        "given_recode_human": human_redable_json(replace_data),
                        "total_recode_json": "",
                        "total_recode_human": "",
                        "value_counts": value_counts_before_str,
                        "value_counts_human": human_redable_json(value_counts_before_str),
                        "value_counts_after_recode": value_counts_before_str,
                        "value_counts_after_recode_human": human_redable_json(value_counts_before_str),
                        "answer_options_nl_json": question_enum_nl_str,
                        "answer_options_nl_human": human_redable_json(question_enum_nl_str),
                        "answer_options_nl_after_recode_json": question_enum_nl_str,
                        "answer_options_nl_after_recode_human": human_redable_json(question_enum_nl_str),
                        "answer_options_en_json": question_enum_nl_str,
                        "answer_options_en_human": human_redable_json(question_enum_nl_str),
                        "answer_options_en_after_recode_json": question_enum_nl_str,
                        "answer_options_en_after_recode_human": human_redable_json(question_enum_nl_str),
                    }
                    total_df_info_per_question.append(dict_data)

                if question_first_time_asked:
                    question_first_time_asked = False

    # recoding is ready
    df_recode_overview = pd.DataFrame(total_df_info_per_question)

    # Create export with the recoded data
    export_df(df_data_filtered_recoded, output_path, "participants_filtered_recoded_answers")

    # create recoding overview document
    export_df(df_recode_overview, output_path, "recode_info_overview")

    print("recoding ready", df_data_filtered_recoded.shape)
    ###
    ### // recode answers
    ###


    ###
    ### Filter timeseries participants
    ###
    df_data_filtered_recoded = pd.read_pickle(sys.argv[9])
    print("df_data_filtered_recoded", df_data_filtered_recoded)

    df_data_filtered_recoded = df_data_filtered_recoded.loc[df_participant_longitudinal_selection.iloc[:, 0], :]
    export_df(df_data_filtered_recoded, output_path, "participants_filtered_recoded_answers_longitudinal_filtered")
    print("filtering longitudinal ready", df_data_filtered_recoded.shape)
    ###
    ### // Filter timeseries participants
    ###


def create_recode_dict(df, model_type, recode_info):
    final_recode_dict = {}
    if pd.isna(recode_info) == False:

        # process the recode information
        replace_info_from_json = json.loads(recode_info)
        for key, replace_key in replace_info_from_json.items():
            # translate "NA" string to an numpy nan
            if replace_key == "NA":
                final_recode_dict[float(key)] = np.nan
            else:
                final_recode_dict[float(key)] = float(replace_key)

    unique_values = df.dropna().unique().astype(float)
    if model_type == "ordinal-ordered-turned":
        # turn the answer direction
        if len(final_recode_dict) > 0:
            unique_values_temp = unique_values
            unique_values = []
            for unique_val in unique_values_temp:
                if unique_val not in final_recode_dict.keys():
                    unique_values.append(unique_val)
        key_range = list(range(int(np.min(unique_values)), int(np.max(unique_values)) + 1))
        recode_table_temp = {}
        for key, replace_key in dict(zip(key_range, key_range[::-1])).items():
            recode_table_temp[float(replace_key)] = float(key)
        print(recode_table_temp)
        final_recode_dict = {**recode_table_temp, **final_recode_dict}
    elif model_type == "binomial":
        # recode the 2 (no answer) to 0, if no other recoding is given
        if len(final_recode_dict) > 0:
            if 1.0 not in final_recode_dict.keys() and 2.0 not in final_recode_dict.keys():
                final_recode_dict[2.0] = 0.0
        else:
            if np.array_equal(np.sort(unique_values), np.sort(np.array([1.0, 2.0]))):
                final_recode_dict[2.0] = 0.0
    return final_recode_dict


def get_enum_of_question(df_enumerations, variable_id, en=False):
    label_id = "ENUMERATION_EN" if en else "ENUMERATION_NL"
    answer_options = df_enumerations.loc[df_enumerations["VARIABLE_NAME"] == variable_id, ["ENUMERATION_CODE", label_id]]
    answer_options = answer_options.set_index("ENUMERATION_CODE")
    answer_options.index = answer_options.index.astype(float)
    answer_options = answer_options.loc[:, label_id]
    return answer_options.to_dict()


def create_json_str(json_dict):
    if pd.isna(json_dict):
        return ""
    return json.dumps(json_dict)


def human_redable_json(json_str):
    if pd.isna(json_str):
        return ""
    return json_str.replace("{", "").replace("}", "").replace("\"", "")


def recode_enum(enum_info, recode_info):
    if len(recode_info) > 0:
        recoded_enum = {}
        for key, value in enum_info.items():
            if key in recode_info.keys():
                if pd.isna(recode_info[key]):
                    recoded_enum[key] = np.nan
                else:
                    recoded_enum[recode_info[key]] = value
            else:
                recoded_enum[key] = value
        return recoded_enum
    return enum_info


def export_df(df, output_path, suffix):
    name = "questionnaire_subset_{suffix}_{date}".format(
        suffix=suffix,
        date=datetime.now().strftime("%d-%m-%Y")
    )
    df.to_pickle(os.path.join(output_path, "{}.pkl".format(name)))
    df.to_csv(os.path.join(output_path, "{}.txt".format(name)), sep="\t")
    print("export files created: {}.pkl/txt".format(name))


def create_dir(dir_name):
    if os.path.isdir(dir_name) == False:
        os.mkdir(dir_name)


if __name__ == '__main__':
    main()




