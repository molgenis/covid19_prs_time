"""
File:         fit_pgs_score_to_phenotypes.py
Created:      03/03/2021
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
from datetime import datetime
import os
import pandas as pd
import numpy as np
import json
import statsmodels.api as sm
from statsmodels.miscmodels.ordinal_model import OrderedModel
from statsmodels.tools.sm_exceptions import ConvergenceWarning


def main():
    ###
    ### input paths
    ###
    input_df_path = sys.argv[1] # dataframe with the combined questionnaire data
    input_question_overview = sys.argv[2] # datafrane with the question ids over time
    input_pgs_path = sys.argv[3] # dataframe with the PGS values of the participants
    output_dir_ori_path = sys.argv[4] # output directory path
    suffix = sys.argv[5] # suffix fot the outputs
    pgs_id = int(sys.argv[7]) # int of pgs number, used for multi node analysis on cluster
    question_prs_selection_path = sys.argv[6] # input path of question and prs selection
    selected_trait_file_path = sys.argv[7] # input path of trait selection file
    model_selection_file = sys.argv[8] # input path of trait selection file
    df_question_prs_selection = pd.read_csv(question_prs_selection_path, sep="\t")

    # create output directory
    create_dir(output_dir_ori_path)

    # create log file
    log_file_path = "log_file_{id}_{suffix}_{date}.txt".format(
        id=pgs_id,
        suffix=suffix,
        date=datetime.now().strftime("%d-%m-%Y")
    )
    logfile = open(os.path.join(output_dir_ori_path, log_file_path), "w")

    # read input files
    df_quest = pd.read_pickle(input_df_path)
    df_question_ids_total = pd.read_csv(input_question_overview, sep="\t", index_col=0, dtype="str")
    df_question_ids = df_question_ids_total.loc[:, df_question_ids_total.columns.difference(["Number of timepoints", "Question answers"])]
    df_question_ids.columns = df_question_ids.columns.astype(float)
    df_question_ids = df_question_ids.T
    df_question_ids = df_question_ids.sort_index()
    trait_subset = pd.read_pickle(input_pgs_path)

    # Filter questionnaire question data on meta information
    start_columns = df_quest.columns[df_quest.columns.str.startswith("covt")]
    date_columns = df_quest.columns[df_quest.columns.str.endswith("DATE")]
    age_columns = df_quest.columns[df_quest.columns.str.endswith("AGE")]
    gender_columns = df_quest.columns[df_quest.columns.str.endswith("GENDER")]
    date_variant_id = df_quest.columns[df_quest.columns.str.endswith("VARIANT_ID")]
    date_zip_code = df_quest.columns[df_quest.columns.str.endswith("ZIP_CODE")]
    date_response_rate = df_quest.columns[df_quest.columns.str.contains("responsedate")]
    skip_columns = ["covt17_COVID172TXT", "covt17_COVID177TXT", "covt17_COVID192A", "covt17_COVID192B", "covt17_PSEUDOIDEXT"]
    selection_columns = start_columns.difference(date_columns).difference(age_columns).difference(gender_columns).difference(date_variant_id).difference(date_zip_code).difference(date_response_rate).difference(skip_columns)
    df_quest = df_quest.loc[df_quest.index.intersection(trait_subset.index), :]

    # trait selection
    trait_subset.columns = trait_subset.columns.str.replace("/", ".").str.replace(" ", ".").str.replace("-", ".").str.replace("(", ".").str.replace(")", ".")
    df_selected_traits = pd.read_csv(selected_trait_file_path, header=None)
    trait_subset = trait_subset.loc[:, trait_subset.columns.intersection(df_selected_traits.iloc[:, 0])]
    trait_subset = trait_subset.iloc[:, [pgs_id]]
    print("Process trade: {}".format(trait_subset.columns[0]))

    #  model covariate columns
    correction_columns = ["age_recent", "age2_recent", "chronic_recent", "household_recent", "have_childs_at_home_recent", "gender_recent"]
    correction_df = df_quest.loc[:, correction_columns]

    # question selection
    df_selected_questions = pd.read_csv(model_selection_file, sep="\t", index_col="Question")
    df_selected_questions = df_selected_questions.dropna(subset=["Type"])
    question_list = df_selected_questions.index.intersection(df_question_ids.columns)

    # create output dataframes
    multiIndex_columns = pd.MultiIndex.from_product([question_list, df_question_ids.index], names=["question", "quest_nr"])
    df_betas_per_question = pd.DataFrame(index=trait_subset.columns, columns=multiIndex_columns)
    df_pvalues_per_question = pd.DataFrame(index=trait_subset.columns, columns=multiIndex_columns)
    df_se_values_per_question = pd.DataFrame(index=trait_subset.columns, columns=multiIndex_columns)
    n_values_per_question = []
    value_counts_per_question = []

    # process all questions
    for column in question_list:
        question_ids_of_question = df_question_ids.loc[:, column]
        question_ids_of_question = question_ids_of_question[~question_ids_of_question.isna()]

        # create temp lists for model outputs
        betas_per_week = []
        pvalues_per_week = []
        se_values_per_week = []
        n_values_per_week = {}
        values_counts_per_week = {}
        # process all datapoints of the question
        for index, quest_id in question_ids_of_question.iteritems():
            if quest_id in selection_columns:

                # create model input data
                df_subset = df_quest.loc[:, quest_id]
                df_subset = df_subset.astype(float)
                df_subset = df_subset.dropna()

                # find the correct model
                model_type = df_selected_questions.loc[column, "Type"]
                if model_type == "ordinal":
                    df_subset = df_subset.sort_index()
                    df_subset = df_subset.astype(int).astype("category")
                    model_type = "ordinal"
                elif model_type == "ordinal-ordered" or model_type == "ordinal-ordered-turned":
                    df_subset = df_subset.sort_index()
                    df_subset = df_subset.astype(int).astype("category").cat.as_ordered()
                    model_type = "ordinal"


                # fit all models per trait
                beta_values_per_model = {}
                pvalues_per_model = {}
                se_values_per_model = {}
                for pgs_column_name in trait_subset.columns:
                    if ((df_question_prs_selection["Question"] == column) & (df_question_prs_selection["prs"] == pgs_column_name)).any():

                        #write logfile info
                        logfile.write("\n\nProcess: {question}, {question_id}, {pgs}, {model}\n".format(
                            question=column,
                            question_id=quest_id,
                            pgs=pgs_column_name,
                            model=model_type
                        ))
                        print("process", column, pgs_column_name)
                        try:
                            # create model and set the model parameters
                            df_model_input = pd.merge(trait_subset.loc[:, [pgs_column_name]], correction_df, left_index=True, right_index=True)
                            fit_args = {}
                            if model_type == "binomial":
                                df_model_input["intercept"] = 1.0
                                mod = sm.Logit(df_subset, df_model_input.loc[df_subset.index, :])
                                fit_args["maxiter"] = 10000
                            elif model_type == "ordinal":
                                mod = OrderedModel(df_subset, df_model_input.loc[df_subset.index, :], distr="logit")
                                fit_args["method"] = 'bfgs'
                                fit_args["maxiter"] = 10000
                            else:
                                df_model_input["intercept"] = 1.0
                                mod = sm.OLS(df_subset, df_model_input.loc[df_subset.index, :])

                            # fit the model
                            res = mod.fit(**fit_args)

                            #write model output
                            print(res.summary())
                            logfile.write(res.summary().as_text())
                            logfile.write("\n")

                            #save output from the model
                            beta_values_per_model[pgs_column_name] = res.params[pgs_column_name]
                            pvalues_per_model[pgs_column_name] = res.pvalues[pgs_column_name]
                            se_values_per_model[pgs_column_name] = res.bse[pgs_column_name]
                        except sm.tools.sm_exceptions.PerfectSeparationError:
                            print("Exception PerfectSeparationError", pgs_column_name, quest_id, column)
                            logfile.write("Error PerfectSeparationError: {question_id}, {pgs}\n".format(
                                question_id=quest_id,
                                pgs=pgs_column_name,
                            ))
                            continue
                        except np.linalg.LinAlgError:
                            print("Exception LinAlgError", pgs_column_name, quest_id, column)
                            logfile.write("Error LinAlgError: {question_id}, {pgs}\n".format(
                                question_id=quest_id,
                                pgs=pgs_column_name,
                            ))
                            continue
                        except UnboundLocalError:
                            print("Exception UnboundLocalError", pgs_column_name, quest_id, column)
                            logfile.write("Error UnboundLocalError: {question_id}, {pgs}\n".format(
                                question_id=quest_id,
                                pgs=pgs_column_name,
                            ))
                            continue
                        except Exception:
                            print("Exception", pgs_column_name, quest_id, column)
                            logfile.write("Error Exception (general): {question_id}, {pgs}\n".format(
                                question_id=quest_id,
                                pgs=pgs_column_name,
                            ))
                            continue

                # save all the model output per time point
                model_betas = pd.Series(beta_values_per_model, name=index)
                betas_per_week.append(model_betas)

                model_pvalues = pd.Series(pvalues_per_model, name=index)
                pvalues_per_week.append(model_pvalues)

                model_se_values = pd.Series(se_values_per_model, name=index)
                se_values_per_week.append(model_se_values)

                n_values_per_week[index] = df_subset.shape[0]

                val_counts = df_subset.value_counts()
                val_counts.index = val_counts.index.astype(int).astype(str)
                values_counts_per_week[index] = json.dumps(val_counts.to_dict())
        # save all the model information per PRS, per time point, per question (multi column table)
        if len(betas_per_week) > 0:
            df_model_betas = pd.concat(betas_per_week, axis=1)
            df_model_betas.columns = pd.MultiIndex.from_product([[column], df_model_betas.columns])
            df_betas_per_question.loc[df_model_betas.index, df_model_betas.columns] = df_model_betas

            df_model_pvalues = pd.concat(pvalues_per_week, axis=1)
            df_model_pvalues.columns = pd.MultiIndex.from_product([[column], df_model_pvalues.columns])
            df_pvalues_per_question.loc[df_model_pvalues.index, df_model_pvalues.columns] = df_model_pvalues

            df_model_se_values = pd.concat(se_values_per_week, axis=1)
            df_model_se_values.columns = pd.MultiIndex.from_product([[column], df_model_se_values.columns])
            df_se_values_per_question.loc[df_model_se_values.index, df_model_se_values.columns] = df_model_se_values

            n_values_per_question.append(pd.Series(n_values_per_week, name=column))
            value_counts_per_question.append(pd.Series(values_counts_per_week, name=column))

    df_nvalues = pd.concat(n_values_per_question, axis=1)
    df_value_counts = pd.concat(value_counts_per_question, axis=1)
    logfile.close()

    print("correlations calculation ready")

    # save the results per PRS in separated folders
    for PRS_name in df_betas_per_question.index:
        print("PRS_name", PRS_name)
        #get the export data per PRS from the dataframes
        df_prs_subset = df_betas_per_question.loc[PRS_name, :]
        df_prs_subset = df_prs_subset.unstack(level=-1)

        df_prs_subset_pvalues = df_pvalues_per_question.loc[PRS_name, :]
        df_prs_subset_pvalues = df_prs_subset_pvalues.unstack(level=-1)

        df_prs_subset_se_values = df_se_values_per_question.loc[PRS_name, :]
        df_prs_subset_se_values = df_prs_subset_se_values.unstack(level=-1)

        time_series_info_corr = df_prs_subset.copy()
        time_series_info_corr.columns = time_series_info_corr.columns.astype(str)

        # create the PRS output dir
        output_prs_dir = os.path.join(output_dir_ori_path, PRS_name)
        create_dir(output_prs_dir)

        # export the files
        df_prs_subset_pvalues.columns = df_prs_subset_pvalues.columns.astype(str)
        export_df(df_prs_subset_pvalues, output_prs_dir, PRS_name, "p_values_{}".format(suffix))

        df_prs_subset_se_values.columns = df_prs_subset_se_values.columns.astype(str)
        export_df(df_prs_subset_se_values, output_prs_dir, PRS_name, "se_values_{}".format(suffix))

        df_nvalues = df_nvalues.T
        print(df_nvalues)
        df_nvalues.columns = df_nvalues.columns.astype(str)
        export_df(df_nvalues, output_prs_dir, PRS_name, "n_values_{}".format(suffix))

        df_value_counts = df_value_counts.T
        df_value_counts.columns = df_value_counts.columns.astype(str)
        export_df(df_value_counts, output_prs_dir, PRS_name, "value_counts_{}".format(suffix))

        time_series_info_corr["Question answers"] = df_selected_questions.loc[time_series_info_corr.index, "Question answers"]
        time_series_info_corr = time_series_info_corr.loc[:, ["Question answers", *list(df_prs_subset.columns.astype(str))]]
        export_df(time_series_info_corr, output_prs_dir, PRS_name, "correlations_{}".format(suffix))


def export_df(df, output_path, trait_name, suffix):
    name = "covid_export_correlations_{trait_name}_{suffix}_{date}".format(
        trait_name=trait_name,
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
