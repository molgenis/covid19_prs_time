"""
File:         create_combined_matrixes.py
Created:      09/03/2021
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
import glob
import os
from datetime import datetime
import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
idx = pd.IndexSlice
from matplotlib.backends.backend_pdf import PdfPages

def main():
    selected_cols = ['1.0', '2.0', '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0', '10.0', '11.0', '12.0', '13.0', '14.0', '15.0', '15.5', '16.0', '16.5', '17.0']
    dates_dict = {
        1: "2020-03-30",
        2: "2020-04-06",
        3: "2020-04-13",
        4: "2020-04-20",
        5: "2020-04-27",
        6: "2020-05-04",
        7: "2020-05-18",
        8: "2020-06-01",
        9: "2020-06-15",
        10: "2020-07-06",
        11: "2020-07-13",
        12: "2020-08-10",
        13: "2020-09-07",
        14: "2020-10-12",
        15: "2020-11-02",
        15.5: "2020-11-16",
        16: "2020-11-30",
        16.5: "2020-12-14",
        17: "2021-01-11"
    }
    series_date = pd.to_datetime(pd.Series(dates_dict), format="%Y-%m-%d")
    days_differ = series_date - pd.to_datetime(dates_dict[1], format="%Y-%m-%d")
    days_differ = days_differ.dt.days

    input_dir = sys.argv[1] # input directory path (ouput directory of the meta analysis
    input_questions_path = sys.argv[2] # questionnaire ids input file path
    df_input_question_info = pd.read_csv(input_questions_path, sep="\t", index_col=0, dtype=str)
    output_dir = sys.argv[3] # output directory path
    question_selection_path = sys.argv[4]

    # Path to the translation of the questionnaire items
    translation_questions_path = sys.argv[4]
    df_translation_info = None
    if translation_questions_path is not None:
        df_translation_info = pd.read_csv(translation_questions_path, sep="\t", index_col="Question")

    # path the the translation of the pgs items
    translation_prs_path = sys.argv[5]
    df_translation_prs_info = None
    if translation_prs_path is not None:
        df_translation_prs_info = pd.read_csv(translation_prs_path, sep="\t", index_col="PRS's")

    # remove weight question, because we add the BMI
    questions_to_remove = ["wat is uw lichaamsgewicht (in kg)?"]

    #create output dir
    create_dir(output_dir)

    # read the quality control question information
    question_selection_path = "/Users/henrywiersma/Desktop/PRS_risk_vs_measures/results-3-5-2021/questions_after_selection_prs_single.pkl"
    df_question_qual_control = pd.read_pickle(question_selection_path)
    df_question_qual_control = df_question_qual_control.loc[:, selected_cols]

    # do not use quality control for the positive tested results
    df_question_qual_control.loc["Positive tested cumsum", :] = ["covt01_positive_tested_cumsum", "covt02_positive_tested_cumsum", "covt03_positive_tested_cumsum", "covt04_positive_tested_cumsum", "covt05_positive_tested_cumsum", "covt06_positive_tested_cumsum", "covt07_positive_tested_cumsum", "covt08_positive_tested_cumsum", "covt09_positive_tested_cumsum", "covt10_positive_tested_cumsum", "covt11_positive_tested_cumsum", "covt12_positive_tested_cumsum", "covt13_positive_tested_cumsum", "covt14_positive_tested_cumsum", "covt15_positive_tested_cumsum", "covt15b_positive_tested_cumsum", "covt16_positive_tested_cumsum", "covt16b_positive_tested_cumsum", "covt17_positive_tested_cumsum"]

    # find the input file paths
    export_paths = glob.glob(os.path.join(input_dir, "*", "covid_meta_analysis_z_scores_*.txt"))

    # create output lists
    columns = []
    columns_pvalues = []
    columns_sevalues = []
    columns_correlations = []
    columns_length = []

    # perform analysis per PGS
    for path in export_paths:
        # Read input data with z-scores
        pgs_name = path.split("/")[-2]
        df = pd.read_csv(path, sep="\t", index_col=0)
        df = df.loc[df.index.difference(questions_to_remove), selected_cols]

        # remove the z score values which does not pass the question quality control
        intersect_index = df.index.intersection(df_question_qual_control.index)
        intersect_col = df.columns.intersection(df_question_qual_control.columns)
        df = df.loc[intersect_index, intersect_col]
        df = df.loc[intersect_index, intersect_col].where(~df_question_qual_control.loc[intersect_index, intersect_col].isna(), np.nan)

        # read p-value data
        df_p_values_path = glob.glob(os.path.join(input_dir, pgs_name, "covid_meta_analysis_pvalues_*.txt"))[0]
        df_p_values = pd.read_csv(df_p_values_path, sep="\t", index_col=0)
        df_p_values = df_p_values.loc[df_p_values.index.difference(questions_to_remove), selected_cols]

        # remove the pvalues values which does not pass the question quality control
        df_p_values = df_p_values.loc[intersect_index, intersect_col]
        df_p_values = df_p_values.loc[intersect_index, intersect_col].where(~df_question_qual_control.loc[intersect_index, intersect_col].isna(), np.nan)

        # read the standard error data
        df_se_values_path = glob.glob(os.path.join(input_dir, pgs_name, "covid_meta_analysis_se_values_*.txt"))[0]
        df_se_values = pd.read_csv(df_se_values_path, sep="\t", index_col=0)

        # remove the standard errors values which does not pass the question quality control
        df_se_values = df_se_values.loc[df_se_values.index.difference(questions_to_remove), selected_cols]
        df_se_values = df_se_values.loc[intersect_index, intersect_col]
        df_se_values = df_se_values.loc[intersect_index, intersect_col].where(~df_question_qual_control.loc[intersect_index, intersect_col].isna(), np.nan)

        # Read the correlations input dataframe
        df_correlations_path = glob.glob(os.path.join(input_dir, pgs_name, "covid_meta_analysis_correlations_*.txt"))[0]
        df_correlations = pd.read_csv(df_correlations_path, sep="\t", index_col=0)

        # remove the correlations which does not pass the question quality control
        df_correlations = df_correlations.loc[intersect_index, intersect_col]
        df_correlations = df_correlations.loc[intersect_index, intersect_col].where(~df_question_qual_control.loc[intersect_index, intersect_col].isna(), np.nan)

        # Process dataframes to create a multi column dataframe
        df_se_values.columns = pd.MultiIndex.from_product([[pgs_name], df_se_values.columns])
        columns_sevalues.append(df_se_values)

        df_p_values.columns = pd.MultiIndex.from_product([[pgs_name], df_p_values.columns])
        columns_pvalues.append(df_p_values)

        df_correlations.columns = pd.MultiIndex.from_product([[pgs_name], df_correlations.columns])
        columns_correlations.append(df_correlations)

        df.columns = pd.MultiIndex.from_product([[pgs_name], df.columns])
        columns.append(df)

        columns_length.append(df.shape[1])

    df_total = pd.concat(columns, axis=1)
    df_total_pvalues = pd.concat(columns_pvalues, axis=1)
    df_total_sevalues = pd.concat(columns_sevalues, axis=1)
    df_total_correlations = pd.concat(columns_correlations, axis=1)

    # export dataframes
    export_df(df_total, output_dir, "meta_analysis_z_scores_all")
    export_df(df_total_correlations, output_dir, "meta_analysis_correlations_all")
    export_df(df_total_pvalues, output_dir, "meta_analysis_pvalues_all")
    export_df(df_total_sevalues, output_dir, "meta_analysis_sevalues_all")


    def func_first_value(x):
        if x.first_valid_index() is None:
            return np.NaN
        else:
            return x[x.first_valid_index()]

    def func_first_val_per_prs(x):
        # method to return the first time when a question is asked
        # except for the number of ever positive tested cases, use for that
        # question the last value
        x_unst = x.unstack()
        x_unst = x_unst.loc[:, selected_cols]

        if x.name == "Positive tested cumsum":
            x_unst = x_unst.loc[:, x_unst.columns[::-1]]
        x_values = x_unst.apply(func_first_value, axis=1)
        x_values.name = x.name
        x_values.index = pd.MultiIndex.from_product([x_values.index, ["first_value"]])
        return x_values

    # find the values of the first time when the quesiton was asked
    df_out_first_z_scores = df_total.apply(func_first_val_per_prs, axis=1)
    df_out_first_p_values = df_total_pvalues.apply(func_first_val_per_prs, axis=1)
    export_df(df_out_first_p_values, output_dir, "p_values_first_values")
    export_df(df_out_first_z_scores, output_dir, "z_scores_first_values")

    # filter on bonf significance for the heatmap
    df_out_first_p_values_bonf_selection = df_out_first_p_values < 0.05/(df_out_first_z_scores.shape[0] * df_out_first_z_scores.shape[1])
    df_out_first_z_scores_bonf_selection = df_out_first_z_scores.where(df_out_first_p_values_bonf_selection, np.nan)
    df_plot_z_score = df_out_first_z_scores_bonf_selection.loc[:, idx[:, "first_value"]]
    df_plot_z_score = df_plot_z_score.T.reset_index(level=1, drop=True).T
    df_plot_z_score = df_plot_z_score.dropna(how="all", axis=0)
    df_plot_z_score = df_plot_z_score.dropna(how="all", axis=1)

    #  Create heatmap
    df_plot_z_score_for_plot = df_plot_z_score
    export_df(df_plot_z_score_for_plot, output_dir, "z_scores_in_heatmap")
    df_plot_z_score_for_plot = df_plot_z_score_for_plot.fillna(0)
    df_plot_z_score_for_plot = df_plot_z_score_for_plot.rename(index=df_translation_info.iloc[:, 0].to_dict())

    sns.clustermap(df_plot_z_score_for_plot,
                   linewidth=0.1,
                   linecolor='gray',
                   figsize=(20,40),
                   standard_scale=None,
                   center=0.00,
                   cmap=sns.diverging_palette(220, 20, as_cmap=True))
    plt.xticks(size=12, rotation=45, ha='right')
    plt.yticks(size=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "heatmap_all_z_scores_bonf_first_values_strong_selection.png"), dpi=300, bbox_inches='tight')
    plt.clf()

    quest_selection = df_input_question_info.index
    jbf_z_score_table = df_plot_z_score.loc[df_plot_z_score.index.intersection(quest_selection)]

    plot_data = []
    jbf_results = []
    prs_per_question_for_long_analysis = {}

    #
    # perform the sensitivity analysis
    #

    # Loop through all the questions
    for question_index, perform_pgs in jbf_z_score_table.iterrows():
        number_of_questions = int(df_input_question_info.loc[question_index, "Number of timepoints"])
        pgs_per_question = []
        perform_pgs = perform_pgs.dropna().index

        # Loop through the PRS information
        for prs_index in perform_pgs:
            time_series_info = df_total_correlations.loc[question_index, idx[prs_index, :]]
            time_series_info = time_series_info.reset_index(level=0, drop=True)
            time_series_info = time_series_info.dropna()

            # two timepoints are nessesary to calculate the pearson correlation over time
            if number_of_questions > 2:
                pgs_per_question.append(prs_index)

                #calculate the pearson correlation over time
                time_series_info.index = time_series_info.index.astype(float)
                time_series_info = time_series_info.rename(index=days_differ)
                corr, pvalue = scipy.stats.pearsonr(time_series_info.index, time_series_info)

                # calculate the s zscore
                z_score = np.abs(scipy.stats.norm.ppf(1-(pvalue/2))) * np.sign(corr)

                # Find the correct translations of all the labels
                label_en = question_index
                if question_index in df_translation_info.index:
                    label_en = df_translation_info.loc[question_index, "label_en"]

                label_prs_en = prs_index
                if prs_index in df_translation_prs_info.index:
                    label_prs_en = df_translation_prs_info.loc[prs_index, "English Label"]

                # save the information
                jbf_result = {
                    "question": question_index,
                    "label_en": label_en,
                    "label_prs_en": label_prs_en,
                    "trait": prs_index,
                    "pvalue": pvalue,
                    "pearson_corr": corr,
                    "r2": corr**2,
                    "zscore": z_score
                }
                jbf_results.append(jbf_result)

                # save information to create the plot
                single_plot_info = {"title": question_index,
                                    "label_en": label_en,
                                    "label_prs_en": label_prs_en,
                                    "trait": prs_index,
                                    "data": time_series_info,
                                    "pvalue": pvalue,
                                    "r2": corr**2,
                                    "zscore": z_score}
                plot_data.append(single_plot_info)
        if len(pgs_per_question) > 0:
            prs_per_question_for_long_analysis[question_index] = ";".join(pgs_per_question)

    # export long results
    df_prs_per_question_for_long_analysis = pd.DataFrame.from_dict(prs_per_question_for_long_analysis, orient="index", columns=["GWAS"])
    df_prs_per_question_for_long_analysis.name = "Question"
    df_prs_per_question_for_long_analysis.to_csv(os.path.join(output_dir, "gwasses_to_perform_filtered.txt"), sep="\t")

    # export jbf results
    df_jbf_results = pd.DataFrame(jbf_results)
    df_jbf_results = df_jbf_results.set_index("question")
    export_df(df_jbf_results, output_dir, "jbf_results")

    # create plots
    plot_file_name = "jbf_correlations_plot_{date}.pdf".format(
        date=datetime.now().strftime("%d-%m-%Y")
    )
    pdf_pages = PdfPages(os.path.join(output_dir, plot_file_name))

    for index, data in enumerate(plot_data):
        fig = plt.figure(figsize=(11.69, 8.27), dpi=100)
        ax = plt.gca()
        sns.regplot(data["data"].index, data["data"],
                    ax=ax,
                    scatter_kws={"color": "dimgray", "edgecolor": "dimgray", 'clip_on':False, "linewidth": 0},
                    line_kws={"color": "black", "alpha": 0.8})

        # Create confidence interval lines
        question_se = df_total_sevalues.loc[data["title"], idx[data["trait"], :]]
        question_se = question_se.reset_index(level=0, drop=True)
        question_se = question_se.dropna()
        question_se.index = question_se.index.astype(float)
        question_se = question_se.rename(index=days_differ)
        for time_pont, se_val in question_se.iteritems():
            start = data["data"][time_pont] - (se_val * 1.96)
            end = data["data"][time_pont] + (se_val * 1.96)
            ax.plot([time_pont, time_pont], [start, end], '-', color ="dimgray")

        # style the plots
        ax.grid(False)
        ax.xaxis.set_ticks(np.arange(0, 301, 50))
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_position(('axes', -0.05))
        ax.yaxis.set_ticks_position('left')
        ax.spines['left'].set_position(('axes', -0.05))
        ax.spines['left'].set_color("dimgray")
        ax.spines['bottom'].set_color("dimgray")

        #create plot title
        pval_value, pval_exp = "{:0.2E}".format(data["pvalue"]).split("E")
        title = "{question}\n{prs}\n$R^2$: {r2:0.02f}, p-value: ${pval_value}x10^{{{pval_exp}}}$".format(
            question=shorten_string(data["label_en"], 50),
            prs=shorten_string(data["label_prs_en"], 50),
            r2=data["r2"],
            pvalue=data["pvalue"],
            pval_value=pval_value,
            pval_exp=pval_exp
        )

        # set the x and y axes
        days_key_list = [0, 98, 196, 287]
        days_name_list = ["30-Mar-2020", "06-Jul-2020", "12-Oct-2020", "11-Jan-2021"]

        ax.set_xticks(days_key_list)
        ax.set_xticklabels(days_name_list, rotation=45, ha='right', fontsize=8)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Date", fontsize=10)
        ax.set_ylabel(shorten_string("Regression coefficients of {question} modeled against {prs} PGS".format(
            question=data["label_en"],
            prs=data["label_prs_en"]), 50), fontsize=10)

        # save the plots
        plt.tight_layout(pad=5.0)
        pdf_pages.savefig()
        plt.clf()


    pdf_pages.close()
    return


def export_df(df, output_path, suffix):
    name = "covid_combined_matrix_{suffix}_{date}".format(
        suffix=suffix,
        date=datetime.now().strftime("%d-%m-%Y")
    )
    df.to_pickle(os.path.join(output_path, "{}.pkl".format(name)))
    df.to_csv(os.path.join(output_path, "{}.txt".format(name)), sep="\t")
    df.to_excel(os.path.join(output_path, "{}.xlsx".format(name)))
    print("export files created: {}.pkl/txt".format(name))


def create_dir(dir_name):
    if os.path.isdir(dir_name) == False:
        os.mkdir(dir_name)

def shorten_string(title, title_len=20):
    if isinstance(title, str):
        if len(title) > title_len:
            title = title.split(" ")
            title_length = np.cumsum([len(x) for x in title])
            title_id = np.argmax(title_length > title_len)
            return "{}\n{}".format(" ".join(title[:title_id]), " ".join(shorten_string(title[title_id:], title_len)))
    return title

if __name__ == '__main__':
    main()