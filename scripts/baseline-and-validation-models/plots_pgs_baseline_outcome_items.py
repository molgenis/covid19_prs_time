"""
File:         plots_values_pgs.py
Created:      12/05/2021
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
import glob
import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
import json
from matplotlib.backends.backend_pdf import PdfPages

def main():
    selected_cols_ids = ['1.0', '2.0', '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0', '10.0', '11.0', '12.0', '13.0', '14.0', '15.0', '15.5', '16.0', '16.5', '17.0']
    remove_questions = ["pokken / kunt u aangeven welke van de onderstaande (kinder)ziektes u gehad hebt?", "op hoeveel momenten van de dag eet u iets?"]
    input_df_path = sys.argv[1] # input questionnaire data
    input_question_overview = sys.argv[2] # input of the questionnaire items
    input_recode_info = sys.argv[3] # input file with the recoding and model information
    input_question_pgs_combinations = sys.argv[4] # selection file of the question x prs combinations
    input_question_pgs_combinations_pvalues = sys.argv[5] #  pvalues from the meta analysis
    input_pgs_path_ugli = sys.argv[6] # path to the directory with the Global Screening Array input data
    input_pgs_path_cyto = sys.argv[7] # path to the directory with the HumanCytoSNP-12 input data
    input_analysis_output_ugli_dir = sys.argv[8] # path to the directory with the Global Screening Array output results
    input_analysis_output_cyto_dir = sys.argv[9] # path to the directory with the HumanCytoSNP-12 output results
    translation_questions_path = sys.argv[10] # input path with the question translations
    translation_prs_path = sys.argv[11]  # input path wit the prs translations
    output_dir = sys.argv[12] # output dit

    # create output dit
    create_dir(output_dir)

    # read the input data
    df = pd.read_pickle(input_df_path)
    df_question_ids_total = pd.read_csv(input_question_overview, sep="\t", index_col=0, dtype="str")
    df_recode_info = pd.read_pickle(input_recode_info)
    df_recode_info = df_recode_info.set_index("question_id")
    df_question_pgs = pd.read_pickle(input_question_pgs_combinations)
    df_question_pgs_pvalues = pd.read_pickle(input_question_pgs_combinations_pvalues)

    #process ugli
    df_pgs_ugli = pd.read_pickle(input_pgs_path_ugli)
    df_pgs_ugli.columns = df_pgs_ugli.columns.str.replace("/", ".").str.replace(" ", ".").str.replace("-", ".").str.replace("(", ".").str.replace(")", ".")
    df_pgs_ugli = df_pgs_ugli.loc[df_pgs_ugli.index.intersection(df.index), :]

    #process cyto
    df_pgs_cyto = pd.read_pickle(input_pgs_path_cyto)
    df_pgs_cyto.columns = df_pgs_cyto.columns.str.replace("/", ".").str.replace(" ", ".").str.replace("-", ".").str.replace("(", ".").str.replace(")", ".")
    df_pgs_cyto = df_pgs_cyto.loc[df_pgs_cyto.index.intersection(df.index), :]

    # process the translations
    df_translation_info = None
    if translation_questions_path is not None:
        df_translation_info = pd.read_csv(translation_questions_path, sep="\t", index_col="Question")

    df_translation_prs_info = None
    if translation_prs_path is not None:
        df_translation_prs_info = pd.read_csv(translation_prs_path, sep="\t", index_col="PRS's")

    # create output lists
    plot_data = []

    # process the questions per question
    for index, row in df_question_pgs.iterrows():
        if index not in remove_questions:
            single_question_ids = df_question_ids_total.loc[index, selected_cols_ids]

            # find the first question id (except for the ever positive tested, we there for use the last question id)
            if index == "Positive tested cumsum":
                single_question_ids = single_question_ids.iloc[::-1]
            single_question_first_id = func_first_value(single_question_ids)

            # get model type to create the correct plot
            model_type = df_recode_info.loc[single_question_first_id, "model_type"]
            answer_options_str = None
            if index in df_translation_info.index:
                answer_options_str = df_translation_info.loc[index, "Answers_options_plots"]

            if isinstance(model_type, pd.core.series.Series):
                model_type = model_type.iloc[0]

            answer_options = dict()
            if answer_options_str is not None and pd.isna(answer_options_str) == False:
                if isinstance(answer_options_str, pd.core.series.Series):
                    answer_options_str = answer_options_str.iloc[0]
                answer_options = json.loads(answer_options_str)

            # get the right plot type
            graph_type = "boxplots"
            if model_type == "gaussian":
                graph_type = "scatter"

            row_filtered = row.dropna()

            # loop through the PRS to create the plot per PRS
            for pgs_id, single_z_score in row_filtered.iteritems():
                single_pvalue = float(df_question_pgs_pvalues.loc[index, pgs_id])
                y_ugli = df_pgs_ugli.loc[:, pgs_id]
                x_ugli = df.loc[y_ugli.index, single_question_first_id]
                x_ugli = x_ugli.dropna()

                # fix plotting issue
                if single_question_first_id == "covt16_christmas_adu_q_1_b":
                    x_ugli = x_ugli.replace({3.0: 2.0})
                y_ugli = y_ugli.loc[x_ugli.index]

                y_cyto = df_pgs_cyto.loc[:, pgs_id]
                x_cyto = df.loc[y_cyto.index, single_question_first_id]
                x_cyto = x_cyto.dropna()

                # fix plotting issue
                if single_question_first_id == "covt16_christmas_adu_q_1_b":
                    x_cyto = x_cyto.replace({3.0: 2.0})
                y_cyto = y_cyto.loc[x_cyto.index]

                # calculate z scores and p values for the
                # individual plots of ugli and cyto
                z_score_ugli, pvalue_ugli = caalculate_z_scrore_from_models(input_analysis_output_ugli_dir, pgs_id, index, single_question_first_id)
                z_score_cyto, pvalue_cyto = caalculate_z_scrore_from_models(input_analysis_output_cyto_dir, pgs_id, index, single_question_first_id)

                # Find the correct translations of question and prs
                label_en = index
                if index in df_translation_info.index:
                    label_en = df_translation_info.loc[index, "label_en"]

                label_prs_en = pgs_id
                if pgs_id in df_translation_prs_info.index:
                    label_prs_en = df_translation_prs_info.loc[pgs_id, "English Label"]

                # Save the plotting information
                single_plot_data = {
                    "x_ugli": x_ugli,
                    "x_cyto": x_cyto,
                    "y_ugli": y_ugli,
                    "y_cyto": y_cyto,
                    "question_title": label_en,
                    "pgs_title": label_prs_en,
                    "graph_type": graph_type,
                    "answers": answer_options,
                    "z_score": single_z_score,
                    "pvalue": single_pvalue,
                    "z_score_ugli": z_score_ugli,
                    "pvalue_ugli": pvalue_ugli,
                    "z_score_cyto": z_score_cyto,
                    "pvalue_cyto": pvalue_cyto,
                }
                plot_data.append(single_plot_data)

    ###
    ### create plots
    ###

    # The PDF document
    plot_file_name = "question_pgs_plots_{date}.pdf".format(
        date=datetime.now().strftime("%d-%m-%Y")
    )
    pdf_pages = PdfPages(os.path.join(output_dir, plot_file_name))

    for index, data in enumerate(plot_data):
        fig = plt.figure(figsize=(11.69, 8.27), dpi=100)
        figure_title = "{question}\nPGS: {pgs}\nMeta analysis Z-score: {zscore:0.02f}, p-value: {pvalue:0.2E}".format(
            question=data["question_title"],
            pgs=data["pgs_title"],
            zscore=data["z_score"],
            pvalue=data["pvalue"]
        )
        fig.suptitle(figure_title, y=0.96, fontsize=12)

        ##
        ## UGLI
        ##
        plt.subplot2grid((1,2), (0,0), fig=fig)
        ax_ugli = plt.gca()
        if data["graph_type"] == "scatter":
            sns.regplot(data["x_ugli"], data["y_ugli"], ax=ax_ugli,
                        scatter_kws={"color": "black",  "s":20, 'alpha':0.3, "rasterized": True, "clip_on": False})
        else:
            sns.boxplot(data["x_ugli"], data["y_ugli"], ax=ax_ugli, color="w", showfliers = False)

            # add n-values above the boxes
            unique_labels = data["x_ugli"].unique()
            n_value_labels = []
            hivalues = []
            lovalues = []
            for unique_label in unique_labels:
                single_data = data["y_ugli"].where(data["x_ugli"] == unique_label).dropna()
                q1, med, q3 = np.percentile(single_data, [25, 50, 75])
                iqr = q3 - q1
                hival = q3 + 1.5 * iqr
                loval = q1 - 1.5 * iqr

                wiskhi = single_data[single_data <= hival]
                if len(wiskhi) == 0 or np.max(wiskhi) < q3:
                    hival = q3
                else:
                    hival = np.max(wiskhi)
                hivalues.append(hival)
                wisklo = single_data[single_data >= loval]
                if len(wisklo) == 0 or np.min(wisklo) > q1:
                    loval = q1
                else:
                    loval = np.min(wisklo)
                lovalues.append(loval)
                if data["question_title"] == "Average time spend sitting per weekend day" or data["question_title"] == "Average time spend sitting per working day":
                    unique_label = unique_label -1
                n_value_labels.append({
                    "x_pos": unique_label -1 if min(unique_labels) != 0 else unique_label,
                    "y_pos": hival,
                    "label": "n={}".format(single_data.shape[0])
                })
            text_distance = (max(hivalues) - min(lovalues)) * 0.005
            for label in n_value_labels:
                if unique_labels.shape[0] < 6:
                    ax_ugli.text(label["x_pos"], label["y_pos"] + text_distance, label["label"] , ha='center', va='bottom', color="gray")
                else:
                    ax_ugli.text(label["x_pos"], label["y_pos"] + text_distance, label["label"] , ha='center', va='bottom', color="gray", fontsize=6)

        # Add title
        title_ugli = "Global Screening Array\nZ-score: {zscore:0.02f}, p-value: {pvalue:0.2E}".format(
            zscore=data["z_score_ugli"],
            pvalue=data["pvalue_ugli"]
        )
        ax_ugli.set_title(title_ugli, fontsize=12)

        # set the labels
        ax_ugli.set_xlabel(short_str(data["question_title"], 50), fontsize=10)
        ax_ugli.set_ylabel("PRS: {}".format(short_str(data["pgs_title"], 50)), fontsize=10)

        # edit the design of the plot
        ax_ugli.grid(False)
        ax_ugli.spines['right'].set_color('none')
        ax_ugli.spines['top'].set_color('none')
        ax_ugli.spines['bottom'].set_position(('axes', -0.05))
        ax_ugli.yaxis.set_ticks_position('left')
        ax_ugli.spines['left'].set_position(('axes', -0.05))
        ax_ugli.spines['left'].set_color("dimgray")
        ax_ugli.spines['bottom'].set_color("dimgray")

        ##
        ## CYTO
        ##
        plt.subplot2grid((1,2), (0,1), fig=fig)
        ax_cyto = plt.gca()

        if data["graph_type"] == "scatter":
            sns.regplot(data["x_cyto"], data["y_cyto"], ax=ax_cyto,
                        scatter_kws={"color": "black",  "s":20, 'alpha':0.3, "rasterized": True, "clip_on": False})
        else:
            sns.boxplot(data["x_cyto"], data["y_cyto"], ax=ax_cyto, color="w", showfliers = False)

            # Add the n-values above the boxes
            unique_labels = data["x_cyto"].unique()
            if True:
                n_value_labels = []
                hivalues = []
                lovalues = []
                for unique_label in unique_labels:
                    single_data = data["y_cyto"].where(data["x_cyto"] == unique_label).dropna()
                    q1, med, q3 = np.percentile(single_data, [25, 50, 75])
                    iqr = q3 - q1
                    hival = q3 + 1.5 * iqr
                    loval = q1 - 1.5 * iqr

                    wiskhi = single_data[single_data <= hival]
                    if len(wiskhi) == 0 or np.max(wiskhi) < q3:
                        hival = q3
                    else:
                        hival = np.max(wiskhi)
                    hivalues.append(hival)
                    wisklo = single_data[single_data >= loval]
                    if len(wisklo) == 0 or np.min(wisklo) > q1:
                        loval = q1
                    else:
                        loval = np.min(wisklo)
                    lovalues.append(loval)
                    if data["question_title"] == "Average time spend sitting per weekend day" or data["question_title"] == "Average time spend sitting per working day":
                        unique_label = unique_label -1
                    n_value_labels.append({
                        "x_pos": unique_label -1 if min(unique_labels) != 0 else unique_label,
                        "y_pos": hival,
                        "label": "n={}".format(single_data.shape[0])
                    })
                text_distance = (max(hivalues) - min(lovalues)) * 0.005
                for label in n_value_labels:
                    if unique_labels.shape[0] < 6:
                        ax_cyto.text(label["x_pos"], label["y_pos"] + text_distance, label["label"] , ha='center', va='bottom', color="gray")
                    else:
                        ax_cyto.text(label["x_pos"], label["y_pos"] + text_distance, label["label"] , ha='center', va='bottom', color="gray", fontsize=6)

        # Add the title
        title_ugli = "HumanCytoSNP-12\nZ-score: {zscore:0.02f}, p-value: {pvalue:0.2E}".format(
            zscore=data["z_score_cyto"],
            pvalue=data["pvalue_cyto"]
        )
        ax_cyto.set_title(title_ugli, fontsize=12)

        # Set the labels
        ax_cyto.set_xlabel(short_str(data["question_title"], 50), fontsize=10)
        ax_cyto.set_ylabel("PRS: {}".format(short_str(data["pgs_title"], 50)), fontsize=10)

        # edit the design of the plot
        ax_cyto.grid(False)
        ax_cyto.spines['right'].set_color('none')
        ax_cyto.spines['top'].set_color('none')
        ax_cyto.spines['bottom'].set_position(('axes', -0.05))
        ax_cyto.yaxis.set_ticks_position('left')
        ax_cyto.spines['left'].set_position(('axes', -0.05))
        ax_cyto.spines['left'].set_color("dimgray")
        ax_cyto.spines['bottom'].set_color("dimgray")

        #change tick labels
        if len(data["answers"]) > 0:
            answers_df = pd.DataFrame.from_dict(data["answers"], orient="index", columns=["answer_str"])
            answers_df.index = answers_df.index.astype("float")
            if np.min(answers_df.index) > 0:
                if data["graph_type"] != "scatter":
                    answers_df.index = answers_df.index - 1
            answers_df = answers_df.sort_index()
            answers_df["answer_str"] = answers_df["answer_str"].replace("NaN", np.nan)
            answers_df["answer_str"] = answers_df["answer_str"].apply(short_str)
            answers_df = answers_df.dropna()

            ax_ugli.set_xticks(list(answers_df.index))
            ax_ugli.set_xticklabels(list(answers_df["answer_str"]), rotation=45, ha='right')
            ax_cyto.set_xticks(list(answers_df.index))
            ax_cyto.set_xticklabels(list(answers_df["answer_str"]), rotation=45, ha='right')

        # save the plots
        plt.tight_layout()
        plt.subplots_adjust(top=0.80, bottom=0.27, left=0.13, right=0.94, wspace=0.3)
        pdf_pages.savefig()
        plt.clf()

    pdf_pages.close()
    print("end script")

def caalculate_z_scrore_from_models(pgs_output_path, pgs, question, question_id):

    #find the pvalues and correlations files
    file_path_corr = os.path.join(pgs_output_path, pgs, "covid_export_correlations_*_correlations_*.pkl")
    file_path_pvalues = os.path.join(pgs_output_path, pgs, "covid_export_correlations_*_p_values_*.pkl")
    file_path_corr = glob.glob(file_path_corr)[0]
    file_path_pvalues = glob.glob(file_path_pvalues)[0]

    # read the files
    df_corr = pd.read_pickle(file_path_corr)
    df_pvalues = pd.read_pickle(file_path_pvalues)

    # fint the corresponded time point from the id
    time_point = question_id.split("_")[0].replace("covt", "")
    time_point = str(float(time_point.replace("b", ".5")))
    corr_value = df_corr.loc[question, time_point]
    pvalue_value = df_pvalues.loc[question, time_point]

    # calculate the z-score
    z_score = np.abs(scipy.stats.norm.ppf((pvalue_value/2))) * np.sign(corr_value)
    return z_score, pvalue_value


def func_first_value(x):
    if x.first_valid_index() is None:
        return np.NaN
    else:
        return x[x.first_valid_index()]

def create_dir(dir_name):
    if os.path.isdir(dir_name) == False:
        os.mkdir(dir_name)

def short_str(title, title_len=20):
    if isinstance(title, str):
        if len(title) > title_len:
            title = title.split(" ")
            title_length = np.cumsum([len(x) for x in title])
            title_id = np.argmax(title_length > title_len)
            return "{}\n{}".format(" ".join(title[:title_id]), " ".join(short_str(title[title_id:], 10)))
    return title

if __name__ == '__main__':
    main()