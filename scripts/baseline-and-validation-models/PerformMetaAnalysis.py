"""
File:         meta_analysis_pgs_scores.py
Created:      08/03/2021
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


def main():
    ugli_dir_path = sys.argv[1] # path to the directory with the Global Screening Array
    cyto_dir_path = sys.argv[2] # path to the directory with the HumanCytoSNP-12
    output_dir = sys.argv[3] # output directory path

    # create output directory
    create_dir(output_dir)

    # find all the input directory paths and extract the PRS names
    pgs_names_ugli = glob.glob(os.path.join(ugli_dir_path, "*", "covid_export_correlations_*_correlations_*.txt"))
    pgs_names_ugli = [x.split("/")[-2] for x in pgs_names_ugli]
    pgs_names_cyto = glob.glob(os.path.join(cyto_dir_path, "*", "covid_export_correlations_*_correlations_*.txt"))
    pgs_names_cyto = [x.split("/")[-2] for x in pgs_names_cyto]
    pgs_names = set([*pgs_names_ugli, *pgs_names_cyto])

    # perform meta analysis per PRS
    for pgs_name in pgs_names:
        print("Process PGS: {}".format(pgs_name))

        # get the input file paths
        # ugli input paths
        perform_meta_analysis = True
        ugli_corr_path = glob.glob(os.path.join(ugli_dir_path, pgs_name, "covid_export_correlations_*_correlations_*.txt"))
        if len(ugli_corr_path) > 0:
            ugli_corr_path = ugli_corr_path[0]
            ugli_se_values_path = glob.glob(os.path.join(ugli_dir_path, pgs_name, "covid_export_correlations_*_se_values_*.txt"))[0]
        else:
            perform_meta_analysis = False
            print("Not in ugli")

        # cyto input paths
        cyto_corr_path = glob.glob(os.path.join(cyto_dir_path, pgs_name, "covid_export_correlations_*_correlations_*.txt"))
        if len(cyto_corr_path) > 0:
            cyto_corr_path = cyto_corr_path[0]
            cyto_se_values_path = glob.glob(os.path.join(cyto_dir_path, pgs_name, "covid_export_correlations_*_se_values_*.txt"))[0]
        else:
            perform_meta_analysis = False
            print("Not in cyto")

        # perform meta analysis if PRS is present in both data sets
        if perform_meta_analysis:
            # select columns
            selected_cols = ['1.0', '2.0', '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0', '10.0', '11.0', '12.0', '13.0', '14.0', '15.0', '15.5', '16.0', '16.5', '17.0']

            # read the input files
            df_ugli_corr = pd.read_csv(ugli_corr_path, sep="\t", index_col=0)
            df_cyto_corr = pd.read_csv(cyto_corr_path, sep="\t", index_col=0)
            index_selection = df_ugli_corr.index.intersection(df_cyto_corr.index)

            df_ugli_se = pd.read_csv(ugli_se_values_path, sep="\t", index_col=0)
            df_cyto_se = pd.read_csv(cyto_se_values_path, sep="\t", index_col=0)

            # remove unnsessesary rows and columns
            df_ugli_corr = df_ugli_corr.loc[index_selection, selected_cols]
            df_cyto_corr = df_cyto_corr.loc[index_selection, selected_cols]
            df_ugli_se = df_ugli_se.loc[index_selection, selected_cols]
            df_cyto_se = df_cyto_se.loc[index_selection, selected_cols]

            # remove values which are not present in both dataframes
            df_ugli_corr_original = df_ugli_corr.copy()
            df_ugli_se_original = df_ugli_se.copy()
            selection = np.logical_and(~df_ugli_corr.loc[index_selection, selected_cols].isna(), ~df_cyto_corr.loc[index_selection, selected_cols].isna())
            df_ugli_corr.loc[index_selection, selected_cols] = df_ugli_corr.loc[index_selection, selected_cols].where(selection, np.nan)
            df_cyto_corr.loc[index_selection, selected_cols] = df_cyto_corr.loc[index_selection, selected_cols].where(selection, np.nan)

            df_ugli_se.loc[index_selection, selected_cols] = df_ugli_se.loc[index_selection, selected_cols].where(selection, np.nan)
            df_cyto_se.loc[index_selection, selected_cols] = df_cyto_se.loc[index_selection, selected_cols].where(selection, np.nan)

            # Perform meta analysis with the inverse-variance weighting method
            df_ugli_se2 = df_ugli_se**2
            df_cyto_se2 = df_cyto_se**2

            yDivededSe2_ugli = df_ugli_corr / df_ugli_se2
            yDivededSe2_cyto = df_cyto_corr / df_cyto_se2
            yDivededSe2_total = yDivededSe2_ugli + yDivededSe2_cyto

            oneDividedBySe2_ugli = 1 / df_ugli_se2
            oneDividedBySe2_cyto = 1 / df_cyto_se2
            oneDividedBySe2_total = oneDividedBySe2_ugli + oneDividedBySe2_cyto

            meta_y = yDivededSe2_total / oneDividedBySe2_total
            meta_se = np.sqrt(1/oneDividedBySe2_total)
            meta_z = meta_y / meta_se
            meta_p = pd.DataFrame(2*scipy.stats.norm.cdf(-np.abs(meta_z.loc[index_selection, selected_cols])), index=index_selection, columns=selected_cols)

            # fill correlations and standard errors with the
            # Global Screening Array data if the value is not present
            # in the meta analysis
            meta_y = meta_y.fillna(df_ugli_corr_original)
            meta_se = meta_se.fillna(df_ugli_se_original)

            #  create pgs output directory
            output_meta_path = os.path.join(output_dir, pgs_name)
            create_dir(output_meta_path)

            #export the files
            output_corr_path = "covid_meta_analysis_correlations_{pgs}_{date}.txt".format(
                pgs=pgs_name,
                date=datetime.now().strftime("%d-%m-%Y")
            )
            output_corr_path = os.path.join(output_meta_path, output_corr_path)
            meta_y.loc[index_selection, selected_cols].to_csv(output_corr_path, sep="\t")

            output_se_values_path = "covid_meta_analysis_se_values_{pgs}_{date}.txt".format(
                pgs=pgs_name,
                date=datetime.now().strftime("%d-%m-%Y")
            )
            output_se_values_path = os.path.join(output_meta_path, output_se_values_path)
            meta_se.loc[index_selection, selected_cols].to_csv(output_se_values_path, sep="\t")

            output_z_score_path = "covid_meta_analysis_z_scores_{pgs}_{date}.txt".format(
                pgs=pgs_name,
                date=datetime.now().strftime("%d-%m-%Y")
            )
            output_z_score_path = os.path.join(output_meta_path, output_z_score_path)
            meta_z.loc[index_selection, selected_cols].to_csv(output_z_score_path, sep="\t")

            output_pvalues_path = "covid_meta_analysis_pvalues_{pgs}_{date}.txt".format(
                pgs=pgs_name,
                date=datetime.now().strftime("%d-%m-%Y")
            )
            output_pvalues_path = os.path.join(output_meta_path, output_pvalues_path)
            meta_p.loc[index_selection, selected_cols].to_csv(output_pvalues_path, sep="\t")


def export_df(df, output_path, trait_name, value_info):
    name = "covid_meta_analysis_{value_info}_{trait_name}_{date}".format(
        value_info=value_info,
        trait_name=trait_name,
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