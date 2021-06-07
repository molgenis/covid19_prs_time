"""
File:         combine_pgs_scores.py
Created:      26/02/2021
Last Changed: 26/02/2021
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
import glob
import os
import pandas as pd


def main():
    # input parameters
    input_pgs_dir = sys.argv[1] # Directory with the PGS results
    mapping_file = sys.argv[2] #Mapping file of the PGS with the name and direcory name information
    output_dir_path = sys.argv[3] # Output directory
    chip_type = sys.argv[4] # The Chip type
    filter_path = sys.argv[5] # Filter file to remove PGS results of participants
    # create output dir
    create_dir(output_dir_path)

    # read PGS mapping file, to get the trait names
    pgs_mapping = read_mapping_file(mapping_file)

    # read the PGS files and merge them
    df_pgs = read_pgs_files(input_pgs_dir, pgs_mapping)

    # connect the PROJECT ids to the PGS traits
    df_pgs = set_index_to_project_id(df_pgs, chip_type)

    df_pgs = df_pgs.loc[df_pgs.index.dropna(), :]

    if chip_type != "cyto":
        filter_df = pd.read_csv(filter_path, sep="\t", index_col=0)
        df_pgs = df_pgs.loc[df_pgs.index.difference(filter_df.index)]

    # export the output file
    export_df(df_pgs, output_dir_path, chip_type)


def read_mapping_file(mapping_file_path):
    df_mapping = pd.read_csv(mapping_file_path, sep="\t",
                             encoding='utf8', index_col=False, engine='python')
    df_mapping = df_mapping.set_index("summaryStatistics")
    return df_mapping.loc[:, "trait"].to_dict()


def read_pgs_files(input_path, mapping_dict):
    pgs_path = os.path.join(input_path, "*", "full.*.pgs.profile")
    total_df = pd.DataFrame()
    # read all pgs paths
    for path in glob.glob(pgs_path):
        # find PGS name
        pgs_name = path.split("/")[-2]
        if pgs_name in mapping_dict:
            pgs_name = mapping_dict[pgs_name]
        else:
            print("PGS not in mapping file:", pgs_name)
        print("process:", pgs_name)

        # read file
        df_pgs = pd.read_csv(path, sep="\t")

        #change values to float and set the index
        df_pgs["SCORESUM"] = df_pgs["SCORESUM"].astype(float)
        df_pgs = df_pgs.set_index("IID")
        df_pgs = df_pgs.loc[:, ["SCORESUM"]]

        # rename the column to the pgs name
        df_pgs = df_pgs.rename(columns=dict(SCORESUM=pgs_name))

        # merge the pgs with the total df
        if total_df.shape[0] > 0:
            total_df = pd.merge(total_df, df_pgs,
                                left_index=True,
                                right_index=True,
                                suffixes=[None, "_2"],
                                how="outer")
        else:
            total_df = df_pgs
    return total_df


def set_index_to_project_id(df, chip_type):
    if chip_type == "cyto":
        # path to the participant-chip id file to the participant file
        chip_id_to_pseudo_id_path = "//groups/umcg-lifelines/tmp01/releases/cytosnp_linkage_files/v4/cytosnp_linkage_file.dat"
        match_column = "cytosnp_ID"
    else:
        # path to the participant-chip id file to the participant file
        chip_id_to_pseudo_id_path = "/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat"
        match_column = "UGLI_ID"

    # path to the old project id to the new project id
    project_id_to_pseudo_id_path = "/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v1/phenotype_linkage_file_project_pseudo_id.txt"
    project_id_to_pseudo_id = pd.read_csv(project_id_to_pseudo_id_path,
                                          sep="\t", dtype=str)

    #add project ids via pseudoids
    chip_id_to_pseudo_id = pd.read_csv(chip_id_to_pseudo_id_path,
                                       sep="\t", dtype=str)
    chip_id_to_pseudo_id = pd.merge(chip_id_to_pseudo_id,
                                    project_id_to_pseudo_id,
                                    on="PSEUDOIDEXT", how="left")

    # add project id to df with PGS scores
    df = pd.merge(df, chip_id_to_pseudo_id, left_index=True,
                  right_on=match_column, how="left")
    df = df.set_index("PROJECT_PSEUDO_ID", drop=True)
    df = df.loc[:, df.columns.difference(["PSEUDOIDEXT",
                                          match_column, "genotyping_name"])]
    return df


def export_df(df, output_path, suffix):
    name = "PGS_combined_{suffix}_{date}".format(
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
