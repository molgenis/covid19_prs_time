import sys
from datetime import datetime
import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    questionniare_id_file_path = sys.argv[1] # input path to the questionnaire id file
    input_questionair_dir_path = sys.argv[2] # input path of the questionnaires directory
    participants_length_file_path = sys.argv[3] # input path to the pseudo id to participant id file
    output_dir_path = sys.argv[4]
    df_quest = pd.read_csv(questionniare_id_file_path, sep="\t", index_col=0, dtype=str)

    # create dir
    create_dir(output_dir_path)
    cols_corona_test_outcome = {
        "t01": "covt01_infection_adu_q_1_a",
        "t02": "covt02_infection_adu_q_1_a",
        "t03": "covt03_infection_adu_q_1_a",
        "t04": "covt04_infection_adu_q_1_a",
        "t05": "covt05_infection_adu_q_1_a",
        "t06": "covt06_infection_adu_q_2_a",
        "t07": "covt07_infection_adu_q_2_a",
        "t08": "covt08_infection_adu_q_2_a",
        "t09": "covt09_infection_adu_q_2_a",
        "t10": "covt10_infection_adu_q_2_a",
        "t11": "covt11_infection_adu_q_2_a",
        "t12": "covt12_infection_adu_q_2_a",
        "t13": "covt13_infection_adu_q_2_a",
        "t14": "covt14_infection_adu_q_2_a",
        "t15": "covt15_infection_adu_q_2_a",
        "t15b": "covt15b_coronatest_adu_q_1_a",
        "t16": "covt16_coronatest_adu_q_1_a",
        "t16b": "covt16b_coronatest_adu_q_1_a",
        "t17": "covt17_coronatest_adu_q_1_a",
    }
    # 1 = positief getest
    # 2 = negatief getest

    # merge questionnaires
    df = merge_new_format(input_questionair_dir_path, cols_corona_test_outcome)

    #  add covariate columns
    df = correct_weight(df)
    df = correct_alcohol(df)
    df = postive_test_cumsum(df)
    df = add_bmi(df, participants_length_file_path)
    df = add_combined_mini_columns(df, df_quest)
    df = add_recent_columns(df)
    df = add_measure_score(df)
    export_df(df, output_dir_path, "inclusive_correction_columns")
    print(df)
    print(list(df.columns))

def correct_time_values(row, min_value=20, max_value=220, max_diff=20, max_diff_thirth_week=30):
    new_row = row.dropna()
    new_row[np.logical_or(new_row < min_value, new_row > max_value)] = np.nan
    last_reported_value = None
    for index, col in new_row.iteritems():
        index_pos = new_row.index.get_loc(index)
        if index_pos > 0:
            last_val = new_row.iloc[index_pos - 1]
            if pd.isna(last_val):
                if last_reported_value != np.nan and np.abs(last_reported_value - col) > max_diff:
                    new_row[index] = np.nan
                else:
                    last_reported_value = col
            else:
                if last_val != np.nan and np.abs(last_val - col) > max_diff:
                    new_row[index] = np.nan
                else:
                    last_reported_value = col
        else:
            if new_row.shape[0] > 2:
                next_val = new_row.iloc[1]
                sec_next_val = new_row.iloc[2]
                if next_val != np.nan and sec_next_val != np.nan and np.abs(next_val - col) > max_diff and np.abs(sec_next_val - col) > max_diff_thirth_week:
                    new_row[index] = np.nan
            elif new_row.shape[0] > 1:
                next_val = new_row.iloc[1]
                if next_val != np.nan and np.abs(next_val - col) > max_diff:
                    new_row[index] = np.nan
                    new_row.iloc[1] = np.nan
            last_reported_value = col

    row[new_row.index] = new_row
    return row


def correct_weight(df):
    weight_cols = ["covt01_bodyweight_adu_q_1", "covt02_bodyweight_adu_q_1",
                   "covt03_bodyweight_adu_q_1", "covt04_bodyweight_adu_q_1",
                   "covt05_bodyweight_adu_q_2", "covt06_bodyweight_adu_q_2",
                   "covt07_bodyweight_adu_q_2", "covt08_bodyweight_adu_q_2",
                   "covt09_bodyweight_adu_q_2", "covt11_bodyweight_adu_q_2",
                   "covt12_bodyweight_adu_q_2", "covt13_bodyweight_adu_q_2",
                   "covt14_bodyweight_adu_q_2", "covt15_bodyweight_adu_q_2",
                   "covt16_bodyweight_adu_q_2", "covt17_bodyweight_adu_q_2"]
    df.loc[:, weight_cols] = df.loc[:, weight_cols].apply(lambda row: correct_time_values(row), axis=1)
    return df


def add_bmi(df, length_df_path, project_id_to_pseudo_id_path):
    project_id_to_pseudo_id = pd.read_csv(project_id_to_pseudo_id_path, sep="\t", dtype=str)
    length_df = pd.read_csv(length_df_path, sep="\t")
    length_df.loc[:, "PSEUDOIDEXT"] = length_df.loc[:, "PSEUDOIDEXT"].astype(str)

    # process dataframe with person length
    length_df.loc[:, "ENCOUNTERCODE_NUMBER"] = length_df.loc[:, "ENCOUNTERCODE"].replace({'Baseline assessment (1A)': 2, 'Second assessment (2A)': 1})
    length_df = length_df.sort_values(["PSEUDOIDEXT", "ENCOUNTERCODE_NUMBER"])
    length_df = length_df.drop_duplicates(keep="first", subset=["PSEUDOIDEXT"])
    length_df = pd.merge(length_df, project_id_to_pseudo_id, on="PSEUDOIDEXT", how="left")
    length_df = length_df.set_index("PROJECT_PSEUDO_ID")
    length_df.loc[:, "LENGTE"] = length_df.loc[:, "LENGTE"].replace({" ": np.nan})
    length_df.loc[:, "LENGTE"] = length_df.loc[:, "LENGTE"].astype(float)
    length_df = length_df.rename(columns={"LENGTE": "body_length"})
    length_df = length_df.loc[:, ["body_length"]]
    df = pd.merge(df.loc[:, df.columns.difference(length_df.columns)], length_df, left_index=True, right_index=True)
    weight_cols = ["covt01_bodyweight_adu_q_1", "covt02_bodyweight_adu_q_1",
                   "covt03_bodyweight_adu_q_1", "covt04_bodyweight_adu_q_1",
                   "covt05_bodyweight_adu_q_2", "covt06_bodyweight_adu_q_2",
                   "covt07_bodyweight_adu_q_2", "covt08_bodyweight_adu_q_2",
                   "covt09_bodyweight_adu_q_2", "covt11_bodyweight_adu_q_2",
                   "covt12_bodyweight_adu_q_2", "covt13_bodyweight_adu_q_2",
                   "covt14_bodyweight_adu_q_2", "covt15_bodyweight_adu_q_2",
                   "covt16_bodyweight_adu_q_2", "covt17_bodyweight_adu_q_2"]
    subset = df.loc[:, [*weight_cols, "body_length"]]
    subset["body_length"] = (subset["body_length"] / 100)**2
    subset.loc[:, weight_cols] = subset.loc[:, weight_cols].divide(subset["body_length"], axis=0)
    subset = subset.loc[:, weight_cols]
    subset.columns = subset.columns.str.split("_").str[0]
    subset = subset.add_suffix("_bmi")
    subset["bmi_recent"] = get_last_values(subset, subset.columns)
    df = pd.merge(df.loc[:, df.columns.difference(subset.columns)], subset, left_index=True, right_index=True, how="left")
    return df


def add_combined_mini_columns(df, df_quest):
    df_quest = df_quest.iloc[:, :-2]
    A1 = df_quest.loc["hebt u zich tijdens de afgelopen 7 dagen voortdurend somber of depressief gevoeld gedurende het grootste gedeelte van de dag, en dit bijna elke dag?" ,:]
    A2 = df_quest.loc["hebt u tijdens de afgelopen 7 dagen voortdurend het gevoel gehad nergens meer zin in te hebben of geen interesse meer te hebben voor dingen die u normaal wel interesseren?" ,:]
    A31 = df_quest.loc["was uw eetlust merkbaar veranderd, of is uw gewicht toegenomen of afgenomen, zonder dat dit de bedoeling was? (in de afgelopen 7 dagen)" ,:]
    A32 = df_quest.loc["hebt u bijna elke nacht slaapproblemen gehad (moeilijk inslapen, wakker worden tijdens de nacht of te vroeg in de ochtend, of juist teveel slapen)? (in de afgelopen 7 dagen)" ,:]
    A33 = df_quest.loc["praatte of bewoog u trager dan gewoonlijk, of voelde u zich juist rusteloos, gejaagd en kon u moeilijk stil blijven zitten? bijna elke dag? (in de afgelopen 7 dagen)" ,:]
    A34 = df_quest.loc["voelde u zich bijna elke dag waardeloos of schuldig? (in de afgelopen 7 dagen)" ,:]
    A35 = df_quest.loc["kon u zich bijna elke dag moeilijk concentreren of moeilijk beslissingen nemen? (in de afgelopen 7 dagen)" ,:]
    A36 = df_quest.loc["voelde u zich bijna elke dag moe of futloos? (in de afgelopen 14 dagen)" ,:]
    A37 = df_quest.loc["hebt u overwogen zichzelf iets aan te doen, wenste u dat u dood was, of had u zelfmoordgedachten? (in de afgelopen 7 dagen)" ,:]
    questionnaire_selection = A1.index[~A1.isna()]

    questionnaire_selection = questionnaire_selection.intersection(A37.index[~A37.isna()])
    A1 = A1.loc[A1.index.intersection(questionnaire_selection)]
    A2 = A2.loc[A2.index.intersection(questionnaire_selection)]
    A31 = A31.loc[A31.index.intersection(questionnaire_selection)]
    A32 = A32.loc[A32.index.intersection(questionnaire_selection)]
    A33 = A33.loc[A33.index.intersection(questionnaire_selection)]
    A34 = A34.loc[A34.index.intersection(questionnaire_selection)]
    A35 = A35.loc[A35.index.intersection(questionnaire_selection)]
    A36 = A36.loc[A36.index.intersection(questionnaire_selection)]
    A37 = A37.loc[A37.index.intersection(questionnaire_selection)]
    A3 = [A31, A32, A33, A34, A35, A37]

    def process_participant(row, A1, A2, A3):
        mini_output_per_week = {}
        for week_id in A1.index:
            mini_outcome = False
            base_score = 0
            mini_score = 0
            A3_ness_answers = 0
            A1_id = A1.loc[week_id]
            A2_id = A2.loc[week_id]
            if pd.isna(row.loc[A1_id]) and pd.isna(row.loc[A2_id]):
                mini_outcome = np.nan
            elif np.float(row.loc[A1_id]) == 2.0 and np.float(row.loc[A2_id]) == 2.0:
                continue
            else:
                base_score = np.int(np.float(row.loc[A1_id]) == 1.0) + np.int(np.float(row.loc[A2_id]) == 1.0)
                if base_score > 0:
                    A3_ness_answers = 3 if base_score == 2 else 4
                    for A3_col in A3:
                        A3_single_col = A3_col.loc[week_id]
                        if pd.isna(A3_single_col) == False:
                            mini_score += np.int(np.float(row.loc[A3_single_col]) == 1.0)
                    mini_outcome = mini_score >= A3_ness_answers
            total_score = base_score + mini_score

            mini_output_per_week["covt{}_mini_combined_depressive".format(str(week_id).replace(".5", "b").replace(".0", ""))] = mini_outcome
            mini_output_per_week["covt{}_mini_depressive_score".format(str(week_id).replace(".5", "b").replace(".0", ""))] = total_score
        return pd.Series(mini_output_per_week)

    df_mini_cols = df.apply(lambda row: process_participant(row, A1, A2, A3), axis=1)
    df_mini_cols = df_mini_cols.iloc[0:, :]
    summed_data = np.sum(df_mini_cols)
    df = pd.merge(df, df_mini_cols, left_index=True, right_index=True, how='left')
    return df


def add_measure_score(df):
    print("calculate measure score")
    questions = [
        "covt14_compliance_adu_q_1_a",
        "covt14_compliance_adu_q_1_b",
        "covt14_compliance_adu_q_1_c",
        "covt14_compliance_adu_q_1_d",
        "covt14_compliance_adu_q_1_e",
        "covt14_compliance_adu_q_1_f",
        "covt14_compliance_adu_q_1_g",
        "covt14_compliance_adu_q_1_h",
        "covt14_compliance_adu_q_1_i",
        "covt14_compliance_adu_q_1_j",
        "covt14_compliance_adu_q_1_k",
        "covt14_compliance_adu_q_1_l"
    ]

    df.loc[:, "measure_score"] = np.nan
    for quest_id in questions:
        df.loc[:, quest_id] = df.loc[:, quest_id].astype(float)
        outcome = np.logical_and(df.loc[:, quest_id].isin([1.0, 2.0]).astype(int), ~df.loc[:, quest_id].isna())
        df.loc[outcome[outcome].index, "measure_score"] = df.loc[outcome[outcome].index, "measure_score"].where(np.logical_and(~df.loc[outcome[outcome].index, "measure_score"].isna(), df.loc[outcome[outcome].index, "measure_score"] > 0), 0)
        df.loc[:, "measure_score"] = df.loc[:, "measure_score"] + outcome
    return df


def correct_alcohol(df):
    alcohol_cols_skipped = ["covt01_alcohol_adu_q_1", "covt02_alcohol_adu_q_1",
                            "covt03_alcohol_adu_q_2", "covt04_alcohol_adu_q_2",
                            "covt05_alcohol_adu_q_3", "covt06_alcohol_adu_q_3"]
    alcohol_cols = ["covt07_alcohol_adu_q_4", "covt08_alcohol_adu_q_4",
                    "covt09_alcohol_adu_q_4", "covt11_alcohol_adu_q_4",
                    "covt14_alcohol_adu_q_4", "covt15_alcohol_adu_q_4",
                    "covt16_alcohol_adu_q_4", "covt17_alcohol_adu_q_4"]
    df.loc[:, alcohol_cols_skipped] = np.nan
    df.loc[:, alcohol_cols] = df.loc[:, alcohol_cols].astype(float)
    df.loc[:, alcohol_cols] = df.loc[:, alcohol_cols].where(df.loc[:, alcohol_cols] > 0.0, np.nan)
    df.loc[:, alcohol_cols] = df.loc[:, alcohol_cols].where(df.loc[:, alcohol_cols] < 200.0, np.nan)
    return df


def postive_test_cumsum(df):
    positive_test_cols = ["covt01_infection_adu_q_1_a", "covt02_infection_adu_q_1_a", "covt03_infection_adu_q_1_a", "covt04_infection_adu_q_1_a", "covt05_infection_adu_q_1_a", "covt06_infection_adu_q_2_a", "covt07_infection_adu_q_2_a", "covt08_infection_adu_q_2_a", "covt09_infection_adu_q_2_a", "covt10_infection_adu_q_2_a", "covt11_infection_adu_q_2_a", "covt12_infection_adu_q_2_a", "covt13_infection_adu_q_2_a", "covt14_infection_adu_q_2_a", "covt15_infection_adu_q_2_a", "covt15b_coronatest_adu_q_1_a", "covt16_coronatest_adu_q_1_a", "covt16b_coronatest_adu_q_1_a", "covt17_coronatest_adu_q_1_a"]
    df_subset = df.loc[:, positive_test_cols].copy()
    df_subset = df_subset == 1
    df_subset = np.cumsum(df_subset, axis=1) > 0
    df_subset = df_subset.astype(float)
    df_subset.columns = df_subset.columns.str.split("_").str[0]
    df_subset = df_subset.add_suffix("_positive_tested_cumsum")
    df = pd.merge(df, df_subset, left_index=True, right_index=True)
    return df


def add_recent_columns(df):
    # code of questionaire columns used for linear correction model
    gender_col = ["covt01_GENDER", "covt02_GENDER", "covt03_GENDER", "covt04_GENDER", "covt05_GENDER", "covt06_GENDER", "covt07_GENDER", "covt08_GENDER", "covt09_GENDER", "covt10_GENDER", "covt11_GENDER", "covt12_GENDER", "covt13_GENDER", "covt14_GENDER", "covt15_GENDER", "covt15b_GENDER", "covt16_GENDER", "covt16b_GENDER", "covt17_GENDER"]
    age_col = ["covt01_AGE", "covt02_AGE", "covt03_AGE", "covt04_AGE", "covt05_AGE", "covt06_AGE", "covt07_AGE", "covt08_AGE", "covt09_AGE", "covt11_AGE", "covt12_AGE", "covt13_AGE", "covt14_AGE", "covt15_AGE", "covt15b_AGE", "covt16_AGE", "covt16b_AGE", "covt17_AGE"]
    household_col = ["covt02_household_adu_q_1", "covt03_household_adu_q_1", "covt04_household_adu_q_1", "covt05_household_adu_q_1", "covt06_household_adu_q_1", "covt07_household_adu_q_1", "covt08_household_adu_q_1", "covt09_household_adu_q_1", "covt11_household_adu_q_1", "covt12_household_adu_q_1", "covt13_household_adu_q_1", "covt14_household_adu_q_1", "covt15_household_adu_q_1", "covt15b_household_adu_q_1", "covt16_household_adu_q_1", "covt16b_household_adu_q_1", "covt17_household_adu_q_1"]
    household_are_childs_1 = ["covt01_children_fam_q_1", "covt02_children_fam_q_1", "covt06_children_fam_q_1"]
    household_are_childs_2 = ["covt10_children_fam_q_2", "covt14_children_fam_q_2", "covt15_children_fam_q_2"]
    housemates_0_12 = ["covt01_household_adu_q_1_a", "covt02_household_adu_q_1_a", "covt06_household_adu_q_1_a", "covt10_household_adu_q_1_a", "covt13_household_adu_q_1_a", "covt17_household_adu_q_1_a"]
    housemates_12_18 = ["covt01_household_adu_q_1_b", "covt02_household_adu_q_1_b", "covt06_household_adu_q_1_b", "covt10_household_adu_q_1_b", "covt13_household_adu_q_1_b", "covt17_household_adu_q_1_b"]
    df.loc[:, gender_col] = df.loc[:, gender_col].replace({'MALE': 1, 'FEMALE': 0})

    # gender
    df["gender_recent"] = get_last_values(df, gender_col)

    # age
    df.loc[:, age_col] = df.loc[:, age_col].astype(float)
    df["age_recent"] = get_last_values(df, age_col)
    df["age2_recent"] = df["age_recent"]**2

    # household with multiple housemates
    df["household_recent"] = get_last_values(df, household_col) == 1
    df["household_recent"] = df["household_recent"].astype(float)

    # household with childs at home
    df["household_are_childs_1_recent"] = get_last_values(df, household_are_childs_1)
    df["household_are_childs_2_recent"] = get_last_values(df, household_are_childs_2)
    df["housemates_0_12_recent"] = get_last_values(df, housemates_0_12)
    df["housemates_12_18_recent"] = get_last_values(df, housemates_12_18)
    have_childs = np.logical_or(df.loc[:, "household_are_childs_1_recent"].isin([1, 3]), df.loc[:, "household_are_childs_2_recent"] == 1)
    have_housemates_below_18 = np.sum(df.loc[:, ["housemates_0_12_recent", "housemates_12_18_recent"]], axis=1) > 0
    df["have_childs_at_home_recent"] = np.logical_and(have_childs, have_housemates_below_18)
    df["have_childs_at_home_recent"] = df["have_childs_at_home_recent"].astype(float)

    # chronic disease
    df["chronic_recent"] = np.logical_or(get_last_values(df, ["covt01_chronic_adu_q_1", "covt02_chronic_adu_q_1"]) == 1, False)
    df["chronic_recent"] = df["chronic_recent"].astype(float)

    # remove all rows which do not contain data for age and gender
    df = df.loc[~df["gender_recent"].isna(), :]
    df = df.loc[~df["age_recent"].isna(), :]

    return df


def merge_new_format(dir_path, corona_test_outcome):
    print("start merging new format")
    df_total = pd.DataFrame()
    quest_numbers = []
    for path in glob.glob(dir_path):
        # process name from the filename
        name = os.path.basename(path).split(".")[0].split("_")[2]
        quest_number = float(name[1:].replace("b", ".5"))
        if quest_number <= 17:
            quest_numbers.append(quest_number)
            total_name = "cov{}".format(name)
            print("process questionaire:", name)

            # read questionaire
            df = read_df(path)

            # print number of positive cases for control
            print("Positive cases: ",
                  np.sum(df.loc[:, corona_test_outcome[name]] == 1))

            # replace column names without prefix
            replace_table = {}
            for col in df.columns[~df.columns.str.startswith("covt")]:
                replace_table[col] = "{}_{}".format(total_name, col)
            df = df.rename(columns=replace_table)

            # merge dataframe
            if df_total.shape[0] > 0:
                df_total = pd.merge(df_total, df,
                                    left_index=True,
                                    right_index=True,
                                    how="outer")
            else:
                df_total = df
    # export the total dataframe
    return df_total

def merge_quest_17(df_total):
    ques_17_path = "/groups/umcg-lifelines/tmp01/projects/ov20_0554/data/raw/covid_questionnaires/week17/rl02/covid19-week17-2-values.dat"
    df_week_17 = pd.read_csv(ques_17_path, sep='\t', dtype=str)
    df_week_17 = df_week_17.replace({"8888": np.nan, "9999": np.nan})

    # convert columns to float
    for column in df_week_17.columns.difference(["DATE_QUESTIONNAIRE", "COVID172TXT", "COVID177TXT", "COVID192A", "COVID192B", "PSEUDOIDEXT"]):
        df_week_17[column] = df_week_17[column].astype(float)

    # merge ids
    project_id_to_pseudo_id_path = "/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v1/phenotype_linkage_file_project_pseudo_id.txt"
    project_id_to_pseudo_id = pd.read_csv(project_id_to_pseudo_id_path, sep="\t", dtype=str)
    df_week_17 = pd.merge(df_week_17, project_id_to_pseudo_id, on="PSEUDOIDEXT", how="left")
    df_week_17 = df_week_17.rename(columns={"AGE_QUESTIONNAIRE": "AGE", "DATE_QUESTIONNAIRE": "DATE"})
    df_week_17 = df_week_17.set_index("PROJECT_PSEUDO_ID").add_prefix("covt17_")
    df_total = pd.merge(df_total, df_week_17, left_index=True, right_index=True, how="outer")
    return df_total


def export_df(df, output_path, suffix):
    numbers = df.columns[df.columns.str.startswith("covt")].str.split("_").str[0].str.replace("covt", "").str.replace("b", ".5").astype(float)
    name = "covid_export_questionaire_{min}-{max}_{suffix}_{date}".format(
        min=str(int(min(numbers))).replace(".5", "b"),
        max=str(int(max(numbers))).replace(".5", "b"),
        suffix=suffix,
        date=datetime.now().strftime("%d-%m-%Y")
    )
    df.to_pickle(os.path.join(output_path, "{}.pkl".format(name)))
    df.to_csv(os.path.join(output_path, "{}.txt".format(name)), sep="\t")
    print("export files created: {}.pkl/txt".format(name))


def read_df(path):
    df = pd.read_csv(path, dtype=str, sep=",", quotechar='"', encoding='utf8', index_col=False, engine='python')
    df = df.set_index("PROJECT_PSEUDO_ID")
    df = df.replace("$7", np.nan)

    column_selection = np.logical_and(df.columns.str.startswith('covt'), ~df.columns.str.contains('responsedate'))
    for column in df.columns[column_selection]:
        print(column)
        df[column] = df[column].astype(float)
    print(df)
    print(df.dtypes)
    return df


def get_last_values(df, subset):
    df_subset = df.loc[:, subset].ffill(axis=1)
    return df_subset.loc[:, subset[-1]]


def create_dir(dir_name):
    if os.path.isdir(dir_name) == False:
        os.mkdir(dir_name)


if __name__ == '__main__':
    main()
