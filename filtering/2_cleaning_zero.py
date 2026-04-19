# -*- coding: utf-8 -*-
"""
Created on Mon May 12 14:51:05 2025

@author: Samantha Crouch
"""


import pandas as pd
import os

# Set working directory
os.chdir("YourFolder")


def clean_and_save(df, filename_out):
    # Delete columns with only zeros (except the first column)
    cols_to_check = df.columns[1:]
    cols_with_only_zeros = [col for col in cols_to_check if (df[col] == 0).all()]
    df = df.drop(columns=cols_with_only_zeros)
    num_cols_deleted = len(cols_with_only_zeros)

    # Delete rows with only zeros (except the first column)
    rows_with_only_zeros = (df[cols_to_check] == 0).all(axis=1)
    num_rows_deleted = rows_with_only_zeros.sum()
    df = df[~rows_with_only_zeros]

    # Save the result
    df.to_csv(filename_out, index=False)

    # Output summary
    print(f"{filename_out}: {num_cols_deleted} column(s) deleted, {num_rows_deleted} row(s) deleted.")


# --- Apply the function to your files ---
df = pd.read_csv("GSM7068361_ETO_day_0_count.csv", sep=',')
df1 = pd.read_csv("GSM7068362_ETO_day_1_count.csv", sep=',')
df2 = pd.read_csv("GSM7068363_ETO_day_2_count.csv", sep=',')
df3 = pd.read_csv("GSM7068364_ETO_day_4_count.csv", sep=',')
df4 = pd.read_csv("GSM7068365_ETO_day_7_count.csv", sep=',')
df5 = pd.read_csv("GSM7068366_ETO_day_10_count.csv", sep=',')
df6 = pd.read_csv("GSM7068354_CTRL_2_count.csv", sep=',')

# Apply cleaning function to the selected dataframes
clean_and_save(df, "GSM7068361_ETO_day_0_count_cleaned.csv")
clean_and_save(df1, "GSM7068362_ETO_day_1_count_cleaned.csv")
clean_and_save(df2, "GSM7068363_ETO_day_2_count_cleaned.csv")
clean_and_save(df3, "GSM7068364_ETO_day_4_count_cleaned.csv")
clean_and_save(df4, "GSM7068365_ETO_day_7_count_cleaned.csv")
clean_and_save(df5, "GSM7068366_ETO_day_10_count_cleaned.csv")
clean_and_save(df6, "GSM7068354_CTRL_2_count_cleaned.csv")
