"""
Heidi Eren
2/12/23
pka-potentials.py: A reuseable library for creating a table of top ten couplings
"""

import pandas as pd
import numpy as np
import re
import os


def read_mom_file(filename):
    """
    read in mom file, convert to string dtype, add new column with only amino acid
    :param file: str, name of file
    :return: df
    """
    mom_df = pd.read_csv(filename, header=None, names=['combined'])
    mom_df.index = np.arange(1, len(mom_df) + 1)
    mom_df = mom_df.astype('string')
    mom_df['amino'] = mom_df['combined'].str[:7]

    return mom_df

def read_pka_file(filename):
    """
    read in pkaS-potentials file, convert dtype
    :param filename: str, name of file
    :return: df
    """
    pka_df = pd.read_csv(filename, header=None)
    pka_df = pka_df[0].astype('string')
    return pka_df

def clean_pka(df):
    """
    clean pka df by removing whitespace and separating into columns
    :param df: df
    :return: new df
    """
    # remove whitespace and add cleaned line to list
    pka_lst = []
    for index, row in df.items():
        if len(row) == 29:
            new = re.sub("\s+", ",", row.strip())
            pka_lst.append(new)
        else:
            df = df.drop(index)

    # create new df, separating each int by column
    new_df = pd.DataFrame(pka_lst, columns=['pka'])
    new_df[['a', 'b', 'c', 'd']] = new_df.pka.str.split(",", expand=True)

    # convert to float data type
    feature = ['a', 'b', 'c', 'd']
    for feat in feature:
        new_df[feat] = new_df[feat].astype(float)

    # convert scientific notation to float
    new_df[["c"]].style.format("{:.0f}")

    # add new colum that calculates pKa value
    new_df['pKa'] = new_df['a'] - new_df['b'] * new_df['c']/1.34

    return new_df

def read_excel(filename):
    """
    read in excel file, data from r-script
    :param filename: str, name of file
    :return: df
    """
    sample_df = pd.read_excel(filename, header=None, names=['numb', 'e'])

    """if len(sample_df) > 10:
        sample_df1 = sample_df.iloc[:10]
        sample_df2 = sample_df.iloc[10:]

        return sample_df1, sample_df2
    else:
        return sample_df"""

    return sample_df

def combine_df(mom_df, pka_df, sample_df):
    """
    combine df to make table with amino acid, potential energy, pka
    :return: df
    """
    final = []
    for i in range(len(sample_df)):
        value = sample_df.loc[i, 'numb']
        data = {'residue': mom_df.iloc[value - 1]['amino'],
                'e': round(sample_df.loc[i, 'e'], 2),
                'pka': round(pka_df.iloc[value - 1]['pKa'], 1),
                'charge': pka_df.iloc[value-1]['b']}
        final.append(data)

    df = pd.DataFrame(final)
    return df



def main():

    # input protein name
    protein = '1b57'

    # read in mom file
    mom_file = open(os.path.expanduser(f'~/Downloads/POOL/{protein}_2/{protein}.mom'))

    mom_df = read_mom_file(mom_file)

    # read in pka file
    pka_file = open(os.path.expanduser(f'~/Downloads/POOL/{protein}_2/pkaS-potentials'))

    pka_df = read_pka_file(pka_file)


    final_pka_df = clean_pka(pka_df)

    # read in excel file
    sample_df = read_excel('sample.xlsx')

    # consolidate into df table
    final_df = combine_df(mom_df, final_pka_df, sample_df)

    print(final_df)



    # Write the dataframe data to XlsxWriter. Turn off the default header and
    # index and skip one row to allow us to insert a user defined header.
    final_df.to_excel(f'{protein}_table.xlsx')



if __name__ == '__main__':
    main()


