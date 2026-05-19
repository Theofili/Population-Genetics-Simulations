import pandas
import numpy
import matplotlib.pyplot as plt


def visualization(df_expected, df_observed):

    """Create a plot to compare expected genotype frequencies and observed due to null alleles
    
        Args:

            df_expected     (dataframe)     : Dataframe with the structure of the output of ado_model function

            df_observed     (dataframe)     : Dataframe with the structure of the output of ado_model function 
"""

    df_expected = df_expected[df_expected['Genotypes'].notna() & (df_expected['Genotypes'] != 'NA')].copy()
    df_observed = df_observed[df_observed['Genotypes'].notna() & (df_observed['Genotypes'] != 'NA')].copy()

    df_expected['Frequency'] = pandas.to_numeric(df_expected['Frequency'], errors='coerce')
    df_observed['Frequency'] = pandas.to_numeric(df_observed['Frequency'], errors='coerce')

    merged = pandas.merge(
        df_expected,
        df_observed,
        on='Genotypes',
        how='outer',
        suffixes=('_expected', '_observed')
    )
 
    x_labels = merged['Genotypes']
    y1 = merged['Frequency_expected']
    y2 = merged['Frequency_observed']

    x = numpy.arange(len(x_labels))  # positions for each genotype
    width = 0.40  # width of the bars


    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman"]


    plt.figure(figsize=(20, 13))
    plt.bar(x - width/2, y1, width, color='tomato', label='Real', edgecolor='black')
    plt.bar(x + width/2, y2, width, color='lightblue', label='Lab', edgecolor='black')

    plt.xticks(x, x_labels, rotation=45, fontsize=15)
    plt.xlabel('Genotypes')
    plt.ylabel('Frequency')
    plt.title('Expected (Real) vs Observed (Lab) Genotype Frequencies', fontsize=18)
    plt.legend(fontsize=18)
    plt.show()