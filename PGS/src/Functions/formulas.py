
import pandas
import numpy
import math
import matplotlib.pyplot as plt

def define_alleles(*arg):

    """ 
    """

    if math.isclose(sum(arg), 1) is False:
        raise ValueError('Add up is not one')
    
    Alleles = []
    for i in range(len(arg)):
        allele_prefix = chr(65)
        al = f'{allele_prefix}{i+1}'
        Alleles.append(al)
    
    frequencies = [item for item in arg]    
    d = {'Alleles': Alleles, 'Freq': frequencies }
    df = pandas.DataFrame(data=d)
    
    return df,frequencies, Alleles

def expected_formulas(frequencies, allele):

    n = len(frequencies)
    gen_num = int(n*(n+1)/2)

    genotypes = []
    gen_freqs = []
    for i in range(gen_num):
        for j in range(i, n):

            gen = f'{allele[i]}/{allele[j]}'

            if allele[i] == allele[j]:
                freq = frequencies[i]*frequencies[i]
            else:
                freq = 2*frequencies[i]*frequencies[j]
            
            genotypes.append(gen)
            gen_freqs.append(freq)   

    d = {'Genotypes': genotypes,
          'Freq': gen_freqs }
    df = pandas.DataFrame(data=d)

    return df

def lab_formulas(frequencies, allele):

    n = len(frequencies)
    gen_num = int(n*(n+1)/2)

    genotypes = []
    gen_freqs = []
    for i in range(gen_num):
        for j in range(i, n):

            gen = f'{allele[i]}/{allele[j]}'

            if allele[i] == f'A{n}' and allele[i] == allele[j]:
                freq = 0
                lost_sample = frequencies[i]*frequencies[i]
            elif allele[j] == f'A{n}':
                added_freq = 2*frequencies[i]*frequencies[j]
                freq = 0
                if f'A{j}/A{j}' in genotypes:
                    idx = genotypes.index(f'A{j}/A{j}')
                    gen_freqs[idx] += added_freq
            elif allele[i] == allele[j]:
                freq = frequencies[i]*frequencies[i]
            else:
                freq = 2*frequencies[i]*frequencies[j]
            
            genotypes.append(gen)
            gen_freqs.append(freq)   

    d = {'Genotypes': genotypes,
          'Freq_Null': gen_freqs }
    
    df = pandas.DataFrame(data=d)

    return df, lost_sample

def visualization_formulas(df_gen, df_null):

    merge_df = pandas.merge(df_gen, df_null, how='outer', on='Genotypes', suffixes=('_Real', '_Lab'))
    x = merge_df[merge_df.columns[0]]
    y1 = merge_df[merge_df.columns[1]]
    y2 = merge_df[merge_df.columns[2]]

    dist = y1 - y2
    merge_df['Difference'] = dist


    w = 0.4
    bar1 = numpy.arange(len(x))
    bar2 = [i+w for i in bar1]

    plt.rcParams['font.family'] = 'DejaVu Serif'
    plt.rcParams['font.serif'] = ['Times New Roman']

    plt.figure(figsize=(20, 13))
    plt.bar(bar1, y1, w, label='Expected', color='tomato', edgecolor='black')
    plt.bar(bar2, y2, w, label='Observed', color='lightblue', edgecolor='black')

    plt.title('Real and Lab Genotype Results', fontsize=20)
    plt.xticks(bar1, x, rotation=45, fontsize=15)
    plt.xlabel('Genotypes', fontsize=15)
    plt.ylabel('Frequencies', fontsize=15)
    plt.legend(fontsize=20)
    plt.show()

    return merge_df