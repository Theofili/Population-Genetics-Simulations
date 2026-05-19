
import pandas
import numpy
import matplotlib.pyplot as plt


# Step 5: Visualize Heterozygous frequencies


def visualization_multiple(df):

    """Create a graph to visualize genotype frequencies with changing allele probabilities

    Args:

        df      (dataframe)      :   Dataframe with the structured of the output of data_collection or ado_model_multiple function
    
"""

    p1 = []
    p2 = []
    p3 = []
    Ho = []
    inbr_ind = []

       

    for j in range(0, len(df), 2):

        if not pandas.isna(df.iloc[j, 3]):
            inbr_ind.append(float(df.iloc[j, 3]))

        start = 4


        al_prob = df.iloc[j, start]
        al_prob = al_prob.split(',')

        prob1 = float(al_prob[0])
        p1.append(prob1)
        prob2 = float(al_prob[1])
        p2.append(prob2)
        prob3 = float(al_prob[-1])
        p3.append(prob3)

        number_of_heterozygous = 0

        for i in range(6):
            
            genotype = df.iloc[j, i+start+1]
                        
            if genotype == str('NA'): 
                continue 
            if pandas.isna(genotype): 
                continue 
            if genotype == 'NULL/NULL': 
                continue
            if genotype[1] != genotype[4]: 
                
                nof = df.iloc[j+1, i+5]
                if nof == 'NA':
                    continue
                else:
                    nof = float(df.iloc[j+1, i+5]) * 1000 
                    
                    number_of_heterozygous += nof  
                    

        H_observed = number_of_heterozygous / 1000
        if H_observed > 1:
            continue
        else:
            Ho.append(H_observed)

    
    x = []
    
    for i in range(len(Ho)):   
        xx = f'Run_{i}'
        x.append(xx)


    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman"]


    if len(inbr_ind) != 0:


        data = list(zip(x, Ho, p1, p2, p3, inbr_ind))
        data_sorted = sorted(data, key=lambda item: item[5])
        x_srt, Ho_srt, p1_srt, p2_srt, p3_srt , inbr_ind_srt = zip(*data_sorted)
        
        plt.figure(figsize=(25, 12))
        #plt.plot(x_srt, p1_srt, color='#1B63A6', label='Major Allele Freq', linewidth=2)
        #plt.plot(x_srt, p2_srt, color='#2A9D8F', label='Freq of Allele 2', linewidth=2)
        #plt.plot(x_srt, p3_srt, color='#D62828', label='Freq of NULL Allele', linewidth=3)
        plt.bar(x_srt, Ho_srt, color="#B8B8B8A7", label='Freq of Heterozygous')
        plt.plot(x_srt, inbr_ind_srt, color='#D62828', label='Inbreeding index', linewidth=3)


        plt.ylabel("Frequency", fontsize='20', font='Times New Roman')
        plt.title("Corelation between Number of Heterozygous Individuals and Inbreeding Index", fontsize='30')
        plt.xticks(rotation=45, fontsize='18')
        plt.yticks(fontsize='20')
        plt.legend(fontsize='25')
        plt.tight_layout() 
        plt.show()

        
    
    else:
        

        data = list(zip(x, Ho, p1, p2, p3))
        data_sorted = sorted(data, key=lambda item: item[3])
        x_srt, Ho_srt, p1_srt, p2_srt, p3_srt = zip(*data_sorted)

        plt.rcParams['font.family'] = 'DejaVu Serif'
        plt.rcParams['font.serif'] = ['Times New Roman']

        plt.figure(figsize=(25, 12))
        plt.bar(x_srt, Ho_srt, color="#B8B8B8A7", label='Freq of Heterozygous')
        plt.plot(x_srt, p1_srt, color='#1B63A6', label='Major Allele Freq', linewidth=2)
        plt.plot(x_srt, p2_srt, color='#2A9D8F', label='Freq of Allele 2', linewidth=2)
        plt.plot(x_srt, p3_srt, color='#D62828', label='Freq of NULL Allele', linewidth=3)

        
        plt.ylabel("Frequency", fontsize='20', font='Times New Roman')
        plt.title("Corelation between Number of Heterozygous Individuals and Allele Frequencies", fontsize='30')
        plt.xticks(rotation=45, fontsize='18')
        plt.yticks(fontsize='20')
        plt.legend(fontsize='25')
        plt.tight_layout() 
        plt.show()
