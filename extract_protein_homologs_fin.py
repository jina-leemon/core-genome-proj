import re
import argparse
import numpy as np
import pandas as pd
from analysis_and_converter import report_high_p_prot

class Homology:
    #qseqid sacc sscinames stitle pident evalue
    #e.g) FDH98_gp150 gb|QTH80280.1| Pseudomonas phage pPa_SNUABM_DT01 putative UvsX protein [Pseudomonas phage pPa_SNUABM_DT01] 39.407 0.0
    def __init__(self, sacc, sscinames, stitle, pident, evalue):
        #gb|QTH80280.1|
        self.sacc = sacc
        #Pseudomonas phage pPa_SNUABM_DT01
        self.sscinames = sscinames
        #putative UvsX protein [Pseudomonas phage pPa_SNUABM_DT01]
        self.stitle = stitle
        #39.407
        self.pident = pident
        #0.0
        self.evalue = evalue
      
    def aboveThreshold(self, parameter, threshold=0):
        if threshold == 0:
            if parameter <= self.pident:
                return True
            else: return False
        else:
            if parameter >= self.evalue:
                return True
            else: return False

# determines if a string is np.NAN
def isNaN(string):
    return string != string

# makes a list of proteins in the phage genome of interest
    # proteins: file with genome/ protein list in .fasta format

def make_prot_list(proteins):
    prot_no, prot_list = 0, []
    for prot_line in proteins:
        if ">" in prot_line:
            prot_no += 1
            prot_list.append(prot_line[1:])
    #print(len(prot_list))
    return  prot_list

# creates a chart of all homologies from the psiblast output if the phage is in the phage_list
    # proteins: file with genome/ protein list in .fasta format
    # phage_list: list of phages that the homologies are in
    # results: the data from the psiblast output
    # output: dataframe with the phage of interest's genes and the homologies

def make_homology_chart(proteins, phage_list, results):
    print("making chart of all homologies")
    prot_list = make_prot_list(proteins)
    homology_chart = [[[] for p in prot_list] for pl in phage_list]


    for result in results:
        if len(result) == 0 or "CONVERGED!" in result:
            continue

        para = result.split(",")

        pl = len(para)

        stitle = ' '.join(para[i] for i in range(3, pl-2))
        sscinames = re.findall(r'(?<=\[).+?(?=\])', stitle)

        if len(sscinames) != 0:
            sscinames = sscinames[0]
        stitle = re.findall(r'(.*?)\s*\[', stitle)
        if len(stitle) != 0:
            stitle = stitle[0]

        new_homology = Homology(para[1], sscinames, stitle, float(para[pl-2]), float(para[pl-1]))

        curr_prot_i = 0
        curr_prot = prot_list[curr_prot_i]
        next_prot = para[0]
        while curr_prot != next_prot:
            curr_prot_i += 1
            #print(curr_prot)
            curr_prot = prot_list[curr_prot_i]

        if new_homology.sscinames in phage_list:
            homology_chart[phage_list.index(new_homology.sscinames)][curr_prot_i].append(new_homology)

    return homology_chart

# counts the number of phages that have a homologous protein to a certain gene, adds a row to the bottom of dataframe
    # input: the output dataframe

def count_row_column(phage_list, output_df):
    print("generating counts per column")

    output_df['% ray protein'] = output_df.apply(lambda x: x.isnull().sum(), axis = 1)

    # number of proteins that are in the row
    output_df.loc[len(output_df.index)] = len(phage_list) - output_df.isnull().sum() -1

    # percentage with that protein
    output_df.loc[len(output_df.index)] = output_df.loc[len(output_df.index)-1] / (len(phage_list)-1) * 100


# goes through the ouptut dataframe and outputs only the ones above threshold
# the psiblast was ran with a threshold of e value 0.01
    # homology_chart: output from make_homology_chart (output_dataframe)
    # proteins: list of phage proteins of interest
    # phage_list: list of phages to compare with
    # parameter: 1 --> e value 0 --> percent identity
    # threshold: desired value

def generate_output(homology_chart, proteins, phage_list, parameter, threshold):
    print("filtering through for ones above threshold \n")
    prot_list = make_prot_list(proteins)


    array = [[np.NaN for p in prot_list] for pl in phage_list]
    output_array = np.array(array)
    output_df = pd.DataFrame(output_array, columns=prot_list, index=phage_list)

    for phage_i in range(len(phage_list)):
        curr_phage = phage_list[phage_i]

        for prot_i in range(len(prot_list)):
            curr_prot = prot_list[prot_i]
            homology = homology_chart[phage_i][prot_i]
            if len(homology) != 0:
                for h in homology:
                    if h.aboveThreshold(parameter, threshold):
                    
                        if isNaN(output_df.loc[curr_phage, curr_prot]):
                            output_df.loc[curr_phage, curr_prot] = ''

                        list = output_df.at[curr_phage, curr_prot]
                        stitle = h.stitle.replace("hypothetical protein", "")
                        text = " " + stitle + " [" + h.sacc + "]"
                        if stitle == "":
                            if h.sacc not in list: output_df.at[curr_phage, curr_prot] += " [" + h.sacc + "]"
                        else:
                            if text not in list: output_df.at[curr_phage, curr_prot] += text
    
    count_row_column(phage_list, output_df)
    return output_df

def generate_csv(output_df, output_file_name):
    output_df.to_csv(output_file_name)

# analysis; reports the proteins that have homologies in more than threshold percent of all the phages in list
    # output_df: dataframe from input
    # threshold: % conserved

def report_high_p_prot(output_df, threshold):
    proteins = []
    protnum = len(output_df.columns)
    for pi in range(protnum):
        p = output_df.loc[len(output_df.index)-1][pi]
        if float(p) >= threshold:
            proteins.append(output_df.columns[pi])

    proteins = proteins[:len(proteins)-1]
    print("the number of proteins that are present in more than " + str(threshold) + " % of the phages are: " + str(len(proteins)))
    print("\nthese are the proteins:")
    print(proteins)
    return proteins

if __name__ == '__main__':
    # argument:
    # python extract_protein_homologs_fin.py ray_genome_format.fasta phage_list.txt ray_psiblast_result_fin -p 1 -t 0.01

    parser = argparse.ArgumentParser("This code organizes bulk psiblast outputs. \n Example command: python extract_protein_homologs_fin.py ray_genome_format.fasta phage_list.txt ray_psiblast_result_fin -p 1 -t 0.01")
    parser.add_argument("proteins", help = "The name of the file with the proteins of interest")
    parser.add_argument("phage_list", help = "The text file with the list of phages, separated by line change")
    parser.add_argument("output", help = "The output of the batch psiblast")
    parser.add_argument("-p", "--peridentity", type=int, default = 1, help = "If you want to compare e-value, enter 1")
    parser.add_argument("-t", "--threshold", type=float, default = 0.01, help = "The threshold for p identity/ e value; default set to e-value of 0.01")
    args = parser.parse_args()

    output_file_name = 'homology_chart_' + args.output + '.csv'

    prot_file = open(args.proteins)
    phages_file = open(args.phage_list)
    input = open(args.output) 
    pident = args.peridentity
    threshold = args.threshold

    print("starting analysis... may take a bit of time")

    results = input.read().splitlines()
    phage_list = phages_file.read().splitlines()
    proteins = prot_file.read().splitlines()

    homology_chart = make_homology_chart(proteins, phage_list, results)

    # you can also test multiple thresholds at once by copying the following command calling generate_output
    # and chanding the paramteres. This is faster than re-running the program
    # since the homology chart will already be available
    output_df = generate_output(homology_chart, proteins, phage_list, threshold, pident)

    report_high_p_prot(output_df, 90)
    # un-commnet the next line if you want a .csv file for analysis (which will be most of the times)
    print(output_file_name)
    generate_csv(output_df, output_file_name)