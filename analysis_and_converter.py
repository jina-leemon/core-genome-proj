import argparse
import re
import pandas as pd
import numpy as np

# takes the csv file and converts it to a dataframe
    # result: csv file into 2D array

def csv_to_df(result):
    output_array = np.array(result)
    data = [d[1:] for d in output_array[1:]]
    output_df = pd.DataFrame(data, columns= output_array[0][1:], index=[item[0] for item in output_array[1:]])
    return output_df

'''analysis'''
# analysis; reports the proteins that have homologies in more than threshold percent of all the phages in list
    # output_df: dataframe from input
    # threshold: % conserved

def report_high_p_prot(output_df, threshold, list=False):
    proteins = []
    protnum = len(output_df.columns)
    for pi in range(protnum):
        p = output_df.loc[str(len(output_df.index)-1)][pi]
        if float(p) >= threshold:
            proteins.append(output_df.columns[pi])
    
    proteins = proteins[:len(proteins)-1]
    print("\nthe number of proteins that are present in more than " + str(threshold) + " % of the phages are: " + str(len(proteins)))
    print("\nthese are the proteins:")
    if list:
        for p in proteins:
            print(p)
    else:
        print(proteins)
    return proteins

# Sort through the proteins of a certain phage and reports the ones that are above the threshold 
    # p_name: name of the phage
    # threshold: the threshold for reporting the protein
    # output_df: input (dataframe with the result from extract_protein_homologs.py)
    # gp_names: the array of strings that denote the signature for gp number

def report_specific_phage(p_name, threshold, output_df):
    proteins = []
    protnum = len(output_df.columns)

    # you may have to add an identifier to this depending on which phage you're looking for
    gp_names = ["gp_"]

    for pi in range(protnum):
        # p is the percent of phages that have a homolog to this ray protein
        p = output_df.loc[str(len(output_df.index)-1)][pi]

        if float(p) >= threshold:
            name = output_df.loc[p_name].iloc[pi]
            print(name)

            gps, gPs, match_gp = [], [], []

            for gp_name in gp_names:
                cmd = r"(?<=" + re.escape(gp_name) + r")(\d+)"
                gps = re.findall(cmd, name, re.IGNORECASE)
                gPs += gps

            for gP in gPs: match_gp.append("gp_" + gP)
            proteins.append(match_gp)
    
    proteins = proteins[:len(proteins)-1]
    print('\nFor ' + p_name + ', the proteins that are in the core genome (threshold: ' + str(threshold) + '% of species)')
    print('\nthese are the proteins:')
    print(proteins)

    return proteins 


'''conversions'''
# Creates a dictionary with the accession and the protein
    # file_name: the file with the protein id (acession) as key and the corresponding gene name

def acession_gp_dict(file_name):
    phage = open(file_name)
    phage_genome = phage.readlines()

    # This creates a dictionary that converts the accession numbers to the gene name
    # Example phiPA3 name:
    # >lcl|NC_028999.1_prot_YP_009217216.1_134 [gene=137] [locus_tag=AVT69_gp135] [db_xref=GeneID:26643665] [protein=hypothetical protein] [protein_id=YP_009217216.1] [location=122733..123266] [gbkey=CDS]
    conversion_list = {}
    print(file_name)
    if "PhiPA3" in file_name or "phiPA3" in file_name:
        for row in phage_genome:
            if ">" in row:
                locus_tag = re.findall(r'(?<=\[gene=).+?(?=\])', row)[0]
                accession = re.search(r"(?<=protein_id=)(\w+)", row)[0]
                conversion_list[accession] = locus_tag
    
    elif "PhiKZ" in file_name or "phiKZ" in file_name:
        for row in phage_genome:
            if ">" in row:
                locus_tag = re.findall(r'(?<=\[protein=PHIKZ).+?(?=\])', row)[0]
                
                accession = re.search(r"(?<=protein_id=)(\w+)", row)[0]
                conversion_list[accession] = locus_tag
    
    #[protein_id=YP_009820689.1], [locus_tag=HOV27_gp004]
    elif "Goslar" in file_name or "goslar" in file_name:
        for row in phage_genome:
            if ">" in row:
                locus_tag = re.findall(r'(?<=\[locus_tag=HOV27_gp).+?(?=\])', row)[0]
                
                accession = re.search(r"(?<=protein_id=)(\w+)", row)[0]
                conversion_list[accession] = locus_tag
    
    else:
        gp_no = 0
        for row in phage_genome:
            if ">" in row:
                gp_no += 1
                locus_tag = gp_no
                
                accession = re.search(r"(?<=protein_id=)(\w+)", row)[0]
                conversion_list[accession] = locus_tag
        
    #assert conversion_list["YP_009217089"] == '007'
    print(conversion_list)
    return conversion_list


# Inserts a gp number (depending on the desired annotation) in the gene name section
# Outputs datatframe with changed names
    # phage_df: the dataframe with the desired phage
    # conversion_list: the dictionary from the conversion format that fits the phage

def converter(phage_df, conversion_list):
    
    for hi in range(len(phage_df)):
        gene = phage_df[hi]
        # if there no gene (np.NaN)
        if gene != gene:
            continue

        gene_txt = gene.split(" ")
        for gi in range(len(gene_txt)):
            ac_des = gene_txt[gi]

            # takes acession and inserts gp number
            if "[" in ac_des:
                accession = ac_des[1:len(ac_des)-1]
                print(accession)
                if accession not in conversion_list:
                    continue
                print("!")
                gp_no = conversion_list[accession]
                if "gp_" + gp_no not in gene_txt[gi-1]:
                    gene_txt.insert(gi, "gp_" + gp_no)

        new_txt = ""
        for g in gene_txt: new_txt += g + " "
        phage_df[hi] = new_txt
        print(new_txt)

    return phage_df   


'''
for Goslar analysis:
    python analysis_and_converter.py homology_chart_ray_psiblast_result_fin_df.csv -t 85 -g "Escherichia phage vB_EcoM_Goslar"
'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser("This code performs some analysis")
    parser.add_argument("output", help = "The name of the file with the output from extract_prot_homologs")
    parser.add_argument("-g", "--genome", help = "The name of the file with the phage genome (for conversions)")
    parser.add_argument("-c", "--convert", help = "The name of the phage (for conversions)")
    parser.add_argument("-n", "--name", help = "The name of the phage; how it would show up on the BLAST output (to output homologies from specific phage); e.g. Pseudomonas Phage PhiPA3")
    parser.add_argument("-t", "--threshold", type=float, help = "% of the phages that have a homolog to the protein")
    args = parser.parse_args()
    
    result = open(args.output)
    result = result.readlines()
    result = [a.strip().split(",") for a in result]

    output_df = csv_to_df(result)
    # This is the converter where the gp number is inserted depending on the acession number
    if args.convert != None:
        for ri in range(len(result)):
            phage_result = result[ri]
            if phage_result[0] == args.convert:
                if args.genome != None:
                    conversion_list = acession_gp_dict(args.genome)
                    result[ri][1:] = converter(phage_result[1:], conversion_list)

        # This converts the already made dataframe(input) to the dataframe
        output_df = csv_to_df(result)
        phage_name = args.convert.split(" ")[::-1][0]
        name = phage_name
        output_df.to_csv("formatted_" + name + ".csv")
            


    if args.convert != None and args.threshold != None:
        report_high_p_prot(output_df, 70)
        report_high_p_prot(output_df, 90)
        report_high_p_prot(output_df, 100)
        report_high_p_prot(output_df, args.threshold, True)

    # This part sorts through the proteins of a certain phage and reports the ones that are above the threshold 
    if args.name != None:
        pseudomonas = report_specific_phage(args.name, args.threshold, output_df)
        for pd in pseudomonas:
            text = ""
            for p in pd:
                text += " " + p
            print(text)
            
    