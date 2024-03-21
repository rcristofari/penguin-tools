import argparse

parser = argparse.ArgumentParser(description='Extract protein ID brom a batch blastp run')
parser.add_argument('--blastp', help='Results of the Blastp run')
parser.add_argument('--out', help='Output file')
args = parser.parse_args()

import re, difflib, operator

inside_block = False
data = None
counter = 0

with open(args.blastp, 'r') as ifile, open(args.out, 'w') as ofile:
    for line in ifile:

        if not inside_block and not line.startswith("Query="):
            next(ifile)

        elif line.startswith("Query="):
            inside_block = True
            gene_name = line.strip("\n").split(" ")[1]
            #print(gene_name, end="")
            next(ifile); next(ifile); next(ifile); next(ifile); next(ifile)
            data = []

        elif inside_block and line.startswith(">") or "Lambda" in line:
            inside_block = False
            counter = 0
            substring_counts={}
            if len(data) > 5:
                for i in range(0, len(data)):
                    for j in range(i+1,len(data)):
                        string1 = data[i]
                        string2 = data[j]
                        match = difflib.SequenceMatcher(None, string1, string2, autojunk=True).find_longest_match(0, len(string1), 0, len(string2))
                        matching_substring=string1[match.a:match.a+match.size]
                        if(matching_substring not in substring_counts):
                            substring_counts[matching_substring]=1
                        else:
                            substring_counts[matching_substring]+=1

               # print(substring_counts) #{'myKey_': 5, 'myKey_apples': 1, 'o': 1, '': 3}
                pruned = {k: v for k, v in substring_counts.items()
                          if ((len(k) > 3 and re.match("^[A-Z]+", k)) or (len(k) > 6 and (k.strip(" ") not in ["REDUCTASE", "DEDHYDROGENASE", "FACTOR", "SYNTHASE", "LIGASE", "KINASE", "KERATIN-", "KERATIN"])))}
                try:
                    max_occurring_substring = max(pruned.items(), key=operator.itemgetter(1))
                    consensus = max_occurring_substring
                except ValueError:
                    consensus = ["***** No hits found *****"]
            else:
                consensus = ["***** No hits found *****"]
            #print("\t" + consensus[0])
            print(gene_name + "\t" + consensus[0])

        elif line != "\n" and counter < 20:
            protein_line = re.sub("[ ]{2,}", "\t", line).strip("\n").split("\t")
            if float(protein_line[2]) < 0.001:
                counter += 1
                protein = " ".join(protein_line[0].split(" ")[1:]).split("[")[0].strip(" ").upper()
                remove_str = ["HYPOTHETICAL", "PARTIAL", "PREDICTED", "ISOFORM", "LOW QUALITY", ":", ".", ","]
                for x in remove_str:
                    protein = protein.replace(x, "")
                if "PROTEIN" in protein:
                    if not re.match("PROTEIN [0-9.]+", protein):
                        protein = protein.replace("PROTEIN", "")
                protein = re.sub("[ ]+", " ", protein).strip(" ")
                data.append(protein)


