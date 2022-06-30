import gzip
genome_file = "/scratch/project_2003907/King/ref/GCA_010087175.1_BGI_Apat.V1_genomic.fna.gz"
position = 0
with gzip.open(genome_file, 'rt') as ifile:
    for line in ifile:
        if line.startswith(">"):
            pass
        else:
            for k, base in enumerate(line[:-1]):
                position += 1

                if position == 183763:
                    print("---------------")
                    print(position)
                    try:
                        print(line[k-5:k].lower() + line[k].upper() + line[k+1:k+5].lower())
                    except:
                        try:
                            print(line[k-5:k].lower() + line[k].upper() + line[k+1:].lower())
                        except:
                            print(line[:k].lower() + line[k].upper() + line[k+1:k+5].lower())
                    print(next(ifile))
                    quit()

        if position % 10000 == 0:
            print(position)