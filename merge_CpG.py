import gzip, os, argparse, time, datetime

parser = argparse.ArgumentParser(description="Merge coverage and methylation level from both strands from Bismark CpG_report files.")
parser.add_argument('--input', help='Input file path')
args = parser.parse_args()
path_in = args.input

if not os.path.isfile(path_in):
    print("Input file could not be found")

path_out = path_in.split(".")[0] + ".CpG_merged.cov"
path_report = path_in.split(".")[0] + ".CpG_merged.report.txt"
path_strands = path_in.split(".")[0] + ".CpG_both_strands.cov"

print(f"\nInput file: {path_in}")
print(f"Output file: {path_out}")
print(f"Report file: {path_report}")

i = 0
start_time = time.time()

with (gzip.open if path_in.endswith("gz") else open)(path_in, 'rt') as ifile, open(path_out, 'w') as ofile, open(path_strands, 'w') as strandfile:

    prev_data = [None, None, None, None, None, None, None]
    nCpG_plus, nCpG_minus, nCpG_both = 0, 0, 0
    covCpG_plus, covCpG_minus, covCpG_both = 0, 0, 0
    covCpG_plus_in_both, covCpG_minus_in_both = 0, 0
    percent_diff_btw_both = 0

    strandfile.write('plus\tminus\n')

    for line in ifile:
        this_data = line.strip("\n").split("\t")
        try:
            if this_data[0] == prev_data[0] and int(this_data[1]) == int(prev_data[1]) + 1 and this_data[2] == '-' and prev_data[2] == '+':
                total_C = int(prev_data[3]) + int(this_data[3])
                total_T = int(prev_data[4]) + int(this_data[4])
                percent = round( (total_C / (total_C + total_T) ) * 100, 6)
                ofile.write("\t".join([prev_dataq[0], str(prev_data[1]), str(prev_data[1]), str(percent), str(total_C), str(total_T)]) + "\n")

                i += 1
                if i % 1000 == 0:
                    print(f"Processed {i} CpG sites")

                if int(prev_data[3]) + int(prev_data[4]) == 0 and int(this_data[3]) + int(this_data[4]) > 0:
                    nCpG_minus += 1
                    covCpG_minus += (total_C + total_T)

                elif int(prev_data[3]) + int(prev_data[4]) > 0 and int(this_data[3]) + int(this_data[4]) == 0:
                    nCpG_plus += 1
                    covCpG_plus += (total_C + total_T)

                elif int(prev_data[3]) + int(prev_data[4]) > 0 and int(this_data[3]) + int(this_data[4]) > 0:
                    nCpG_both += 1
                    covCpG_both += (total_C + total_T)
                    covCpG_plus_in_both += int(prev_data[3]) + int(prev_data[4])
                    covCpG_minus_in_both += int(this_data[3]) + int(this_data[4])
                    percent_diff_btw_both += abs(covCpG_plus_in_both - covCpG_minus_in_both) / (covCpG_plus_in_both + covCpG_minus_in_both)
                    strandfile.write(f"{int(prev_data[3]) + int(prev_data[4])}\t{int(this_data[3]) + int(this_data[4])}")

            else:
                prev_data = this_data[:]

        except:

            pass

with open(path_report, 'w') as ofile:
    nCpG = [nCpG_plus, nCpG_minus, nCpG_both]
    covCpG = [covCpG_plus, covCpG_minus, covCpG_both]
    cov = [round(covCpG[i]/nCpG[i],4) for i, x in enumerate(nCpG)]
    ofile.write("# Count and mean coverage per CpG site, by strand:\n\n")
    ofile.write("\tplus\tminus\tboth\n")
    ofile.write(f"n_CpG\t{nCpG_plus}\t{nCpG_minus}\t{nCpG_both}\n")
    ofile.write(f"cov\t{cov[0]}\t{cov[1]}\t{cov[2]}\n")
    ofile.write("\n# Strand-specific coverage for CpGs sequenced on both strand:\n\n")
    ofile.write(f"plus-strand:\t{round(covCpG_plus_in_both/nCpG_both, 2)} X\n")
    ofile.write(f"minus-strand:\t{round(covCpG_minus_in_both/nCpG_both, 2)} X\n")

    ofile.write(f"\nAverage variation between strands:\t{round(percent_diff_btw_both/nCpG_both, 2)}\n")

    ofile.write(f"\n# {round((nCpG_both / sum(nCpG))*100, 2)}% of all CpG sites are sequenced on both strands")

print(f"Completed in {str(datetime.timedelta(seconds=round(time.time() - start_time)))}")

