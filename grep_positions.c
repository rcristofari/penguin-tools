#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>

int main(int argc, char *argv[]) {

    int opt; // options in the argument parser

    FILE *ifile; // input file, and optional genome coordinate file
    char *sp; // string pointer for parsing strings
    char chrom[16]; chrom[0] = '\0'; // the chromosome name
    int range_start=0, range_end=0, abs_start=0;
    char *line = NULL;
    size_t len = 0;
    ssize_t read = 0;
    char this_chrom[16], this_sample[32], this_feature[32], this_cpgi[1]; this_chrom[0] = '\0'; // buffers to read the GLM file lines
    int this_pos, this_fullpos, this_C, this_T;
    bool reached_region = false, keep_searching = true;
    char c;

    // Command line argument parsing
    while ((opt = getopt(argc, argv, "i:r:")) != -1) {
        switch (opt) {
            case 'i':
                // Input GLM file
                ifile = fopen(optarg, "r");
                break;
            case 'r':
                // Range of positions chrom:start-end.
                sp = strtok(optarg, ":-");
                strcpy(chrom, sp);
                sp = strtok(NULL, ":-");
                range_start = atoi(sp);
                sp = strtok(NULL, ":-");
                range_end = atoi(sp);
                break;
            default:
                fprintf(stderr, "Usage: %s -f input matrix -r genomic range\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // Check that input files were correctly opened
    if (ifile == NULL){
        printf("Could not open the GLM data file.\nTerminating.\n");
        exit(0);
    }

// Skip the header line
while (c != '\n') {
  c = fgetc(ifile);
};

while (fscanf(ifile, "%s\t%i\t%i\t%s\t%i\t%i\t%s\t%s\n", this_chrom, &this_pos, &this_fullpos, this_sample, &this_C, &this_T, this_feature, this_cpgi) == 8) {
  if(strcmp(this_chrom, chrom)==0 & this_pos >= range_start & this_pos <= range_end){
      printf("%s\t%i\t%i\t%s\t%i\t%i\t%s\t%s\n", this_chrom, this_pos, this_fullpos, this_sample, this_C, this_T, this_feature, this_cpgi);
      reached_region = true;
  } else if(reached_region == true){
      exit(0);
  }
}
if (feof(ifile))
{
  exit(0);
}
else
{
  // some other error interrupted the read
    exit(0);
}
}