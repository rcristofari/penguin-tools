#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv[]) {

    int opt; // options in the argument parser
    FILE *ifile = fopen(argv[1], "r");
    char chrBase[32], chr[32]; // chromosome
    char buff[512];
    int base, coverage; // position and coverage
    float freqC, freqT; // base frequencies
    int nC, nT;
    char strand[8]; // F or R
    char c = '\0';

    // Skip the header line
    while (c != '\n') {
      c = fgetc(ifile);
      //printf("%c", c);
    };

    // parse and transform the file, print to stdout
    while (fscanf(ifile, "%s\t%s\t%i\t%s\t%i\t%f\t%f\n", chrBase, chr, &base, strand, &coverage, &freqC, &freqT) == 7) {
       nC = round((freqC*coverage)/100);
       nT = round((freqT*coverage)/100);
       //printf("%s\t%s\t%i\t%s\t%i\t%02.02f\t%02.02f\n", chrBase, chr, base, strand, coverage, freqC, freqT);
       printf("%s\t%i\t%i\t%02.02f\t%i\t%i\n", chr, base, base, freqC, nC, nT);
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