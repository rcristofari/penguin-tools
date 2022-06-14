#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

int main(int argc, char *argv[]) {

    int opt; // options in the argument parser
    float upper_quantile, lower_quantile;
    bool do_upper_quantile = false, do_lower_quantile = false, do_auto = false;
    int quantile_index, lower_quantile_id, upper_quantile_id; // indices of the quantiles in the array
    FILE *ifile; // input file pointer
    int nlines = 0; // number of records in the input file
    char c; // characters in the input file
    float *depth; // pointer to the depth float array
    int line_length = 512; // the maximum length of one line
    char line[line_length]; // the line array
    char *sp; // the string pointer to put elements of the line
    int i = 0; // loop counter

    // Command line argument parsing
    while ((opt = getopt(argc, argv, "i:u:l:a")) != -1) {
        switch (opt) {
            case 'i':
                ifile = fopen(optarg, "r");
                printf("Input file : %s\n", optarg);
                break;
            case 'u':
                upper_quantile = atof(optarg);
                do_upper_quantile = true;
                printf("Upper quantile = %f\n", upper_quantile);
                break;
            case 'l':
                lower_quantile = atof(optarg);
                do_lower_quantile = true;
                printf("Lower quantile = %f\n", lower_quantile);
                break;
            case 'a':
                do_auto = true;
                printf("Evaluating standard quantiles\n");
                break;
            default:
                fprintf(stderr, "Usage: %s [-ul] [file...]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }


    // Check that file was correctly opened
    if (ifile == NULL){
        printf("Could not open the specified file.\nTerminating.\n");
        exit(0);
    }

    // Counting lines in the input file:
    for (c = getc(ifile); c != EOF; c = getc(ifile)){
        if(c == '\n'){
            nlines += 1;
        }
    }
    printf("The input file contains %d records.\n", nlines);

    // Allocating memory for the depth values:
    depth = (float*)malloc(nlines*sizeof(float));
    //printf("Address of array: %x\n", &depth);

    // Scanning the file:
    fseek(ifile, 0, SEEK_SET); // go back to the start of the file (we already took it to the end while counting lines)
    for (i = 0; i < nlines; i++){
        fgets(line, line_length, ifile);
        // Extract chromosome name
        sp = strtok(line, "\t");
        // Extract position
        sp = strtok(NULL, "\t");
        // Extract depth
        sp = strtok(NULL, "\t");
        depth[i] = atof(sp);
        // No need to parse var_depth.
    }
    fclose(ifile);

    // The function to compare values
    int cmpfunc (const void * a, const void * b) {
       return ( *(int*)a - *(int*)b );
    }

    // Sort the depth list
    qsort(depth, nlines, sizeof(int), cmpfunc);

    // Get the quantiles:


    if (do_lower_quantile == true) {
        quantile_index = (int)round(lower_quantile * nlines);
        printf("Lower %f percentile at depth = %f X\n", lower_quantile*100, depth[quantile_index]);
    }
    if (do_upper_quantile == true) {
        quantile_index = (int)round(upper_quantile * nlines);
        printf("Upper %f percentile at depth = %f X\n", upper_quantile*100, depth[quantile_index]);
    }
    if (do_auto == true) {
        float levels[8] = {0.001, 0.01, 0.05, 0.1, 0.9, 0.95, 0.99, 0.999};
        printf("PROB\tDEPTH\n");
        for (i = 0; i < 8; i++){
            quantile_index = (int)round(levels[i] * nlines);
            printf("%f\t%3.1f\n", levels[i], depth[quantile_index]);
        }
    }
   return(0);

}