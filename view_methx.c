#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>

int main(int argc, char *argv[]) {

    int opt; // options in the argument parser

    FILE *ifile;
    int includes = 0, excludes = 0, n_samples=0;
    char keep_string[1024]; // unlikely to be exceeded
    int filesize, n_rows, n_items = 0;
    int array_size = 1000;
    int array_index = 0;
    int pos;
    char flag;
    int n_samples_kept;
    int sample_order[512];
    int i, m, index = 0;
    float *data, *new_data; // pointer to the data array
    bool do_subset = false;
    char *sp;

    // Command line argument parsing
    while ((opt = getopt(argc, argv, "i:s:n:x:k:")) != -1) {
        switch (opt) {
            case 'i':
                ifile = fopen(optarg, "r");
                printf("Input file : %s\n", optarg);
                break;
            case 's':
                n_samples = atoi(optarg);
                printf("Matrix file contains %d samples\n", n_samples);
                break;
            case 'n':
                includes = atoi(optarg);
                printf("Including %d\n", includes);
                break;
            case 'x':
                excludes = atoi(optarg);
                printf("Excluding %d\n", excludes);
                break;
            case 'k':
                strcpy(keep_string, optarg);
                do_subset = true;
                n_samples_kept = (strlen(optarg)+1)/2;
                printf("Number of individuals to retain: %d\n", (strlen(optarg)+1)/2);
                printf("Keeping individuals %s\n", keep_string);
                break;
            default:
                fprintf(stderr, "Usage: %s -i input_matrix -s number_of_samples [-n include_flag] [-x exclude_flag]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    float meth[n_samples];
    if (n_samples == 0){
        printf("You need to define the number of samples in the matrix.");
        exit(0);
    }

    // Check that file was correctly opened
    if (ifile == NULL){
        printf("Could not open the specified file.\nTerminating.\n");
        exit(0);
    }


    // Handling the subsetting: if -k is set, parse, otherwise make a list for all samples
    if (do_subset == true){
        sp = strtok(keep_string, ",");
        sample_order[0] = atoi(sp);
        for (i = 1; i < n_samples_kept; i++) {
            sp = strtok(NULL, ",");
            printf("%d\n", sp);
            sample_order[i] = atoi(sp);
        }
    } else {
        n_samples_kept = n_samples;
        for (i=0; i<n_samples; i++){
            printf("%d ", i);
            sample_order[i] = i;
        }
        printf("\n");
    }

    // Find out he number of lines in the file
    fseek(ifile, 0, SEEK_END);
    filesize = ftell(ifile);
    fseek(ifile, 0, SEEK_SET);
    n_rows = filesize / (sizeof(int) + sizeof(char) + sizeof(float) * n_samples);

    printf("Number of records in the input file: %d\n", n_rows);

    data = (float*)malloc(array_size * sizeof(float) * n_samples_kept);

    for (i = 0; i < n_rows; i++){

        fread(&pos, sizeof(int), 1, ifile);
        fread(&flag, sizeof(char), 1, ifile);
        //printf("Iteration %d - position: %d, flag: %d\n", i, pos, flag);
        fread(&meth, sizeof(float), n_samples, ifile);
        //printf("%f\n", meth[0]);

        if (((flag & includes) == includes) && ((flag & excludes) == 0)) {
            //printf("%d\t", index);
            for (m = 0; m < n_samples_kept; m++){
                *(data+index) = meth[sample_order[m]];
                printf("%f ", meth[sample_order[m]]);
                index++;
            }
            printf("\n");
            n_items += 1;
            //printf("Here - %d\n", n_items);
        }

        if (n_items >= array_size) {
            //printf("----------Here----------\n");
            array_size *= 2;
            printf("Array new size: %d\n", array_size);
            new_data = (float*)realloc(data, array_size*sizeof(float)*n_samples_kept);
            //printf("There\n");
            if (new_data == NULL) {
                printf("Memory could not be allocated");
            } else {
                data = new_data;
            }
        }

    }
    printf("Matching rows: %d out of %d\n", n_items, n_rows);

free(data);

}