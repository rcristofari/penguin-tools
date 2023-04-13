#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>

int main(int argc, char *argv[]) {

    int opt; // options in the argument parser

    FILE *ifile, *gfile; // input file, and optional genome coordinate file
    int includes = 0, excludes = 0, n_samples=0;
    char keep_string[1024]; // unlikely to be exceeded
    int filesize = 0, n_rows = 0, n_items = 0;
    int array_size = 1000;
    int array_index = 0;
    int pos = 0; // the current position
    int range_start = 0, range_end = 0; // start and end positions if we extract from a range
    int abs_start = 0;
    char chrom[64]; // the chromosome name if we extract from a range
    char flag = 0;
    int n_samples_kept = 0;
    int sample_order[512];
    int index = 0;
    float *data, *new_data; // pointer to the data array
    bool do_subset = false, do_range = false;
    char *sp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read = 0;

    printf("+--------------------------------------------------------------------------------\n");

    // Command line argument parsing
    while ((opt = getopt(argc, argv, "m:s:n:i:x:k:r:f:")) != -1) {
        switch (opt) {
            case 'm':
                // Input matrix file
                ifile = fopen(optarg, "r");
                printf("+ Input file : %s\n", optarg);
                break;
            case 's':
                // Total number of Samples in the matrix file
                n_samples = atoi(optarg);
                break;
            case 'i':
                // Including types of positions (in binary code format)
                includes = atoi(optarg);
                break;
            case 'x':
                // eXcluding types of positions (in binary code format)
                excludes = atoi(optarg);
                break;
            case 'k':
                // We Keep only a list of samples (comma-separated)
                strcpy(keep_string, optarg);
                do_subset = true;
                n_samples_kept = (strlen(optarg)+1)/2;
                //printf("+ Number of individuals to retain: %d\n", (strlen(optarg)+1)/2);
                printf("+ Keeping %d individuals: %s\n", (strlen(optarg)+1)/2, keep_string);
                break;
            case 'r':
                // We are only exporting data for a particular range. This is provided as a standard string: chrom:start-end.
                do_range = true;
                sp = strtok(optarg, ":-");
                strcpy(chrom, sp);
                sp = strtok(NULL, ":-");
                range_start = atoi(sp);
                sp = strtok(NULL, ":-");
                range_end = atoi(sp);
                break;
            case 'f':
                // In case we want to export a range, we need a genome index file
                gfile = fopen(optarg, "r");
                printf("+ Genome index file: %s\n", optarg);
                break;
            default:
                fprintf(stderr, "Usage: %s -m input matrix -s number of samples [-i include flag] [-x exclude flag] [-k keep list] [-r genomic range] [-f fai index]\n", argv[0]);
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

    // Handling the subsetting option:
    // if -k is set, parse, otherwise make a list for all samples
    if (do_subset == true){
        sp = strtok(keep_string, ",");
        sample_order[0] = atoi(sp);
        for (int i = 1; i < n_samples_kept; i++) {
            sp = strtok(NULL, ",");
            //printf("%d\n", sp);
            sample_order[i] = atoi(sp);
        }
        printf("\n");
    } else {
        n_samples_kept = n_samples;
        for (int i = 0; i<n_samples; i++){
            //printf("%d ", i);
            sample_order[i] = i;
        }
    }

    // Handling of the ranging options:
    // if -r is set, check that the genome index is given and that the range is parsed correctly
    if (do_range == true){
        if (gfile == NULL) {
            printf("Could not open the genome index file.\nTerminating.\n");
            exit(0);
        } else {
            // Parse the index file until the target scaffold is found:

            while ((read = getline(&line, &len, gfile)) != -1) {
                //printf("Retrieved line of length %zu:\n", read);
                //printf("%s", line);
                sp = strtok(line, "\t");
                //printf("First chunk %s\n", sp);
                if (strcmp(sp, chrom) == 0) {
                    //printf("%s == %s\n", chrom, line);
                    range_start += abs_start;
                    range_end += abs_start;
                    //printf("Absolute range: %d to %d\n", range_start, range_end);
                    break;
                }
                sp = strtok(NULL, "\t");
                abs_start += atoi(sp);
            }
            fclose(gfile);

        }
        //printf("Absolute start: %d\n", abs_start); // USE THIS WITH RANGE COORDS TO TEST THE POSITIONS FOR OUTPUTTING. WHEN REACHING END, BREAK.
    } else {
        // If we are not selecting a range, we set the range to the full file
        range_end = 2147483647; // max int
    }

    // Find out he number of lines in the file
    fseek(ifile, 0, SEEK_END);
    filesize = ftell(ifile);
    fseek(ifile, 0, SEEK_SET);
    n_rows = filesize / (sizeof(int) + sizeof(char) + sizeof(float) * n_samples);

    printf("+ Matrix file contains %d samples\n", n_samples);
    if (do_range == true) {
        printf("+ Exporting data for scaffold %s at positions %d to %d\n", chrom, range_start-abs_start, range_end-abs_start);
    } else {
        printf("+ Exporting data for all available positions\n");
    }
    printf("+ Including %d\n", includes);
    printf("+ Excluding %d\n", excludes);
    printf("+ Number of records in the methylation matrix: %d\n", n_rows);
    printf("+--------------------------------------------------------------------------------\n");

    data = (float*)malloc(array_size * sizeof(float) * n_samples_kept);

    for (int i = 0; i < n_rows; i++){

        fread(&pos, sizeof(int), 1, ifile);
        fread(&flag, sizeof(char), 1, ifile);
        //printf("Iteration %d - position: %d, flag: %d\n", i, pos, flag);
        fread(&meth, sizeof(float), n_samples, ifile);
        //printf("%f\n", meth[0]);

        if (((flag & includes) == includes) && ((flag & excludes) == 0) && (pos >= range_start) && (pos <= range_end)) {
            printf("%d ", flag);
            for (int m = 0; m < n_samples_kept; m++){
                *(data+index) = meth[sample_order[m]];
                printf("%f ", meth[sample_order[m]]);
                index++;
            }
            printf("\n");
            n_items += 1;
            //printf("Here - %d\n", n_items);
        } else if (pos > range_end) {
            //printf("Reached final position: %d > %d\n", pos, range_end);
            break;
        }

        if (n_items >= array_size) {
            array_size *= 2;
            //printf("Array new size: %d\n", array_size);
            new_data = (float*)realloc(data, array_size*sizeof(float)*n_samples_kept);
            //printf("There\n");
            if (new_data == NULL) {
                printf("Memory could not be allocated");
            } else {
                data = new_data;
            }
        }

    }
    printf("+ Matching rows: %d out of %d\n", n_items, n_rows);

free(data);

}