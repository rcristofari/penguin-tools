#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

float view_methx(const char *ifile_path, const int n_samples, const int includes, const int excludes, const char keep_string) {

    FILE *ifile;
    //char keep_string_local[1024]; // unlikely to be exceeded
    int filesize, n_rows, n_items = 0;
    int array_size = 1000;
    int pos;
    char flag;
    int n_samples_kept;
    int sample_order[512];
    int i, m, index = 0;
    float *data, *new_data; // pointer to the data array
    char *sp;

    ifile = fopen(ifile_path, "r");

    n_samples_kept = (strlen(keep_string)+1)/2;

    float meth[n_samples];
    if (n_samples == 0){
        printf("You need to define the number of samples in the matrix.");
        exit(0);
    }

    // Check that file was correctly opened
    if (ifile == NULL){
<<<<<<< Updated upstream
        printf("Could not open the specified file.\nTerminating.\n");
=======
        printf("Could not open the specified file.\n");
>>>>>>> Stashed changes
        exit(0);
    }

    sp = strtok(keep_string, ",");
    sample_order[0] = atoi(sp);
    for (i = 1; i < n_samples_kept; i++) {
        sp = strtok(NULL, ",");
        printf("%d\n", sp);
        sample_order[i] = atoi(sp);
    }

    // Find out he number of lines in the file
    fseek(ifile, 0, SEEK_END);
    filesize = ftell(ifile);
    fseek(ifile, 0, SEEK_SET);
    n_rows = filesize / (sizeof(int) + sizeof(char) + sizeof(float) * n_samples);

    // Allocate a starting amount of memory to store methylation data
    data = (float*)malloc(array_size * sizeof(float) * n_samples_kept);

    // Start iterating through the methylation matrix:
    for (i = 0; i < n_rows; i++){

        fread(&pos, sizeof(int), 1, ifile);
        fread(&flag, sizeof(char), 1, ifile);
        fread(&meth, sizeof(float), n_samples, ifile);

        if (((flag & includes) == includes) && ((flag & excludes) == 0)) {
            for (m = 0; m < n_samples_kept; m++){
                *(data+index) = meth[sample_order[m]];
                printf("%f ", meth[sample_order[m]]);
                index++;
            }
            printf("\n");
            n_items += 1;
        }

        if (n_items >= array_size) {
            array_size *= 2;
            printf("Array new size: %d\n", array_size);
            new_data = (float*)realloc(data, array_size*sizeof(float)*n_samples_kept);
            if (new_data == NULL) {
                printf("Memory could not be allocated");
            } else {
                data = new_data;
            }
        }

    }

return(*data);

}