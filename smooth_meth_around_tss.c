#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

/* NOTE: for scaffold boundaries, windows spilling out are simply clipped off (but the framing around TSS remains the same) */

// Min and max functions
int min(int x, int y){
  return y ^ ((x ^ y) & -(x < y));
}

int max(int x, int y){
  return x ^ ((x ^ y) & -(x < y));
}

// Main

int main(int argc, char *argv[]) {

    int opt;

    // I/O files:
    FILE *tssfile; // output from findTSS
    FILE *methfile; // binary methylation matrix
    FILE *ofile; // output file
    char * line = NULL;
    size_t len = 0, read;

    int n_samples;
    int pos_in_matrix = -1;
    int keep_scanning = 1;

    // Sliding window parameters:
    int range = 6000, window_size = 500, window_step = 50;
    int start_pos, end_pos;
    int w_start, w_end, n_pos_in_window=0, n_windows=0;
    int abs_window_center;

    // Loop variables
    int i, j, k;

    // String parsing:
    char *sp;
    int pos, abs_pos;
    char *chr, *strand;

    // Methylation matrix parameters
    int meth_pos = 0, meth_flag;
    float meth_level;
    float *meth_array;

    // Averaging variables
    int sum = 0, n = 0;

    // parse command line
    while ((opt = getopt(argc, argv, "i:m:o:n:r:w:s:")) != -1) {
        switch (opt) {
            case 'i':
                tssfile = fopen(optarg, "r");
                break;
            case 'm':
                methfile = fopen(optarg, "rb");
                break;
            case 'o':
                ofile = fopen(optarg, "w");
                break;
            case 'n':
                n_samples = atoi(optarg);
                break;
            case 'r':
                range = atoi(optarg);
                break;
            case 'w':
                window_size = atoi(optarg);
                break;
            case 's':
                window_step = atoi(optarg);
                break;
            default:
                fprintf(stderr, "Usage: %s -g annotation.gtf -d depth.txt\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // Allocate memory for the sliding window array (maximum possible size)
    meth_array = (float*)malloc(n_samples * window_size * sizeof(float));

    // Initialise the output file
    fprintf(ofile, "#site\twindow\twindow_center\tabs_window_center\tstrand\tcpg_dens");
    for (i = 0; i < n_samples; i++){
        fprintf(ofile, "\tsample_%d", i);
    }
    fprintf(ofile, "\n");

    // Start to loop through the TSS file. At each line, do the averaging on the meth matrix and output the results.
    while ((read = getline(&line, &len, tssfile)) != -1) {

        if (*line != '#'){

            chr = strtok(line, "\t");
            sp = strtok(NULL, "\t");
            pos = atoi(sp);
            sp = strtok(NULL, "\t");
            abs_pos = atoi(sp);
            sp = strtok(NULL, "\t");
            strand = strtok(NULL, "\n");

            printf("Processing TSS on %s at %d (%d)                  \r", chr, pos, abs_pos);

            if (abs_pos < range/2) {
                start_pos = 0;
            } else {
                start_pos = abs_pos - range/2;
            }
            end_pos = abs_pos + range/2;

            fseek(methfile, pos_in_matrix * (sizeof(int) + sizeof(char) + n_samples * sizeof(float)), SEEK_SET);

            // Scan the matrix until we reach positions in the focal range
            while(meth_pos < start_pos) {
                pos_in_matrix++;
                fread(&meth_pos, sizeof(int), 1, methfile);
                fseek(methfile, sizeof(char) + n_samples * sizeof(float), SEEK_CUR);
            }

            // We have reached the focal region and start windowing:
            w_start = start_pos;
            w_end = start_pos + window_size;

            while(w_end <= end_pos) {

                // Reset the file to the start of the focal range
                fseek(methfile, pos_in_matrix * (sizeof(int) + sizeof(char) + n_samples * sizeof(float)), SEEK_SET);
                fread(&meth_pos, sizeof(int), 1, methfile);
                fseek(methfile, pos_in_matrix * (sizeof(int) + sizeof(char) + n_samples * sizeof(float)), SEEK_SET);

                // Skip the positions that are too early in the range (before the current window):
                while(meth_pos < w_start){
                    fread(&meth_pos, sizeof(int), 1, methfile);
                    fseek(methfile, sizeof(char) + n_samples * sizeof(float), SEEK_CUR);
                }

                // When we are in the window area, we register the data:
                while(meth_pos >= w_start && meth_pos < w_end){
                    fread(&meth_pos, sizeof(int), 1, methfile);
                    fread(&meth_flag, sizeof(char), 1, methfile);
                    fread(meth_array + n_pos_in_window * n_samples, sizeof(float), n_samples, methfile);
                    n_pos_in_window++;
                }

                // Provided the window doesn't spill over the beginning of the scaffold, we write out the data

                abs_window_center = (pos - (range / 2)) + n_windows * window_step + (window_size / 2) ;

                if(abs_window_center >= (window_size / 2)) {

                    if(strcmp(strand, "+") == 0){
                        fprintf(ofile, "%s:%d\t%d\t%d\t%d\t%s\t%f", chr, pos, n_windows, abs_window_center, w_start - window_size / 2, strand, (float)n_pos_in_window / window_size);
                    } else {
                        fprintf(ofile, "%s:%d\t%d\t%d\t%d\t%s\t%f", chr, pos, ((range - window_size) / window_step) - n_windows, abs_window_center, w_start - window_size / 2, strand, (float)n_pos_in_window / window_size);
                    }

                    for(i = 0; i < n_samples; i++){
                        n = n_pos_in_window; // to start with, we assume no missing data
                        for(j = 0; j < n_pos_in_window; j++){
                            if(meth_array[i + j * n_samples] <= 1){
                                sum += meth_array[i + j * n_samples];
                            } else {
                                n--;
                            }
                        }
                        fprintf(ofile, "\t%f", (float)sum / n);
                        sum = 0;
                    }
                    fprintf(ofile, "\n");

                }


                w_start += window_step;
                w_end += window_step;
                n_windows++;
                n_pos_in_window = 0;

            }

        meth_pos = 0;
        n_windows = 0;
        pos_in_matrix = max(pos_in_matrix - window_size, 0);
        }
    }

printf("\nFinished.\n");
return 0;

}