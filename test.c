#include <stdio.h>

int matching_sites(const int includes, const int excludes, const int n_samples){
    FILE *ifile;
    int filesize, n_rows, n_items = 0;
    int pos;
    char flag;
    float meth[n_samples];
    int i;
    float *data; // pointer to the data array

    ifile = fopen("/scratch/project_2003907/King/features/methylation_features.matrix", "r");

    fseek(ifile, 0, SEEK_END); filesize = ftell(ifile); fseek(ifile, 0, SEEK_SET);
    n_rows = filesize / (sizeof(int) + sizeof(char) + sizeof(float) * n_samples);

    for (i = 0; i < n_rows; i++){
        fread(&pos, sizeof(int), 1, ifile);
        fread(&flag, sizeof(char), 1, ifile);
        fread(&meth, sizeof(float), n_samples, ifile);
        if (((flag & includes) == includes) && ((flag & excludes) == 0)) {
            n_items += 1;
        }
    }

    return(n_items);
}


int squre(const int i){
    return (i*i);
}