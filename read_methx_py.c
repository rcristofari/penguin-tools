#include <stdio.h>

//float subset_methx(FILE *, const int includes, const int excludes, const int n_samples) {
//
//    const int INIT_SIZE = 10000;
//    float data[INIT_SIZE][n_samples];
//    return data;
//}


int main(){

    FILE *ifile;
    const int includes = 5, excludes = 32, n_samples = 16;
    //float data[][];
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
    printf("Matching rows: %d out of %d\n", n_items, n_rows);

    data = (float*)malloc(n_items*sizeof(float)*n_samples);



}