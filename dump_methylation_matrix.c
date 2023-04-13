#include <stdio.h>

int main(int argc, char *argv[]){
    FILE *ifile = fopen(argv[1], "r");
    int n_samples = atoi(argv[2]);
    int filesize;
    int pos;
    char flag;
    float meth[n_samples];

    int i;
    int j;

    fseek(ifile, 0, SEEK_END);
    filesize = ftell(ifile);
    fseek(ifile, 0, SEEK_SET);

    int n_rows = filesize / (sizeof(int) + sizeof(char) + sizeof(float) * n_samples);

    for(i = 0; i < n_rows; i++) {
        fread(&pos, sizeof(int), 1, ifile);
        fread(&flag, sizeof(char), 1, ifile);
        fread(&meth, sizeof(float), n_samples, ifile);
        printf("%d\t", pos);
        printf("%d\t", flag);
        for(j = 0; j < n_samples; j++){
            printf("%f\t", meth[j]);
        }
        printf("\n");
    }

}