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

    int n_rows = filesize / (4*(n_samples+1)+1);

    for(i = 0; i < n_rows; i++) {
        fread(&pos, 4, 1, ifile);
        fread(&flag, 1, 1, ifile);
        fread(&meth, 4, 16, ifile);
        printf("%d\t", pos);
        printf("%d\t", flag);
        for(j = 0; j < 16; j++){
            printf("%f\t", meth[j]);
        }
        printf("\n");
    }

}