#include <stdio.h>

int main(int argc, char *argv[]){
    FILE *ifile = fopen(argv[1], "r");
    int filesize, nrows;
    int prev_pos = -1;
    char flag;
    int i;

    fseek(ifile, 0, SEEK_END); filesize = ftell(ifile); fseek(ifile, 0, SEEK_SET);
    nrows = filesize / sizeof(char);

    for(i=0; i<nrows; i++){
        fread(&flag, sizeof(char), 1, ifile);
        if ((flag & 1) == 1) {
            if ((i - prev_pos)==1){
                printf("%d\n", prev_pos);
            }
            prev_pos = i;
        }
    }
}