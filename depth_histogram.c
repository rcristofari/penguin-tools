#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

int main(int argc, char *argv[]){

    FILE *ifile = fopen(argv[1], "r");
    int depth[1024] = { 0 };
    char *line = NULL;
    size_t len, read;
    char *sp;
    int this_depth = 0;
    int k=0, i;


    while ((read = getline(&line, &len, ifile)) != -1) {
        k++;
        sp = strtok(line, "\t");
        sp = strtok(NULL, "\t");
        sp = strtok(NULL, "\t");
        this_depth = atoi(sp);
        if (this_depth < 1024){
            depth[this_depth] ++;
        } else {
            depth[1023] ++;
        }
//        if (k % 1000000 == 0){
//            printf("Processed %d positions\r", k);
//        }
    }
//    printf("Processed %d positions\n", k);

    printf("\ndepth\tN\n");
    for(i=0; i<1024; i++){
        printf("%d\t%d\n", i, depth[i]);
    }

  }