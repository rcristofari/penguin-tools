#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

int main(int argc, char *argv[]) {

    int opt; // options in the argument parser
    FILE *gfile;
    FILE *pfile;
    int i = 1;
    int chr_size = 0;
    int abs_start = 0;
    char * line = NULL;
    size_t len = 0;
    size_t read;
    char *sp;
    char *chr;
    int fullpos;
    int n_scaffolds = 0;

    // a struct to contain scaffolds and their coordinates:
    struct Scaffold {char chr[32]; int abs;};

    // Allocate starting memory for the first 1000 scaffolds:
    int initial_size = 1000;
    struct Scaffold *scaffolds, *new_scaffolds;
    scaffolds = malloc(initial_size * sizeof(struct Scaffold));

    // Command line argument parsing
    while ((opt = getopt(argc, argv, "g:p:")) != -1) {
        switch (opt) {
            case 'g':
                gfile = fopen(optarg, "r");
                break;
            case 'p':
                pfile = fopen(optarg, "r");
                break;
            default:
                fprintf(stderr, "Usage: %s -g genome.fasta -p position.list\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    while ((read = getline(&line, &len, gfile)) != -1) {
        if (*line == '>') {

            if (n_scaffolds >= initial_size) {
                initial_size *= 2;
                new_scaffolds = (struct Scaffold*)realloc(scaffolds, initial_size*sizeof(struct Scaffold));
                if (new_scaffolds == NULL) {
                    printf("Memory could not be allocated");
                } else {
                    scaffolds = new_scaffolds;
                }
            }

            abs_start += chr_size;
            sp = strtok(line, " ");
            chr = sp + 1;

            strcpy(scaffolds[n_scaffolds].chr, chr);
            scaffolds[n_scaffolds].abs = abs_start;

            chr_size = 0;
            n_scaffolds += 1;

        } else {
            chr_size += (read-1);
        }
    }

    fclose(gfile);
    if (line)
        free(line);

    while ((read = getline(&line, &len, pfile)) != -1) {
        fullpos = atoi(line);

        // To avoid re-scanning the list everytime, we restart just ahead of where we left off:
        i -= 1;
        while (scaffolds[i].abs < fullpos) {
            i++;
            if (i >= n_scaffolds){
                break;
            }
        }
        printf("%s\t%d\t%d\n", scaffolds[i-1].chr, fullpos - scaffolds[i-1].abs, fullpos - scaffolds[i-1].abs);
    }


    exit(EXIT_SUCCESS);



}