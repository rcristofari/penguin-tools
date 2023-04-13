#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

/* REMAINS TO DO: Implement various methods for cutoff: threshold, stdev, sliding window...*/

// Simple math functions (max, min, mean, stdev):
int min(int x, int y){
  return y ^ ((x ^ y) & -(x < y));
}

int max(int x, int y){
  return x ^ ((x ^ y) & -(x < y));
}

float mean(int *array, int n){
    int i, sum = 0;
    for (i=0; i<n; i++) {
        sum += array[i];
    }
    return (float)sum / n ;
}

float stdev(int *array, float mean, int n){
    int i, sumsq=0;
    for (i=0; i<n; i++){
        sumsq += pow((mean - array[i]), 2);
    }
    return pow((float)sumsq / n, 0.5) ;
}


int main(int argc, char *argv[]) {

    int opt;

    // Input files:
    FILE *reffile; // reference genome fasta file
    FILE *gtffile; // annotation GTF file
    FILE *depthfile; // per-base depth file from samtools

    // Bin file for depth
    FILE *bfile;
    int do_bin_depth = 0; // whether we need to recompute the depth bin file
    char bin_path[1024];

    // Output file
    FILE *ofile;
    int ofile_set = 0;

    // File-reading variables:
    char * line = NULL;
    size_t len = 0;
    size_t read;
    char *sp=NULL;

    // Genome coordinate variables:
    int i = 1;
    int chr_size = 0;
    int abs_start = 0;
    int genome_size = 0;
    int n_scaffolds = 0;
    char *chr;
    // a struct to contain scaffolds and their coordinates:
    struct Scaffold {char chr[32]; int abs; int len;};
    // Allocate starting memory for the first 1000 scaffolds:
    int initial_size = 1000;
    struct Scaffold *scaffolds, *new_scaffolds;
    scaffolds = malloc(initial_size * sizeof(struct Scaffold));

    // Depth file variables:
    int depth, max_depth = 0;

    // GTF file scanning variables:
    char *scaf, *feature, *strand, *start, *end;
    int search_start, search_len;
    int j=1;

    // TSS finding variables:
    int depth_threshold = 3, transcript_is_missing = 0;
    float this_mean = 0.0, this_stdev = 0.0;
    int begin_around_start, start_codon, TSS;
    int size_around_start = 50;
    int around_start[size_around_start]; // an array to store depths right after the ATG (as a proxy of expected expression of transcript)
    int size_utr_search = 1000;
    int utr_search[size_utr_search]; //an array to store positions upstream of ATG until TSS is found (we will explore 1kb upstream of the ATG)
    int window_size = 3, window_sum = 0, k = 0;
    float rolling_mean = 0.0;

    printf("-------------------------------------------------------------------------------------------------------\n");

    // parse command line
    while ((opt = getopt(argc, argv, "r:g:d:b:t:o:")) != -1) {
        switch (opt) {
            case 'r':
                if (access(optarg, F_OK) == 0) {
                    reffile = fopen(optarg, "r");
                    printf("Reference genome file: %s\n", optarg);
                } else {
                    printf("ERROR: reference genome file %s does not exist.\n", optarg);
                    return(1);
                }
                break;
            case 'g':
                if (access(optarg, F_OK) == 0) {
                    gtffile = fopen(optarg, "r");
                    printf("Annotation GTF file: %s\n", optarg);
                } else {
                    printf("ERROR: annotation gtf file %s does not exist.\n", optarg);
                    return(1);
                }
                break;
            case 'd':
                if (access(optarg, F_OK) == 0) {
                    depthfile = fopen(optarg, "r");
                    do_bin_depth = 1;
                    strcpy(bin_path, "./depth.bin");
                    printf("Alignment depth file: %s\n", optarg);
                } else {
                    printf("ERROR: alignment depth file %s does not exist.\n", optarg);
                    return(1);
                }
                break;
            case 'b':
                if (access(optarg, F_OK) == 0) {
                    strcpy(bin_path, optarg);
                    do_bin_depth = 0;
                    printf("Pre-computed binary depth file: %s\n", bin_path);
                } else {
                    printf("ERROR: pre-computed binary depth file %s does not exist.\n", optarg);
                    return(1);
                }
                break;
            case 't':
                depth_threshold = atoi(optarg);
                printf("Discarding genes with first exon coverage < %d\n", depth_threshold);
                break;
            case 'o':
                ofile_set = 1;
                ofile = fopen(optarg, "w");
                printf("Writing TSS coordinates to %s\n", optarg);
                break;
            default:
                fprintf(stderr, "Usage: %s -r reference.fasta -g annotation.gtf -d depth.txt -t depth_threshold\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // Output file path:
    if (ofile_set == 0){
        ofile = fopen("./tss.coord", "w");
        printf("Writing TSS coordinates to tss.coord\n");
        ofile_set = 1;
    }

    fprintf(ofile, "#chr\ttss\tfull_pos\t5utr_len\tstrand\n");

    //----------------------------------------------------------------------------------------------------------------//
    // Start with indexing the genome (scaffold names and absolute start):
    printf("Indexing scaffolds...");
    fflush(stdout); // force output even without a newline
    while ((read = getline(&line, &len, reffile)) != -1) {
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
            if (n_scaffolds > 0) {
                scaffolds[n_scaffolds-1].len = chr_size;
            }
            chr_size = 0;
            n_scaffolds += 1;

        } else {

            chr_size += (read-1);
            genome_size += (read-1);
        }
    }

    scaffolds[n_scaffolds-1].len = chr_size;

    printf("done.\nFound %d scaffolds.\nGenome size: %d bp.\n", n_scaffolds, genome_size);

    fclose(reffile);
    if (line) {
        free(line);
    }
    // End of the genome coordinate indexing

    //----------------------------------------------------------------------------------------------------------------//
    // Change the depth file into a binary integer file for quick access
    if (do_bin_depth == 1) {
        printf("Indexing the depth file...");
        fflush(stdout);
        bfile = fopen(bin_path, "wb");
        while ((read = getline(&line, &len, depthfile)) != -1) {
            sp = strtok(line, "\t"); // first field, scaffold
            sp = strtok(NULL, "\t"); // second field, position
            sp = strtok(NULL, "\t"); // third field, depth (the only one we are interested in here)
            depth = atoi(sp);
            //max_depth = max(depth, max_depth);
            //printf("%d\r", max_depth);
            fwrite(&depth, sizeof(int), 1, bfile);
        }
        fclose(depthfile);
        fclose(bfile);
        printf("done.\n");
    }

    //----------------------------------------------------------------------------------------------------------------//
    // Start opening the GTF file
    printf("Interating through the annotations:\n");
    bfile = fopen(bin_path, "rb");

    while ((read = getline(&line, &len, gtffile)) != -1) {
        if (*line != '#'){
            scaf = strtok(line, "\t");
            sp = strtok(NULL, "\t"); // origin of the annotation
            feature = strtok(NULL, "\t");
            start = strtok(NULL, "\t");
            end = strtok(NULL, "\t");
            sp = strtok(NULL, "\t"); // score (null)
            strand = strtok(NULL, "\t");

            // This is a start codon, we will look up the coordinates in absolute values
            if (strcmp(feature, "start_codon") == 0) {
                printf("%s\t%s\t%s\t%s\t%s\n", scaf, feature, start, end, strand);

                j--;
                while(strcmp(scaf, scaffolds[j].chr) != 0) {
                    j++;
                }

                // Search coordinates, 0-based:
                if (strcmp(strand, "+") == 0) {
                    search_start = max(0, atoi(start) - size_utr_search - 1) + scaffolds[j].abs;
                    search_len = (atoi(start) + scaffolds[j].abs - 1) - search_start;
                } else {
                    search_start = atoi(end) + scaffolds[j].abs - 1;
                    search_len = (min(scaffolds[j].len - 1, atoi(end) + size_utr_search - 1) + scaffolds[j].abs) - search_start;
                }

                //printf("Searching from %d to %d\n", search_start, search_end);

                // Retrieve the depth data around the start codon:
                if(strcmp(strand, "+") == 0){
                    begin_around_start = atoi(start) + scaffolds[j].abs - 1 ;
                } else {
                    begin_around_start = atoi(start) - size_around_start + scaffolds[j].abs - 1 ;
                }



                fseek(bfile, begin_around_start*sizeof(int), SEEK_SET);
                fread(&around_start, sizeof(int), size_around_start, bfile);
//                for (i=0; i<size_around_start; i++){
//                    printf("%d ", around_start[i]);
//                }

                this_mean = mean(around_start, size_around_start);
                this_stdev = stdev(around_start, this_mean, size_around_start);
                //printf("Mean depth: %f +- %f\n", this_mean, this_stdev);

                // If the mean in the first exon is too low, we consider we have no data for the gene:
                if (this_mean < depth_threshold) {
                    transcript_is_missing = 1;
                    //printf("Transcript is missing / too low coverage\n---\n");
                } else {

                // Retrieve the depth data in the search area:
                fseek(bfile, search_start*sizeof(int), SEEK_SET);
                fread(&utr_search, sizeof(int), size_utr_search, bfile);
//                this_mean = mean(utr_search, size_utr_search);
//                this_stdev = stdev(utr_search, this_mean, size_utr_search);
//                printf("Mean depth: %f +- %f\n", this_mean, this_stdev);

//                if (strcmp(strand, "+")==0){
//                    i = size_utr_search-3;
//                    rolling_mean = this_mean;
//                    //while(i > 0 && utr_search[i] > (this_mean - 2*this_stdev)) {
//                    while(i > ((window_size-1)/2) && rolling_mean > depth_threshold) {
//                    //while(i > 2 && rolling_mean > 5) {
//                        window_sum = 0;
//                        for (k=(1-window_size)/2; k<(window_size-1)/2; k++){
//                            window_sum += utr_search[i+k];
//                        }
//                        rolling_mean = (float)window_sum / window_size ;
//                        i--;
//                    }
//                    printf("Boundary %d bp before start codon, coverage drops to %f\n---\n", (size_utr_search-i), rolling_mean);
//
//                } else {
//                    i = 2;
//                    rolling_mean = this_mean;
//                    while(i < (size_utr_search-(((window_size-1)/2)+1)) && rolling_mean > depth_threshold) {
//                    //while(i < (size_utr_search-3) && rolling_mean > 5) {
//                        window_sum = 0;
//                        for (k=(1-window_size)/2; k<(window_size-1)/2; k++){
//                            window_sum += utr_search[i+k];
//                        }
//                        rolling_mean = (float)window_sum / window_size ;
//                        i++;
//                    }
//                    printf("Boundary %d bp before start codon, coverage drops to %f\n---\n", i, rolling_mean);
//
//                }

                if (strcmp(strand, "+")==0){
                    i = search_len - 1;
                    rolling_mean = this_mean;
                    //while(i > 0 && utr_search[i] > (this_mean - 2*this_stdev)) {
                    while(i > 0 && utr_search[i] > depth_threshold) {
                        i--;
                    }
                    //printf("Boundary %d bp before start codon, coverage drops to %d\n---\n", (size_utr_search-i), utr_search[i]);
                    //printf("%d\n", size_utr_search-i);
                    fprintf(ofile, "%s\t%d\t%d\t%d\t%s\n", scaf, atoi(start)-(search_len-i), atoi(start)-(search_len-i)+scaffolds[j].abs, size_utr_search-i, strand);

                } else {
                    i = 0;
                    //while(i < size_utr_search && utr_search[i] > (this_mean - 2*this_stdev)) {
                    while(i < search_len && utr_search[i] > depth_threshold) {
                        i++;
                    }
                    //printf("Boundary %d bp before start codon, coverage drops to %d\n---\n", i, utr_search[i]);
                    //printf("%d\n", i);
                    fprintf(ofile, "%s\t%d\t%d\t%d\t%s\n", scaf, atoi(end)+i, atoi(end)+i+scaffolds[j].abs, i, strand);
                }



                }

            }
        }
    }
    fclose(bfile);
    fclose(gtffile);

}