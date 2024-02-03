#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "debug.h"

int main(int argc, char **argv)
{
    if(validargs(argc, argv))
        USAGE(*argv, EXIT_FAILURE);
    if(global_options == HELP_OPTION)
        USAGE(*argv, EXIT_SUCCESS);

    FILE *infile = stdin;

    if(global_options == MATRIX_OPTION){
        if(read_distance_data(infile) == 0){
            return EXIT_SUCCESS;
        } else {
            fprintf(stderr, "Error with -m option\n");
            return EXIT_FAILURE;
        }
    } else if(global_options == NEWICK_OPTION){
        if(read_distance_data(infile) == 0){
            return EXIT_SUCCESS;
        } else {
            fprintf(stderr, "Error with -n option\n");
            return EXIT_FAILURE;
        }
    } else{
         if(read_distance_data(infile) == 0){
            return EXIT_SUCCESS;
        } else {
            fprintf(stderr, "Error with no option\n");
            return EXIT_FAILURE;
        }
    }

    return EXIT_FAILURE;
}

/*
 * Just a reminder: All non-main functions should
 * be in another file not named main.c
 */
