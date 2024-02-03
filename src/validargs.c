#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "debug.h"

int stringcomp(char *str1, char *str2);

/**
 * @brief Validates command line arguments passed to the program.
 * @details This function will validate all the arguments passed to the
 * program, returning 0 if validation succeeds and -1 if validation fails.
 * Upon successful return, the various options that were specified will be
 * encoded in the global variable 'global_options', where it will be
 * accessible elsewhere in the program.  For details of the required
 * encoding, see the assignment handout.
 *
 * @param argc The number of arguments passed to the program from the CLI.
 * @param argv The argument strings passed to the program from the CLI.
 * @return 0 if validation succeeds and -1 if validation fails.
 * @modifies global variable "global_options" to contain an encoded representation
 * of the selected program options.
 */
int validargs(int argc, char **argv)
{

    if(argc == 1){
        global_options = 0;
        return 0;
    }

    char *m = "-m";
    char *n = "-n";
    char *o = "-o";
    char *h = "-h";

    if(stringcomp(*(argv + 1), h)){
        global_options = HELP_OPTION;
        return 0;
    }

    if(stdin == NULL){
        return -1;
    }

//m and any other
    for(int i = 1; i < argc; i++){
        if(stringcomp(*(argv + i), m)){
            for(int j = 1; j < argc; j++){
                if(j != i)
                    if(stringcomp(*(argv+j), m) == 1 || stringcomp(*(argv+j), n) == 1  || stringcomp(*(argv+j), o) == 1  || stringcomp(*(argv+j), h) == 1 )
                        return -1;
            }
        }
    }

//h not at the front
    if(argc > 2){
        for(int i = 2; i < argc; i++){
            if(stringcomp(*(argv+i), h) == 1){
                return -1;
            }
        }
    }

//o in front of n
    for(int i = 1; i < argc; i++){
        if(stringcomp(*(argv + i), n)){
            for(int j = 1; j < i; j++){
                if(j != i)
                    if(stringcomp(*(argv+j), o) == 1)
                        if(j < i)
                            return -1;
            }
        }
    }

//no arg for o
    for(int i = 1; i < argc; i++){
        if(stringcomp(*(argv + i), o)){
            if(*(argv + i+1) == NULL || *(*(argv + i+1)+0) == '-'){
                return -1;
            }
        }
    }

//m and any other
    for(int i = 1; i < argc; i++){
        if(*(*(argv + i)+0) == '-')
            if(*(*(argv + i)+1) != 'm' && *(*(argv + i)+1) != 'n' && *(*(argv + i)+1) != 'o' && *(*(argv + i)+1) != 'h'){
                printf("%c\n",*(*(argv + i)+1));
                return -1;
            }
    }

//no args
    int tax = 0;
    for(int i = 1; i < argc; i++){
        if(stringcomp(*(argv + i), m) == 1 || stringcomp(*(argv + i), n) == 1 ||  stringcomp(*(argv + i), o) == 1 || stringcomp(*(argv + i), h) == 1){
            tax = 1;
        }
    }
    if(tax == 0)
        return -1;


    for(int i = 1; i < argc; i++){
        if(stringcomp(*(argv + i), m) == 1){
            global_options = MATRIX_OPTION;
        }
    }
    int temp = 0;
    for(int i = 1; i < argc; i++){
        if(stringcomp(*(argv + i), n) == 1){
            global_options = NEWICK_OPTION;
        }
        if(stringcomp(*(argv + i), o) == 1){
            outlier_name = *(argv + i+1);
            temp = 1;
        }
    }
    if(temp == 0){
        outlier_name = NULL;
    }
    return 0;
}

int stringcomp(char *str1, char *str2){
    while(*str1 != '\0' && *str2 != '\0'){
        if(*str1 != *str2)
            return 0;
        str1++;
        str2++;
    }
    return 1;
}
