#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "debug.h"

#define DBL_MAX 1.7976931348623158e+308
void calculateQ1(int *min_i, int *min_j);
void findMin1(int *min_i, int *min_j);
double sumRow(int row);
void printInOrder(struct node* root, struct node* node2, FILE *out);
int stringcomp1(char *str1, char *str2);
int callerHelp(struct node *root,struct node *current);
int countDigits(int b);
int tempcount = 0;
int taxacount = 0;
int track1 = 0;
int numtemp = 0;
int numtemp2 = 0;
int numtemp3 = 0;

int countDigits(int b) {
    int numDigits = 0;
    while (b != 0) {
        b /= 10;
        numDigits++;
    }
    return numDigits;
}

/**
 * @brief  Read genetic distance data and initialize data structures.
 * @details  This function reads genetic distance data from a specified
 * input stream, parses and validates it, and initializes internal data
 * structures.
 *
 * The input format is a simplified version of Comma Separated Values
 * (CSV).  Each line consists of text characters, terminated by a newline.
 * Lines that start with '#' are considered comments and are ignored.
 * Each non-comment line consists of a nonempty sequence of data fields;
 * each field iis terminated either by ',' or else newline for the last
 * field on a line.  The constant INPUT_MAX specifies the maximum number
 * of data characters that may be in an input field; fields with more than
 * that many characters are regarded as invalid input and cause an error
 * return.  The first field of the first data line is empty;
 * the subsequent fields on that line specify names of "taxa", which comprise
 * the leaf nodes of a phylogenetic tree.  The total number N of taxa is
 * equal to the number of fields on the first data line, minus one (for the
 * blank first field).  Following the first data line are N additional lines.
 * Each of these lines has N+1 fields.  The first field is a taxon name,
 * which must match the name in the corresponding num_all_nodes of the first line.
 * The subsequent fields are numeric fields that specify N "distances"
 * between this taxon and the others.  Any additional lines of input following
 * the last data line are ignored.  The distance data must form a symmetric
 * matrix (i.e. D[i][j] == D[j][i]) with zeroes on the main diagonal
 * (i.e. D[i][i] == 0).
 *
 * If 0 is returned, indicating data successfully read, then upon return
 * the following global variables and data structures have been set:
 *   num_all_nodes - set to the number N of taxa, determined from the first data line
 *   num_all_nodes - initialized to be equal to num_all_nodes
 *   num_active_nodes - initialized to be equal to num_all_nodes
 *   node_names - the first N entries contain the N taxa names, as C strings
 *   distances - initialized to an NxN matrix of distance values, where each
 *     row of the matrix contains the distance data from one of the data lines
 *   nodes - the "name" fields of the first N entries have been initialized
 *     with pointers to the corresponding taxa names stored in the node_names
 *     array.
 *   active_node_map - initialized to the identity mapping on [0..N);
 *     that is, active_node_map[i] == i for 0 <= i < N.
 *
 * @param in  The input stream from which to read the data.
 * @return 0 in case the data was successfully read, otherwise -1
 * if there was any error.  Premature termination of the input data,
 * failure of each line to have the same number of fields, and distance
 * fields that are not in numeric format should cause a one-line error
 * message to be printed to stderr and -1 to be returned.
 */

int read_distance_data(FILE *in) {
    if(in == NULL) {
        fprintf(stderr, "Error: File is NULL.\n");
        return -1;
    }

    int c;
    int skipLine = 0;
    int x = 0;
    //matrix with only the information
    while((c = fgetc(in)) != EOF){

        //ignore # comment
        if(c == '#') {
            skipLine = 1;
        }
        while(c != EOF && c != '\n'){
            if(skipLine == 0){
                *(input_buffer+x) = c;
                x++;
            }
            c = fgetc(in);
        }

        if(skipLine == 0){
            *(input_buffer+x) = c;
            x++;
        }
        skipLine = 0;
    }
    *(input_buffer+x)= '\0';
    num_all_nodes = 0;
    num_active_nodes = 0;
    num_taxa = 0;

    double a = *(input_buffer+0);
    x = 0;
    int y = 0;
    while (a != '\0' && a != '\n') {
        if(a != ','){
            *(*(node_names + y)+0) = a;
            y++;
            num_active_nodes++;
            num_all_nodes++;
            (nodes + num_all_nodes)->name = *(node_names + num_all_nodes);
        }
        x++;
        a = *(input_buffer+x);
    }

    //sets all the names in the nodes and then moves the matrix up by 1 to get rid of names
    a = *(input_buffer+0);
    int counter = 0;
    x++;
    while(*(input_buffer+x) != '\0'){
        *(input_buffer+counter) = *(input_buffer+x);
        counter++;
        x++;
    }
    *(input_buffer+counter) = '\0';
    x = counter;
    int a1 = 0;
    int b1 = 0;
    int nL = 0;
    for(int i = 0; i < x; i++){
        if(nL == 0){
            while(a != ','){
                i++;
                a = *(input_buffer+i);
            }
            nL = 1;
        } else if (nL == 1){
            a = *(input_buffer+i);
        }
        if(a == '\n'){
            nL = 0;
        }
        if(a == '0' || a == '1' || a == '2'|| a == '3' || a == '4' || a == '5' ||
                         a == '6' || a == '7' || a == '8' || a == '9' || a == '.'){
            if(a != '.'){
                a = a-48;
                double b = 0;
                b = *(input_buffer+i+1);
                while(b == '0' || b == '1' || b == '2'|| b == '3' || b == '4' || b == '5' ||
                          b == '6' || b == '7' || b == '8' || b == '9'){
                        b = b-48;
                        a = ((a*10)+ b);
                        i++;
                        b = *(input_buffer+i+1);
                }
                if(b == '.'){
                    i++;
                    b = *(input_buffer+i+1);
                    double ta = 0;
                    while(b == '0' || b == '1' || b == '2'|| b == '3' || b == '4' || b == '5' ||
                          b == '6' || b == '7' || b == '8' || b == '9'){
                        b = b-48;
                        ta = ((ta*10)+ b);
                        i++;
                        b = *(input_buffer+i+1);
                    }
                        int numDigits = countDigits(ta);
                        double divisor = 1.0;
                        for (int k = 0; k < numDigits; k++) {
                            divisor *= 10.0;
                        }
                        double decimalb = (double) ta / divisor;
                        a+= decimalb;
                }
            } else{
                a = 0;
                double b = 0;
                b = *(input_buffer+i+1);
                double ta = 0;
                while(b == '0' || b == '1' || b == '2'|| b == '3' || b == '4' || b == '5' ||
                      b == '6' || b == '7' || b == '8' || b == '9'){
                    b = b-48;
                    ta = ((ta*10)+ b);
                    i++;
                    b = *(input_buffer+i+1);
                }
                int numDigits = countDigits(ta);
                double divisor = 1.0;
                for (int k = 0; k < numDigits; k++) {
                    divisor *= 10.0;
                }
                double decimalb = (double) ta / divisor;
                a+= decimalb;
            }
            if(b1 < num_all_nodes){
                *(*(distances+a1)+b1) = a;
            } else{
                a1++;
                b1 = 0;
                *(*(distances+a1)+b1) = a;
            }
            b1++;
        }
    }
    if(global_options == 0){
        FILE *file = stdout;
        if(emit_distance_matrix(file) == -1){
            fprintf(stderr, "An error occurred: Printing distances failed.\n");
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;

    } else if(global_options == MATRIX_OPTION){
        FILE *file = stdout;
        if(emit_distance_matrix(file) == -1){
            fprintf(stderr, "An error occurred: Distance matrix failed.\n");
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    } else if(global_options == NEWICK_OPTION){
        FILE *file = stdout;
        if(emit_newick_format(file) == -1){
            fprintf(stderr, "An error occurred: Newick format failed.\n");
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    }
    return EXIT_SUCCESS;
}

void calculateQ1(int *min_i, int *min_j){
    for(int i = 0; i < MAX_NODES; i++){
        *(row_sums+i) = 0;
    }

    for(int i = 0; i < num_taxa; i++){
        for(int j = num_active_nodes; j < num_taxa+num_active_nodes; j++){
            *(row_sums + i) += *(*(distances + *(active_node_map + i)) + j);
        }
    }
    int r = 100;
    int q = 0;
    for(int i = 0; i < num_taxa; i++){
        for(int j = num_active_nodes; j< num_taxa+num_active_nodes; j++){
            if(i != j){
                q = (num_taxa-2) * *(*(distances + *(active_node_map + i)) + j) - *(row_sums+i) - *(row_sums+j);
                if(q < r){
                    r = q;
                    *min_i = i;
                    *min_j = *(active_node_map+j);
                }
            }
        }
    }
}

void findMin1(int *min_i, int *min_j) {
    double min_dist = *(row_sums);

    for (int i = 0; i < num_taxa; i++) {
        for (int j = i + 1; j < num_taxa; j++) {
            if (*(row_sums + i + j) < min_dist) {
                min_dist = *(row_sums + i + j);
                *min_i = i;
                *min_j = j;
            }
        }
    }
}

double sumRow(int row) {
    double sum = 0.0;
    for (int i = 0; i < num_taxa; i++) {
        if (i != row) {
            sum += *(*(distances + *(active_node_map+row)) + *(active_node_map+i));
        }
    }
    return sum;
}

/**
 * @brief  Emit a representation of the phylogenetic tree in Newick
 * format to a specified output stream.
 * @details  This function emits a representation in Newick format
 * of a synthesized phylogenetic tree to a specified output stream.
 * See (https://en.wikipedia.org/wiki/Newick_format) for a description
 * of Newick format.  The tree that is output will include for each
 * node the name of that node and the edge distance from that node
 * its parent.  Note that Newick format basically is only applicable
 * to rooted trees, whereas the trees constructed by the neighbor
 * joining method are unrooted.  In order to turn an unrooted tree
 * into a rooted one, a root will be identified according by the
 * following method: one of the original leaf nodes will be designated
 * as the "outlier" and the unique node adjacent to the outlier
 * will serve as the root of the tree.  Then for any other two nodes
 * adjacent in the tree, the node closer to the root will be regarded
 * as the "parent" and the node farther from the root as a "child".
 * The outlier node itself will not be included as part of the rooted
 * tree that is output.  The node to be used as the outlier will be
 * determined as follows:  If the global variable "outlier_name" is
 * non-NULL, then the leaf node having that name will be used as
 * the outlier.  If the value of "outlier_name" is NULL, then the
 * leaf node having the greatest total distance to the other leaves
 * will be used as the outlier.
 *
 * @param out  Stream to which to output a rooted tree represented in
 * Newick format.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.  If the global variable "outlier_name" is
 * non-NULL, then it is an error if no leaf node with that name exists
 * in the tree.
 */
int emit_newick_format(FILE *out) {
    numtemp = num_all_nodes;
    numtemp2 = num_active_nodes;
    numtemp3 = num_taxa;
    if(build_taxonomy(out) == EXIT_SUCCESS){
        if(out != NULL){
            num_all_nodes--;

            //finding greatest total distance outlier.
            int sum = 0;
            int sum1 = 0;
            int index = 0;
            for (int i = 0; i < numtemp; i++) {
                sum = 0;
                for (int j = 0; j < numtemp; j++) {
                    sum += *(*(distances + i) + j);
                }
                if(sum1 < sum){
                    index = i;
                    sum1 = sum;
                }
            }
            if(outlier_name != NULL){
                for(int i = 0; i < numtemp; i++){
                    if(stringcomp1(outlier_name, *(node_names + i)) == 1){
                        index = i;
                    }
                }
            }
            printInOrder(*((nodes + index)->neighbors + 0), (nodes+index), out);
            fprintf(out, "\n");
        }

        num_all_nodes = numtemp;
        num_active_nodes = numtemp2;
        num_taxa = numtemp3;
        return EXIT_SUCCESS;
    } else{
        return EXIT_FAILURE;
    }
}

void printInOrder(struct node* node1, struct node* node2, FILE *out) {
    int comma = 0;
    if(*(node1 -> neighbors + 1) != NULL || *(node1 -> neighbors + 2) != NULL ){
        fprintf(out, "(");
    }
    // Recursively print the subtree
    if(*(node1 -> neighbors + 0) != NULL && *(node1 -> neighbors + 0) != node2){
        printInOrder(*(node1 -> neighbors + 0), node1, out);

        if(*(node1 -> neighbors + 1) != NULL || *(node1 -> neighbors + 2) != NULL ) {
            fprintf(out,",");
            comma = 1;
        }
    }

    if(*(node1 -> neighbors + 2) != NULL && *(node1 -> neighbors + 2) != node2){
        printInOrder(*(node1 -> neighbors + 2), node1, out);


    }
    if(*(node1 -> neighbors + 1) != NULL && *(node1 -> neighbors + 1) != node2){
        if (!comma)
            fprintf(out,",");
        printInOrder(*(node1 -> neighbors + 1), node1, out);
    }
    if(*(node1 -> neighbors + 1) != NULL || *(node1 -> neighbors + 2) != NULL ){
        fprintf(out,")");
        fflush(stdout);
    }
    if(node1 != NULL){
        for(int i = 0; i < num_all_nodes+1; i++){
            if(stringcomp1(*(node_names + i), node1 -> name) == 1){
                if(i < numtemp){
                    fprintf(out,"%s:", node1 -> name);
                } else{
                    fprintf(out, "%c", *(*(node_names + i)));
                    fprintf(out, "%d:", (*(*(node_names + i)+1)-48));
                }
            }
        }
        fflush(stdout);
        int x = 0;
        int y = 0;
        for(int i = 0; i< num_all_nodes+1; i++){
            if(stringcomp1((node1 -> name), *(node_names+i))){
                x = i;
            }
            if(stringcomp1((node2 -> name), *(node_names+i))){
                y = i;
            }
        }
        fprintf(out,"%.2f",*(*(distances+x)+y));
    }
}


int stringcomp1(char *str1, char *str2){
    while(*str1 != '\0' && *str2 != '\0'){
        if(*str1 != *str2)
            return 0;
        str1++;
        str2++;
    }
    return 1;
}

/**
 * @brief  Emit the synthesized distance matrix as CSV.
 * @details  This function emits to a specified output stream a representation
 * of the synthesized distance matrix resulting from the neighbor joining
 * algorithm.  The output is in the same CSV form as the program input.
 * The number of rows and num_all_nodes of the matrix is equal to the value
 * of num_all_nodes at the end of execution of the algorithm.
 * The submatrix that consists of the first num_leaves rows and num_all_nodess
 * is identical to the matrix given as input.  The remaining rows and num_all_nodess
 * contain estimated distances to internal nodes that were synthesized during
 * the execution of the algorithm.
 *
 * @param out  Stream to which to output a CSV representation of the
 * synthesized distance matrix.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.
 */
int emit_distance_matrix(FILE *out) {
    numtemp = num_all_nodes;
    numtemp2 = num_active_nodes;
    numtemp3 = num_taxa;

    if(build_taxonomy(out) == EXIT_SUCCESS){
        num_all_nodes = numtemp;
        num_active_nodes = numtemp2;
        num_taxa = numtemp3;
        return EXIT_SUCCESS;
    } else{
        return EXIT_FAILURE;
    }
}

/**
 * @brief  Build a phylogenetic tree using the distance data read by
 * a prior successful invocation of read_distance_data().
 * @details  This function assumes that global variables and data
 * structures have been initialized by a prior successful call to
 * read_distance_data(), in accordance with the specification for
 * that function.  The "neighbor joining" method is used to reconstruct
 * phylogenetic tree from the distance data.  The resulting tree is
 * an unrooted binary tree having the N taxa from the original input
 * as its leaf nodes, and if (N > 2) having in addition N-2 synthesized
 * internal nodes, each of which is adjacent to exactly three other
 * nodes (leaf or internal) in the tree.  As each internal node is
 * synthesized, information about the edges connecting it to other
 * nodes is output.  Each line of output describes one edge and
 * consists of three comma-separated fields.  The first two fields
 * give the names of the nodes that are connected by the edge.
 * The third field gives the distance that has been estimated for
 * this edge by the neighbor-joining method.  After N-2 internal
 * nodes have been synthesized and 2*(N-2) corresponding edges have
 * been output, one final edge is output that connects the two
 * internal nodes that still have only two neighbors at the end of
 * the algorithm.  In the degenerate case of N=1 leaf, the tree
 * consists of a single leaf node and no edges are output.  In the
 * case of N=2 leaves, then no internal nodes are synthesized and
 * just one edge is output that connects the two leaves.
 *
 * Besides emitting edge data (unless it has been suppressed),
 * as the tree is built a representation of it is constructed using
 * the NODE structures in the nodes array.  By the time this function
 * returns, the "neighbors" array for each node will have been
 * initialized with pointers to the NODE structure(s) for each of
 * its adjacent nodes.  Entries with indices less than N correspond
 * to leaf nodes and for these only the neighbors[0] entry will be
 * non-NULL.  Entries with indices greater than or equal to N
 * correspond to internal nodes and each of these will have non-NULL
 * pointers in all three entries of its neighbors array.
 * In addition, the "name" field each NODE structure will contain a
 * pointer to the name of that node (which is stored in the corresponding
 * entry of the node_names array).
 *
 * @param out  If non-NULL, an output stream to which to emit the edge data.
 * If NULL, then no edge data is output.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.
 */
int build_taxonomy(FILE *out) {
    num_active_nodes = 0;
    num_taxa = num_all_nodes - num_active_nodes;
    int min_i, min_j;
    for(int i =0; i < num_taxa; i++){
        *(active_node_map+i) = i;
    }
    for(int i = 0; i< num_taxa; i++){
        (nodes+i) -> name = *(node_names+i);
    }

    if(global_options == MATRIX_OPTION || global_options == NEWICK_OPTION || global_options == HELP_OPTION){
        while(num_taxa>3){
            num_taxa = num_all_nodes - num_active_nodes;
            calculateQ1(&min_i, &min_j);
            *(*(node_names + num_all_nodes)) = '#';
            *(*(node_names + num_all_nodes)+1) = (char)(num_all_nodes+48);
            (nodes + num_all_nodes)->name = *(node_names + num_all_nodes);
            // Add a new node
            for (int k = 0; k < num_taxa; k++) {
                double limb_i = 0.5 * *(*(distances + *(active_node_map+min_i)) +
                *(active_node_map+min_j)) +  1.0 / (2.0*(num_taxa - 2.0)) * (sumRow(min_i) - sumRow(min_j));
                double limb_j = *(*(distances + *(active_node_map+min_i)) + *(active_node_map+min_j)) - limb_i;

                if(k == num_all_nodes){
                    *(*(distances+num_all_nodes)+num_all_nodes) = 0;
                    *(*(distances+*(active_node_map+k))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+k));
                } else if(*(active_node_map+k) == *(active_node_map+min_i) && num_taxa != 3){

                    //sets parent node
                    *((nodes+*(active_node_map+k))->neighbors + 0) = (nodes + (num_all_nodes));

                    //sets children node
                    *((nodes+(num_all_nodes))->neighbors + 1) = (nodes + *(active_node_map+k));

                    *(*(distances+num_all_nodes)+*(active_node_map+k)) = limb_i;
                    *(*(distances+*(active_node_map+k))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+k));

                } else if(*(active_node_map+k) == *(active_node_map+min_j) && num_taxa != 3){

                    //sets parent node
                    *((nodes+*(active_node_map+k))->neighbors + 0) = (nodes + (num_all_nodes));

                    //sets children node
                    *((nodes+(num_all_nodes))->neighbors + 2) = (nodes + *(active_node_map+k));

                    *(*(distances+num_all_nodes)+*(active_node_map+k)) = limb_j;
                    *(*(distances+*(active_node_map+k))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+k));

                } else {
                    if(num_taxa != 3){
                        *(*(distances+num_all_nodes)+*(active_node_map+k)) = (*(*(distances+*(active_node_map+min_i))+*(active_node_map+k)) +
                            *(*(distances+*(active_node_map+min_j))+*(active_node_map+k)) -
                            *(*(distances+*(active_node_map+min_i))+*(active_node_map+min_j))) / 2;

                        *(*(distances+*(active_node_map+k))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+k));
                    } else if(num_taxa == 3){
                        double first = 0.5 * *(*(distances + *(active_node_map+0)) + *(active_node_map+1)) + 1.0 / (2.0*(num_taxa - 2.0)) * (sumRow(0) - sumRow(1));
                        //sets parent node
                        *((nodes+*(active_node_map+0))->neighbors + 0) = (nodes + (num_all_nodes));
                        *((nodes+*(active_node_map+1))->neighbors + 0) = (nodes + (num_all_nodes));
                        *((nodes+*(active_node_map+2))->neighbors + 0) = (nodes + (num_all_nodes));

                        //sets children node
                        *((nodes+(num_all_nodes))->neighbors + 0) = (nodes + *(active_node_map+0));
                        *((nodes+(num_all_nodes))->neighbors + 1) = (nodes + *(active_node_map+1));
                        *((nodes+(num_all_nodes))->neighbors + 2) = (nodes + *(active_node_map+2));

                        *(*(distances+num_all_nodes)+*(active_node_map+0)) = first;
                        *(*(distances+num_all_nodes)+*(active_node_map+1)) = *(*(distances + *(active_node_map+0)) + *(active_node_map+1)) - first;
                        *(*(distances+num_all_nodes)+*(active_node_map+2)) = *(*(distances + *(active_node_map+0)) + *(active_node_map+2)) - first;
                        *(*(distances+*(active_node_map+0))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+0));
                        *(*(distances+*(active_node_map+1))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+1));
                        *(*(distances+*(active_node_map+2))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+2));
                        break;
                    }
                }
            }

            *(active_node_map+min_i) = num_all_nodes;
            if(min_j != num_taxa-1){
                *(active_node_map+min_j) = *(active_node_map+num_taxa-1);
            } else{
                *(active_node_map+min_j) = *(active_node_map+num_taxa-2);
            }
            num_all_nodes++;
            num_active_nodes += 2;
        }

        if(global_options == MATRIX_OPTION){
            if(out != NULL){
                for(int i = 0; i < numtemp; i++){
                    fprintf(out, ",%s", *(node_names + i));
                }
                for(int i = numtemp; i < num_all_nodes; i++){
                    fprintf(out, ",%c", *(*(node_names + i)));
                    fprintf(out, "%d", (*(*(node_names + i)+1)-48));
                }

                fprintf(out, "\n");
                for(int i = 0; i < num_all_nodes; i++){
                    if(stringcomp1(*(node_names + i), *(node_names + i)) == 1){
                        if(i < numtemp){
                            fprintf(out,"%s", *(node_names + i));
                        } else{
                            fprintf(out, "%c", *(*(node_names + i)));
                            fprintf(out, "%d", (*(*(node_names + i)+1)-48));
                        }
                    }
                    for(int j = 0; j < num_all_nodes; j++){
                        fprintf(out, ",%.2f", *(*(distances+i)+j));
                    }
                    fprintf(out, "\n");
                }
                fprintf(out, "\n");
            }
        }
        return EXIT_SUCCESS;
    } else {
        while(num_taxa>3){
            num_taxa = num_all_nodes - num_active_nodes;
            calculateQ1(&min_i, &min_j);
            *(*(node_names + num_all_nodes)) = '#';
            *(*(node_names + num_all_nodes)+1) = (char)(num_all_nodes+48);
            (nodes + num_all_nodes)->name = *(node_names + num_all_nodes);

            // Add a new node
            for (int k = 0; k < num_taxa; k++) {
                double limb_i = 0.5 * *(*(distances + *(active_node_map+min_i)) +
                *(active_node_map+min_j)) +  1.0 / (2.0*(num_taxa - 2.0)) * (sumRow(min_i) - sumRow(min_j));
                double limb_j = *(*(distances + *(active_node_map+min_i)) + *(active_node_map+min_j)) - limb_i;

                if(k == num_all_nodes){
                    *(*(distances+num_all_nodes)+num_all_nodes) = 0;
                    *(*(distances+*(active_node_map+k))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+k));
                } else if(*(active_node_map+k) == *(active_node_map+min_i) && num_taxa != 3 ){
                    fprintf(out, "%d,%d,%.2f\n",*(active_node_map+k),num_all_nodes,limb_i);
                    *(*(distances+num_all_nodes)+*(active_node_map+k)) = limb_i;
                    *(*(distances+*(active_node_map+k))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+k));

                } else if(*(active_node_map+k) == *(active_node_map+min_j) && num_taxa != 3){
                    fprintf(out, "%d,%d,%.2f\n",*(active_node_map+k),num_all_nodes,limb_j);
                    *(*(distances+num_all_nodes)+*(active_node_map+k)) = limb_j;
                    *(*(distances+*(active_node_map+k))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+k));

                } else {
                    if(num_taxa != 3){
                        *(*(distances+num_all_nodes)+*(active_node_map+k)) = (*(*(distances+*(active_node_map+min_i))+*(active_node_map+k)) +
                            *(*(distances+*(active_node_map+min_j))+*(active_node_map+k)) -
                            *(*(distances+*(active_node_map+min_i))+*(active_node_map+min_j))) / 2;

                        *(*(distances+*(active_node_map+k))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+k));
                    } else if(num_taxa == 3){
                        double first = 0.5 * *(*(distances + *(active_node_map+0)) + *(active_node_map+1)) + 1.0 / (2.0*(num_taxa - 2.0)) * (sumRow(0) - sumRow(1));
                        fprintf(out, "%d,%d,%.2f\n",*(active_node_map+0),num_all_nodes,first);
                        fprintf(out, "%d,%d,%.2f\n",*(active_node_map+1),num_all_nodes,*(*(distances + *(active_node_map+0)) + *(active_node_map+1)) - first);
                        fprintf(out, "%d,%d,%.2f\n",*(active_node_map+2),num_all_nodes,*(*(distances + *(active_node_map+0)) + *(active_node_map+2)) - first);
                        *(*(distances+num_all_nodes)+*(active_node_map+0)) = first;
                        *(*(distances+num_all_nodes)+*(active_node_map+1)) = *(*(distances + *(active_node_map+0)) + *(active_node_map+1)) - first;
                        *(*(distances+num_all_nodes)+*(active_node_map+2)) = *(*(distances + *(active_node_map+0)) + *(active_node_map+2)) - first;
                        *(*(distances+*(active_node_map+0))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+0));
                        *(*(distances+*(active_node_map+1))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+1));
                        *(*(distances+*(active_node_map+2))+num_all_nodes) = *(*(distances+num_all_nodes)+*(active_node_map+2));
                        break;
                    }
                }
            }

            *(active_node_map+min_i) = num_all_nodes;
            if(min_j != num_taxa-1){
                *(active_node_map+min_j) = *(active_node_map+num_taxa-1);
            } else{
                *(active_node_map+min_j) = *(active_node_map+num_taxa-2);
            }

            num_all_nodes++;
            num_active_nodes += 2;
        }
    }
    return 0;
}