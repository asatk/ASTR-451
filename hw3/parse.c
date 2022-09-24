#include <math.h>
#include <regex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAXLINE 1 << 16

typedef struct AList {
    double mcol;
    double temp;
    struct AList *tail;
} *LList;

static void substring(char *s, char *buf, int b, int e);
double parsenum(char *exp);
static LList create_list(double mcol, double temp, LList tail);
static LList insert_list(LList head, double mcol, double temp);
static void delete_list(LList list);
static void print_list(LList list);
int parsekurucz(char *f, double **arrayptr);

/**
 * Determine the substring of a string from an initial position to final
 * position in the target string. The result is stored in 'buf'.
 * 
 * s    target string
 * buf  destination buffer for substring
 * b    beginning position/byte offset of substring in target string
 * e    ending position/byte offset of substring in target string
 */
static void substring(char *s, char *buf, int b, int e) {
    int i;
    char c, *ptr;

    for(i = 0; b + i < e; i++) {
        c = *(s + b + i * sizeof(char));
        ptr = buf + i;
        *ptr = c;
    }
    buf[i] = '\0';
}

/**
 * Parse a string in decimal format or scientific notation for its double
 * representation. Returns '0.0' for a string that does no match.
 * 
 * exp  string containing numeric expression
 */
double parsenum(char *exp) {

    double val;
    char *buf1, *buf2;
    regex_t re;
    regmatch_t match[4];

    regcomp(&re, "^(-?[0-9]+[.]?[0-9]*)(e(-*[0-9]+))?$", REG_EXTENDED);

    if (regexec(&re, exp, 4, match, 0) == REG_NOMATCH) {
        // No match for the regex exp: return 0.0
        printf("value does not match this pattern:\n" \
            "^([0-9]+[.]?[0-9]*)(e(-*[0-9]+))?$\n");
        regfree(&re);
        return 0.;
    } else {
        // The string matches the regex exp: derive decimal value from subexps
        
        // Obtain coefficient from the scietific notation representation
        buf1 = (char *) malloc(match[1].rm_eo - match[1].rm_so);
        substring(exp, buf1, match[1].rm_so, match[1].rm_eo);
        
        // Obtain exponent from the scientific notation representation if exists
        if (match[2].rm_eo > match[2].rm_so) {
            buf2 = (char *) malloc(match[3].rm_eo - match[3].rm_so);
            substring(exp, buf2, match[3].rm_so, match[3].rm_eo);
        } else {
            buf2 = (char *) malloc(2);
            strcpy(buf2, "0\0");
        }

        // Calculate value from coefficient and exponent fields
        val = atof(buf1) * pow(10.0, (double) atof(buf2));
        
        free(buf1);
        free(buf2);
    }

    regfree(&re);
    return val;
}

static LList create_list(double mcol, double temp, LList tail) {
    LList head = malloc(sizeof(struct AList));
    head->mcol = mcol;
    head->temp = temp;
    head->tail = tail;
    return head;
}

/**
 * Inserts list at the head of the input list. returns the current jead of the list.
 */
static LList insert_list(LList head, double mcol, double temp) {
    return create_list(mcol, temp, head);
}

/**
 * Deletes the entire list from the provided node onwards. Tail recursion.
 */
static void delete_list(LList list) {
    if (list == NULL)
        return;  
    if (list->tail != NULL)
        delete_list(list->tail);

    free(list);
}

/**
 * Print the data of the linked list starting with list as the head.
 */
static void print_list(LList list) {
    LList p = list;
    int count = 0;
    while(p != NULL) {
        printf("(%lf, %lf) -> ",p->mcol,p->temp);
        p = p->tail;
        if(++count%5 == 0) {
            printf("\n");
        }
    }
    printf("|\n");
    fflush(stdout);
}

/**
 * Parses kurucz file for data. Returns a double array of values, the first
 * index for the row and the second for the column (0: mass col, 1: temp).
 * 
 * f - string name for file being parsed
 */
int parsekurucz(char *f, double **arrayptr) {
    
    regex_t regp;
    regmatch_t regm[3];

    int nrows, i;
    char line[MAXLINE], *buf1, *buf2;
    double mcoli, tempi, *p;
    LList head;

    FILE *file = fopen(f, "r");
    if (file == NULL) {
        printf("cannot open %s\n",f);
    }

    regcomp(&regp, "[[:space:]]*([0-9.]+)[[:space:]]+([0-9.]+)", REG_EXTENDED);

    nrows = 0;
    head = NULL;
    // Iterate through each line in the file
    while(fgets(line, MAXLINE, file) != NULL) {
        // Determine if the line has data in the desired representation
        if(!regexec(&regp, line, 3, regm, 0)) {
            nrows += 1;

            // Pare the numeric data substrings from the line string
            buf1 = (char *) malloc(regm[1].rm_eo - regm[1].rm_so);
            buf2 = (char *) malloc(regm[2].rm_eo - regm[2].rm_so);
            substring(line, buf1, regm[1].rm_so, regm[1].rm_eo);
            substring(line, buf2, regm[2].rm_so, regm[2].rm_eo);
            
            // Convert the substrings into floating-point numbers
            mcoli = parsenum(buf1);
            tempi = parsenum(buf2);

            // Append the data into a linked list, aggregating the file's data
            head = insert_list(head, mcoli, tempi);
        }
    }

    // Check that the end of the file was reached while reading it.
    if (!feof(file))
        printf("reading (first pass) error in %s\n",f);

    // Create array that will store the data
    *arrayptr = (double *) malloc(2 * nrows * sizeof(double));

    // Convert the linked list into a 2D array (pointer)
    p = *arrayptr;
    i = nrows - 1;
    for(i = nrows - 1; i >= 0; i--) {
        // Assign entry's 1st column - mcol
        p = *arrayptr + 2 * i;
        *p = head->mcol;
        // Assign entry's 2nd column - temp
        p += 1;
        *p = head->temp;

        // Move to the next element in the reverse-ordered linked list
        head = head->tail;
    }

    fclose(file);
    delete_list(head);

    return nrows;
}

// Prevent "unused function" and "unused variable" warnings.
static const void *dummy_ref[] = {print_list, dummy_ref};