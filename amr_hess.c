/* This program takes input a .ele file or a .node files with switches -e and -n resp*/
/* There is a small bug, the name of the file without extension should be 3 or more characters length*/
// Here we are marking at the beginning. So the newly formed elements are not marked and hence may cause to exit the program
#ifdef SINGLE
#define REAL float
#else
#define REAL double
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include "mmap_lib.h"
/* Maximum number of characters in a file name (including the null).*/
#define FILENAMESIZE 2048
#define INPUTLINESIZE 1024

/* structure to represent triangles*/
struct element{
	int ele_no,nodes[3],mark_flag;
	int neighbour[3];
	struct element *next;
}
