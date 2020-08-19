
#ifndef make_aln_h
#define make_aln_h

#include <stdio.h>
#include <curses.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>

#include "sam.h"
#include "fastq.h"
#include "nuc.h"
#include "qual.h"
#include "uthash.h"
#include "io.h"
#include "array.h"
#include "pick_reads.h"
#include "options.h"
#include "myfun.h"

int make_alignment(options opt);
#endif /* make_aln_h */
