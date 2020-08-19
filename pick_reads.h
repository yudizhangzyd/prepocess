//
//  pick_reads.h
//  project
//
//  Created by Yudi Zhang on 3/24/20.
//  Copyright © 2020 Yudi Zhang. All rights reserved.
//

#ifndef pick_reads_h
#define pick_reads_h

#include <stdint.h>
#include "sam.h"
#include "options.h"
#include "array.h"

typedef struct _ref_info ref_info;
typedef struct _ref_entry ref_entry;
typedef struct _ref_options options_rf;

struct _ref_options {
	char const * sam_file;			// reference sam file
	const char * rsam_files[N_FILES]; 	// read sam files
	const char *samtools_command;
	const char * extracted_rf[N_FILES]; 	// targeted reference fsa files
	const char * fsa_files[N_FILES];	// reference fsa files
	char filter_unmapped;
	char *delim_ref;
	char *delim_len;
	const char *fastq_file;
};

struct _ref_entry {
	char *name_A;
	char *name_B;
	size_t start_A;
	size_t start_B;
	size_t end_A;
	size_t end_B;
	unsigned int strand_A;  	/* 0:forward, 1:reverse */
	unsigned int strand_B;
	int *idx_map;		 	/* which base in B is aligned to which in A  */
};


struct _ref_info {
	ref_entry *info;
	sam *ref_sam;
};

void ll_all_align(sam *sds, const char *ref_file);
void make_options(options_rf *opt);
int make_targets_info(options_rf opt, ref_info **ref_info);
int pickreads(ref_info *ref_info, options_rf *opt, sam **sds);
int parse_rf_options(options_rf *opt, int argc, char *argv[]);
void extract_ref(const char *samtools_command, char *region, const char *ref_file, const char *ext_rf);
void output_selected_reads(const char *f, sam **sds, merge_hash *mh);
void output_data(FILE *fp, sam_entry *se, unsigned int id);
//void match_pair(ref_info *rf_info, size_t my_refs);
double sub_prob_given_q_with_encoding(data_t s, data_t r, int es, int er, data_t q, int logged, void *vptr);
void adjust_alignment(sam_entry *se, data_t *ref, unsigned int strand, int *id_uni, size_t uni_aln_len, long *real_id);
#endif /* pick_reads_h */
