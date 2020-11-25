#ifndef __OPTION_H__
#define __OPTION_H__

#include <stdio.h>
#include "uthash.h"
#include "sam.h"
#include "fastq.h"
#include "nuc.h"
#include "qual.h"

#define N_FILES 2

typedef struct _options options;


struct _options {
	unsigned char drop_unmapped;	/*<! drop reads that are unmapped
					 * in either alignment
					 */
	unsigned char drop_secondary;	/*<! drop secondary alignments */
	unsigned char display_alignment;/*<! display alignments to stderr */
	int drop_soft_clipped;		/*<! drop reads soft-clipped by this
					 * length of more in either alignment
					 */
	int drop_indel;			/*<! drop reads whose alignments
					 * contain indel this long or more
					 * in either alignment
					 */
//	int proptest_screen;		/*<! drop reads according to abundance test */
//	double coverage_screen;		/*<! drop reads whose coverage */
//	double min_expected_coverage;	/*<! abort if subgenome coverage < this */
//	int posthoc_coverage_test;	/*<! perform post hoc coverage test */
//	double min_genotype_post_prob;	/*<! min. posterior probability to call heterozygote */
//	double min_alignment_post_prob;	/*<! min. posterior probability of alignment */
	int min_length;			/*<! drop reads shorter than this */
	int max_length;			/*<! drop reads longer than this */
	double max_eerr;		/*<! drop reads with more expected
					 * errors than this
					 */
//	double min_log_likelihood;	/*<! drop reads smaller log likelihood */
//	int max_quality_score;		/*<! censory quality scores here */
	unsigned int n_sample;		/*<! number of Monte Carlo samples */
	const char *out_file;		/*<! out_file */
	const char *uni_geno_file;	/*<! uni_geno_file */
	FILE *error_file;		/*<! error data file */
	unsigned int use_bam;		/*<! files are in bam format */
	const char *sbam_files[N_FILES];/*<! sam/bam files */
	const char *fsa_files[N_FILES];	/*<! fsa files */
	const char *ref_names[N_FILES];	/*<! name of references */
	char const * sam_file;		/*<! reference sam file */
	const char *uni_genome;		/*<! uni_genome */
    const char *selected_fq;        /*<! selected reads in a fastq */
    const char *splited_fq[N_FILES];        /*<! selected reads in a splited fastq */
};

typedef struct nuc_state {
	UT_hash_handle hh;		/* 56 bytes w/ pointer alignment */
	double prob;			/* 8 bytes */
	unsigned int pos;		/* 4 */
	data_t true_nuc, read_nuc;	/* 1, 1 bytes */
	data_t quality;			/* 1 */
	/* total: 71 bytes round up to
	 * 72 with 1 byte slop
	 */
} nuc_state;

typedef struct mlogit_stuff {
	nuc_state *ns;
	int pos;
} mlogit_stuff;

int default_options(options *opt);
int parse_options(options *opt, int argc, const char **argv);
void fprint_usage(FILE *fp, const char *cmdname, void *obj);
#endif
