#include <stdio.h>
#include <curses.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>

#include "make_aln.h"

#ifndef STANDALONE
#include <Rinternals.h>
#define PRINTF(str, ...) Rprintf((str), __VA_ARGS__)
#define EPRINTF(str, ...) REprintf((str), __VA_ARGS__)
#else
#define PRINTF(str, ...) fprintf(stdout, (str), __VA_ARGS__)
#define EPRINTF(str, ...) fprintf(stderr, (str), __VA_ARGS__)
#endif

#ifdef STANDALONE
int main(int argc, const char *argv[])
{
	options opt;
	default_options(&opt);
	if (parse_options(&opt, argc, argv))
		exit(mmessage(ERROR_MSG, INVALID_CMDLINE, ""));
//	int debug_level = QUIET;//ABSOLUTE_SILENCE;//MINIMAL;//DEBUG_I;//
	make_alignment(opt);
}
#endif

// if the alignment has hard/soft clips, then the current code my not work, need to discrad the reads that align to the hard/soft clips regions
int make_alignment(options opt) {
	int err = NO_ERROR;
	sam *sds[N_FILES] = {NULL, NULL};
	fastq_options fop = {.read_encoding = IUPAC_ENCODING, .read_names = 1};
	sam_hash *by_name[N_FILES] = {NULL, NULL};
	size_t my_refs[N_FILES] = {0, 0};
	FILE *fp = NULL;
	unsigned int j, i;
	options_rf opt_rf;
	make_options(&opt_rf);
	opt_rf.sam_file = opt.sam_file;
	ref_info *rf_info = NULL;
	fastq_data *fdr = NULL;
	
    if ((err = read_fastq(opt.uni_genome, &fdr, &fop))) {
#ifdef STANDALONE
        exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Reading '%s' "
        "failed with error '%s' (%d).\n",
        opt.uni_genome, fastq_error_message(err),
        err));
#else
        error("Reading '%s' failed with error '%s' (%d).\n",
        opt.uni_genome, fastq_error_message(err),
        err);
#endif
    }
		
	// get the selected reference (match the name)
	unsigned int A_id = 0, B_id = 0;
	unsigned int rptr = 0, rptr_b = 0;
	char *A_name = NULL;
    char *B_name = NULL;
	char *names = NULL;
//		printf("%lu\n", strlen(opt.ref_names[0]));
	//	strlen has to be the same when comparing
//
	if (fdr->n_max_length == fdr->n_min_length) {
		fdr->n_lengths = malloc(fdr->n_reads * sizeof(*fdr->n_lengths));
		for (j = 0; j < fdr->n_reads; ++j)
			fdr->n_lengths[j] = fdr->n_min_length;
	}
	
	for (j = 0; j < fdr->n_reads; ++j) {
		names = fdr->names;
		for (i = 0; i < j; ++i)
			names += fdr->name_lengths[i];
		A_name = malloc((fdr->name_lengths[j] + 1) * sizeof(*A_name));
        B_name = malloc((fdr->name_lengths[j + 1] + 1) * sizeof(*B_name));
		//		for (i = 0; i < fdr->name_lengths[j]; ++i)
		//			A_name[i] = names[i];
		strncpy(A_name, names, fdr->name_lengths[j]);
		A_name[fdr->name_lengths[j]] = '\0';
        strncpy(B_name, names + fdr->name_lengths[j], fdr->name_lengths[j + 1]);
        B_name[fdr->name_lengths[j + 1]] = '\0';
//				printf("%lu %s\n", strlen(A_name), A_name);
		if (!strcmp(opt.ref_names[0], A_name) && !strcmp(opt.ref_names[1], B_name)) {
			A_id = j;
            B_id = j+1; // note here, it is possible that A aligns to multiple B, better to match by B name
			rptr_b = rptr + fdr->n_lengths[j];
			break;
		}
		rptr += fdr->n_lengths[j];
	}
    free(A_name);
    free(B_name);
	// make universal genome, gap from A as I(4), gap from B as J(5), mismatch mas M (6)
	data_t to_xy[NUM_IUPAC_SYMBOLS] = {
		0, XY_A, XY_C, 0, XY_G, 0, 0, 0, XY_T, 0, 0, 0, 0, 0, 0, 0
	};
	data_t *uni_genome = malloc(fdr->n_lengths[A_id] * sizeof(*uni_genome));
	
	for (j = 0; j < fdr->n_lengths[A_id]; ++j) {
		if (fdr->reads[rptr + j] == fdr->reads[rptr_b + j]) {
			uni_genome[j] = to_xy[fdr->reads[rptr + j]];
		} else if (fdr->reads[rptr + j] != 0 && fdr->reads[rptr_b + j] != 0 &&
			   fdr->reads[rptr + j] != fdr->reads[rptr_b + j]) {
			uni_genome[j] = 6;
		} else if (fdr->reads[rptr + j] != 0 && fdr->reads[rptr_b + j] == 0) {
			uni_genome[j] = 5;
		} else {
			uni_genome[j] = 4;
		}
	}
	FILE *uni_fa = NULL;
	if (opt.uni_geno_file) {
		uni_fa = fopen(opt.uni_geno_file, "w");
        if (!uni_fa) {
#ifdef STANDALONE
            exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt.uni_geno_file));
#else
            error("%s open error", opt.uni_geno_file);
#endif
        }
			
		// output the universal alignment
		fprintf(uni_fa, ">uni_genome\n");
		for (j = 0; j < fdr->n_lengths[A_id]; ++j) {
			if (uni_genome[j] == 5) {
				fprintf(uni_fa, "J"); // GAP IN B
			} else if (uni_genome[j] == 4) {
				fprintf(uni_fa, "I"); // GAP IN A
			} else
				fprintf(uni_fa, "M");
		}
		fprintf(uni_fa, "\n");
		fclose(uni_fa);
	}
//		PRINT_VECTOR(uni_genome, fdr->n_lengths[A_id]);
	// store the index after alignment, if gap use -1, this alignment does not contain softclips, when map the reads back, the start and position should consider its length
	unsigned int gap_a = 0;
	unsigned int gap_b = 0;
	int *id_A = malloc(fdr->n_lengths[A_id] * sizeof(*id_A));
	int *id_B = malloc(fdr->n_lengths[A_id] * sizeof(*id_B));
	for (j = 0; j < fdr->n_lengths[A_id]; ++j) {
		id_A[j] = j - gap_a;
		id_B[j] = j - gap_b;
		if (uni_genome[j] == 4) {
			gap_a++;
			id_A[j] = -1;
		} else if (uni_genome[j] == 5) {
			gap_b++;
			id_B[j] = -1;
		}
	}
    free(uni_genome);
//    PRINT_VECTOR(id_A, fdr->n_lengths[A_id]);
//    printf("B\n");
//    PRINT_VECTOR(id_B, fdr->n_lengths[A_id]);
    mmessage(INFO_MSG, NO_ERROR, "make_targets_info\n");
	/* store information to the reference targeted sam file */
	make_targets_info(opt_rf, &rf_info);
	
	// read in the alignments to A B
	for (j = 0; j < N_FILES; ++j) {
		fp = fopen(opt.sbam_files[j], "r");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR,
				      opt.sbam_files[j]));
		read_sam(fp, &sds[j]);
		fclose(fp);
	}
	fp = NULL;
    mmessage(INFO_MSG, NO_ERROR, "pickreads\n");
	// pick the reads, more reads could be picked
	pickreads(rf_info, &opt_rf, sds);
	
	char *strand;
	for (j = 0; j < N_FILES; ++j) {
		/* find selected references index in sam files */
		size_t rchar = 0;
		unsigned char found = 0;
		for (size_t i = 0; i < sds[j]->n_se; ++i) {
			sam_entry *se = &sds[j]->se[i];
			
			/* skip unmapped */
			if ((se->flag & (1 << 2)))
				continue;
			if (se->which_ref == -1)
				continue;
			
			if (!strcmp(opt.ref_names[j], se->ref_name)) {
				se->name_s = NULL;
				/* strand for hashing on strand and name */
				if ((se->flag & 16) == 0) {
					strand = "+";
					// flip the strand if A is aligned to reverse complement of B
					if (j == 1 && rf_info->info[se->which_ref].strand_B == 1)
						strand = "-";
				} else {
					strand = "-";
					if (j == 1 && rf_info->info[se->which_ref].strand_B == 1)
						strand = "+";
				}
				
				size_t length = strlen(se->name) + strlen(strand) + 1;
				se->name_s = malloc(length);
				sprintf(se->name_s, "%s%s", se->name, strand);
				my_refs[j] = se->which_ref; // this my_refs index should be adjusted
				found = 1;
//				printf("%s %lu\t", se->name_s, se->pos);
			}
			rchar += strlen(se->ref_name) + 1;
		}
//		printf("\n");
        if (!found) {
#ifdef STANDALONE
            exit(mmessage(ERROR_MSG, INVALID_USER_INPUT, "no "
            "reference '%s' in fasta file '%s'",
            opt.ref_names[j], opt.sbam_files[j]));
#else
            error("no reference '%s' in fasta file '%s'",
            opt.ref_names[j], opt.sbam_files[j]);
#endif
        }
			
		/* hash sam file to reference (use n_se since some references are repeated in the targted sam file) */
		hash_sam(sds[j], &by_name[j], HASH_REFERENCE, my_refs[j], rf_info->ref_sam->n_se,
			 opt.drop_unmapped, opt.drop_secondary,
			 opt.drop_soft_clipped, opt.drop_indel,
			 opt.min_length, opt.max_length, opt.max_eerr);
		
		mmessage(INFO_MSG, NO_ERROR, "Number of %u alignments: %zu\n",
			 j, sds[j]->n_per_ref[my_refs[j]]);
	}
	
	//merge read pair, if read only align to one genome, keep the alignment information, but do not pick the best alignment
	merge_hash *mh = NULL;
	size_t total_read = hash_merge(&mh, N_FILES, sds, my_refs);
#ifdef STANDALONE
	printf("total picked reads %lu\n", total_read);
#endif
    
	// 0-based, so plus 1 to match with what the rest related code designed for
	sam_entry *fse = &rf_info->ref_sam->se[my_refs[0]];
	ref_entry *re = &rf_info->info[my_refs[0]];
	re->start_B++;
	//	re->end_B;
	re->start_A++;
	//	re->end_A;
	// store aligned index used in the real genome
	
	long *real_id_A = malloc(fdr->n_lengths[A_id] * sizeof(*real_id_A));
	long *real_id_B = malloc(fdr->n_lengths[A_id] * sizeof(*real_id_B));
	
	real_id_A[0] = re->start_A + fse->pos - 1; // A start index, 1 based
	if (re->strand_B) { // if B reversed
#ifdef STANDALONE
		printf("Genome B is reverse complemented\n");
#endif
		real_id_B[0] = re->end_B - 1;
		if (fse->cig->ashes[0].type == CIGAR_SOFT_CLIP || fse->cig->ashes[0].type == CIGAR_HARD_CLIP)
			real_id_B[0] -= fse->cig->ashes[0].len;
	} else {
		real_id_B[0] = re->start_B;
		if (fse->cig->ashes[0].type == CIGAR_SOFT_CLIP || fse->cig->ashes[0].type == CIGAR_HARD_CLIP)
			real_id_B[0] += fse->cig->ashes[0].len; // length of unaligned in B
	}
	
	for (j = 1; j < fdr->n_lengths[A_id]; ++j) {
		if (id_A[j] == -1) {
			real_id_A[j] = -1;
			if (re->strand_B) // if B is reversed
				real_id_B[j] = real_id_B[0] - id_B[j];
			else
				real_id_B[j] = real_id_B[0] + id_B[j];
			
		} else if (id_B[j] == -1) {
			real_id_B[j] = -1;
			real_id_A[j] = real_id_A[0] + id_A[j];
		} else {
			real_id_A[j] = real_id_A[0] + id_A[j];
			if (re->strand_B) // if B is reversed
				real_id_B[j] = real_id_B[0] - id_B[j];
			else
				real_id_B[j] = real_id_B[0] + id_B[j];
		}
	}
//#ifdef STANDALONE
//	for (j = 1; j < fdr->n_lengths[A_id]; ++j)
//		printf("%ld: %c\t\t\t", real_id_A[j], iupac_to_char[fdr->reads[rptr + j]]);
//	for (j = 1; j < fdr->n_lengths[A_id]; ++j)
//		printf("%ld: %c\t", real_id_B[j], iupac_to_char[fdr->reads[rptr_b + j]]);
//	printf("\nreal start-end in genome A %ld-%ld B %ld-%ld\n", real_id_A[0], real_id_A[fdr->n_lengths[A_id] - 1], real_id_B[0], real_id_B[fdr->n_lengths[A_id] - 1]);
//	PRINT_VECTOR(real_id_A, fdr->n_lengths[A_id]);
//	printf("B\n");
//	PRINT_VECTOR(real_id_B, fdr->n_lengths[A_id]);
//#endif
	unsigned int strand_genome;
	int *id_uni = NULL;
	//	size_t uni_len;
	long *real_id = NULL;
	unsigned int rf_id = 0;
	
	//	for (i = 0; i < fdr->n_lengths[A_id]; ++i) {
	//		printf("%d:%d:%c:%ld\t", i, id_A[i], iupac_to_char[fdr->reads[i + rptr]], real_id_A[i]);
	//	}
	//	printf("\n");
	//	for (i = 0; i < fdr->n_lengths[B_id]; ++i) {
	//		printf("%d:%d:%c:%ld\t", i, id_B[i], iupac_to_char[fdr->reads[i + rptr_b]], real_id_B[i]);
	//	}
	//	printf("\n");
//	printf("start!\n");
	if (opt.out_file) {
		fp = fopen(opt.out_file, "w");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt.out_file));
	}
    FILE *splitted[2] = {NULL, NULL};
    if (opt.splited_fq[0])
        for (i = 0; i < N_FILES; ++i) {
            splitted[i] = fopen(opt.splited_fq[i], "w");
            if (!splitted[i])
                exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt.splited_fq[i]));
        }
    
	unsigned int n_read = 1;
	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
		sam_entry *se;
		
		// only mapped to one reference, adjust the alignment to universal
		if (me->nfiles != N_FILES) {
			me->exclude = 1;
//			if (me->indices[0]) {
//				se = &sds[0]->se[me->indices[0][0]]; // indices[0] represents A
//				strand_genome = re->strand_A;
//				id_uni = id_A;
//				rf_id = rptr;
//				//				uni_len = fdr->n_lengths[A_id] - gap_a;
//				real_id = real_id_A;
//			} else {
//				se = &sds[1]->se[me->indices[1][0]];
//				strand_genome = re->strand_B;
//				id_uni = id_B;
//				rf_id = rptr_b;
//				//				uni_len = fdr->n_lengths[A_id] - gap_b;
//				real_id = real_id_B;
//			}
#ifdef STANDALONE
            mmessage(WARNING_MSG, NO_ERROR, "Read %u aligns once.\n",
                                      sds[j]->se[me->indices[j][0]].name_s);
#endif
//			adjust_alignment(se, &fdr->reads[rf_id], strand_genome, id_uni, fdr->n_lengths[A_id], real_id);
//			// output the final alignment
//			output_data(fp, se, n_read);
//			n_read++;
			continue;
		}
		
		/* force one alignment per sub-genome */
		double max_ll = -INFINITY;
		unsigned int max_id = 0;
		for (j = 0; j < N_FILES; ++j) {
			//			printf("genome %d\n", j);

            if (me->count[j] > 1) { // if reads align to multiple places, skip this read
#ifdef STANDALONE
                mmessage(WARNING_MSG, NO_ERROR,
                                      "Read %u aligns twice in genome %s.\n",
                                      sds[j]->se[me->indices[j][0]].name_s, j);
#else
                warning("Read %u aligns twice in genome %s.\n",
                                      sds[j]->se[me->indices[j][0]].name_s, j);
#endif
                me->exclude = 1;
                break;
            }
            
			se = &sds[j]->se[me->indices[j][0]];
			if (j == 0) {
				strand_genome = re->strand_A;
				id_uni = id_A;
				rf_id = rptr;
				//				uni_len = fdr->n_lengths[A_id] - gap_a;
				real_id = real_id_A;
			} else {
				strand_genome = re->strand_B;
				id_uni = id_B;
				rf_id = rptr_b;
				//				uni_len = fdr->n_lengths[A_id] - gap_b;
				real_id = real_id_B;
			}
#ifdef STANDALONE
			printf("%s is adjusting\n", se->name_s);
#endif
			int exclude = adjust_alignment(se, &fdr->reads[rf_id], strand_genome, id_uni, fdr->n_lengths[A_id], real_id);
            if (exclude) {
                me->exclude = 1;
#ifdef STANDALONE
                mmessage(WARNING_MSG, NO_ERROR, "Read %u aligns wrong.\n",
                         se->name_s);
#endif
                break;
            }
            // if the same, then ramdom sample
			if (max_ll < se->ll_aln) {
				max_ll = se->ll_aln;
				max_id = j;
			}
            else if (max_ll == se->ll_aln) {
                double which_ref = rand() / (RAND_MAX + 1.);
                if (which_ref <= 0.5)
                    max_id = j;
            }
		}
#ifdef STANDALONE
		printf("choose %d\n", max_id);
#endif
        if (me->exclude == 1)
            continue;
        
		se = &sds[max_id]->se[me->indices[max_id][0]];
//        if(se->gap_in) { // if there is gap on one alignment but not in another, then make up this info
//            int left = 0;
//            if (!max_id)
//                left = 1;
//            sam_entry *se2 = &sds[left]->se[me->indices[left][0]];
//            printf("\nfill gap: %zu %zu\n", se->aln_len, se2->aln_len);
//            for (size_t m = 0; m < se->aln_len; ++m) {
//                if (se->uni_aln[m] == 4) {
//                    se->uni_aln[m] = se2->uni_aln[m];
//                    printf("%d:%d ",  se2->uni_aln[m], se->uni_aln[m]);
//                }
//            }
//            fprintf(stderr, "\n");
//        }
		// output the needed format
		output_data(fp, se, n_read, &opt_rf);
        // output the sam file, split A and B
        if (opt.splited_fq[0]) {
            fprintf(splitted[max_id], "@%s\n", se->name);
            fwrite_nuc_segment(splitted[max_id], se->read, XY_ENCODING, 0,
                       se->read->len);
            fprintf(splitted[max_id], "\n+\n");
            fwrite_qual_sequence(splitted[max_id], se->qual);
            fprintf(splitted[max_id], "\n");
        }
		n_read++;
	}
    
    // output selected reads
    if (opt.selected_fq)
        output_selected_reads(opt.selected_fq, sds, mh);
    
#ifdef STANDALONE
	printf("read: %d\n", n_read);
#endif
	fclose(fp);
    if (opt.splited_fq[0])
        for (i = 0; i < N_FILES; ++i)
            fclose(splitted[i]);
    free(real_id_A);
    free(real_id_B);
	return err;
}


void output_sam(FILE *fp, sam_entry *se, long ref_len) {
	fprintf(fp, "@SQ	SN:");
	fprintf(fp, "%s	LN:%lu\n", se->ref_name, ref_len);
	fprintf(fp, "%s	%d	%s	%zu	60	", se->name, se->flag, se->ref_name, se->pos); // give mapping quality 60 (GATK does not depend on MAPQ)
	cigar *cig = se->cig;
	for (unsigned int l = 0; l < cig->n_ashes; ++l) {
		fprintf(fp, "%d", cig->ashes[l].len);
		if (cig->ashes[l].type == CIGAR_DELETION)
			fprintf(fp, "D");
		else if (cig->ashes[l].type == CIGAR_MATCH ||
			 cig->ashes[l].type == CIGAR_MISMATCH ||
			 cig->ashes[l].type == CIGAR_MMATCH)
			fprintf(fp, "M");
		else if (cig->ashes[l].type == CIGAR_SKIP)
			fprintf(fp, "N");
		else if (cig->ashes[l].type == CIGAR_PAD)
			fprintf(fp, "P");
		else if (cig->ashes[l].type == CIGAR_INSERTION)
			fprintf(fp, "I");
		else if (cig->ashes[l].type == CIGAR_SOFT_CLIP)
			fprintf(fp, "S");
		else if (cig->ashes[l].type == CIGAR_HARD_CLIP)
			fprintf(fp, "H");
	}
	fprintf(fp, "	*	0	0	");
	fwrite_nuc_sequence(fp, se->read, XY_ENCODING);
	fprintf(fp, "	");
	fwrite_qual_sequence(fp, se->qual);
	fprintf(fp, "\n");
}
