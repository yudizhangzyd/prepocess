//
//  pick_reads.c
//  pick reads aligned to each targeted region
//
//  Created by Yudi Zhang on 3/24/20.
//  Copyright © 2020 Yudi Zhang. All rights reserved.
//

#include <stdlib.h>
#include <string.h>

#include "fastq.h"
#include "nuc.h"
#include "qual.h"
#include "uthash.h"
#include "io.h"
#include "cmdline.h"
#include "error.h"
#include "pick_reads.h"
#include "myfun.h"

void make_options(options_rf *opt)
{
	opt->sam_file = NULL;
	opt->filter_unmapped = 1;
	opt->delim_ref = ":";
	opt->delim_len = "-";
	opt->samtools_command = "samtools";
	for (int i = 0; i < N_FILES; ++i) {
		opt->rsam_files[i] = NULL;
		opt->fsa_files[i] = NULL;
	}
	opt->extracted_rf[0] = "ref_A.fasta";
	opt->extracted_rf[1] = "ref_B.fasta";
	opt->fastq_file = NULL;
	
} /* make_options */

/**
* Read sam file containing alignment(s) of homoeologous region(s) of
* subgenomic references. Record names of aligned references, start
* and end locations of regions along genome (0-based), start and end
* locations of aligned region (1-based). These are stored in ref
* entry objects. Original sam_entry objects are also kept. Both in
* ref_info object.
*
*
* @param opt        reference options object
* @param ref_in    reference information object, to be created
* @return        error status
*/
int make_targets_info(options_rf opt, ref_info **ref_in)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	sam *sd_ref  = NULL;
	FILE *fp = fopen(opt.sam_file, "r");
	
	if (!fp)
		exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt.sam_file));
	
	read_sam(fp, &sd_ref);	/* assumes XY_ENCODING */
	fclose(fp);
	fp = NULL;
	
	ref_info *rf_info;
	
	*ref_in = malloc(sizeof **ref_in);
	if (!*ref_in)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "ref_in");
	rf_info = *ref_in;
	
	//get the length of each reference name
	size_t *rchar = NULL;
	rchar = malloc(sd_ref->n_ref * sizeof(*rchar));
	size_t rchar_in = 0;
	for (unsigned int i = 0; i < sd_ref->n_ref; ++i) {
		rchar[i] = rchar_in;
		rchar_in += strlen(&sd_ref->ref_names[rchar_in]) + 1;
	}
	//
	//	for (size_t i = 0; i < sd_ref->n_se; ++i) {
	//		sam_entry *se = &sd_ref->se[i];
	//		printf("%d\t", se->ref);
	//		printf("%s %s\n", se->name, &sd_ref->ref_names[rchar[se->ref]]);
	//	}
	
	//store the name[include chrosome] and starting and ending positions in the targeted A, B genome
	// i.e. split aradu.V14167.gnm2.chr02:696301-697301
	ref_entry *re = malloc(sd_ref->n_se * sizeof *re);
	
	for (size_t i = 0; i < sd_ref->n_se; ++i) {
		re[i].name_A = NULL;
		re[i].name_B = NULL;
		sam_entry *se = &sd_ref->se[i];
		/* filter unmapped */
		if (opt.filter_unmapped && se->flag >> 2 & 1)
			continue;
		
		/* strand of A, B, in order to find the reads */
		re[i].strand_A = 0;
		if ((se->flag & 16) == 0)
			re[i].strand_B = 0;
		else
			re[i].strand_B = 1;
		/* find names, the start and end in A, B subgenomes */
		char temp_B[strlen(se->name) + 1];
		strcpy(temp_B, se->name);
		char *ptr_B = strtok(temp_B, opt.delim_ref);
		re[i].name_B = malloc(strlen(ptr_B) + 1);
		strcpy(re[i].name_B, ptr_B);
		ptr_B = strtok(NULL, opt.delim_ref);
		char *ptr_B_pos = strtok(ptr_B, opt.delim_len);
		re[i].start_B = atoi(ptr_B_pos);
		ptr_B_pos = strtok(NULL, opt.delim_len);
		re[i].end_B = atoi(ptr_B_pos);
        
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			  "name: %s start %zu end %zu\n", re[i].name_B, re[i].start_B, re[i].end_B);
		
		char temp_A[strlen(&sd_ref->ref_names[rchar[se->ref]]) + 1];
		strcpy(temp_A, &sd_ref->ref_names[rchar[se->ref]]);
		char *ptr_A = strtok(temp_A, opt.delim_ref);
		re[i].name_A = malloc(strlen(ptr_A) + 1);
		strcpy(re[i].name_A, ptr_A);
		ptr_A = strtok(NULL, opt.delim_ref);
		char *ptr_A_pos = strtok(ptr_A, opt.delim_len);
		re[i].start_A = atoi(ptr_A_pos);
		ptr_A_pos = strtok(NULL, opt.delim_len);
		re[i].end_A = atoi(ptr_A_pos);
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			  "name: %s start %zu end %zu\n", re[i].name_A, re[i].start_A, re[i].end_A);
		
		size_t length = 0;
		for (unsigned int m = 0; m < se->cig->n_ashes; ++m)
			if (se->cig->ashes[m].type == CIGAR_INSERTION
			    || se->cig->ashes[m].type == CIGAR_MATCH
			    || se->cig->ashes[m].type == CIGAR_MMATCH
			    || se->cig->ashes[m].type == CIGAR_MISMATCH)
				length += se->cig->ashes[m].len;
	}
	rf_info->info = re;
	rf_info->ref_sam = sd_ref;
	free(rchar);
	return NO_ERROR;
}/* make_targets_info */

/**
* Find reads aligned to homeologous reference regions. The reads are aligned to
* whole subgenome A and separately to whole subgenome B. We are focused on
* subsets of homeologous target regions. Here, we seek those reads that align
* to any of these homoeologous regions.
*
* @param ref_info    information about homoeologous aligned reference regions
* @param opt        options about homeologous reference regions
* @param sds        sam file objects of reads aligned to subgenomes
* @return        error status
*/
int pickreads(ref_info *ref_info, options_rf *opt, sam **sds) {
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	unsigned int i, j, m;
    size_t *rchar[N_FILES];

    /* get start index of each reference name */
    for (j = 0; j < N_FILES; ++j) {
        size_t rchar_in = 0;

        rchar[j] = malloc(sds[j]->n_ref * sizeof **rchar);
        
        for (i = 0; i < sds[j]->n_ref; ++i) {
            rchar[j][i] = rchar_in;

            debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "%zu %s\n",
                rchar[j][i], &sds[j]->ref_names[rchar[j][i]]);
            
            rchar_in += strlen(&sds[j]->ref_names[rchar_in]) + 1;
        }
//        for (m = 0; m < sds[j]->n_se; ++m) {
//            sam_entry *se = &sds[j]->se[m];
//            if ((se->flag & (1 << 2)))
//                continue;
//            cigar *cig = se->cig;
//            printf("%d: %zu ", m, cig->length_rf);
//            printf("%s\t",  &sds[j]->ref_names[rchar[j][se->ref]]);
//        }
//        printf("\n");
    }
    /* find the reference region the read aligns to, if any */
	for (j = 0; j < N_FILES; ++j) {
		for (m = 0; m < sds[j]->n_se; ++m) {
			sam_entry *se = &sds[j]->se[m];
			se->which_ref = -1;
			se->ref_name = NULL;
			if ((se->flag & (1 << 2)))
				continue;
			cigar *cig = se->cig;
			for (i = 0; i < ref_info->ref_sam->n_se; ++i) {
				ref_entry *re = &ref_info->info[i];
				sam_entry *rse = &ref_info->ref_sam->se[i];
				if (rse->flag >> 11 & 1) //skip secondary
					continue;
				rse->cig->length_rf;
				char *ref_names[N_FILES] = {re->name_A, re->name_B};
				size_t start_pos[N_FILES] = {re->start_A, re->start_B}; // 0 based
				size_t end_pos[N_FILES] = {re->end_A, re->end_B};// 1 based
//				printf("%s %s\n ", ref_names[j], &sds[j]->ref_names[rchar[j][se->ref]]);
				if (!strcmp(&sds[j]->ref_names[rchar[j][se->ref]], ref_names[j])) {
					size_t len_ref = end_pos[j] - start_pos[j];
					size_t rf_index_s = se->pos - 1;
					size_t rf_index_e = rf_index_s + cig->length_rf;
//					printf("%zu %zu || %zu %zu\n", rf_index_s, rf_index_e, start_pos[j], end_pos[j]);
					if(cig->length_rf < len_ref) {
						if ((rf_index_s >= start_pos[j] && rf_index_e <= end_pos[j]) ||
						    (rf_index_s < start_pos[j] && rf_index_e >= start_pos[j] + 1) ||
						    (rf_index_s + 1 < end_pos[j] && rf_index_e >= end_pos[j])) {
							size_t length = strlen(ref_names[j]) + strlen(opt->delim_len) + strlen(opt->delim_ref) + (int)(log10(end_pos[j]) + 1) + 1;
							if (start_pos[j] != 0) {
								length += (int)(log10(start_pos[j]) + 1);
							} else
								length += 1;
							se->ref_name = malloc(length);
							sprintf(se->ref_name, "%s%s%zu%s%zu", ref_names[j], opt->delim_ref, start_pos[j], opt->delim_len, end_pos[j]);
							se->which_ref = i;
							debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "REF_ID: %d REF_NAME: %s \n", se->which_ref, se->ref_name);
							break;
						}
					} else {
						if(start_pos[j] >= rf_index_s && end_pos[j] <= rf_index_e) {
							size_t length = strlen(ref_names[j]) + strlen(opt->delim_len) + strlen(opt->delim_ref) + (int)(log10(start_pos[j]) + 1) + (int)(log10(end_pos[j]) + 1) + 1;
							se->ref_name = malloc(length);
							sprintf(se->ref_name, "%s%s%zu%s%zu", ref_names[j], opt->delim_ref, start_pos[j], opt->delim_len, end_pos[j]);
							se->which_ref = i;
							debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "REF_ID: %d REF_NAME: %s \n", se->which_ref, se->ref_name);
							break;
						}
					}
				}
			}
		}
	}
	for (j = 0; j < N_FILES; ++j)
		free(rchar[j]);
	
	return NO_ERROR;
}/* pickreads */

/**
* Extract selected regions from fasta file.
*
* @param samtools_command    samtools executable
* @param region        samtools region specification: chr:from-to
* @param ref_file        fasta file to extract regions from
* @param ext_rf        fasta file to write with chosen region
*/
void extract_ref(const char *samtools_command, char *region,
                const char *ref_file, const char *ext_rf)
{
    // index the whole reference genome file
    if (!ref_file)
        mmessage(ERROR_MSG, FILE_NOT_FOUND, "Reference file '%s' not found\n", ref_file);
    
    unsigned int cmd_len = strlen(samtools_command) + strlen(ref_file)
        + strlen(" faidx  -o ") + strlen(region) + strlen(ext_rf) + 8;
    char *command = NULL;

    mmessage(INFO_MSG, NO_ERROR, "Length of command: %u\n", cmd_len);

    command = malloc(cmd_len * sizeof *command);
    if (!command)
        mmessage(ERROR_MSG, MEMORY_ALLOCATION, "samtools command");
    
    sprintf(command, "%s faidx %s %s -o %s",
        samtools_command, ref_file, region, ext_rf);
    
    mmessage(INFO_MSG, NO_ERROR, "Running samtools: '%s'\n", command);

    system(command);

    free(command);
    
}/* extract_ref */

int parse_rf_options(options_rf *opt, int argc, char *argv[])
{
	int argv_idx = 1;
	
	while (argv_idx < argc) {
		int j = 0;
		int len = strlen(argv[argv_idx]);
		while (j < len && argv[argv_idx][j] == '-')
			++j;
		char c = argv[argv_idx][j];
		switch (c) {
			case 'u':
				opt->filter_unmapped = 0;
				break;
			case 'l':
				opt->delim_len = argv[++argv_idx];
				break;
			case 'r':
				opt->delim_ref = argv[++argv_idx];
				break;
			case 'f':
				opt->sam_file = argv[++argv_idx];
				break;
			case 'b':
				for (j = 0; j < N_FILES; ++j) {
					opt->rsam_files[j] = argv[++argv_idx];
					fprintf(stderr, " %s", opt->rsam_files[j]);
				}
				fprintf(stderr, "\n");
				break;
			case 'd':
				for (j = 0; j < N_FILES; ++j) {
					opt->fsa_files[j] = argv[++argv_idx];
					fprintf(stderr, " %s", opt->fsa_files[j]);
				}
				fprintf(stderr, "\n");
				break;
			case 's':
				opt->samtools_command = argv[++argv_idx];
				mmessage(INFO_MSG, NO_ERROR, "Samtools "
					 "command: '%s'\n",
					 opt->samtools_command);
			case 'h':
			default:
				if (c != 'h')
					mmessage(ERROR_MSG, INVALID_CMD_OPTION,
						 argv[argv_idx]);
				usage_error((const char **)argv, argv_idx, opt);
				return 1;
		}
		++argv_idx;
	}
	return 0;
} /* parse_rf_options */

// output the reads aligned to the target including the one not aligned to both
void output_selected_reads(const char *f, sam **sds, merge_hash *mh) {
	FILE *fpp = NULL;
	fpp = fopen(f, "w");
	if (!fpp)
		exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, f));
	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
		sam_entry *se;
		if (me->nfiles != N_FILES) {
			if (me->indices[0])
				se = &sds[0]->se[me->indices[0][0]];
			else
				se = &sds[1]->se[me->indices[1][0]];
		} else {
			se = &sds[0]->se[me->indices[0][0]];
		}
		
		fprintf(fpp, "@%s\n", se->name);
		fwrite_nuc_segment(fpp, se->read, XY_ENCODING, 0,
				   se->read->len);
		fprintf(fpp, "\n+\n");
		fwrite_qual_sequence(fpp, se->qual);
		fprintf(fpp, "\n");
		continue;
	}
	fclose(fpp);
}

// NOTE: NUCMER PRODUCES SECONDARY ALIGNMENTS, FILTER SECONDARY, need to consider the RC in reads
// uni_len is the consumed length in genome
// uni_aln_len is the alignment length
// s_id is universal genome starting aligned position relative to the whole genome
// e_id is uni alignment end position relative to the whole genome
//read:     GT-CT
//A:      AGGTCCT
//U:      ANGTCNT
//B:    CCACGTC-T
//read:     GTC-T
void adjust_alignment(sam_entry *se, data_t *ref, unsigned int strand, int *id_uni, size_t uni_aln_len, long *real_id) {
	
	size_t length = 0;
	int rd_idx = 0;
	
	//	 get the length of consumed reference(need to be adjusted if the starting position is before the selected region)
	size_t s_id = real_id[0];
	size_t e_id = real_id[uni_aln_len - 1];
	for (unsigned int i = 0; i < se->cig->n_ashes; ++i)
		if (se->cig->ashes[i].type == CIGAR_DELETION
		    ||se->cig->ashes[i].type == CIGAR_SKIP
		    || se->cig->ashes[i].type == CIGAR_MATCH
		    || se->cig->ashes[i].type == CIGAR_MMATCH
		    || se->cig->ashes[i].type == CIGAR_MISMATCH)
			length += se->cig->ashes[i].len;
#ifdef STANDALONE
	printf("reference length comsumed: %lu\n", length);
#endif
	// length of reference != length summed by using CIGAR
	// length of reference = position - 1 (if alignment start with S/H) + M + D
	int *rid_map = NULL; // read alignment index in the original alignment, -1 represent deletion
	rid_map = malloc(length * sizeof(*rid_map));
	for (size_t j = 0; j < length; ++j)
		rid_map[j] = -1;
	
	// mapped read index(-1 is deletion) relative to the aligned reference
	// 3M2D1I3M: 0 1 2 -1 -1 4 5 6
	size_t rf_idx = se->pos - 1;
	for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
		if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP) {
			rd_idx += se->cig->ashes[i].len;
			continue;
		} else if (se->cig->ashes[i].type == CIGAR_DELETION
			   || se->cig->ashes[i].type == CIGAR_SKIP) {
			rf_idx += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_INSERTION) {
			rd_idx += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_MATCH
			   || se->cig->ashes[i].type == CIGAR_MMATCH
			   || se->cig->ashes[i].type == CIGAR_MISMATCH) {
			for (size_t m = rf_idx; m < rf_idx + se->cig->ashes[i].len; ++m)
				rid_map[m - se->pos + 1] = rd_idx + m - rf_idx;
			rf_idx += se->cig->ashes[i].len;
			rd_idx += se->cig->ashes[i].len;
		}
	}
#ifdef STANDALONE
	printf("old alignment index, start from %zu\n", se->pos);
#endif
	// adjust and trim the alignment (need to consider genome reverse c and read RC)
	int len = 0;
	se->rd_map = NULL;
	unsigned int *remain = NULL;
	unsigned int gaps = 0;
	int location = 0;
	size_t new_len = 0;
	if (!strand) { // if genome is forward
#ifdef STANDALONE
		printf("genome forward\n");
#endif
		if (se->pos < s_id) { // then rd_map is the read (covers the begining) index when aligned to picked region
#ifdef STANDALONE
			printf("align start before the real homo\n");
#endif
			len = s_id - se->pos; // length outside targeted genome
			remain = malloc(uni_aln_len * sizeof(*remain));
			int tmp = length - (s_id - se->pos);
			for (unsigned int i = 0; i < uni_aln_len; ++i) {
				if (id_uni[i] <= tmp) { //uni alignment within the read's aligned region
					if (id_uni[i] == -1)
						remain[gaps++] = i; // record the position that are gaps in uni
				}
				else {
					break;}
			}
			new_len = length - len + gaps;
//			printf("%d gaps meet in the uni genome, new aligned length %lu\n", gaps, length - len + gaps);
			int consumed = 0;
			se->rd_map = malloc(new_len * sizeof(*se->rd_map));
			int flag = 0;
			for (unsigned int i = 0; i < new_len; ++i) {
				
				for (unsigned int j = 0; j < gaps; ++j) {
					if (i == remain[j]) {
						flag = 1;
						consumed++;
						se->rd_map[i] = -1;
						break;
					}
				}
				if (!flag) {
					se->rd_map[i] = rid_map[i + len - consumed];
				}
				flag = 0;
			}
		} else if (se->pos + length - 1 > e_id) {
#ifdef STANDALONE
			printf("align over the real homo\n");
#endif
			// alignment start location relative to id_uni
			for (unsigned int i = 0; i < uni_aln_len; ++i)
				if (se->pos == real_id[i]) {
					location = i;
					break;
				}
//			len = se->pos + length - e_id; //length outside targeted genome
			remain = malloc((uni_aln_len - location) * sizeof(*remain));
			for (unsigned int i = location; i < uni_aln_len; ++i) {
				if (id_uni[i] == -1)
					remain[gaps++] = i - location; // record the position that are gaps in uni
			}
//			printf("%d gaps meet in the uni genome, new aligned length %lu\n", gaps, e_id - se->pos + 1 + gaps);
			int consumed = 0;
			new_len = e_id - se->pos + 1 + gaps;
			se->rd_map = malloc((new_len) * sizeof(*se->rd_map));
			int flag = 0;
			for (unsigned int i = 0; i < new_len; ++i) {
				
				for (unsigned int j = 0; j < gaps; ++j) {
					if (i == remain[j]) {
						flag = 1;
						consumed++;
						se->rd_map[i] = -1;
						break;
					}
				}
				if (!flag)
					se->rd_map[i] = rid_map[i - consumed];
				flag = 0;
			}
		} else { // alignment in the region
			// alignment start location relative to id_uni
			for (unsigned int i = 0; i < uni_aln_len; ++i)
				if (se->pos == real_id[i]) {
					location = i;
					break;
				}
//			printf(" start align location %d\t\t", location);
			remain = malloc((uni_aln_len - location) * sizeof(*remain));
			long tmp = se->pos + length;
			for (unsigned int i = location; i < uni_aln_len; ++i) {
				if (tmp > real_id[i]) {// if the consumed ref is within this
//					printf("%ld\t", real_id[i]);
					if (id_uni[i] == -1)
						remain[gaps++] = i - location;} // record the position that are gaps in uni
				else {
					break;
				}
			}
//			printf("%d gaps meet in the uni genome, new aligned length %lu\n", gaps, length + gaps);
			int consumed = 0;
			new_len = length + gaps;
			se->rd_map = malloc(new_len * sizeof(*se->rd_map));
			int flag = 0;
			for (unsigned int i = 0; i < new_len; ++i) {
				
				for (unsigned int j = 0; j < gaps; ++j) {
					if (i == remain[j]) {
						flag = 1;
						consumed++;
						se->rd_map[i] = -1;
						break;
					}
				}
				if (!flag) {
					se->rd_map[i] = rid_map[i - consumed];
				}
				flag = 0;
			}
		}
	} else {// if genome is reversed
#ifdef STANDALONE
		printf("genome reversed\n");
#endif
		if (se->pos < e_id) { // then se->rd_map is the read (covers the begining) index when aligned to picked region
#ifdef STANDALONE
			printf("align over the homo\t\t");
#endif
			len = e_id - se->pos; // length outside targeted genome
			remain = malloc(uni_aln_len * sizeof(*remain));
			int tmp = length - len;
			for (unsigned int i = 0; i < uni_aln_len; ++i) {
				if (id_uni[i] < tmp) { //uni alignment within the read's aligned region
					if (id_uni[i] == -1)
						remain[gaps++] = i; // record the position that are gaps in uni
				}
				else {
					break;}
			}
			int consumed = 0;
//			printf("%d gaps meet in the uni genome, %d cut length, new aligned length %lu\n", gaps, len, length - len + gaps);
			new_len = length - len + gaps - 1; // temporary - 1
			se->rd_map = malloc(new_len * sizeof(*se->rd_map));
			int flag = 0;
			for (unsigned int i = 0; i < new_len; ++i) {
				
				for (unsigned int j = 0; j < gaps; ++j) {
					if (i == remain[j]) {
						flag = 1;
						consumed++;
						se->rd_map[i] = -1;
						break;
					}
				}
				if (!flag) {
//					printf("%d %lu\t", i, length - (i - consumed) - 1);
					se->rd_map[i] = rid_map[length - (i - consumed) - 1];
				}
				flag = 0;
			}
			location = uni_aln_len - new_len;
		}  else if (se->pos + length - 1 > s_id) {
#ifdef STANDALONE
			printf("align before the homo\n");
#endif
			gaps = 0;
			// alignment start location relative to id_uni
			for (unsigned int i = 0; i < uni_aln_len; ++i)
				if (se->pos == real_id[i]) {
					location = i;
					break;
				}
			remain = malloc((location) * sizeof(*remain));
//			len = se->pos + length - s_id; //length outside targeted genome
			for (unsigned int i = location; i <= 0; --i) {
				if (id_uni[i] == -1)
					remain[gaps++] = location; // record the position that are gaps in uni
			}
			int consumed = 0;
			new_len = s_id - (se->pos - 1) + gaps + 1;
			se->rd_map = malloc(new_len * sizeof(*se->rd_map));
//			printf("%d gaps meet in the uni genome, new aligned length %lu\n", gaps, s_id - se->pos + gaps + 1);
			
			int flag = 0;
			for (unsigned int i = 0; i < new_len; ++i) {
				
				for (unsigned int j = 0; j < gaps; ++j) {
					if (i == remain[j]) {
						flag = 1;
						consumed++;
						se->rd_map[i] = -1;
						break;
					}
				}
				if (!flag) {
//					printf("%d %lu\t", i, s_id - se->pos - (i - consumed));
					se->rd_map[i] = rid_map[s_id - se->pos - (i - consumed) + 1];
				}
				flag = 0;
			}
			location = 0;
		} else { // alignment in the region
			gaps = 0;
			// alignment start location relative to id_uni
			for (unsigned int i = 0; i < uni_aln_len; ++i)
				if (se->pos == real_id[i]) {
					location = i;
					break;
				}
			remain = malloc((location) * sizeof(*remain));
			long tmp = se->pos + length;
			for (unsigned int i = location; i <= 0; --i) {
				if(tmp > real_id[i]){ // if the consumed ref is within this
					if (id_uni[i] == -1)
						remain[gaps++] = location; }// record the position that are gaps in uni
				else {
					break;
				}
			}
			int consumed = 0;
			new_len = length + gaps;
			se->rd_map = malloc(new_len * sizeof(*se->rd_map));
//			printf("%d gaps meet in the uni genome, new aligned length %lu\n", gaps, length + gaps);
			int flag = 0;
			for (unsigned int i = 0; i < new_len; ++i) {
				
				for (unsigned int j = 0; j < gaps; ++j) {
					if (i == remain[j]) {
						flag = 1;
						consumed++;
						se->rd_map[i] = -1;
						break;
					}
				}
				if (!flag) {
//					printf("%d %lu\t", i, length - (i - consumed) - 1);
					se->rd_map[i] = rid_map[length - (i - consumed) - 1];
				}
				flag = 0;
			}
			location = location - new_len + 1 + 1;
		}
	}
	
	// show the alignment relative to the uni genome
	se->uni_aln = NULL;
	se->uni_aln = malloc(new_len * sizeof(*se->uni_aln));
	data_t *read = NULL;
	read = malloc(se->read->len * sizeof(*read));
	for (size_t i = 0; i < se->read->len; ++i)
		read[i] = read_char(se->read, &_xy_sequence_opt, i);
	
//	fwrite_nuc_segment(stdout, se->read, XY_ENCODING, 0,
//			   se->read->len);
	se->new_pos = location;
	data_t const to_compliment[NUM_NUCLEOTIDES] = {
		XY_T,
		XY_G,
		XY_A,
		XY_C
	};
	
	for (size_t j = 0; j < new_len; ++j) {
		// to get the quality score: get_qual(se->qual, se->rd_map[j]);
		if (se->rd_map[j] == -1)
			se->uni_aln[j] = 4; // here use 4 to represent gaps
		else {
			if (strand)
				se->uni_aln[j] = to_compliment[read[se->rd_map[j]]];
			else
				se->uni_aln[j] = read[se->rd_map[j]];
		}
	}
	
	se->aln_len = new_len;
#ifdef STANDALONE
//	PRINT_VECTOR(se->rd_map, length);
	printf("alignment of read to universal, start from %zu, length %zu\n", se->new_pos, se->aln_len);
	PRINT_VECTOR(se->uni_aln, se->aln_len);
	printf("alignment of ref\n");
	for (size_t j = 0; j < se->aln_len; ++j) {
		printf("%d\t", iupac_to_xy[ref[se->new_pos + j]]); // se->new_pos is 0 based
	}
#endif
	// compute the alignment probability, if, ther is a gap in the read alignment but not in uni genome, then add a penalty
	se->ll_aln = 0;
	unsigned int gap_in = 0;
	for (size_t j = 0; j < new_len; ++j) {
		if (se->uni_aln[j] == 4) {
			gap_in++;
			continue;
		}
		double llt = sub_prob_given_q_with_encoding(ref[se->new_pos + j], se->uni_aln[j],
							    IUPAC_ENCODING, XY_ENCODING,
							    get_qual(se->qual, se->rd_map[j]), 1, (void *) se);
		se->ll_aln += llt;
	}
	if (!gap_in) {
		double penalty = gap_in * log(1e-5 * new_len) - 1e-5 * new_len;
		for (unsigned int i = 0; i < gap_in; ++i)
			penalty -= log(i+1);
//		printf("penalty %lf\n", penalty);
		se->ll_aln += penalty;
	}
#ifdef STANDALONE
	printf("%lf\n", se->ll_aln);
#endif
	free(remain);
	free(read);
}

void output_data(FILE *fp, sam_entry *se, unsigned int id) {
	for (unsigned int i = 0; i < se->aln_len; ++i) {
		fprintf(fp, "%d ", id);
		if (se->uni_aln[i] == 4 && i != 0) {
			unsigned int id = i - 1;
			while (se->rd_map[id] == -1) {
				id--;
			}
			fprintf(fp, "%d ", se->rd_map[id]);
		} else
			fprintf(fp, "%d ", se->rd_map[i]);
		fprintf(fp, "%lu ", se->new_pos + i);
		if (se->uni_aln[i] == 4) {
			fprintf(fp, "-1 -\n");
		} else
			fprintf(fp, "%d %c\n", get_qual(se->qual, se->rd_map[i]), xy_to_char[se->uni_aln[i]]);
	}
}


