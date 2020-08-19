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

SEXP r_read_sam (SEXP samfile_r, SEXP ref_name_r, SEXP fastq_file_r, SEXP datafile_r)
{
	sam *sd = NULL;
	char const * ref_name = NULL;
	char const * fastq_file = NULL;
	char const * sam_file = NULL;
	char const * datafile = NULL;
	char filter_unmapped = 1;
	unsigned int min_length = 0;
	unsigned int max_length = 0;
	
	sam_file = CHAR(STRING_ELT(samfile_r, 0));
	if (ref_name)
		ref_name = CHAR(STRING_ELT(ref_name_r, 0));
	
	if (fastq_file)
		fastq_file = CHAR(STRING_ELT(fastq_file_r, 0));
	
	datafile = CHAR(STRING_ELT(datafile_r, 0));
	
	uint32_t rf_idx = -1;
	FILE *fp = fopen(sam_file, "r");
	
	if (!fp)
		exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, sam_file));
	
	read_sam(fp, &sd);	/* assumes XY_ENCODING */
	fclose(fp);
	fp = NULL;
	
	if (fastq_file) {
		fp = fopen(fastq_file, "w");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, fastq_file));
	}
	
	/* find index of selected reference */
	if (ref_name) {
		int rchar = 0;
		int found = 0;
		for (unsigned int i = 0; i < sd->n_ref; ++i) {
			if (!strcmp(ref_name, &sd->ref_names[rchar])) {
				rf_idx = i;
				found = 1;
				break;
			}
			rchar += strlen(&sd->ref_names[rchar]) + 1;
		}
		if (!found)
			error("Bad reference name: '%s'\n", ref_name);
//			exit(mmessage(ERROR_MSG, INVALID_USER_INPUT, "Bad reference name: '%s'\n", ref_name));
	}

	FILE *fp_dat = fopen(datafile, "w");
	if (!fp_dat)
		exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, datafile));
	unsigned int count = 0;
	for (size_t i = 0; i < sd->n_se; ++i) {
		sam_entry *se = &sd->se[i];

		/* filter unmapped */
		if (filter_unmapped && se->flag >> 2 & 1)
			continue;

		/* filter too short or long */
		if ((min_length && se->read->len < min_length)
		    || (max_length && se->read->len > max_length))
			continue;

		/* filter mapped to wrong reference */
		if (ref_name && se->ref != rf_idx)
			continue;
		
		/* output to fastq file */
		if (fastq_file) {
			fprintf(fp, "@%s\n", se->name);
			fwrite_nuc_segment(fp, se->read, XY_ENCODING, 0,
					   se->read->len);
			fprintf(fp, "\n+\n");
			fwrite_qual_sequence(fp, se->qual);
			fprintf(fp, "\n");
		}
		
		se->index = ++count;
		output_error_data(fp_dat, se, NULL, NAN);
	}
	fclose(fp);
	fp = NULL;
	fclose(fp_dat);
	fp_dat = NULL;
	
	return R_NilValue;
}

SEXP r_ampliclust_init(SEXP ampliclust_command_r, SEXP fastq_file_r, SEXP ac_outfile_r)
{
	int err = NO_ERROR;
	const char *ampliclust_command = NULL;	/*<! screen reads w/ ampliCI command */
	const char *fastq_file = NULL;	/*<! ampliCI fastq filename */
	const char *fsa_file = NULL;
	const char *ac_outfile = NULL; /*<! ampliCI output filename */
	
	ampliclust_command = CHAR(STRING_ELT(ampliclust_command_r, 0));
	fastq_file = CHAR(STRING_ELT(fastq_file_r, 0));
	ac_outfile = CHAR(STRING_ELT(ac_outfile_r, 0));
	
	unsigned int cmd_len = strlen(ampliclust_command) + strlen(fastq_file) + strlen(" -f  -n --run -i amplici -lb  -ll -Inf -o ") + strlen(ac_outfile) + 8;
	//fprintf(stderr, "Length of command: %u\n", cmd_len);
	char *command = malloc(cmd_len * sizeof *command);
	sprintf(command, "%s -f %s -n --run -i amplici -lb 1.5 -ll -Inf -o %s",
		ampliclust_command, fastq_file, ac_outfile);
	
	PRINTF("Running ampliclust: '%s'\n", command);
	system(command);
	free(command);
	
	char *outfile_hap = malloc((strlen(ac_outfile)+5) * sizeof (char));
	if (!outfile_hap)
		error("ampliclust output file allocation wrong");
	strcpy(outfile_hap, ac_outfile);
	strcat(outfile_hap, ".fa");
	fsa_file = outfile_hap;
	
	fastq_options *fqo = NULL;
	fastq_data *fdata = NULL;
	
	if ((err = make_fastq_options(&fqo)))
		return R_NilValue;
	fqo->read_encoding = XY_ENCODING;
	read_fastq(fsa_file, &fdata, fqo);
	
	SEXP r_list;
	SEXP r_data_dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(r_data_dim)[0] = fdata->n_reads;
	INTEGER(r_data_dim)[1] = fdata->n_max_length;
	
	SEXP r_reads = PROTECT(allocVector(INTSXP, fdata->n_reads * fdata->n_max_length));
	
	for (unsigned int i = 0; i < fdata->n_reads * fdata->n_max_length; ++i) {
		INTEGER(r_reads)[i] = fdata->reads[i];
	}
	
	PROTECT(r_list = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(r_list, 0, r_reads);
	SET_VECTOR_ELT(r_list, 1, r_data_dim);
	UNPROTECT(3);
	
	if (fdata)
		free_fastq(fdata);
	
	return r_list;
}

SEXP r_read_fasta(SEXP datafile_r)
{
	int err = NO_ERROR;
	unsigned int i;
	fastq_options *fqo = NULL;
	fastq_data *fdata = NULL;
	char const *datafile = NULL;
	SEXP r_list;
	PROTECT(r_list = allocVector(VECSXP, 2));
	datafile = CHAR(STRING_ELT(datafile_r, 0));
	
	if ((err = make_fastq_options(&fqo)))
		return R_NilValue;
	fqo->read_encoding = IUPAC_ENCODING;
	fqo->read_names = 1;
	read_fastq(datafile, &fdata, fqo);
	//SEXP r_name_lengths = PROTECT(allocVector(STRSXP, fdata->n_reads));
	//for (i = 0; i < fdata->n_reads; ++i)
	  //SET_STRING_ELT(r_names, i, mkChar(&fdata->names[i]));
	
	if (fdata->n_max_length == fdata->n_min_length)
	{
		SEXP r_data_dim = PROTECT(allocVector(INTSXP, 2));
		INTEGER(r_data_dim)[0] = fdata->n_reads;
		INTEGER(r_data_dim)[1] = fdata->n_max_length;
		
		SEXP r_reads = PROTECT(allocVector(INTSXP, fdata->n_reads * fdata->n_max_length));
		
		for (i = 0; i < fdata->n_reads * fdata->n_max_length; ++i)
			INTEGER(r_reads)[i] = fdata->reads[i];
		
		SET_VECTOR_ELT(r_list, 0, r_reads);
		SET_VECTOR_ELT(r_list, 1, r_data_dim);
	} else {
		unsigned int total_length = 0;
		SEXP r_data_lengths = PROTECT(allocVector(INTSXP, fdata->n_reads));
		for (i = 0; i < fdata->n_reads; ++i) {
			total_length += fdata->n_lengths[i];
			INTEGER(r_data_lengths)[i] = fdata->n_lengths[i];
		}
		
		SEXP r_reads = PROTECT(allocVector(INTSXP, total_length));
		for (i = 0; i < total_length; ++i)
			INTEGER(r_reads)[i] = fdata->reads[i];
		
		SET_VECTOR_ELT(r_list, 0, r_reads);
		SET_VECTOR_ELT(r_list, 1, r_data_lengths);
	}
	//SET_VECTOR_ELT(r_list, 2, r_names);
	UNPROTECT(3);
	if (fdata)
		free_fastq(fdata);
	
	return r_list;
}

SEXP r_read_fastq(SEXP datafile_r)
{
  int err = NO_ERROR;
  fastq_options *fqo = NULL;
  fastq_data *fdata = NULL;
  char const *datafile = NULL;
  SEXP r_list;
  
  
  datafile = CHAR(STRING_ELT(datafile_r, 0));
  
  if ((err = make_fastq_options(&fqo)))
    return R_NilValue;
  fqo->read_encoding = XY_ENCODING;
  read_fastq(datafile, &fdata, fqo);
  
  SEXP r_data_dim = PROTECT(allocVector(INTSXP, 2));
  INTEGER(r_data_dim)[0] = fdata->n_reads;
  INTEGER(r_data_dim)[1] = fdata->n_max_length;
  
  SEXP r_reads = PROTECT(allocVector(INTSXP, fdata->n_reads * fdata->n_max_length));
  SEXP r_quals = PROTECT(allocVector(INTSXP, fdata->n_reads * fdata->n_max_length));
  
  for (unsigned int i = 0; i < fdata->n_reads * fdata->n_max_length; ++i) {
    INTEGER(r_reads)[i] = fdata->reads[i];
    INTEGER(r_quals)[i] = fdata->quals[i] + fdata->min_quality;
  }
  
  PROTECT(r_list = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(r_list, 0, r_reads);
  SET_VECTOR_ELT(r_list, 1, r_quals);
  SET_VECTOR_ELT(r_list, 2, r_data_dim);
  UNPROTECT(4);
  
  if (fdata)
    free_fastq(fdata);
  
  return r_list;
}


SEXP r_make_aln (SEXP ref_nameA,
		   SEXP ref_nameB,
		   SEXP ref_fsa,
		   SEXP ref_sam,
		   SEXP alnA,
		   SEXP alnB,
		   SEXP out_file,
		   SEXP uni_geno_file
		   )
{
	options opt;
	default_options(&opt);
	opt.sbam_files[0] = CHAR(STRING_ELT(alnA, 0));
	opt.sbam_files[1] = CHAR(STRING_ELT(alnB, 0));
	opt.ref_names[0] = CHAR(STRING_ELT(ref_nameA, 0));
	opt.ref_names[1] = CHAR(STRING_ELT(ref_nameB, 0));
	opt.out_file = CHAR(STRING_ELT(out_file, 0));
	opt.uni_geno_file = CHAR(STRING_ELT(uni_geno_file, 0));
	opt.uni_genome = CHAR(STRING_ELT(ref_fsa, 0));
	opt.sam_file = CHAR(STRING_ELT(ref_sam, 0));
	make_alignment(opt);
	return R_NilValue;
}




