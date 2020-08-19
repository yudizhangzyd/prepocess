//
//  options.c
//  Test
//
//  Created by Yudi Zhang on 4/16/20.
//  Copyright © 2020 Yudi Zhang. All rights reserved.
//

#include "options.h"
#include "cmdline.h"

int default_options(options *opt)
{
	opt->display_alignment = 0;
	opt->drop_unmapped = 1;
	opt->drop_secondary = 1;
	opt->drop_soft_clipped = INFINITY;
	opt->drop_indel = INFINITY;
//	opt->proptest_screen = 0;
//	opt->coverage_screen = 0.5;
	opt->min_length = 0;
	opt->max_length = INFINITY;
//	opt->min_log_likelihood = -INFINITY;
	opt->n_sample = 100;
	opt->max_eerr = INFINITY;
	opt->out_file = NULL;
	opt->uni_geno_file = NULL;
	opt->error_file = NULL;
	//	opt->ampliclust_file = NULL;
//	opt->max_quality_score = INT_MAX;
	opt->use_bam = 0;
	for (int i = 0; i < N_FILES; ++i) {
		opt->sbam_files[i] = NULL;
		opt->fsa_files[i] = NULL;
		opt->ref_names[i] = NULL;
	}
	opt->uni_genome = NULL;
	opt->sam_file = NULL;
    opt->selected_fq = NULL;
    
	return NO_ERROR;
} /* default_options */

int parse_options(options *opt, int argc, const char **argv)
{
	int i, j;
	int err = NO_ERROR;
	char a;
	
	for (i = 1; i < argc; ++i) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < (int) strlen(argv[i]))
			a = argv[i][j];
		switch(a) {
				
				/* cases within switch not indented */
			case 'b':
				if (!strncmp(&argv[i][j], "bam", 3)) {
					if (i + N_FILES >= argc) {
						err = mmessage(ERROR_MSG,
							       INVALID_CMD_ARGUMENT, "Too few "
							       "arguments to --bam_files "
							       "command-line option.\n");
						goto CMDLINE_ERROR;
					}
					opt->use_bam = 1;
					mmessage(INFO_MSG, NO_ERROR, "BAM files:");
					for (j = 0; j < N_FILES; ++j) {
						opt->sbam_files[j] = argv[++i];
						fprintf(stderr, " %s",
							opt->sbam_files[j]);
					}
					fprintf(stderr, "\n");
				} else {
					goto CMDLINE_ERROR;
				}
				break;
			case 'f':
				if (i + N_FILES >= argc) {
					err = mmessage(ERROR_MSG, INVALID_CMD_ARGUMENT,
						       "Too few arguments to --fsa_files "
						       "command-line option.\n");
					goto CMDLINE_ERROR;
				}
				mmessage(INFO_MSG, NO_ERROR, "Fasta files:");
				for (j = 0; j < N_FILES; ++j) {
					opt->fsa_files[j] = argv[++i];
					fprintf(stderr, " %s",
						opt->fsa_files[j]);
				}
				fprintf(stderr, "\n");
				break;
			case 'g':
				mmessage(INFO_MSG, NO_ERROR, "Sam file of aligned targets:");
				opt->sam_file = argv[++i];
				fprintf(stderr, " %s",
					opt->sam_file);
				
				fprintf(stderr, "\n");
				break;
			case 'o':
				mmessage(INFO_MSG, NO_ERROR, "Output file:");
				opt->out_file = argv[++i];
				fprintf(stderr, " %s",
					opt->out_file);
				
				fprintf(stderr, "\n");
				break;
            case 'q':
            mmessage(INFO_MSG, NO_ERROR, "Selected reads file:");
            opt->selected_fq = argv[++i];
            fprintf(stderr, " %s",
                opt->selected_fq);
            fprintf(stderr, "\n");
            break;
			case 'n':
				mmessage(INFO_MSG, NO_ERROR, "Uni genome file:");
				opt->uni_geno_file = argv[++i];
				fprintf(stderr, " %s",
					opt->uni_geno_file);
				
				fprintf(stderr, "\n");
				break;
			case 'h':
				fprint_usage(stderr, argv[0], (void *)opt);
				exit(EXIT_SUCCESS);
			case 'j':
				mmessage(INFO_MSG, NO_ERROR, "Fsa file of aligned targets:");
				opt->uni_genome = argv[++i];
				fprintf(stderr, " %s",
					opt->uni_genome);
				
				fprintf(stderr, "\n");
				break;
			case 'r':
				if (i + N_FILES >= argc) {
					err = mmessage(ERROR_MSG, INVALID_CMD_ARGUMENT,
						       "Too few arguments to --ref_names "
						       "command-line option.\n");
					goto CMDLINE_ERROR;
				}
				mmessage(INFO_MSG, NO_ERROR, "Reference names:");
				for (j = 0; j < N_FILES; ++j) {
					opt->ref_names[j] = argv[++i];
					fprintf(stderr, " %s",
						opt->ref_names[j]);
				}
				fprintf(stderr, "\n");
				break;
			case 's':
				if (!strncmp(&argv[i][j], "se", 2)) {
					opt->drop_secondary = !opt->drop_secondary;
					mmessage(INFO_MSG, NO_ERROR, "%s secondary "
						 "alignments\n", opt->drop_secondary
						 ? "Dropping" : "Keeping");
				} else if (!strncmp(&argv[i][j], "so", 2)) {
					opt->drop_soft_clipped
					= read_int(argc, argv, ++i, opt);
					mmessage(INFO_MSG, NO_ERROR, "Dropping reads "
						 "with soft clip >= %d in either "
						 "alignment.\n", opt->drop_soft_clipped);
				} else if (!strncmp(&argv[i][j], "samp", 4)) {
					opt->n_sample = read_uint(argc, argv, ++i, opt);
					mmessage(INFO_MSG, NO_ERROR, "%u Monte Carlo "
						 "samples.\n", opt->n_sample);
				} else if (!strncmp(&argv[i][j], "sam", 3)) {
					if (i + N_FILES >= argc) {
						err = mmessage(ERROR_MSG,
							       INVALID_CMD_ARGUMENT, "Too few "
							       "arguments to --sam_files "
							       "command-line option.\n");
						goto CMDLINE_ERROR;
					}
					mmessage(INFO_MSG, NO_ERROR, "Sam files:");
					for (j = 0; j < N_FILES; ++j) {
						opt->sbam_files[j] = argv[++i];
						fprintf(stderr, " %s",
							opt->sbam_files[j]);
					}
					fprintf(stderr, "\n");
				}
				break;
				
			case 'u':
				opt->drop_unmapped = !opt->drop_unmapped;
				mmessage(INFO_MSG, NO_ERROR, "%s unmapped reads\n",
					 opt->drop_unmapped ? "Dropping" : "Keeping");
				break;
				
			default:
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
		}
	}
	
	return err;
	
CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */

void fprint_usage(FILE *fp, const char *cmdname, void *obj) {
	options *opt = (options *) obj;
	size_t start = strlen(cmdname) - 1;
	
	while (cmdname[start] != '/' && start) start--;
	if (cmdname[start] == '/') start++;
	
	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "\nNAME\n\t%s - genotype tetraploids preprocessing\n", &cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s --sam_files <fsam1> <fsam2> --fsa_files "
		"<fsa1> <fsa2> --j <fsat>\n --g <ftsam>\t\n",
		&cmdname[start]);
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t--o <outfile> \n\t\tOut file \n");
    fprintf(fp, "\t--q <selected> \n\t\tSelected reads(fastq) file \n");
	fprintf(fp, "\t--n <uni_geno> \n\t\tUni_geno file \n");
	fprintf(fp, "\t--sam_files <fsam1> <fsam2>\n\t\tSpecify sam files "
		"containing alignments (Default: none)\n");
	fprintf(fp, "\t--ref_names <sref1> <sref2>\n\t\tSpecify names of "
		"subgenomic references for target region; must exist in sam "
		" files (Default: none)\n");
	fprintf(fp, "\t--g <ftsam>\n\t\tSpecify name of targeted regions sam file\n");
	fprintf(fp, "\t--j <fsat>\n\t\tSpecify name of targeted regions fsa file\n");
} /* fprint_usage */
