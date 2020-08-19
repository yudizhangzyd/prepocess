CC = gcc

CFLAGS = -Wall -Wextra -pedantic -O3 -std=c11 -D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=500 -I/usr/include/R/ -FPIC		# c11 with POSIX 20 112L and optimization

# -e DBG=1
ifdef DBG
	CFLAGS = -Wall -Wextra -pedantic -g -std=c11 -D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=500 -DDEBUGGING_CODE		# debugging
endif

# -e DSTD=1
ifdef DSTD
	CFLAGS += -DSTANDALONE
endif

IFLAGS = 

LDFLAGS = -lm -lncurses

# -e YUDI=1
ifdef YUDI
	CFLAGS += -DYUDI -I/usr/local/Cellar/r/4.0.2_1/include/ -I/usr/local/Cellar/openblas/0.3.5/include/
	LDFLAGS += -L/usr/local/Cellar/r/4.0.2_1/lib/ -L/usr/local/opt/openblas/lib
endif

# Local variables
srcs = $(wildcard *.c)
hds = $(wildcard *.h)
objs = $(srcs:.c=.o)
deps = $(srcs:.c=.d)

sam_objs = sync_data_r.o sam.o sequence.o nuc.o qual.o error.o fastq.o align.o io.o cmdline.o make_aln.o options.o pick_reads.o

aln_objs = make_aln.o options.o sam.o sequence.o nuc.o qual.o error.o fastq.o align.o io.o cmdline.o pick_reads.o

simu_objs = simulation.o sequence.o nuc.o qual.o error.o fastq.o align.o io.o cmdline.o

Rlib: $(sam_objs)
	R CMD SHLIB $(sam_objs) $(LDFLAGS)

aln:	$(aln_objs)
	$(CC) -o aln $(aln_objs) $(CFLAGS) $(LDFLAGS) -lz

simu:	$(simu_objs)
	$(CC) -o run_simu $(simu_objs) $(CFLAGS) $(LDFLAGS) -lRmath

include $(deps)

%.d : %.c
	-@$(SHELL) -ec '$(CC) -MM $(CFLAGS) $(IFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@'

.PHONY : archive clean clear dummy

archive: 
	tar czvf haplotype.tar.gz $(srcs) $(hds) makefile

clean:	
	-rm $(objs) $(deps) Rlib aln simu 2>/dev/null

clear:
	rm $(objs)

dummy:
	@echo $(objs)
