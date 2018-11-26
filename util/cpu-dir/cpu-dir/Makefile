CC=			gcc
#CC=			clang --analyze
#CFLAGS=		-g -Wall -Wno-unused-function  
#CFLAGS=		-g -Wall -Wno-unused-function -O2
CFLAGS=		-g -Wall -Wno-unused-function -O3
WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
AR=			ar
DFLAGS=		-DHAVE_PTHREAD $(WRAP_MALLOC)
LOBJS=		utils.o kthread.o kstring.o ksw.o bwt.o bntseq.o bwa.o bwamem.o bwamem_pair.o bwamem_extra.o malloc_wrap.o
AOBJS=		QSufSort.o bwt_gen.o bwase.o bwaseqio.o bwtgap.o bwtaln.o bamlite.o \
			is.o bwtindex.o bwape.o kopen.o pemerge.o \
			bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o \
			bwtsw2_chain.o fastmap.o bwtsw2_pair.o
PROG=		cpu
INCLUDES=	-I../zlib-1.2.8/
#INCLUDES=	
#LIBS=		-lm -lz -lpthread
#LIBS=		-lm -lpthread -L../zlib-1.2.8/ -lz
LIBS=		-lm -lpthread 
SUBDIRS=	.

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bwa:libbwa.a $(AOBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) main.o -o $@ -L. -lbwa $(LIBS)

bwamem-lite:libbwa.a example.o
		$(CC) $(CFLAGS) $(DFLAGS) example.o -o $@ -L. -lbwa $(LIBS)

cpu:libbwa.a cpu.o cpuadapter.o cpulinker.o cputag.o cpustag.o cpuread.o cpufq.o cpustat.o cpulmp.o cpumemaln.o ssw.o kopen.o bwase.o bwtaln.o bwtgap.o bwaseqio.o bam.o sam.o sam_header.o bgzf.o razf.o cpucluster.o cpumtdedup.o cpudedup.o cpuspan.o cpupair.o bam_import.o bam_pileup.o faidx.o bam_aux.o libz.a libbwa.a
		$(CC) $(CFLAGS) $(DFLAGS) cpu.o cpuadapter.o cpulinker.o cputag.o cpustag.o cpuread.o cpufq.o cpustat.o cpulmp.o cpumemaln.o ssw.o kopen.o bwase.o bwtaln.o bwtgap.o bwaseqio.o bam.o sam.o sam_header.o bgzf.o razf.o cpucluster.o cpumtdedup.o cpudedup.o cpuspan.o cpupair.o bam_import.o bam_pileup.o faidx.o bam_aux.o libz.a libbwa.a -o $@ -L. $(LIBS)

libbwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )

# DO NOT DELETE THIS LINE -- make depend depends on it.

faidx.o : faidx.h
sam_header.o : sam_header.h
faidx.o : faidx.h
sam.o : sam.h
razf.o : razf.h
bgzf.o : bgzf.h
bam.o : bam.h
cpucluster.o: bwa.h kstring.h kseq.h utils.h bwa.h sam_header.h sam.h cputagcommon.h cputagbam.h cpuspan.h
cpumtdedup.o: bwa.h kstring.h kseq.h utils.h bwa.h sam_header.h sam.h cputagcommon.h cputagbam.h
cpudedup.o: bwa.h kstring.h kseq.h utils.h bwa.h sam_header.h sam.h cputagcommon.h cputagbam.h
cpuspan.o: bwa.h kstring.h kseq.h utils.h bwa.h sam_header.h sam.h cputagcommon.h cputagbam.h
cpupair.o: bwamem.h kvec.h malloc_wrap.h utils.h kseq.h cpuspan.h
cpumemaln.o: bwa.h bntseq.h bwt.h bwamem.h kvec.h malloc_wrap.h utils.h kseq.h
cpulmp.o: utils.h kseq.h kstring.h bwa.h bwamem.h cputag.h cpuadapter.h cpulinker.h
cpustat.o: utils.h kseq.h kstring.h bwa.h bwamem.h cputag.h
cpufq.o: utils.h kseq.h kstring.h bwa.h bwamem.h cputag.h
cpuread.o: utils.h kseq.h kstring.h bwa.h bwamem.h
cpustag.o: utils.h kseq.h kstring.h bwa.h bwamem.h cpustag.h cpuadapter.h cpulinker.h
cputag.o: utils.h kseq.h kstring.h bwa.h bwamem.h cputag.h cpuadapter.h cpulinker.h
cpulinker.o: utils.h kseq.h kstring.h bwa.h bwamem.h cputag.h cpulinker.h
cpuadapter.o: utils.h kseq.h kstring.h bwa.h bwamem.h cputag.h cpuadapter.h
cpu.o: utils.h
QSufSort.o: QSufSort.h
bamlite.o: bamlite.h malloc_wrap.h
bntseq.o: bntseq.h utils.h kseq.h malloc_wrap.h
bwa.o: bntseq.h bwa.h bwt.h ksw.h utils.h kstring.h malloc_wrap.h kseq.h
bwamem.o: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h ksw.h kvec.h
bwamem.o: ksort.h utils.h kbtree.h
bwamem_extra.o: bwa.h bntseq.h bwt.h bwamem.h kstring.h malloc_wrap.h
bwamem_pair.o: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h kvec.h
bwamem_pair.o: utils.h ksw.h
bwape.o: bwtaln.h bwt.h kvec.h malloc_wrap.h bntseq.h utils.h bwase.h bwa.h
bwape.o: ksw.h khash.h
bwase.o: bwase.h bntseq.h bwt.h bwtaln.h utils.h kstring.h malloc_wrap.h
bwase.o: bwa.h ksw.h
bwaseqio.o: bwtaln.h bwt.h utils.h bamlite.h malloc_wrap.h kseq.h
bwt.o: utils.h bwt.h kvec.h malloc_wrap.h
bwt_gen.o: QSufSort.h malloc_wrap.h
bwt_lite.o: bwt_lite.h malloc_wrap.h
bwtaln.o: bwtaln.h bwt.h bwtgap.h utils.h bwa.h bntseq.h malloc_wrap.h
bwtgap.o: bwtgap.h bwt.h bwtaln.h malloc_wrap.h
bwtindex.o: bntseq.h bwt.h utils.h malloc_wrap.h
bwtsw2_aux.o: bntseq.h bwt_lite.h utils.h bwtsw2.h bwt.h kstring.h
bwtsw2_aux.o: malloc_wrap.h bwa.h ksw.h kseq.h ksort.h
bwtsw2_chain.o: bwtsw2.h bntseq.h bwt_lite.h bwt.h malloc_wrap.h ksort.h
bwtsw2_core.o: bwt_lite.h bwtsw2.h bntseq.h bwt.h kvec.h malloc_wrap.h
bwtsw2_core.o: khash.h ksort.h
bwtsw2_main.o: bwt.h bwtsw2.h bntseq.h bwt_lite.h utils.h bwa.h
bwtsw2_pair.o: utils.h bwt.h bntseq.h bwtsw2.h bwt_lite.h kstring.h
bwtsw2_pair.o: malloc_wrap.h ksw.h
example.o: bwamem.h bwt.h bntseq.h bwa.h kseq.h malloc_wrap.h
fastmap.o: bwa.h bntseq.h bwt.h bwamem.h kvec.h malloc_wrap.h utils.h kseq.h
is.o: malloc_wrap.h
kopen.o: malloc_wrap.h
kstring.o: kstring.h malloc_wrap.h
ksw.o: ksw.h malloc_wrap.h
main.o: kstring.h malloc_wrap.h utils.h
malloc_wrap.o: malloc_wrap.h
pemerge.o: ksw.h kseq.h malloc_wrap.h kstring.h bwa.h bntseq.h bwt.h utils.h
test-extend.o: ksw.h
utils.o: utils.h ksort.h malloc_wrap.h kseq.h
ssw.o: ssw.h
