
all:  matcmp veccmp mapMerge rhscmp

rhscmp:  rhscmp.C  rhscmp.h basecmp.o matrixcmp.o vectorcmp.o 
	${CC} ${CARGS} -o rhscmp_$(SYSTEM) matrixcmp.o vectorcmp.o basecmp.o rhscmp.C ${LARGS_BLAS}

matcmp:  matrixcmp.C  matrixcmp.h basecmp.o
	${CC} ${CARGS} -o matcmp_$(SYSTEM) basecmp.o matrixcmp.C ${LARGS}

veccmp:  vectorcmp.C vectorcmp.h basecmp.o
	${CC} ${CARGS} -o veccmp_$(SYSTEM) basecmp.o vectorcmp.C ${LARGS}

mapMerge:  mapMerge.C 
	${CC} ${CARGS} -o mapMerge_$(SYSTEM) mapMerge.C ${LARGS}

basecmp.o: basecmp.C basecmp.h
	${CC} ${CARGS} -c  basecmp.C 

matrixcmp.o:  matrixcmp.C matrixcmp.h
	${CC} ${CARGS} -DRHSCMP_BUILD  -c  matrixcmp.C 

vectorcmp.o:  vectorcmp.C vectorcmp.h
	${CC} ${CARGS} -DRHSCMP_BUILD  -c  vectorcmp.C 

clobber:
	@echo Clobber for $(SYSTEM)
	rm basecmp.o; \
	rm veccmp_$(SYSTEM); \
	rm matcmp_$(SYSTEM); \
	rm mapMerge_$(SYSTEM); \
	rm rhscmp_$(SYSTEM);

CC    = $(CC_$(SYSTEM))
CARGS = $(CARGS_$(SYSTEM))
LARGS = $(LARGS_$(SYSTEM))
LARGS_BLAS = $(LARGS_BLAS_$(SYSTEM))

CC_sgi10k    = CC
CARGS_sgi10k = -LANG:std 
LARGS_sgi10k = -lm

CC_linux     = g++
CARGS_linux  = -g -w -pedantic -D_LINUX_
LARGS_linux  = -lm

CC_freebsd     = g++
CARGS_freebsd  = -g -w -pedantic -D__FreeBSD__
LARGS_freebsd  = -lm

CC_OSX     = g++
CARGS_OSX  = -g -w -pg
LARGS_OSX  = -lm
LARGS_BLAS_OSX  = -Wl,-framework,Accelerate -lm

