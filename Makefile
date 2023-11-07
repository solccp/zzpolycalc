SRCDIR := src/
OBJDIR := obj/
MODDIR := obj/
BINDIR := bin/

#FC := ifort
#CC := icc
FC := gfortran
CC := gcc

#ifort recommended flags
#FFLAGS := -O3 -no-prec-div -static -fp-model fast=2 -axSSE4.2,AVX,CORE-AVX2 -funroll-loops -module ${MODDIR} 

FFLAGS := -g -O2  -J ${MODDIR}

#CFLAGS := -O3 -no-prec-div -static -fp-model fast=2 -axSSE4.2,AVX,CORE-AVX2

CFLAGS := -g -O2

LDFLAGS := -static 


COMMON_FILES :=  hash_modules
PREPROC_FILES := zhang
ZHANG_FILES := types_module input functions daughters findZZ polynomial decompose operators print hexagon schlegel options getopt sort
C_FILES := md5

COMMON_SRCS := $(addsuffix .F90, ${COMMON_FILES})
COMMON_SRCS := $(addprefix ${SRCDIR}, ${COMMON_SRCS})

PREPROC_SRCS := $(addsuffix .F90, ${PREPROC_FILES})
PREPROC_SRCS := $(addprefix ${SRCDIR}, ${PREPROC_SRCS})

PREPROC_OBJS := $(addsuffix .o, ${PREPROC_FILES})
PREPROC_OBJS := $(addprefix ${OBJDIR}, ${PREPROC_OBJS})

COMMON_OBJS := $(addsuffix .o, ${COMMON_FILES})
COMMON_OBJS := $(addprefix ${OBJDIR}, ${COMMON_OBJS})


C_SRCS       := $(addsuffix .c, ${C_FILES})
C_SRCS       := $(addprefix ${SRCDIR}, ${C_SRCS})

C_OBJS       := $(addsuffix .o, ${C_FILES})
C_OBJS       := $(addprefix ${OBJDIR}, ${C_OBJS})


ZHANG_SRCS := $(addsuffix .f90, ${ZHANG_FILES})
ZHANG_SRCS := $(addprefix ${SRCDIR}, ${ZHANG_SRCS})

ZHANG_OBJS := $(addsuffix .o, ${ZHANG_FILES})
ZHANG_OBJS := $(addprefix ${OBJDIR}, ${ZHANG_OBJS})

SRCS := ${COMMON_SRCS} ${PREPROC_SRCS} ${ZHANG_SRCS} ${C_SRCS}

.PHONY: all
all: zhang

.PHONY: zhang
zhang: ${BINDIR}ZZPolyCalc

.PHONY: clean
clean:
	rm -f ${OBJDIR}*.o
	rm -f ${MODDIR}*.mod

.PHONY: mrproper
mrproper: clean
	rm -f ${BINDIR}*
	rm -f ${SRCDIR}*~
	rm -f *~

${BINDIR}ZZPolyCalc: ${ZHANG_OBJS} ${PREPROC_OBJS} ${COMMON_OBJS} ${C_OBJS}
	@echo [LD] $@
	@mkdir -p ${BINDIR}
	@${FC} ${LDFLAGS} -o $@ $^ 

${OBJDIR}%.o ${MODDIR}%.mod:${SRCDIR}%.f90
	@mkdir -p $(@D)
	@echo [FC] $<
	@${FC} ${FFLAGS} -c -o $@ $<

${OBJDIR}%.o ${MODDIR}%.mod:${SRCDIR}%.F90
	mkdir -p $(@D)
	@echo [FC] $<
	@${FC} ${FFLAGS} -c -o $@ $< 


${OBJDIR}%.o :${SRCDIR}%.c
	mkdir -p $(@D)
	@echo [CC] $<
	@${CC} ${CFLAGS} -c -o $@ $< 


include .depend

.PHONY:	depend
depend .depend: ${SRCS}
	@echo Finding dependencies
#	@bin/makedepf90 -b ${OBJDIR} ${SRCS} > .depend

