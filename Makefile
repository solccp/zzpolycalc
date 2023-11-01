SRCDIR := src/
OBJDIR := obj/
MODDIR := obj/
BINDIR := bin/

FC := ifort
CC := icc
#FC := gfortran
#CC := gcc
#FFLAGS := -profile-functions -profile-loops=all -profile-loops-report=2  -ipo -O3 -funroll-loops -no-prec-div -module ${MODDIR} 
#FFLAGS := -ipo -O3 -funroll-loops -no-prec-div -module ${MODDIR} 
#FFLAGS := -g -debug -trace -ipo -O3 -funroll-loops -no-prec-div -module ${MODDIR} 
#FFLAGS :=-ipo  -Ofast  -axSSE4.2,AVX,CORE-AVX2 -module ${MODDIR} 
#FFLAGS :=-ipo  -Ofast  -axSSE4.2,AVX,CORE-AVX2 -funroll-loops -no-prec-div -module ${MODDIR} 
#FFLAGS := -parallel -ipo -O3 -no-prec-div -static -fp-model fast=2 -axSSE4.2,AVX,CORE-AVX2 -funroll-loops -module ${MODDIR} 
FFLAGS := -O3 -no-prec-div -static -fp-model fast=2 -axSSE4.2,AVX,CORE-AVX2 -funroll-loops -module ${MODDIR} 
#FFLAGS :=-profile-functions -profile-loops=all -profile-loops-report=2 -Ofast  -axSSE4.2,AVX,CORE-AVX2 -module ${MODDIR} 


#-check all -debug all
#FFLAGS := -g3 -O0 -module ${MODDIR} -debug -trace
#FFLAGS := -g3 -O0 -module ${MODDIR} -debug -trace -check all -debug all

#CFLAGS = -ipo -O3 -funroll-loops -no-prec-div
#CFLAGS :=-ipo  -Ofast -axSSE4.2,AVX,CORE-AVX2 
#CFLAGS := -parallel -ipo -O3 -no-prec-div -static -fp-model fast=2 -axSSE4.2,AVX,CORE-AVX2
CFLAGS := -O3 -no-prec-div -static -fp-model fast=2 -axSSE4.2,AVX,CORE-AVX2

#CFLAGS :=-profile-functions -profile-loops=all -profile-loops-report=2 -Ofast -axSSE4.2,AVX,CORE-AVX2 
#CFLAGS := -g3 -O0 
#-debug -trace -check=stack,uninit -debug all


LDFLAGS := -static  -lpthread -qopenmp
# -debug all
#LDFLAGS := -static -lmkl_em64t -lguide -lpthread

COMMON_FILES := zhang hash_modules
ZHANG_FILES := types_module input functions daughters findZZ polynomial decompose operators print hexagon schlegel
C_FILES := md5 xxhashwrapper

COMMON_SRCS := $(addsuffix .F90, ${COMMON_FILES})
COMMON_SRCS := $(addprefix ${SRCDIR}, ${COMMON_SRCS})

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

SRCS := ${COMMON_SRCS} ${ZHANG_SRCS} ${C_SRCS}

.PHONY: all
all: zhang

.PHONY: zhang
zhang: ${BINDIR}zhangadjtmp

.PHONY: clean
clean:
	rm -f ${OBJDIR}*.o
	rm -f ${MODDIR}*.mod
	rm .depend

.PHONY: mrproper
mrproper: clean
	rm -f ${BINDIR}*
	rm -f ${SRCDIR}*~
	rm -f *~

${BINDIR}zhangadjtmp: ${ZHANG_OBJS} ${COMMON_OBJS} ${C_OBJS}
	@echo [LD] $@
	@${FC} ${LDFLAGS} -o $@ $^  -L ~/ZZ/xxHash/ -l xxhash

${OBJDIR}%.o ${MODDIR}%.mod:${SRCDIR}%.f90
	@echo [FC] $<
	@${FC} ${FFLAGS} -c -o $@ $<

${OBJDIR}%.o ${MODDIR}%.mod:${SRCDIR}%.F90
	@echo [FC] $<
	@${FC} ${FFLAGS} -c -o $@ $< -DUSE_XXHASH


${OBJDIR}%.o :${SRCDIR}%.c
	@echo [CC] $<
	@${CC} ${CFLAGS} -c -o $@ $< -I ~/ZZ/xxHash


include .depend

.PHONY:	depend
depend .depend: ${SRCS}
	@echo Finding dependencies
	@bin/makedepf90 -b ${OBJDIR} ${SRCS} > .depend

