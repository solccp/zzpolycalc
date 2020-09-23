SRCDIR := src/
OBJDIR := obj/
MODDIR := obj/
BINDIR := bin/

FC := ifort
CC := icc
#FFLAGS := -profile-functions -profile-loops=all -profile-loops-report=2  -ipo -O3 -funroll-loops -no-prec-div -module ${MODDIR} 
#FFLAGS := -ipo -O3 -funroll-loops -no-prec-div -module ${MODDIR} 
FFLAGS :=-Ofast  -axSSE4.2,AVX,CORE-AVX2 -module ${MODDIR} 

#-check all -debug all
#FFLAGS := -g3 -O0 -module ${MODDIR} -debug -trace
#FFLAGS := -g3 -O0 -module ${MODDIR} -debug -trace -check all -debug all

#CFLAGS = -ipo -O3 -funroll-loops -no-prec-div
CFLAGS :=-Ofast -axSSE4.2,AVX,CORE-AVX2 


LDFLAGS := -static -lpthread
# -debug all
#LDFLAGS := -static -lmkl_em64t -lguide -lpthread

COMMON_FILES := types_module
ZHANG_FILES := crc zhang input functions daughters findZZ polynomial decompose operators print hexagon schlegel
C_FILES := md5

COMMON_SRCS := $(addsuffix .f90, ${COMMON_FILES})
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
	@${FC} ${LDFLAGS} -o $@ $^

${OBJDIR}%.o ${MODDIR}%.mod:${SRCDIR}%.f90
	@echo [FC] $<
	@${FC} ${FFLAGS} -c -o $@ $<

${OBJDIR}%.o :${SRCDIR}%.c
	@echo [CC] $<
	@${CC} ${CFLAGS} -c -o $@ $<


include .depend

.PHONY:	depend
depend .depend: ${SRCS}
	@echo Finding dependencies
	@bin/makedepf90 -b ${OBJDIR} ${SRCS} > .depend

