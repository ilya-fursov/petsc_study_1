
CFLAGS           =
FFLAGS		 =
CPPFLAGS         =
FPPFLAGS         =
MANSEC           = SNES

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


all: ex3.o  chkopts
	-${CLINKER} -o ex3 ex3.o ${PETSC_SNES_LIB}
	${RM} ex3.o



