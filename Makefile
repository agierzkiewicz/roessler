# a list of all the programs in your package
PROGS = 01-Roessler_a525 \
02-Roessler_a47

# a list of all your units to be linked with your programs
OTHERS = utils

# directory where capd scripts are (e.g. capd-config)
CAPDBINDIR = /home/zorra/capd-build/bin/

# setting compiler and linker flags
CAPDFLAGS = `${CAPDBINDIR}capd-config --cflags`
CAPDLIBS = `${CAPDBINDIR}capd-config --libs`
CXXFLAGS += ${CAPDFLAGS} -O2 -std=c++11

# directory where object and dependancy files will be created
OBJDIR = .obj/


#============ the following should be used for graphics later =========

#LIBS = `$(CAPD)/capd-config --cflags --libs` `$(CAPD)/capd-gui-config --cflags --libs` -D__USE_GRAPHICS__

#GINCLUDE = `$(CAPD)/capd-gui-config --cflags` -D__USE_GRAPHICS__
#GLIBS = `$(CAPD)/capd-gui-config --libs`

#============ the following should not be changed =========

OTHERS_OBJ = ${OTHERS:%=${OBJDIR}%.o}
OBJ_FILES = ${OTHERS_OBJ} ${PROGS:%=${OBJDIR}%.o}

.PHONY: all
all: ${PROGS}

# rule to link executables
${PROGS}: % : ${OBJDIR}%.o ${OTHERS_OBJ}
	${CXX} -o $@ $< ${OTHERS_OBJ} ${CAPDLIBS}

# include files with dependencies
-include ${OBJ_FILES:%=%.d}

#rule to compile .cpp files and generate corresponding files with dependencies
${OBJ_FILES}: ${OBJDIR}%.o : %.cpp
	@mkdir -p ${OBJDIR}
	$(CXX) ${CXXFLAGS} -MT $@ -MD -MP -MF ${@:%=%.d} -c -o $@ $<

# rule to clean all object files, dependencies and executables
.PHONY: clean
clean:
	rm -f ${OBJDIR}*.o ${OBJDIR}*.o.d ${PROGS}cd
