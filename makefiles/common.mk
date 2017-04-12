SRC_DIR = ../../src
SRC = $(wildcard $(SRC_DIR)/*.f90 $(SRC_DIR)/*.F90)
PROG_SRC = $(wildcard $(SRC_DIR)/programs/*.f90 $(SRC_DIR)/programs/*.F90)
TEST_SRC = $(wildcard $(SRC_DIR)/test/*.f90 $(SRC_DIR)/test/*.F90)

OBJ = $(notdir $(SRC:.f90=.o))
OBJ := $(notdir $(OBJ:.F90=.o))

PROG_OBJ = $(notdir $(PROG_SRC:.f90=.o))
PROG_OBJ := $(notdir $(PROG_OBJ:.F90=.o))

PROGS = $(basename $(PROG_OBJ))

TEST_OBJ = $(notdir $(TEST_SRC:.f90=.o))
TEST_OBJ := $(notdir $(TEST_SRC:.F90=.o))

TESTS = $(addprefix test/,$(basename $(TEST_OBJ)))

FCFLAGS += -L../../lib
FCFLAGS += -I../../lib


$(SRC_DIR)/depend: $(SRC) $(PROG_SRC) $(TEST_SRC)
	makedepf90 -W -m"%m.mod" -b"." $^ > $@

-include $(SRC_DIR)/depend

#elastic2d: $(SRC_DIR)/depend $(OBJ) elastic2d.o
#	$(FC) $(FCFLAGS) $(OBJ) elastic2d.o -o elastic2d $(LDFLAGS)

$(PROGS) : % : $(SRC_DIR)/depend $(OBJ) %.o
	$(FC) $(FCFLAGS) $(OBJ) $@.o -o $@ $(LDFLAGS)

$(TESTS) : test/% : $(SRC_DIR)/depend $(OBJ) %.o
	$(FC) $(FCFLAGS) $(OBJ) $(notdir $@).o -o $@ $(LDFLAGS)

# these rules gain more dependencies from 'depend', some of which are
# .mod files (hence the filter on the dependencies)
%.mod:
	$(FC) $(FCFLAGS) $(FCSYNTAXONLY) $(filter %.f90 %.F90, $^)

%.o:
	$(FC) $(FCFLAGS) $(filter %.f90 %.F90, $^) -c -o $@

