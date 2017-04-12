# This makefile doesn't know how to make anything not specified here
.SUFFICES: 

# Use the elastic4_arch environmental variable. Give it a short name.
ifndef elastic4_arch
$(error elastic4_arch is not defined)
endif
arch = $(elastic4_arch)

# Should not treat the build directory like an ordinary target: always
# try to invoke the rule if requested.
.PHONY: build/$(arch)

# Any target depends on the build directory, and are made by that
# rule.  Needs to specify a recipe or it thinks it can't build the
# target, but an empty command (with the ';') is enough.  It is a
# double-colon rule so it is terminal for efficiency reasons.
% :: build/$(arch) ;

# Copy the makefile for the appropriate architecture, and hand over to
# make run in the build directory. MAKECMDGOALS is a make built-in
# variable.
build/$(arch):
	mkdir -p $@/test
	cp makefiles/Makefile.$(arch) $@/Makefile
	$(MAKE) -C $@ $(MAKECMDGOALS)
