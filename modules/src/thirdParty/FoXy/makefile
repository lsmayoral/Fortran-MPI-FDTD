#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = exe/obj/
DMOD    = exe/mod/
DEXE    = exe/
LIBS    =
FC      = mpifort
#OPTSL   =  -module exe/mod
OPTSL   =  -J exe/mod
OPTSC   = -DGNU -cpp -c -frealloc-lhs -O2 $(OPTSL)
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)ADD_ATTRIBUTES: $(MKDIRS) $(DOBJ)add_attributes.o
	@rm -f $(filter-out $(DOBJ)add_attributes.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) ADD_ATTRIBUTES
$(DEXE)PARSE_STRING_SIMPLE: $(MKDIRS) $(DOBJ)parse_string_simple.o
	@rm -f $(filter-out $(DOBJ)parse_string_simple.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) PARSE_STRING_SIMPLE
$(DEXE)DELETE_CONTENT: $(MKDIRS) $(DOBJ)delete_content.o
	@rm -f $(filter-out $(DOBJ)delete_content.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) DELETE_CONTENT
$(DEXE)CREATE_TAG: $(MKDIRS) $(DOBJ)create_tag.o
	@rm -f $(filter-out $(DOBJ)create_tag.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) CREATE_TAG
$(DEXE)ADD_TAG: $(MKDIRS) $(DOBJ)add_tag.o
	@rm -f $(filter-out $(DOBJ)add_tag.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) ADD_TAG
$(DEXE)DELETE_ATTRIBUTES: $(MKDIRS) $(DOBJ)delete_attributes.o
	@rm -f $(filter-out $(DOBJ)delete_attributes.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) DELETE_ATTRIBUTES
$(DEXE)WRITE_TAG: $(MKDIRS) $(DOBJ)write_tag.o
	@rm -f $(filter-out $(DOBJ)write_tag.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) WRITE_TAG
$(DEXE)PARSE_FILE_SIMPLE: $(MKDIRS) $(DOBJ)parse_file_simple.o
	@rm -f $(filter-out $(DOBJ)parse_file_simple.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) PARSE_FILE_SIMPLE
$(DEXE)DELETE_TAG: $(MKDIRS) $(DOBJ)delete_tag.o
	@rm -f $(filter-out $(DOBJ)delete_tag.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) DELETE_TAG
$(DEXE)INDENT_TAG: $(MKDIRS) $(DOBJ)indent_tag.o
	@rm -f $(filter-out $(DOBJ)indent_tag.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) INDENT_TAG

#compiling rules
$(DOBJ)add_attributes.o: src/tests/add_attributes.f90 \
	$(DOBJ)foxy.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parse_string_simple.o: src/tests/parse_string_simple.f90 \
	$(DOBJ)foxy.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)delete_content.o: src/tests/delete_content.f90 \
	$(DOBJ)foxy.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)create_tag.o: src/tests/create_tag.f90 \
	$(DOBJ)foxy.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)add_tag.o: src/tests/add_tag.f90 \
	$(DOBJ)foxy.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)delete_attributes.o: src/tests/delete_attributes.f90 \
	$(DOBJ)foxy.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)write_tag.o: src/tests/write_tag.f90 \
	$(DOBJ)foxy.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)parse_file_simple.o: src/tests/parse_file_simple.f90 \
	$(DOBJ)foxy.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)delete_tag.o: src/tests/delete_tag.f90 \
	$(DOBJ)foxy.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)indent_tag.o: src/tests/indent_tag.f90 \
	$(DOBJ)foxy.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)foxy_xml_file.o: src/lib/foxy_xml_file.f90 \
	$(DOBJ)foxy_xml_tag.o \
	$(DOBJ)penf.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)foxy_xml_tag.o: src/lib/foxy_xml_tag.F90 \
	$(DOBJ)penf.o \
	$(DOBJ)stringifor.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)foxy.o: src/lib/foxy.f90 \
	$(DOBJ)foxy_xml_file.o \
	$(DOBJ)foxy_xml_tag.o \
	$(DOBJ)penf.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stringifor_string_t.o: src/third_party/StringiFor/src/lib/stringifor_string_t.F90 \
	$(DOBJ)befor64.o \
	$(DOBJ)penf.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stringifor.o: src/third_party/StringiFor/src/lib/stringifor.F90 \
	$(DOBJ)penf.o \
	$(DOBJ)stringifor_string_t.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)befor64_pack_data_m.o: src/third_party/BeFoR64/src/lib/befor64_pack_data_m.F90 \
	$(DOBJ)penf.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)befor64.o: src/third_party/BeFoR64/src/lib/befor64.F90 \
	$(DOBJ)penf.o \
	$(DOBJ)befor64_pack_data_m.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)penf_stringify.o: src/third_party/PENF/src/lib/penf_stringify.F90 \
	$(DOBJ)penf_b_size.o \
	$(DOBJ)penf_global_parameters_variables.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)penf.o: src/third_party/PENF/src/lib/penf.F90 \
	$(DOBJ)penf_global_parameters_variables.o \
	$(DOBJ)penf_b_size.o \
	$(DOBJ)penf_stringify.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)penf_global_parameters_variables.o: src/third_party/PENF/src/lib/penf_global_parameters_variables.F90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)penf_b_size.o: src/third_party/PENF/src/lib/penf_b_size.F90 \
	$(DOBJ)penf_global_parameters_variables.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
