include ../../Makefile.inc

CXXFLAGS += -MMD -O2
LD_FLAGS+=-L../../core/lib -lcore
INCLUDE+=-I./inc -I../ -I../../core/ -I../../core/inc/


CPP_FILES1=src/matrix_inverse.cpp
OBJ_FILES1=$(CPP_FILES1:.cpp=.o)
CPP_FILES2=src/expr_parser.cpp
OBJ_FILES2=$(CPP_FILES2:.cpp=.o)
EXE:= bin/matrix_inverse
EXE2:=bin/expr_parser

default: $(EXE) $(EXE2)

.PHONY: directories

$(EXE): directories $(OBJ_FILES1)
	$(CXX) $(LD_FLAGS)  -o $@ $(OBJ_FILES1)

$(EXE2): directories $(OBJ_FILES2)
	$(CXX) $(LD_FLAGS)  -o $@ $(OBJ_FILES2)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@



directories: bin obj

bin:
	mkdir -p bin

obj:
	mkdir -p obj

clean:
	rm -rf bin	\
	rm -rf obj	\
	rm -rf .depend

include $(OBJ_FILES1:.o=.d)
include $(OBJ_FILES2:.o=.d)
