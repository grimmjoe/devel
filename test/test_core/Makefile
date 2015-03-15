include ../../Makefile.inc

CXXFLAGS += -MMD -O2
LD_FLAGS+=-L../../core/lib -lcore
INCLUDE+=-I./inc -I../ -I../../core/ -I../../core/inc/


CPP_FILES1=src/matrix_inverse.cpp
OBJ_FILES1=$(CPP_FILES1:.cpp=.o)

CPP_FILES2=src/expr_parser.cpp
OBJ_FILES2=$(CPP_FILES2:.cpp=.o)

CPP_FILES3=src/deriv_matrix.cpp
OBJ_FILES3=$(CPP_FILES3:.cpp=.o)

CPP_FILES4=src/b_inverse.cpp
OBJ_FILES4=$(CPP_FILES4:.cpp=.o)

CPP_FILES5=src/q_inverse.cpp
OBJ_FILES5=$(CPP_FILES5:.cpp=.o)

CPP_FILES6=src/rank.cpp
OBJ_FILES6=$(CPP_FILES6:.cpp=.o)

EXE:= bin/matrix_inverse
EXE2:=bin/expr_parser
EXE3:=bin/deriv_matrix
EXE4:=bin/b_inverse
EXE5:=bin/q_inverse
EXE6:=bin/rank

default: $(EXE) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6)

.PHONY: directories

$(EXE): directories $(OBJ_FILES1)
	$(CXX) $(LD_FLAGS)  -o $@ $(OBJ_FILES1)

$(EXE2): directories $(OBJ_FILES2)
	$(CXX) $(LD_FLAGS)  -o $@ $(OBJ_FILES2)

$(EXE3): directories $(OBJ_FILES3)
	$(CXX) $(LD_FLAGS)  -o $@ $(OBJ_FILES3)

$(EXE4): directories $(OBJ_FILES4)
	$(CXX) $(LD_FLAGS)  -o $@ $(OBJ_FILES4)

$(EXE5): directories $(OBJ_FILES5)
	$(CXX) $(LD_FLAGS)  -o $@ $(OBJ_FILES5)

$(EXE6): directories $(OBJ_FILES6)
	$(CXX) $(LD_FLAGS)  -o $@ $(OBJ_FILES6)

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

-include $(OBJ_FILES1:.o=.d)
-include $(OBJ_FILES2:.o=.d)
-include $(OBJ_FILES3:.o=.d)
-include $(OBJ_FILES4:.o=.d)
-include $(OBJ_FILES5:.o=.d)
-include $(OBJ_FILES6:.o=.d)
