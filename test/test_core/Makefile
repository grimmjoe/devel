include ../../Makefile.inc

CXXFLAGS += -O2
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

CPP_FILES7=src/bq_inverse.cpp
OBJ_FILES7=$(CPP_FILES7:.cpp=.o)

CPP_FILES8=src/drazin_recursive.cpp
OBJ_FILES8=$(CPP_FILES8:.cpp=.o)

CPP_FILES9=src/drazin_skeleton.cpp
OBJ_FILES9=$(CPP_FILES9:.cpp=.o)

CPP_FILES10=src/drazin_canonical.cpp
OBJ_FILES10=$(CPP_FILES10:.cpp=.o)

CPP_FILES11=src/qreal_inverse.cpp
OBJ_FILES11=$(CPP_FILES11:.cpp=.o)

CPP_FILES12=src/breal_inverse.cpp
OBJ_FILES12=$(CPP_FILES12:.cpp=.o)

CPP_FILES13=src/qmulti_inverse.cpp
OBJ_FILES13=$(CPP_FILES13:.cpp=.o)

CPP_FILES14=src/make_test_matrix.cpp
OBJ_FILES14=$(CPP_FILES14:.cpp=.o)

CPP_FILES15=src/run_algo.cpp
OBJ_FILES15=$(CPP_FILES15:.cpp=.o)

EXE:= bin/matrix_inverse
EXE2:=bin/expr_parser
EXE3:=bin/deriv_matrix
EXE4:=bin/b_inverse
EXE5:=bin/q_inverse
EXE6:=bin/rank
EXE7:=bin/bq_inverse
EXE8:=bin/drazin_recursive
EXE9:=bin/drazin_skeleton
EXE10:=bin/drazin_canonical
EXE11:=bin/qreal_inverse
EXE12:=bin/breal_inverse
EXE13:=bin/qmulti_inverse
EXE14:=bin/make_test_matrix
EXE15:=bin/run_algo

##debug: CXXFLAGS += -DDEBUG -g
##debug: default
##debug: CXXFLAGS+=-DDEBUG -g
##debug: CXXFLAGS:=$(filter-out -O2,$(CXXFLAGS))
##debug: $(EXE10)

default: $(EXE) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(EXE13) $(EXE14) $(EXE15)

.PHONY: directories

$(EXE): directories $(OBJ_FILES1)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES1) $(LIBS)

$(EXE2): directories $(OBJ_FILES2)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES2) $(LIBS)

$(EXE3): directories $(OBJ_FILES3)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES3) $(LIBS)

$(EXE4): directories $(OBJ_FILES4)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES4) $(LIBS)

$(EXE5): directories $(OBJ_FILES5)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES5) $(LIBS)

$(EXE6): directories $(OBJ_FILES6)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES6) $(LIBS)

$(EXE7): directories $(OBJ_FILES7)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES7) $(LIBS)

$(EXE8): directories $(OBJ_FILES8)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES8) $(LIBS)

$(EXE9): directories $(OBJ_FILES9)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES9) $(LIBS)

$(EXE10): directories $(OBJ_FILES10)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES10) $(LIBS)

$(EXE11): directories $(OBJ_FILES11)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES11) $(LIBS)

$(EXE12): directories $(OBJ_FILES12)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES12) $(LIBS)

$(EXE13): directories $(OBJ_FILES13)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES13) $(LIBS)

$(EXE14): directories $(OBJ_FILES14)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES14) $(LIBS)

$(EXE15): directories $(OBJ_FILES15)
	$(CXX) $(LD_FLAGS) -o $@ $(OBJ_FILES15) $(LIBS)

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
	rm -rf .depend \
	rm -rf src/*.o

-include $(OBJ_FILES1:.o=.d)
-include $(OBJ_FILES2:.o=.d)
-include $(OBJ_FILES3:.o=.d)
-include $(OBJ_FILES4:.o=.d)
-include $(OBJ_FILES5:.o=.d)
-include $(OBJ_FILES6:.o=.d)
-include $(OBJ_FILES7:.o=.d)
-include $(OBJ_FILES8:.o=.d)
-include $(OBJ_FILES9:.o=.d)
-include $(OBJ_FILES10:.o=.d)
-include $(OBJ_FILES11:.o=.d)
-include $(OBJ_FILES12:.o=.d)
-include $(OBJ_FILES13:.o=.d)
-include $(OBJ_FILES14:.o=.d)
-include $(OBJ_FILES15:.o=.d)
