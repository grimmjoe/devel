include ../../Makefile.inc

CXXFLAGS+= 
LD_FLAGS+=-L../../core/lib -lcore
INCLUDE=-I./inc -I../ -I../../core/ -I../../core/inc/


CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
EXE:= bin/matrix_inverse

default: $(EXE)

.PHONY: directories

$(EXE): directories $(OBJ_FILES)
	$(CXX) $(LD_FLAGS)  -o $@ $(OBJ_FILES)

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

directories: bin obj

bin:
	mkdir -p bin

obj:
	mkdir -p obj

clean:
	rm -rf bin	\
	rm -rf obj

CXXFLAGS += -MMD
-include $(OBJ_FILES:.o=.d)
