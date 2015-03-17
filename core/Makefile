##CXX=g++
##CXXFLAGS=-Wall -fPIC -std=c++11
##LD_FLAGS=-shared
include ../Makefile.inc
CXXFLAGS+= -fPIC -O2
LD_FLAGS+=-shared
INCLUDE=-Iinc/

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
LIBRARY:= lib/libcore.so

##debug: CXXFLAGS += -DDEBUG -g
##debug: default

default: $(LIBRARY)

.PHONY: directories

$(LIBRARY): directories $(OBJ_FILES)
	$(CXX) $(LD_FLAGS)  -o $@ $(OBJ_FILES)

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

directories: lib obj

lib:
	mkdir -p lib

obj:
	mkdir -p obj

clean:
	rm -rf lib	\
	rm -rf obj

CXXFLAGS += -MMD
-include $(OBJ_FILES:.o=.d)
