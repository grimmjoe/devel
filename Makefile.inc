CXX=g++
CXXFLAGS=-Wall -std=c++11
AR=ar
AR_FLAGS=
RANLIB=ranlib
RM=rm
ECHO=echo
##TBB_PATH=/home/hamlet/thesis/devel/3rd_party/tbb43_20150316oss
TBB_PATH=/usr
LD_FLAGS=
LIBS+=-ltbb -ltbbmalloc_proxy -ltbbmalloc
##LIBS+=-L$(TBB_PATH)/lib -ltbb -ltbbmalloc_proxy -ltbbmalloc
##INCLUDE+=-I$(TBB_PATH)/include

