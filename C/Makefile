CC=nvcc
CFLAGS=-std=c++11 `pkg-config --cflags igraph`
INC=
LIB=`pkg-config --libs igraph`
BIN=
SRC=QuantifyingInherentGraphAnonymity.c
OBJ=$(SRC:.cpp=.o)
EXEC=pj3cpu
NVFLAGS=-arch=sm_30


all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CLFAGS) $(NVFLAGS) $(INC) -o $@ $^ $(LIB)
	
# .cpp.o:
# 	$(CC) $(CFLAGS) $(INC) -c -o $@ $^

clean:
	rm -f $(EXEC)
