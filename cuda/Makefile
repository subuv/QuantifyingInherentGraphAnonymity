CC=nvcc
CFLAGS=-std=c++11 `pkg-config --cflags igraph`
INC=
LIB=`pkg-config --libs igraph`
BIN=
SRC=QuantifyingInherentGraphAnonymity.cu
OBJ=$(SRC:.cpp=.o)
EXEC=LCMGpu
NVFLAGS=-arch=sm_30


all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CLFAGS) $(NVFLAGS) $(INC) -o $@ $^ $(LIB)
	
# .cpp.o:
# 	$(CC) $(CFLAGS) $(INC) -c -o $@ $^

clean:
	rm -f $(EXEC)
