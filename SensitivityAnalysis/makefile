CC = gcc
CFLAGS = -c -Wall
INC_FOLDERS = -I ../Librairies/include
LIB_FOLDERS = -L ../Librairies/lib
LIB_LINK =  ../Librairies/lib/libsundials_ida.a ../Librairies/lib/libsundials_nvecserial.a  ../Librairies/lib/libm.a ../Librairies/lib/liblapack.a

default: main

main: 	Model_Sensitivity.c
	$(CC) $(INC_FOLDERS) $(LIB_FOLDERS) -c Model_Sensitivity.c
	$(CC) $(LIB_FOLDERS) Model_Sensitivity.o $(LIB_LINK) -o output.o
 

clean:
	rm -rf *o output Model_Sensitivity


