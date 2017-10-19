CC = gcc
CFLAGS = -c -Wall
INC_FOLDERS = -I ./Librairies/include
LIB_FOLDERS = -L ./Librairies/lib 
LIB_LINK = ./Librairies/lib/libsundials_ida.a  ./Librairies/lib/libsundials_nvecserial.a

default: main

main: 	Model.c
	$(CC) $(INC_FOLDERS) -c Model.c
	$(CC) $(LIB_FOLDERS) Model.o $(LIB_LINK) -lm -llapack -o output.exe
clean:
	rm -rf *o output.exe Model data.txt


