all: lz-ani-0.1

ANI_ROOT_DIR = .
ANI_MAIN_DIR = src
ANI_LIBS_DIR = .

CC 	= /usr/local/gcc62/bin/g++
CFLAGS	= -Wall -O3 -m64 -std=c++14 -pthread -mavx -I $(ANI_LIBS_DIR)
CLINK	= -lm -O3 -std=c++14 -static -pthread -mavx -fabi-version=6 

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

ani-0.1: $(ANI_MAIN_DIR)/lz-ani.o \
	$(ANI_MAIN_DIR)/worker.o \
	$(ANI_MAIN_DIR)/s_worker.o
	$(CC) $(CLINK) -o $(ANI_ROOT_DIR)/$@  \
	$(ANI_MAIN_DIR)/lz-ani.o \
	$(ANI_MAIN_DIR)/worker.o \
	$(ANI_MAIN_DIR)/s_worker.o


clean:
	-rm $(ANI_MAIN_DIR)/*.o
	-rm $(ANI_LIBS_DIR)/*.o
	-rm lz-ani-0.1
	
