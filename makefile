all: lz-ani-0.1

ANI_ROOT_DIR = .
ANI_MAIN_DIR = src
ANI_LIBS_DIR = .

CC 	= g++
CFLAGS	= -fPIC -Wall -O3 -m64 -std=c++17 -pthread -mavx -I $(ANI_LIBS_DIR) -fpermissive
CLINK	= -lm -lpthread -O3 -std=c++17 -static-libgcc

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

lz-ani-0.1: $(ANI_MAIN_DIR)/lz-ani.o \
	$(ANI_MAIN_DIR)/worker.o \
	$(ANI_MAIN_DIR)/s_worker.o
	$(CC) -o $(ANI_ROOT_DIR)/$@  \
	$(ANI_MAIN_DIR)/lz-ani.o \
	$(ANI_MAIN_DIR)/worker.o \
	$(ANI_MAIN_DIR)/s_worker.o \
	$(CLINK)


clean:
	-rm $(ANI_MAIN_DIR)/*.o
	-rm $(ANI_LIBS_DIR)/*.o
	-rm lz-ani-0.1
	
