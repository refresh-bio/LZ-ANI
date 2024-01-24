all: lz-ani-0.1

ANI_ROOT_DIR = .
ANI_MAIN_DIR = src
ANI_LIBS_DIR = .

MIMALLOC_INLUCDE_DIR = libs/mimalloc/include

CC 	= g++
CFLAGS	= -fPIC -Wall -O3 -m64 -std=c++20 -pthread -mavx -I $(ANI_LIBS_DIR) -I $(MIMALLOC_INLUCDE_DIR) -fpermissive
CLINK	= -lm -lpthread -O3 -std=c++20 -static-libgcc -static-libstdc++

MIMALLOC_OBJ=libs/mimalloc/mimalloc.o

$(MIMALLOC_OBJ):
	$(CC) -DMI_MALLOC_OVERRIDE -O3 -DNDEBUG -fPIC -Wall -Wextra -Wno-unknown-pragmas -fvisibility=hidden -Wstrict-prototypes -ftls-model=initial-exec -fno-builtin-malloc -std=gnu11 -c -I libs/mimalloc/include libs/mimalloc/src/static.c -o $(MIMALLOC_OBJ)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

lz-ani-0.1: $(ANI_MAIN_DIR)/lz-ani.o \
	$(ANI_MAIN_DIR)/filter.o \
	$(ANI_MAIN_DIR)/lz_matcher.o \
	$(ANI_MAIN_DIR)/parser.o \
	$(ANI_MAIN_DIR)/seq_reservoir.o \
	$(ANI_MAIN_DIR)/utils.o \
	$(MIMALLOC_OBJ)
	$(CC) -o $(ANI_ROOT_DIR)/$@  \
	$(MIMALLOC_OBJ) \
	$(ANI_MAIN_DIR)/lz-ani.o \
	$(ANI_MAIN_DIR)/filter.o \
	$(ANI_MAIN_DIR)/lz_matcher.o \
	$(ANI_MAIN_DIR)/parser.o \
	$(ANI_MAIN_DIR)/seq_reservoir.o \
	$(ANI_MAIN_DIR)/utils.o \
	$(CLINK)


clean:
	-rm $(ANI_MAIN_DIR)/*.o
	-rm $(ANI_LIBS_DIR)/*.o
	-rm $(MIMALLOC_OBJ)
	-rm lz-ani-0.1
	
