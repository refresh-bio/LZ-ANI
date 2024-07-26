all: lz-ani

LZANI_ROOT_DIR = .
LZANI_MAIN_DIR = src
LZANI_LIBS_DIR = libs
ISAL_DIR = libs/isa-l
ZLIB_DIR = libs/zlib-ng

INC_DIRS =. libs/mimalloc/include libs/zlib-ng/ libs/isa-l/include
INCLUDE_DIR=$(foreach d, $(INC_DIRS), -I$d)

MIMALLOC_INLUCDE_DIR = libs/mimalloc/include

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
	uname_M := "x86_64"
else                          # If uname not available => 'not'
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
	uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
endif

NASM_V := $(shell nasm --version 2>/dev/null)

ifeq ($(PLATFORM), arm8)
$(info *** ARMv8 with NEON extensions ***)
	ARCH_FLAGS := -march=armv8-a  -DARCH_ARM
else ifeq ($(PLATFORM), m1)
$(info *** Apple M1(or never) with NEON extensions ***)
	ARCH_FLAGS := -march=armv8.4-a  -DARCH_ARM
else ifeq ($(PLATFORM), sse2)
$(info *** x86-64 with SSE2 extensions ***)
	ARCH_FLAGS := -msse2 -m64 -DARCH_X64 
else ifeq ($(PLATFORM), avx)
$(info *** x86-64 with AVX extensions ***)
	ARCH_FLAGS := -mavx -m64  -DARCH_X64
else ifeq ($(PLATFORM), avx2)
$(info *** x86-64 with AVX2 extensions ***)
	ARCH_FLAGS := -mavx2 -m64  -DARCH_X64
else
$(info *** Unspecified platform - use native compilation)
	ifeq ($(uname_M),x86_64)
		ARCH_FLAGS := -march=native -DARCH_X64
	else
		ARCH_FLAGS := -march=native -DARCH_ARM
	endif	
endif

CFLAGS	= -fPIC -static -pthread -Wall -O3 -std=c++20 $(ARCH_FLAGS) $(INCLUDE_DIR) -fpermissive
CLINK	= -lm


ifeq ($(uname_S),Linux)
	CLINK+=-fabi-version=6
	CLINK+=-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
endif

ifeq ($(uname_S),Darwin)
	CLINK += -lpthread -static-libgcc
endif

ifeq ($(uname_M),x86_64)
	ifdef NASM_V
		GZ_LIB:=isa-l.a
		gz_target:=isa-l
		CFLAGS+=-DREFRESH_USE_IGZIP
	else
		GZ_LIB:=libz.a
		gz_target:=ng_zlib
		CFLAGS+=-DREFRESH_USE_ZLIB
	endif
else
	GZ_LIB:=libz.a
	gz_target:=ng_zlib
	CFLAGS+=-DREFRESH_USE_ZLIB
endif

MIMALLOC_OBJ=libs/mimalloc/mimalloc.o



$(MIMALLOC_OBJ):
	$(CC) -DMI_MALLOC_OVERRIDE -O3 -DNDEBUG -fPIC -Wall -Wextra -Wno-unknown-pragmas -fvisibility=hidden -Wstrict-prototypes -ftls-model=initial-exec -fno-builtin-malloc -std=gnu11 -c -I libs/mimalloc/include libs/mimalloc/src/static.c -o $(MIMALLOC_OBJ)

%.o: %.cpp $(gz_target)
	$(CXX) $(CFLAGS) -c $< -o $@

ng_zlib:
	cd $(ZLIB_DIR) && ./configure --zlib-compat && $(MAKE) libz.a
	cp $(ZLIB_DIR)/libz.* $(LZANI_LIBS_DIR)

isa-l:
	cd $(ISAL_DIR) && $(MAKE) -f Makefile.unx
	cp $(ISAL_DIR)/bin/isa-l.a $(LZANI_LIBS_DIR)
	cp $(ISAL_DIR)/bin/libisal.* $(LZANI_LIBS_DIR)

lz-ani: $(gz_target) \
	$(LZANI_MAIN_DIR)/lz-ani.o \
	$(LZANI_MAIN_DIR)/filter.o \
	$(LZANI_MAIN_DIR)/lz_matcher.o \
	$(LZANI_MAIN_DIR)/parser.o \
	$(LZANI_MAIN_DIR)/seq_reservoir.o \
	$(LZANI_MAIN_DIR)/utils.o \
	$(MIMALLOC_OBJ)
	$(CXX) -o $(LZANI_ROOT_DIR)/$@  \
	$(MIMALLOC_OBJ) \
	$(LZANI_MAIN_DIR)/lz-ani.o \
	$(LZANI_MAIN_DIR)/filter.o \
	$(LZANI_MAIN_DIR)/lz_matcher.o \
	$(LZANI_MAIN_DIR)/parser.o \
	$(LZANI_MAIN_DIR)/seq_reservoir.o \
	$(LZANI_MAIN_DIR)/utils.o \
	$(LZANI_LIBS_DIR)/$(GZ_LIB) \
	$(CLINK)


clean:
	-rm $(LZANI_MAIN_DIR)/*.o
	-rm $(LZANI_LIBS_DIR)/*.o
	-rm $(MIMALLOC_OBJ)
	cd $(ZLIB_DIR) && $(MAKE) -f Makefile.in clean
	cd $(ISAL_DIR) && $(MAKE) -f Makefile.unx clean
	-rm lz-ani
	-rm $(LZANI_LIBS_DIR)/libz.*
	-rm $(LZANI_LIBS_DIR)/isa-l.*
	-rm $(LZANI_LIBS_DIR)/libisal.*
