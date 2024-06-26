rem %1 - $(SolutionDir)
rem %2 - $(Configuration)

cd %1\..\libs

@echo "zlib-ng"
cd zlib-ng 
cmake -B build-vs -S . -DZLIB_COMPAT=ON 
cmake --build build-vs --config %2
cd ..

@echo "isa-l"
cd isa-l
nmake -f Makefile.nmake
