## Download of FF parameters from Tinker repository
file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/release/params/amoebabio18.prm
     ${CMAKE_BINARY_DIR}/amoebabio18.prm
     EXPECTED_HASH MD5=18b942176d18f77e5c10d3ed13490f7b) 
file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/release/params/amoeba09.prm
     ${CMAKE_BINARY_DIR}/amoeba09.prm
     EXPECTED_HASH MD5=b28c0e0b33b5bfe790c16847622a5e0b) 
file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/release/params/amber99sb.prm
     ${CMAKE_BINARY_DIR}/amber99sb.prm
     EXPECTED_HASH MD5=e9330e340a3f4045cc8025dc9ce42e33) 

file(CREATE_LINK ${CMAKE_SOURCE_DIR}/tests ${CMAKE_BINARY_DIR}/tests SYMBOLIC)
set(VALGRINDOPT --track-origins=yes --leak-check=full --show-leak-kinds=all) 
set(TESTLANG "C" CACHE STRING "Language to be used for the test suite, the very same test programs are implemented in C99 (C) and Fortran2003 (F03).")


include(tests/TestsCmake.txt)
