## Download of FF parameters from Tinker repository
# Commit of Tinker Tools used to retrive FF (they change sometimes, so it's
# better to take them always from the same commit and not from release for 
# test purpose.
set(PRM_COMMIT 8a8098d10864348ecfcd062561e4ec4405bac3b1)

file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/${PRM_COMMIT}/params/amoebabio18.prm
     ${CMAKE_BINARY_DIR}/amoebabio18.prm
     EXPECTED_HASH MD5=18b942176d18f77e5c10d3ed13490f7b) 
file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/${PRM_COMMIT}/params/amoeba09.prm
     ${CMAKE_BINARY_DIR}/amoeba09.prm
     EXPECTED_HASH MD5=69bab8d2c8b74f1488465541b7811b2f) 
file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/${PRM_COMMIT}/params/amber99sb.prm
     ${CMAKE_BINARY_DIR}/amber99sb.prm
     EXPECTED_HASH MD5=e9330e340a3f4045cc8025dc9ce42e33) 

file(CREATE_LINK ${CMAKE_SOURCE_DIR}/tests ${CMAKE_BINARY_DIR}/tests SYMBOLIC)
set(VALGRINDOPT --track-origins=yes --leak-check=full --show-leak-kinds=all) 

include(tests/test_programs/CMakeLists.txt)
include(tests/TestsCmake.txt)
