cmake_minimum_required(VERSION 3.5)
project(sidh_search)

set(CMAKE_C_STANDARD 99)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_C_FLAGS_RELEASE -O3)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

include_directories(include)

file(GLOB HEADER_FILES CONFIGURE_DEPENDS "include/*.h")
set(SOURCE_FILES src/mont_utils.c src/printing.c)

find_library(GMP gmp)
find_library(SIKE_REF sike_ref lib)

add_executable(find_basis_p33 src/find_basis_p33.c ${SOURCE_FILES} ${HEADER_FILES})
add_executable(find_order src/find_order.c ${SOURCE_FILES} ${HEADER_FILES})
add_executable(explore_graph src/explore_graph.c ${SOURCE_FILES} ${HEADER_FILES})

target_link_libraries(find_order "${GMP}" "${SIKE_REF}")
target_link_libraries(find_basis_p33 "${GMP}" "${SIKE_REF}")
target_link_libraries(explore_graph "${GMP}" "${SIKE_REF}")

