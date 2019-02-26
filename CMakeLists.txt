cmake_minimum_required( VERSION 3.0 )

project( luscher_zeta )

set( CMAKE_C_STANDARD 99 )

add_executable( luscher_zeta.x main.c )
target_link_libraries( luscher_zeta.x -lcerf -lcubature -lm )