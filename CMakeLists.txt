cmake_minimum_required(VERSION 3.12)
project(
        Mirage
        VERSION 2.0.0
        HOMEPAGE_URL https://github.com/TravisWheelerLab/Mirage
)

set(CMAKE_C_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY build)

add_executable(
        MultiSeqNW

        src/MultiSeqNW.c
)
target_link_libraries(MultiSeqNW m)
install(TARGETS MultiSeqNW RUNTIME)

add_executable(
        ExonWeaver

        src/ExonWeaver.c
        src/BasicBio.c
)
target_link_libraries(ExonWeaver m)
install(TARGETS ExonWeaver RUNTIME)

add_executable(
        FastMap2

        src/FastMap2.c
        src/BasicBio.c
)
target_link_libraries(FastMap2 m)
install(TARGETS FastMap2 RUNTIME)

# Custom target for handling dependencies, Perl scripts, and shell scripts. We
# just move all these files to the build directory to aid in distribution and
# installation. These will be automatically run by the build.

add_custom_target(
        build-perl ALL
        COMMAND cp
                ${CMAKE_SOURCE_DIR}/src/*.pl
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
        COMMAND cp
                ${CMAKE_SOURCE_DIR}/src/*.pm
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
)

add_custom_target(
        build-shell ALL
	COMMAND cp
		${CMAKE_SOURCE_DIR}/src/run_mirage2.sh
		${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/mirage2
 	COMMAND ln
 		-s ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/mirage2
 		   ${CMAKE_SOURCE_DIR}/mirage2
	COMMAND cp
		${CMAKE_SOURCE_DIR}/src/run_diviner.sh
		${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/diviner
 	COMMAND ln
 		-s ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/diviner
 		   ${CMAKE_SOURCE_DIR}/diviner
	COMMAND cp
		${CMAKE_SOURCE_DIR}/src/run_bazaar.sh
		${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bazaar
)

add_custom_target(
        build-dependencies ALL
        COMMAND ${CMAKE_SOURCE_DIR}/src/install-dependencies.sh
        VERBATIM
)
