cmake_minimum_required(VERSION 3.21)
project(
        Mirage
        VERSION 1.0.0
        HOMEPAGE_URL https://github.com/TravisWheelerLab/Mirage
)

set(CMAKE_C_STANDARD 17)
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
                ${CMAKE_SOURCE_DIR}/src/*.sh
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
)

add_custom_target(
        build-dependencies ALL
        COMMAND ./install-dependencies.sh
        VERBATIM
)

# Docker targets so that we don't have to type out the tag every time.

add_custom_target(
        build-image
        COMMAND docker build -t traviswheelerlab/mirage_build:latest -f Dockerfile_build .
        VERBATIM
)

add_custom_target(
        push-image
        COMMAND docker push traviswheelerlab/mirage_build:latest
        VERBATIM
)
