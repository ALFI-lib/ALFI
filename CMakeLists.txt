cmake_minimum_required(VERSION 3.22)

project(ALFI
	VERSION 0.0.0.0
	DESCRIPTION "Advanced Library for Function Interpolation"
	HOMEPAGE_URL "https://alfi-lib.github.io/ALFI/"
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -march=native ${ALFI_EXTRA_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -Og ${ALFI_EXTRA_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_SANITIZE "-g3 -Og -fsanitize=address,leak,undefined ${ALFI_EXTRA_CXX_FLAGS_SANITIZE}")
set(CMAKE_CXX_FLAGS_FASTBUILD "-g0 -O0 -DNDEBUG ${ALFI_EXTRA_CXX_FLAGS_FASTBUILD}")
set(CMAKE_CXX_FLAGS_RELEASE "-g0 -s -O3 -flto=auto -DNDEBUG ${ALFI_EXTRA_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g3 -O2 -flto=auto -DNDEBUG ${ALFI_EXTRA_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-g0 -s -Oz -flto=auto -DNDEBUG ${ALFI_EXTRA_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_FAST "-g0 -s -Ofast -fno-finite-math-only -flto=auto -DNDEBUG ${ALFI_EXTRA_CXX_FLAGS_FAST}")
set(CMAKE_CXX_FLAGS_FASTPARALLEL "-g0 -s -Ofast -fno-finite-math-only -flto=auto -fopenmp -DNDEBUG ${ALFI_EXTRA_CXX_FLAGS_FASTPARALLEL}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ALFI_ENABLE_EXAMPLES "Enable examples subdirectory for the ALFI library" ${PROJECT_IS_TOP_LEVEL})
option(ALFI_ENABLE_TESTS "Enable tests subdirectory for the ALFI library" ${PROJECT_IS_TOP_LEVEL})
option(ALFI_ENABLE_BENCHES "Enable benches subdirectory for the ALFI library" ${PROJECT_IS_TOP_LEVEL})

add_subdirectory(ALFI)

if(ALFI_ENABLE_EXAMPLES)
	add_subdirectory(examples)
endif()

if(ALFI_ENABLE_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()

if(ALFI_ENABLE_BENCHES)
	add_subdirectory(benches)
endif()