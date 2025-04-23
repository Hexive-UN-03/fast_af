# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/hexive/code/fast_af/out/build/fast_pre/_deps/eagle-src"
  "/home/hexive/code/fast_af/out/build/fast_pre/_deps/eagle-build"
  "/home/hexive/code/fast_af/out/build/fast_pre/_deps/eagle-subbuild/eagle-populate-prefix"
  "/home/hexive/code/fast_af/out/build/fast_pre/_deps/eagle-subbuild/eagle-populate-prefix/tmp"
  "/home/hexive/code/fast_af/out/build/fast_pre/_deps/eagle-subbuild/eagle-populate-prefix/src/eagle-populate-stamp"
  "/home/hexive/code/fast_af/out/build/fast_pre/_deps/eagle-subbuild/eagle-populate-prefix/src"
  "/home/hexive/code/fast_af/out/build/fast_pre/_deps/eagle-subbuild/eagle-populate-prefix/src/eagle-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/hexive/code/fast_af/out/build/fast_pre/_deps/eagle-subbuild/eagle-populate-prefix/src/eagle-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/hexive/code/fast_af/out/build/fast_pre/_deps/eagle-subbuild/eagle-populate-prefix/src/eagle-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
