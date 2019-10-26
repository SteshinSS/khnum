cmake_minimum_required(VERSION 3.0)
project(glpk)
file(GLOB_RECURSE glpk_sources
                  glpk/amd/*
                  glpk/api/*
                  glpk/bflib/*
                  glpk/colamd/*
                  glpk/draft/*
                  glpk/env/*
                  glpk/intopt/*
                  glpk/minisat/*
                  glpk/misc/*
                  glpk/mpl/*
                  glpk/npp/*
                  glpk/proxy/*
                  glpk/simplex/*
                  glpk/zlib/*)
add_library(glpk STATIC ${glpk_sources})
set_target_properties(glpk PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(glpk PUBLIC ${PROJECT_SOURCE_DIR} glpk glpk/amd/
                                                             glpk/api
                                                             glpk/bflib
                                                             glpk/colamd
                                                             glpk/draft
                                                             glpk/env
                                                             glpk/intopt
                                                             glpk/minisat
                                                             glpk/misc
                                                             glpk/mpl
                                                             glpk/npp
                                                             glpk/proxy
                                                             glpk/simplex
                                                             glpk/zlib)
