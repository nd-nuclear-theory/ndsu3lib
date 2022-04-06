###############################################################################
# Defintions of FetchContent and function to fetch if package not found.
#
# Anna E. McCoy
# Institute for Nuclear Theory 
#
# SPDX-License-Identifier: MIT
#
# 3/29/22 (aem) : Created
###############################################################################
cmake_minimum_required(VERSION 3.20)

include(FetchContent)

###########################################################################
# From where do you fetch
###########################################################################
FetchContent_Declare(
  wigxjpf
  GIT_REPOSITORY git@github.com:nd-nuclear-theory/wigxjpf.git
  # Working commit as of (3/29/22) 15a6330dfe48c4fc85d6592fb892ef1468617fa5 
  GIT_TAG        main 
  GIT_SHALLOW    TRUE
)

FetchContent_Declare(
  mpfun20
  GIT_REPOSITORY git@github.com:nd-nuclear-theory/mpfun20-fort.git
  GIT_TAG        master #0628af82d02c2e07e1e24f8b6005aef3c58b6cd1
  GIT_SHALLOW    TRUE
)
