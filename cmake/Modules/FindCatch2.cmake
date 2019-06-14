# Copyright (c) 2018, K. Kumar (me@kartikkumar.com)

if (CATCH2_INCLUDE_DIRS)
  # in cache already
  set(CATCH2_FOUND TRUE)
else (CATCH2_INCLUDE_DIRS)

  find_path(CATCH2_INCLUDE_DIR
    NAMES
      catch.hpp
      catch.h
    PATHS
      /usr/include
      /usr/local/include
      /opt/local/include
      /sw/include
      C:
      ${PROJECT_PATH}/..
      ${EXTERNAL_PATH}
    PATH_SUFFIXES
      Catch2 Catch2/single_include Catch2/single_include/catch2 Catch2/src/catch2/single_include
  )

  if(CATCH2_INCLUDE_DIR)
    set(CATCH2_INCLUDE_DIRS
        ${CATCH2_INCLUDE_DIR}/..
    )
  endif(CATCH2_INCLUDE_DIR)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Catch2 DEFAULT_MSG CATCH2_INCLUDE_DIRS)

  # show the CATCH2_INCLUDE_DIRS variables only in the advanced view
  mark_as_advanced(CATCH2_INCLUDE_DIRS)

endif (CATCH2_INCLUDE_DIRS)