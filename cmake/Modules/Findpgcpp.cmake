# Copyright (c) <year> <author> (<email>)
# Distributed under the MIT License.
# See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT

if(PGCPP_INCLUDE_DIRS)
  # in cache already
  set(PGCPP_FOUND TRUE)
else(PGCPP_INCLUDE_DIRS)

  find_path(PGCPP_INCLUDE_DIR
    NAMES
      pgcpp.hpp
    PATHS
      /usr/include
      /usr/local/include
      /opt/local/include
      /sw/include
      ${PROJECT_PATH}/..
      ${EXTERNAL_PATH}
    PATH_SUFFIXES
      pgcpp include/pgcpp
  )

  if(PGCPP_INCLUDE_DIR)
    set(PGCPP_INCLUDE_DIRS
      ${PGCPP_INCLUDE_DIR}/..
    )
  endif(PGCPP_INCLUDE_DIR)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(pgcpp DEFAULT_MSG PGCPP_INCLUDE_DIRS)

  # show the PGCPP_INCLUDE_DIRS variables only in the advanced view
  mark_as_advanced(PGCPP_INCLUDE_DIRS)

endif(PGCPP_INCLUDE_DIRS)