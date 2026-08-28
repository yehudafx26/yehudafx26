find_package(Matlab REQUIRED)

if (Matlab_FOUND)
  matlab_add_mex(
      NAME yehuMex
      SRC yehu_mex.cpp
      LINK_TO yehu
  )
else()
  message(WARNING "MATLAB root directory not found. Failed to build MATLAB bindings.")
endif ()