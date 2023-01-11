add_executable(run_p0
  main.cpp 
  draw.cpp)
target_link_libraries(run_p0 PUBLIC 
  util-compatible 
  ${gl_libs})
set_target_properties(run_p0 PROPERTIES 
  CXX_STANDARD 11) 
