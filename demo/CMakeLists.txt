# ---> demo 1 (a rotating triangle)
add_executable(d0_hello_triangle d0_hello_triangle/main.cpp)
target_link_libraries(d0_hello_triangle PUBLIC util-compatible)
set_target_properties(d0_hello_triangle PROPERTIES CXX_STANDARD 11) 

# ---> demo 2 (famous glxgears)
add_executable(d1_glgears d1_glgears/main.cpp d1_glgears/draw.cpp)
target_link_libraries(d1_glgears PUBLIC util-compatible)
set_target_properties(d1_glgears PROPERTIES CXX_STANDARD 11) 

# ---> demo 3 (usage of imgui + opengl2)
add_executable(d2_imgui d2_imgui/main.cpp)
target_link_libraries(d2_imgui PUBLIC util-compatible)
set_target_properties(d2_imgui PROPERTIES CXX_STANDARD 11) 

# ---> demo 3 (usage of imgui + opengl3)
add_executable(d3_imgui d3_imgui/main.cpp)
target_link_libraries(d3_imgui PUBLIC util-core-3.3)
set_target_properties(d3_imgui PROPERTIES CXX_STANDARD 11) 
