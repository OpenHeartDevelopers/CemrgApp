set(CMAKE_AUTOMOC ON)

set(module_dirs
  MyCemrgLib
)

foreach(module_dir ${module_dirs})
  add_subdirectory(${module_dir})
endforeach()
