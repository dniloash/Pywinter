cmake_minimum_required(VERSION 3.18)

project(pywinter
    LANGUAGES C Fortran
)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c
      "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(f2py_module_name "KreatE_inter_m_f")
set(f_src_files 
    "${CMAKE_SOURCE_DIR}/src/fortran/KreatE_inter_m_f.f90")
set(f2py_module_c "${f2py_module_name}module.c")
set(f2py_wrap_f "${f2py_module_name}-f2pywrappers.f")

set(py_src_files
    "${CMAKE_SOURCE_DIR}/src/pywinter/__init__.py"
    "${CMAKE_SOURCE_DIR}/src/pywinter/winter.py")

add_custom_command(
    OUTPUT "${f2py_module_c}" "${f2py_wrap_f}"
    COMMAND Python::Interpreter -m numpy.f2py --lower
                -m "${f2py_module_name}"
                "${f_src_files}"
    COMMAND cmake -E touch "${f2py_wrap_f}"
    DEPENDS "${f_src_files}" # Fortran source
)

python_add_library(${f2py_module_name} MODULE WITH_SOABI
    "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
    "${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrap_f}"
    "${F2PY_INCLUDE_DIR}/fortranobject.c"
    "${f_src_files}")

target_include_directories(${f2py_module_name} PUBLIC ${F2PY_INCLUDE_DIR})
target_link_libraries(${f2py_module_name} PUBLIC Python::NumPy)

install(TARGETS ${f2py_module_name} DESTINATION pywinter)
