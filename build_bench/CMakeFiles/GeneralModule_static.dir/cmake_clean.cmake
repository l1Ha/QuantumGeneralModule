file(REMOVE_RECURSE
  "lib/libgeneral_module.a"
  "lib/libgeneral_module.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/GeneralModule_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
