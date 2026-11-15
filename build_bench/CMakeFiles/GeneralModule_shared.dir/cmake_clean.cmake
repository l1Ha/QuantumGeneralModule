file(REMOVE_RECURSE
  "lib/libgeneral_module.pdb"
  "lib/libgeneral_module.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/GeneralModule_shared.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
