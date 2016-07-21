FILE(REMOVE_RECURSE
  "CMakeFiles/make_clean"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/make_clean.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
