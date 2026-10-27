# Basic Funcs
# Thomams Johnson III
# 07-23-2023


produce_csv = function(created_object, base_str, file_name_str)
{
  write_csv(created_object, paste0(base_str, file_name_str))
  return(file_name_str)
}