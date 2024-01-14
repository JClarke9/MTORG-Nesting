# Function to run all scripts in a directory
run_scripts_in_directory <- function(directory_path) {
  # Get a list of all files in the directory
  files <- list.files(path = directory_path, pattern = "\\.R$", full.names = TRUE)
  
  # Run each script in the directory
  for (script_file in files) {
    cat(paste("Running script:", script_file, "\n"))
    
    # Use source() to execute the script
    source(script_file, local = TRUE)
    
    cat(paste("Script execution completed for:", script_file, "\n\n"))
  }
}

# Replace 'your_directory_path' with the actual path to your subdirectory
directory_path <- 'your_directory_path'
run_scripts_in_directory(directory_path)
