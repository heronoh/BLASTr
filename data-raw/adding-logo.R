fs::file_copy("~/Pictures/logos/BLASTr/BLASTr-logo.svg", "data-raw/")

fs::file_copy("~/Pictures/logos/BLASTr/BLASTr-logo.png", "data-raw/")

usethis::use_logo("data-raw/BLASTr-logo.png")
