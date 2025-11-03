# ## BLASTr 0.1.8 (Development Version)

# Release Date: Unreleased

# Development Changelog: [dev](https://github.com/heronoh/BLASTr/compare/v0.1.7...HEAD)

# ### Added

# * Support for Windows dependencies.

# * New `install_blast_windows()` function to install BLAST+ on Windows.

### Changed

# * `parallel_blast()` will keep a temporary output directory for each thread
#   to avoid file writing conflicts.

# * `parallel_blast()` will automatically clean up temporary output directories
#   after the function completes.
