# The logo source lives in the metabarverse repository
# (logos/final/BLASTr/: SVG and PNGs of 240, 600, 1200 and 3000 px).
# use_logo() does not overwrite in non-interactive sessions; resize to 240 px wide:
img <- magick::image_read("~/prjcts/metabarverse/logos/final/BLASTr/BLASTr-logo-1200px.png")
magick::image_write(magick::image_resize(img, "240x278"), "man/figures/logo.png")
