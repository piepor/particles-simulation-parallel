#!/bin/csh
# This needs Image Magick tools
foreach pp (`ls *.ppm`)
   convert $pp ${pp:r}.jpg
end
