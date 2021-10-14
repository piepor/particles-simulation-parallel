#!/bin/bash

ffmpeg -framerate 10 -i stage%03d.jpg -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4;
