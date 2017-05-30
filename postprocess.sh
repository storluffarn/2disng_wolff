## convert .pbm files to png 
mogrify -format png ~/phd-stuff/courses/comp_phys/final_assignment/*.pbm

## .gif script, not recommended, inefficient
##convert -delay 1 -loop 0 *.png anim.gif

## .mp4 scripts
## use this script to search for specific name pattern
ffmpeg -framerate 50 -i image-%05d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
## use this script to search for all .ng files
##ffmpeg -framerate 50 -pattern_type glob -i '*.png' -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
## use this script if the number of pixels are odd (x264 for all it's greatness assumes even number of pixels)
##ffmpeg -framerate 50 -pattern_type glob -i '*.png' -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" output.mp4

mkdir p256256_T300_J0005_b05_s7d5k

mv *.mp4 *.dat *.csv image7499.png p256256_T300_J0005_b05_s7d5k/

rm *.png
rm *.pbm

