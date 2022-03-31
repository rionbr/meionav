# Remove build
echo "-- Removing old build --"
rm -f *.class
rm -f images/debug.pdf
rm -f images/debug.jpg

# Build
echo "-- Build Main.java --"
javac -cp gephi-toolkit-0.9.2-all.jar PlotNetSingleMod.java

# Run
echo "-- Run Main.class --"
java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input debug/net-debug.graphml --coords debug/net-debug-forceatlas2.txt --output images/debug.pdf --rotate 120 --highlight d9 --highlight-color '#00ff00'

# Pdf to Jpg
echo "-- Converts pdf to jpg --"
#convert -density 300 -resize 3300x2550 output.pdf output.jpg
convert -density 72 -resize 1024x768 images/debug.pdf images/debug.jpg