# HS
# <key id="d6" for="node" attr.name="conserved" attr.type="boolean"/>
# MM
# <key id="d6" for="node" attr.name="conserved" attr.type="boolean"/>
# DM
# <key id="d6" for="node" attr.name="conserved" attr.type="boolean"/>

#echo "-- Build Main.java --"
#javac -cp gephi-toolkit-0.9.2-all.jar PlotNetSingleMod.java

# Define paths
gp="/home/casci/rionbr/spermnet/05-data-analysis/gephi-plotting/results/graphml"
cp="/home/casci/rionbr/spermnet/05-data-analysis/gephi-plotting/results/forceatlas2"
op="/home/casci/rionbr/spermnet/05-data-analysis/gephi-plotting/images/conserved"


# celltype
celltype="enterocyte"
# Network
network="thr"
# Threshold
threshold="0p5"



# Layer
layer="HS"
# Rotation
rotate=345
#
echo "-- Plotting $celltype-$network-$threshold-$layer-conserved --"
#
highlight="d6"
#
input="$gp/$celltype/net-$celltype-$network-$threshold-$layer.graphml"
coords="$cp/$celltype/net-$celltype-$network-$threshold-$layer-forceatlas2.txt"
#
pdfpath="$op/$celltype/$layer"
jpgpath="$op/$celltype/$layer"
#
pdf="$pdfpath/net-$celltype-$network-$threshold-$layer-conserved.pdf"
jpg="$jpgpath/net-$celltype-$network-$threshold-$layer-conserved.jpg"

mkdir -p $pdfpath
mkdir -p $jpgpath

# Run GephiLoadAndPlot
java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input $input --coords $coords --output $pdf --rotate $rotate --highlight $highlight

# Convert PDF to JPG
convert -density 300 -resize 3300x2550 $pdf $jpg



# Layer
layer="MM"
# Rotation
rotate=215
#
echo "-- Plotting $celltype-$network-$threshold-$layer-conserved --"
#
highlight="d6"
#
input="$gp/$celltype/net-$celltype-$network-$threshold-$layer.graphml"
coords="$cp/$celltype/net-$celltype-$network-$threshold-$layer-forceatlas2.txt"
#
pdfpath="$op/$celltype/$layer"
jpgpath="$op/$celltype/$layer"
#
pdf="$pdfpath/net-$celltype-$network-$threshold-$layer-conserved.pdf"
jpg="$jpgpath/net-$celltype-$network-$threshold-$layer-conserved.jpg"

mkdir -p $pdfpath
mkdir -p $jpgpath

# Run GephiLoadAndPlot
java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input $input --coords $coords --output $pdf --rotate $rotate --highlight $highlight

# Convert PDF to JPG
convert -density 300 -resize 3300x2550 $pdf $jpg



# Layer
layer="DM"
# Rotation
rotate=195
#
echo "-- Plotting $celltype-$network-$threshold-$layer-conserved --"
#
highlight="d6"
#
input="$gp/$celltype/net-$celltype-$network-$threshold-$layer.graphml"
coords="$cp/$celltype/net-$celltype-$network-$threshold-$layer-forceatlas2.txt"
#
pdfpath="$op/$celltype/$layer"
jpgpath="$op/$celltype/$layer"
#
pdf="$pdfpath/net-$celltype-$network-$threshold-$layer-conserved.pdf"
jpg="$jpgpath/net-$celltype-$network-$threshold-$layer-conserved.jpg"

mkdir -p $pdfpath
mkdir -p $jpgpath

# Run GephiLoadAndPlot
java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input $input --coords $coords --output $pdf --rotate $rotate --highlight $highlight

# Convert PDF to JPG
convert -density 300 -resize 3300x2550 $pdf $jpg