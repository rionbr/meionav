# Remove build
echo "-- Removing old build and output --"
rm -f *.class

# Build
echo "-- Build Main.java --"
javac -cp gephi-toolkit-0.9.2-all.jar PlotNetSingleMod.java

# Run
echo "-- Running to define layout --"

# Define paths
gp="/home/casci/rionbr/spermnet/4-network/results/graphml"
cp="/home/casci/rionbr/spermnet/4-network/results/forceatlas2"
op="/home/casci/rionbr/spermnet/4-network/images/graphs"


# Network
network="thr"
# Threshold
threshold="0p5"
# Module
module="simple"

for celltype in "enterocyte" "spermatocyte"; do

	for layer in "DM" "MM" "HS"; do
	
		echo "-- Plotting $celltype-$network-$threshold-$layer-mod-$module --"
		#
		input="$gp/$celltype/net-$celltype-$network-$threshold-$layer.graphml"
		coords="$cp/$celltype/net-$celltype-$network-$threshold-$layer-forceatlas2.txt"
		#
		pdfpath="$op/$celltype/pdf"
		jpgpath="$op/$celltype/jpg"
		#
		pdf="$pdfpath/net-$celltype-$network-$threshold-$layer-mod-$module.pdf"
		jpg="$jpgpath/net-$celltype-$network-$threshold-$layer-mod-$module.jpg"
	
		mkdir -p $pdfpath
		mkdir -p $jpgpath

		# Run GephiLoadAndPlot
		java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input $input --coords $coords --output $pdf --rotate 0 --highlight none
	
		# Convert PDF to JPG
		convert -density 300 -resize 3300x2550 $pdf $jpg
		
		# Delete pdf file
		rm $pdf
	done
done