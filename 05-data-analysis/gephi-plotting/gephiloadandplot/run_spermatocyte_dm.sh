#<key id="d6" for="node" attr.name="conserved" attr.type="boolean"/>

#<key id="d13" for="node" attr.name="module-pca-DM-1" attr.type="boolean"/>
#<key id="d12" for="node" attr.name="module-pca-DM-2" attr.type="boolean"/>
#<key id="d10" for="node" attr.name="module-pca-DM-3" attr.type="boolean"/>
#<key id="d11" for="node" attr.name="module-pca-DM-4" attr.type="boolean"/>
#<key id="d14" for="node" attr.name="module-pca-DM-5" attr.type="boolean"/>
#<key id="d8" for="node" attr.name="module-pca-DM-6" attr.type="boolean"/>
#<key id="d9" for="node" attr.name="module-pca-DM-11" attr.type="boolean"/>
#<key id="d7" for="node" attr.name="module-pca-DM-12" attr.type="boolean"/>


# Define paths
gp="/home/casci/rionbr/spermnet/4-network/results/graphml"
cp="/home/casci/rionbr/spermnet/4-network/results/forceatlas2"
op="/home/casci/rionbr/spermnet/4-network/images/graphs"

# celltype
celltype="spermatocyte"
# Network
network="thr"
# Threshold
threshold="0p5"
# Layer
layer="DM"
# Rotation
rotate=15

# Map
declare -A pairs=(
	[simple]=none
	[conserved]=d6
	[1]=d13
	[2]=d12
	[3]=d10
	[4]=d11
	[5]=d14
	#[6]=
	#[7]=d10
	#[8]=d26
	#[9]=d18
	#[10]=d28
	[11]=d9
	[12]=d7
)
for module in "${!pairs[@]}";
do
	echo "-- Plotting $celltype-$network-$threshold-$layer-mod-$module --"
	#
	highlight="${pairs[$module]}"
	#
	input="$gp/$celltype/net-$celltype-$network-$threshold-$layer.graphml"
	coords="$cp/$celltype/net-$celltype-$network-$threshold-$layer-forceatlas2.txt"
	#
	pdfpath="$op/$celltype/$layer/pdf"
	jpgpath="$op/$celltype/$layer/jpg"
	#
	pdf="$pdfpath/net-$celltype-$network-$threshold-$layer-mod-$module.pdf"
	jpg="$jpgpath/net-$celltype-$network-$threshold-$layer-mod-$module.jpg"
	
	mkdir -p $pdfpath
	mkdir -p $jpgpath

	# Run GephiLoadAndPlot
	java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input $input --coords $coords --output $pdf --rotate $rotate --highlight $highlight
	
	# Convert PDF to JPG
	convert -density 300 -resize 3300x2550 $pdf $jpg
done
