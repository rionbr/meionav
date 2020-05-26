#<key id="d10" for="node" attr.name="core" attr.type="boolean"/>
#<key id="d9" for="node" attr.name="module-svd-MM-1" attr.type="boolean"/>
#<key id="d14" for="node" attr.name="module-svd-MM-2" attr.type="boolean"/>
#<key id="d22" for="node" attr.name="module-svd-MM-3" attr.type="boolean"/>
#<key id="d13" for="node" attr.name="module-svd-MM-4" attr.type="boolean"/>
#<key id="d15" for="node" attr.name="module-svd-MM-5" attr.type="boolean"/>
#<key id="d21" for="node" attr.name="module-svd-MM-6" attr.type="boolean"/>
#<key id="d18" for="node" attr.name="module-svd-MM-7" attr.type="boolean"/>
#<key id="d16" for="node" attr.name="module-svd-MM-8" attr.type="boolean"/>
#<key id="d19" for="node" attr.name="module-svd-MM-9" attr.type="boolean"/>


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
layer="MM"
# Rotation
rotate=270

# Map
declare -A pairs=(
	[simple]=none
	[core]=d10
	[1]=d9
	[2]=d14
	[3]=d22
	[4]=d13
	[5]=d15
	[6]=d21
	[7]=d18
	[8]=d16
	[9]=d19
	#[10]=
	#[11]=
	#[12]=
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
