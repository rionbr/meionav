#<key id="d12" for="node" attr.name="core" attr.type="boolean"/>
#<key id="d18" for="node" attr.name="module-svd-HS-1" attr.type="boolean"/>
#<key id="d17" for="node" attr.name="module-svd-HS-2" attr.type="boolean"/>
#<key id="d15" for="node" attr.name="module-svd-HS-3" attr.type="boolean"/>
#<key id="d16" for="node" attr.name="module-svd-HS-4" attr.type="boolean"/>
#<key id="d21" for="node" attr.name="module-svd-HS-5" attr.type="boolean"/>
#<key id="d9" for="node" attr.name="module-svd-HS-6" attr.type="boolean"/>
#<key id="d23" for="node" attr.name="module-svd-HS-7" attr.type="boolean"/>
#<key id="d14" for="node" attr.name="module-svd-HS-8" attr.type="boolean"/>
#<key id="d22" for="node" attr.name="module-svd-HS-9" attr.type="boolean"/>
#<key id="d20" for="node" attr.name="module-svd-HS-12" attr.type="boolean"/>


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
layer="HS"
# Rotation
rotate=90

# Map
declare -A pairs=(
	[simple]=none
	[core]=d12
	[1]=d18
	[2]=d17
	[3]=d15
	[4]=d16
	[5]=d21
	[6]=d9
	[7]=d23
	[8]=d14
	[9]=d22
	#[10]=
	#[11]=
	[12]=20
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
