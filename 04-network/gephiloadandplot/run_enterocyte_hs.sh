#<key id="d15" for="node" attr.name="core" attr.type="boolean"/>
#<key id="d14" for="node" attr.name="module-svd-HS-1" attr.type="boolean"/>
#<key id="d11" for="node" attr.name="module-svd-HS-2" attr.type="boolean"/>
#<key id="d10" for="node" attr.name="module-svd-HS-3" attr.type="boolean"/>
#<key id="d16" for="node" attr.name="module-svd-HS-4" attr.type="boolean"/>
#<key id="d17" for="node" attr.name="module-svd-HS-7" attr.type="boolean"/>
#<key id="d9" for="node" attr.name="module-svd-HS-8" attr.type="boolean"/>
#<key id="d12" for="node" attr.name="module-svd-HS-9" attr.type="boolean"/>


# Define paths
gp="/home/casci/rionbr/spermnet/4-network/results/graphml"
cp="/home/casci/rionbr/spermnet/4-network/results/forceatlas2"
op="/home/casci/rionbr/spermnet/4-network/images/graphs"

# celltype
celltype="enterocyte"
# Network
network="thr"
# Threshold
threshold="0p5"
# Layer
layer="HS"
# Rotation
rotate=0

# Map
declare -A pairs=(
	[simple]=none
	[core]=d15
	[1]=d14
	[2]=d11
	[3]=d10
	[4]=d16
	#[5]=
	#[6]=d9
	[7]=d17
	[8]=d9
	[9]=d12
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
