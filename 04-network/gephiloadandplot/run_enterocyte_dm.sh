#<key id="d13" for="node" attr.name="core" attr.type="boolean"/>
#<key id="d9" for="node" attr.name="module-svd-DM-1" attr.type="boolean"/>
#<key id="d12" for="node" attr.name="module-svd-DM-2" attr.type="boolean"/>
#<key id="d20" for="node" attr.name="module-svd-DM-3" attr.type="boolean"/>
#<key id="d10" for="node" attr.name="module-svd-DM-7" attr.type="boolean"/>
#<key id="d11" for="node" attr.name="module-svd-DM-8" attr.type="boolean"/>
#<key id="d21" for="node" attr.name="module-svd-DM-9" attr.type="boolean"/>


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
layer="DM"
# Rotation
rotate=90

# Map
declare -A pairs=(
	[simple]=none
	[core]=d13
	[1]=d9
	[2]=d12
	[3]=d20
	#[4]=
	#[5]=
	#[6]=
	[7]=d10
	[8]=d11
	[9]=d21
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
	java -cp gephi-toolkit-0.9.2-all.jar:. Main --input $input --coords $coords --output $pdf --rotate $rotate --highlight $highlight
	
	# Convert PDF to JPG
	convert -density 300 -resize 3300x2550 $pdf $jpg
done
