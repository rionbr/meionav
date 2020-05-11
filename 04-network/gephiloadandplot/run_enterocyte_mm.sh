#<key id="d7" for="node" attr.name="core" attr.type="boolean"/>
#<key id="d18" for="node" attr.name="module-svd-MM-1" attr.type="boolean"/>
#<key id="d14" for="node" attr.name="module-svd-MM-2" attr.type="boolean"/>
#<key id="d15" for="node" attr.name="module-svd-MM-3" attr.type="boolean"/>
#<key id="d11" for="node" attr.name="module-svd-MM-4" attr.type="boolean"/>
#<key id="d13" for="node" attr.name="module-svd-MM-5" attr.type="boolean"/>
#<key id="d16" for="node" attr.name="module-svd-MM-7" attr.type="boolean"/>
#<key id="d12" for="node" attr.name="module-svd-MM-8" attr.type="boolean"/>
#<key id="d17" for="node" attr.name="module-svd-MM-9" attr.type="boolean"/>


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
layer="MM"
# Rotation
rotate=270

# Map
declare -A pairs=(
	[simple]=none
	[core]=d7
	[1]=d18
	[2]=d14
	[3]=d15
	[4]=d11
	[5]=d14
	#[6]=
	[7]=d16
	[8]=d12
	[9]=d17
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
