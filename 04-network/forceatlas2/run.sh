# Commands to run ForceAtlas2 via command-line

# java -Djava.awt.headless=true -Xmx8g -cp forceatlas2.jar:gephi-toolkit-0.9.2-all.jar kco.forceatlas2.Main --input $gp/net_threshold-0p5-HS.graphml --output $gp/net_threshold-0p5-HS-forceatlas2.txt --nsteps 2000 --2d --nthreads 12 --format txt --seed 123 --gravity 1.0 --scalingRatio 1.0 --outboundAttractionDistribution true

# Input Path
ip="/home/casci/rionbr/spermnet/4-network/results/graphml"
# Output Path
op="/home/casci/rionbr/spermnet/4-network/results/forceatlas2"

# Network
network="thr"
# Threshold
threshold="0p5"

## Loop celltypes
for celltype in "enterocyte" "spermatocyte"
do
	## Loop Layers
	for layer in HS MM DM
	do
		echo "Computing ForceAtlas2 on $celltype - $layer"

		# input/output file
		input="$ip/$celltype/net-$celltype-$network-$threshold-$layer.graphml"
		output="$op/$celltype/net-$celltype-$network-$threshold-$layer-forceatlas2.txt"
		

		# Run ForceAtlas2
		java -Djava.awt.headless=true -Xmx8g -cp forceatlas2.jar:gephi-toolkit-0.9.2-all.jar kco.forceatlas2.Main --input $input --output $output --nsteps 2000 --2d --nthreads 12 --format txt --seed 123 --gravity 1.0 --scalingRatio 1.0 --outboundAttractionDistribution true
	done	
done
