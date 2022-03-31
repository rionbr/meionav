# coding=utf-8
# Author: Rion B Correia
# Date: Oct 09, 2020
#
# Description: Loads gephi files and request java to plot the networks
# COMPILE: javac -cp gephi-toolkit-0.9.2-all.jar PlotNetSingleMod.java
# COMPILE: javac -cp gephi-toolkit-0.9.2-all.jar PlotNetCategoricalMod.java
#
import os
import subprocess
import xml.etree.ElementTree as ET


def ensurePathExists(path):
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)


if __name__ == '__main__':

    celltypes = ['spermatocyte', 'enterocyte']
    layers = ['HS', 'MM', 'DM']
    #network = 'thr'
    network = 'backbone'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    data = {
        'spermatocyte': {
            'HS': {
                'rotate': 270
            },
            'MM': {
                'rotate': 90
            },
            'DM': {
                'rotate': 15
            },
        },
        'enterocyte': {
            'HS': {
                'rotate': 345
            },
            'MM': {
                'rotate': 215
            },
            'DM': {
                'rotate': 195
            },
        }
    }
    colors = {
        1: '#e377c2',
        2: '#2ca02c',
        3: '#ff7f0e',
        4: '#8c564b',
        5: '#9467bd',
        6: '#bcbd22',
        7: '#17becf',
        8: '#e7ba52',
        9: '#d6616b',
        10: '#ce6dbd',
        11: '#6b6ecf',
        12: '#8ca252',
    }
    #
    ns = {'ns': 'http://graphml.graphdrawing.org/xmlns'}
    #
    gp = '/home/casci/rionbr/spermnet/05-data-analysis/gephi-plotting/results/graphml'
    cp = '/home/casci/rionbr/spermnet/05-data-analysis/gephi-plotting/results/forceatlas2'
    op = '/home/casci/rionbr/spermnet/05-data-analysis/gephi-plotting/images'

    #
    for celltype in celltypes:
        print('Celltype: {celltype:s}'.format(celltype=celltype))

        for layer in layers:
            print('Layer: {layer:s}'.format(layer=layer))

            # Load GraphML
            print('parsing xml')
            rGraphmlFile = '../results/graphml/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            xml = ET.parse(rGraphmlFile)
            root = xml.getroot()

            # [(i, key.attrib.get('attr.name')) for i, key in enumerate(keys)]
            # Extract module identification
            keys = root.findall("ns:key[@for='node']", ns)
            for key in keys:
                
                if 'module-pca' in key.attrib.get('attr.name', ''):
                    continue

                    mid = int(key.attrib.get('attr.name').split('-')[-1])
                    gid = key.attrib.get('id')

                    print("Plotting network module: M{mid:d} ({gid:s})".format(mid=mid, gid=gid))

                    input = '{gp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(gp=gp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    coords = '{cp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-forceatlas2.txt'.format(cp=cp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    #
                    pdfp = '{op:s}/{celltype:s}/{layer:s}/pdf'.format(op=op, celltype=celltype, layer=layer)
                    jpgp = '{op:s}/{celltype:s}/{layer:s}/jpg'.format(op=op, celltype=celltype, layer=layer)
                    #
                    pdf = '{pdfp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-M{mid:d}.pdf'.format(pdfp=pdfp, celltype=celltype, network=network, threshold=threshold_str, layer=layer, mid=mid)
                    jpg = '{jpgp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-M{mid:d}.jpg'.format(jpgp=jpgp, celltype=celltype, network=network, threshold=threshold_str, layer=layer, mid=mid)
                    #
                    rotate = data[celltype][layer]['rotate']
                    node_highlight_color = colors[mid]
                    node_highlight_alt_color = '#c7c7c7'

                    ensurePathExists(pdf)
                    ensurePathExists(jpg)
                    # Trigger PDF subprocess
                    cmd = "java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input {input:s} --coords {coords:s} --output {pdf:s} --rotate {rotate:d} --node-size 20 --edge-thickness 0.4 --edge-opacity 10.0 --node-highlight {node_highlight:s} --node-highlight-color '{node_highlight_color:s}' --node-highlight-alt-color '{node_highlight_alt_color:s}'".format(input=input, coords=coords, pdf=pdf, rotate=rotate, node_highlight=gid, node_highlight_color=node_highlight_color, node_highlight_alt_color=node_highlight_alt_color)
                    subprocess.run(cmd, shell=True)

                    # Trigger JPG subprocess (150dpi = 1754x1240 ; 300dpi = 3508x2480)
                    cmd = "convert -density 150 -resize 1754x1240 -quality 80 {pdf:s} {jpg:s}".format(pdf=pdf, jpg=jpg)
                    subprocess.run(cmd, shell=True)

                # Plot conserved networks
                elif 'conserved' in key.attrib.get('attr.name'):
                    continue

                    gid = key.attrib.get('id')

                    print("Plotting conserved genes ({gid:s})".format(gid=gid))

                    input = '{gp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(gp=gp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    coords = '{cp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-forceatlas2.txt'.format(cp=cp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    #
                    pdfp = '{op:s}/{celltype:s}/{layer:s}/pdf'.format(op=op, celltype=celltype, layer=layer)
                    jpgp = '{op:s}/{celltype:s}/{layer:s}/jpg'.format(op=op, celltype=celltype, layer=layer)
                    #
                    pdf = '{pdfp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-conserved.pdf'.format(pdfp=pdfp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    jpg = '{jpgp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-conserved.jpg'.format(jpgp=jpgp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    #
                    rotate = data[celltype][layer]['rotate']
                    node_highlight_color = '#d62728'
                    node_highlight_alt_color = '#1f77b4'

                    ensurePathExists(pdf)
                    ensurePathExists(jpg)
                    # Trigger PDF subprocess
                    cmd = "java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input {input:s} --coords {coords:s} --output {pdf:s} --rotate {rotate:d} --node-size 20 --edge-thickness 0.4 --edge-opacity 10.0 --node-highlight {node_highlight:s} --node-highlight-color '{node_highlight_color:s}' --node-highlight-alt-color '{node_highlight_alt_color:s}'".format(input=input, coords=coords, pdf=pdf, rotate=rotate, node_highlight=gid, node_highlight_color=node_highlight_color, node_highlight_alt_color=node_highlight_alt_color)
                    subprocess.run(cmd, shell=True)

                    # Trigger JPG subprocess (150dpi = 1754x1240 ; 300dpi = 3508x2480)
                    cmd = "convert -density 150 -resize 1754x1240 -quality 80 {pdf:s} {jpg:s}".format(pdf=pdf, jpg=jpg)
                    subprocess.run(cmd, shell=True)
                
                # Plot core networks
                if 'core' in key.attrib.get('attr.name'):

                    gid = key.attrib.get('id')

                    print("Plotting core genes ({gid:s})".format(gid=gid))

                    input = '{gp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(gp=gp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    coords = '{cp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-forceatlas2.txt'.format(cp=cp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    #
                    pdfp = '{op:s}/{celltype:s}/{layer:s}/pdf'.format(op=op, celltype=celltype, layer=layer)
                    jpgp = '{op:s}/{celltype:s}/{layer:s}/jpg'.format(op=op, celltype=celltype, layer=layer)
                    #
                    pdf = '{pdfp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-core.pdf'.format(pdfp=pdfp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    jpg = '{jpgp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-core.jpg'.format(jpgp=jpgp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    #
                    rotate = data[celltype][layer]['rotate']
                    node_highlight_color = '#2ca02c'
                    node_highlight_alt_color = '#c7c7c7'

                    ensurePathExists(pdf)
                    ensurePathExists(jpg)
                    # Trigger PDF subprocess
                    cmd = "java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input {input:s} --coords {coords:s} --output {pdf:s} --rotate {rotate:d} --node-size 20 --edge-thickness 0.4 --edge-opacity 10.0 --node-highlight {node_highlight:s} --node-highlight-color '{node_highlight_color:s}' --node-highlight-alt-color '{node_highlight_alt_color:s}'".format(input=input, coords=coords, pdf=pdf, rotate=rotate, node_highlight=gid, node_highlight_color=node_highlight_color, node_highlight_alt_color=node_highlight_alt_color)
                    subprocess.run(cmd, shell=True)

                    # Trigger JPG subprocess (150dpi = 1754x1240 ; 300dpi = 3508x2480)
                    cmd = "convert -density 150 -resize 1754x1240 -quality 80 {pdf:s} {jpg:s}".format(pdf=pdf, jpg=jpg)
                    subprocess.run(cmd, shell=True)

                if 'mdlc_dge' in key.attrib.get('attr.name'):

                    gid = key.attrib.get('id')

                    print("Plotting mdlc-DGE genes ({gid:s})".format(gid=gid))

                    input = '{gp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(gp=gp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    coords = '{cp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-forceatlas2.txt'.format(cp=cp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    #
                    pdfp = '{op:s}/{celltype:s}/{layer:s}/pdf'.format(op=op, celltype=celltype, layer=layer)
                    jpgp = '{op:s}/{celltype:s}/{layer:s}/jpg'.format(op=op, celltype=celltype, layer=layer)
                    #
                    pdf = '{pdfp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-mdlc-dge.pdf'.format(pdfp=pdfp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    jpg = '{jpgp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-mdlc-dge.jpg'.format(jpgp=jpgp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    #
                    rotate = data[celltype][layer]['rotate']
                    #
                    node_highlight_v1_value = 'up'
                    node_highlight_v1_color = '#d62728'
                    #
                    node_highlight_v2_value = 'down'
                    node_highlight_v2_color = '#1f77b4'
                    #
                    node_highlight_alt_color = '#c7c7c7'

                    ensurePathExists(pdf)
                    ensurePathExists(jpg)
                    # Trigger PDF subprocess
                    cmd = "java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetCategoricalMod --input {input:s} --coords {coords:s} --output {pdf:s} --rotate {rotate:d} --node-size 20 --edge-thickness 0.4 --edge-opacity 10.0 --node-highlight {node_highlight:s} --node-highlight-v1-value '{node_highlight_v1_value:s}' --node-highlight-v1-color '{node_highlight_v1_color:s}' --node-highlight-v2-value '{node_highlight_v2_value:s}' --node-highlight-v2-color '{node_highlight_v2_color:s}' --node-highlight-alt-color '{node_highlight_alt_color:s}'".format(input=input, coords=coords, pdf=pdf, rotate=rotate, node_highlight=gid, node_highlight_v1_value=node_highlight_v1_value, node_highlight_v1_color=node_highlight_v1_color, node_highlight_v2_value=node_highlight_v2_value, node_highlight_v2_color=node_highlight_v2_color, node_highlight_alt_color=node_highlight_alt_color)
                    subprocess.run(cmd, shell=True)

                    # Trigger JPG subprocess (150dpi = 1754x1240 ; 300dpi = 3508x2480)
                    cmd = "convert -density 150 -resize 1754x1240 -quality 80 {pdf:s} {jpg:s}".format(pdf=pdf, jpg=jpg)
                    subprocess.run(cmd, shell=True)

                if 'mdlc_intron' in key.attrib.get('attr.name'):

                    gid = key.attrib.get('id')

                    print("Plotting mdlc-Introns genes ({gid:s})".format(gid=gid))

                    input = '{gp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(gp=gp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    coords = '{cp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-forceatlas2.txt'.format(cp=cp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    #
                    pdfp = '{op:s}/{celltype:s}/{layer:s}/pdf'.format(op=op, celltype=celltype, layer=layer)
                    jpgp = '{op:s}/{celltype:s}/{layer:s}/jpg'.format(op=op, celltype=celltype, layer=layer)
                    #
                    pdf = '{pdfp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-mdlc-intron.pdf'.format(pdfp=pdfp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    jpg = '{jpgp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-mdlc-intron.jpg'.format(jpgp=jpgp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
                    #
                    rotate = data[celltype][layer]['rotate']
                    node_highlight_color = '#d62728'
                    node_highlight_alt_color = '#c7c7c7'

                    ensurePathExists(pdf)
                    ensurePathExists(jpg)
                    # Trigger PDF subprocess
                    cmd = "java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input {input:s} --coords {coords:s} --output {pdf:s} --rotate {rotate:d} --node-size 20 --edge-thickness 0.4 --edge-opacity 10.0 --node-highlight {node_highlight:s} --node-highlight-color '{node_highlight_color:s}' --node-highlight-alt-color '{node_highlight_alt_color:s}'".format(input=input, coords=coords, pdf=pdf, rotate=rotate, node_highlight=gid, node_highlight_color=node_highlight_color, node_highlight_alt_color=node_highlight_alt_color)
                    subprocess.run(cmd, shell=True)

                    # Trigger JPG subprocess (150dpi = 1754x1240 ; 300dpi = 3508x2480)
                    cmd = "convert -density 150 -resize 1754x1240 -quality 80 {pdf:s} {jpg:s}".format(pdf=pdf, jpg=jpg)
                    subprocess.run(cmd, shell=True)

            # Plot an all gray network
            print("Plotting all-gray network")
            input = '{gp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(gp=gp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            coords = '{cp:s}/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-forceatlas2.txt'.format(cp=cp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            #
            pdfp = '{op:s}/{celltype:s}/{layer:s}/pdf'.format(op=op, celltype=celltype, layer=layer)
            jpgp = '{op:s}/{celltype:s}/{layer:s}/jpg'.format(op=op, celltype=celltype, layer=layer)
            #
            pdf = '{pdfp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-simple.pdf'.format(pdfp=pdfp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            jpg = '{jpgp:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-simple.jpg'.format(jpgp=jpgp, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            #
            gid = 'd2'
            rotate = data[celltype][layer]['rotate']
            node_highlight_color = '#c7c7c7'
            node_highlight_alt_color = '#c7c7c7'

            ensurePathExists(pdf)
            ensurePathExists(jpg)
            # Trigger PDF subprocess
            cmd = "java -cp gephi-toolkit-0.9.2-all.jar:. PlotNetSingleMod --input {input:s} --coords {coords:s} --output {pdf:s} --rotate {rotate:d} --node-size 20 --edge-thickness 0.4 --edge-opacity 10.0 --node-highlight {node_highlight:s} --node-highlight-color '{node_highlight_color:s}' --node-highlight-alt-color '{node_highlight_alt_color:s}'".format(input=input, coords=coords, pdf=pdf, rotate=rotate, node_highlight=gid, node_highlight_color=node_highlight_color, node_highlight_alt_color=node_highlight_alt_color)
            subprocess.run(cmd, shell=True)

            # Trigger JPG subprocess (150dpi = 1754x1240 ; 300dpi = 3508x2480)
            cmd = "convert -density 150 -resize 1754x1240 -quality 80 {pdf:s} {jpg:s}".format(pdf=pdf, jpg=jpg)
            subprocess.run(cmd, shell=True)

