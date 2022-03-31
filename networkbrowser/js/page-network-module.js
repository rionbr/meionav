// HTMLs
var progressbars = document.getElementById('progressbars');
var progressbar_file = document.getElementById('progressbar-file');
var progressbar_stab = document.getElementById('progressbar-stab');
//
var select_find = document.getElementById('select-find');
var select_variable = document.getElementById('select-variable');
//var select_query = document.getElementById('select-query');
var canvas = document.getElementById('canvas');
var file = document.getElementById('canvas').getAttribute("data-file")
var layer = canvas.getAttribute("data-layer");
var algorithm = canvas.getAttribute("data-algorithm");
var module = canvas.getAttribute("data-module");
var status_span = document.getElementById('status-msg');
//
var selected_nodes = document.getElementById('selected-nodes');
var btn_select = document.getElementById('btn-select');
var btn_simulation = document.getElementById('btn-simulation');
var btn_reset = document.getElementById('btn-reset');
var btn_orphans = document.getElementById('btn-orphans');

// Variables
var json;
var network;

// DataSet
var nodes = new vis.DataSet();
var edges = new vis.DataSet();


// DataViews
var nodes_view_selected = new vis.DataView(nodes, {
	filter: function(item) {
		return (item.selected == true);
	}
});
var nodes_view_display = new vis.DataView(nodes, {
	filter: function(item) {
		return (item.display == true);
	}
});

//
var edges_view_display = new vis.DataView(edges, {
	filter: function(item) {
		return (item.display == true);
	}
});


/*
var edges.get({
	filter: function(item) {
		return (
				(item.type == 'cross') && ((item.from == 'A-DM') || (item.to == 'A-DM'))
			);
		}
	});
*/


function returnNetworkObjectsFromString() {
	return {
		'network': network,
		'nodes': nodes_view_display,
		'edges': edges_view_display,
		'select-find': select_find,
		'select-variable': select_variable
	};
}

function colorNetworkFromVariable() {
	objs = returnNetworkObjectsFromString();
	network = objs['network'];
	view = objs['nodes'];
	select = objs['select-variable'];
	var_name = select.options[select.selectedIndex].value;
	var_type = select.options[select.selectedIndex].getAttribute('data-type')
	var_cmap = select.options[select.selectedIndex].getAttribute('data-cmap')
	//
	colorNetworkCmap(network, view, var_name, var_type, eval(var_cmap));
}


function colorNetworkCmap(network, view, var_name, var_type, var_cmap) {
	updates = [];
	view.get({'fields':['id', var_name]}).forEach(function(node) {
		variable = node[var_name];
		if (var_type == 'discrete') {
			color = var_cmap[variable] ? var_cmap[variable] : var_cmap['other']
		} else if (var_type == 'continuous') {
			color = var_cmap(variable).hex();
		} else {
			throw new TypeError('Variable name/type/cmap invalid');
		}
		updates.push({'id':node['id'], 'color':{'background':color}});
	});
	nodes.update(updates)
}

var cmap_Biotype = {
	'protein_coding': '#d62728',
	'other':'#ffffff'
}
var cmap_RNAi = {
	'Yes': '#2ca02c',
	'No': '#d62728', //'darkgray',
	'other':'#ffffff'
}

var cmap_Blues_r = chroma.scale(chroma.brewer.Blues).nodata('white').domain([1,0]);

var cmap_Reds = chroma.scale(chroma.brewer.Reds).nodata('white').domain([0,10]);

var cmap_Boolean = {
	true: '#2ca02c',
	false: '#d62728', //'darkgray',
	'other': '#ffffff', 
}

var cmap_pheno_code = {
	'A': '#d62728',
	'B': '#ce6dbd',
	'C': '#756bb1',
	'D': '#c7e9c0',
	'E': '#9edae5',
	'F': '#fdd0a2',
	'G': '#dadaeb',
	'H': '#bdbdbd',
	'other': '#ffffff'
}

function cmap_modules_func(n_colors) {
	return chroma.scale(chroma.brewer.Paired).nodata('white').colors(n_colors);
}


// Network data
var network_data;

// Network Options
var network_options = {
	nodes: {
		fixed: false,
		shape: 'dot',
		size: 10,
		font: {
			color: 'black',
		},
		borderWidth: 1,
		color: {
			background: '#d9d9d9',
			border: '#636363',
			highlight: '#00ff00',
		},
	},	
	edges: {
		labelHighlightBold: false,
		smooth: {
			enabled: false,
			type: 'discrete',
			roundness: 0.5,
		},
		color: {
			color:'#969696',
			highlight:'#00ff00',
		},
		scaling: {
			min: 1,
			max: 30,
		},
	},
	interaction: {
		dragNodes: true,
		dragView: true,
		hideEdgesOnZoom: true,
		navigationButtons: false,
		tooltipDelay: 200,
		multiselect: false,
	},
	layout: {
		improvedLayout: false,
	},
	physics: {
		maxVelocity: 150,
		minVelocity: 1,
		solver: 'forceAtlas2Based',
		timestep: 0.9,
		stabilization: {
			enabled: true,
			iterations: 100,
			updateInterval: 50
		},
	},
};


function startStopSimulation(button) {
	isOn = network.physics.options.enabled;
	if (isOn) {
		network.setOptions( { physics: false } );
		button.innerHTML = 'Run';
		button.className = 'btn btn-sm btn-success';
	} else {
		network.setOptions( { physics: true } );
		button.innerHTML = 'Stop';
		button.className = 'btn btn-sm btn-danger';
	}
}


function zoomInNode(network_name) {
	objs = returnNetworkObjectsFromString(network_name);
	network = objs['network'];
	select = objs['select-find'];
	
	var options = {
			scale:1,
			animation: {
				duration:650,
				easingFunction:'easeInOutCubic',
			}
		}
	nid = select.options[select.selectedIndex].value;
	network.focus(nid, options);
}

function displayStatusMessage(msg) {
	status_span.innerHTML = msg;
}


$(document).ready(function() {
	console.log('Loading network JSON file.');
	displayStatusMessage('Loading network JSON file.'); 
	
	var loader = $.ajax({
		type: "GET",
		url: file,
		dataType: "text/json",
		xhr: function() {
			var xhr = $.ajaxSettings.xhr();
			xhr.onprogress = function(e) {
				if (e.lengthComputable) {
					progressbar_file.style.width = (e.loaded/e.total*100)+'%';
				} else {
					console.log('Length not computable. Setting to 10% and waiting completion.');
					progressbar_file.style.width = '10%';
				}
			};
			xhr.onload = function(e) {
				progressbarfile.style.width = '100%';
			};
			return xhr;
		}
	})
	loader.always(function(e) {

		try {
			json = JSON.parse(e.responseText);
		} catch (err) {
			displayStatusMessage('ERROR: bad format JSON.');
			throw new TypeError('bad JSON format');
		}
		console.log('Parsing nodes and edges.');
		displayStatusMessage('Parsing nodes and edges.');
		
		// Nodes
		list_nodes = [];
		$.each(json.nodes, function(k,node) {
			layer = node['layer'];
			node['display'] = true;
			node['title'] = node['label']
			node['size'] = 25
			//delete node['label'];
			//node['color'] = '#ff9896';//node['color-fert-rate'];
			//delete node['color'];
			list_nodes.push(node);
		});
		nodes.add(list_nodes);
		
		// Populate Selects
		nodes.get({'order':'title'}).forEach( function(n) {
			var opt = document.createElement("option");
			opt.value = n['id'];
			opt.innerHTML = n['title'];
			select_find.appendChild(opt);
			//
		});
		
		// Edges
		list_edges = [];
		$.each(json.edges, function(k,edge) {
			edge['display'] = true;
			list_edges.push(edge);
		});
		edges.add(list_edges);
		
		if ((nodes.length<=0) || (edges.length<=0)) {
			displayStatusMessage('ERROR: no node/edge to display.');
			throw new Error('no node/edge to display');
		}
		
		console.log('Drawing Network.');
		displayStatusMessage('Drawing network.');
		// create a network

		network_data = {nodes: nodes_view_display, edges: edges_view_display};
		network = new vis.Network(canvas, network_data, network_options);
		
		network.on("stabilizationProgress", function(p) {
			btn_simulation.innerHTML = 'Stop';
			btn_simulation.className = 'btn btn-danger';
			progressbar_stab.style.width = (p.iterations/p.total*100)+'%';
		});
		network.on("stabilizationIterationsDone", function () {
			network.setOptions( { physics: false } );
			progressbar_stab.style.width = '100%';
			setTimeout(function () { progressbars.style.display = 'none'; }, 500);
			btn_simulation.innerHTML = 'Run';
			btn_simulation.className = 'btn btn-success';
		});
		network.on("selectNode", function(params) {
			node = nodes.get(params.nodes[0]);
		});
		network.on("deselectNode", function(params) {
		});
	});
	
	// Tooltips
	$(document).ready(function(){
		$('[data-toggle="tooltip"]').tooltip(); 
	});
	
});