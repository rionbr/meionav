// HTMLs
var progressbars_dm = document.getElementById('progressbars-dm');
var progressbar_file = document.getElementById('progressbar-file');
var progressbar_dm_stab = document.getElementById('progressbar-dm-stab');
//
var progressbars_mm = document.getElementById('progressbars-mm');
var progressbar_mm_stab = document.getElementById('progressbar-mm-stab');
//
var progressbars_hs = document.getElementById('progressbars-hs');
var progressbar_hs_stab = document.getElementById('progressbar-hs-stab');
//
var select_dm_find = document.getElementById('select-dm-find');
var select_mm_find = document.getElementById('select-mm-find');
var select_hs_find = document.getElementById('select-hs-find');
var select_dm_variable = document.getElementById('select-dm-variable');
var select_mm_variable = document.getElementById('select-mm-variable');
var select_hs_variable = document.getElementById('select-hs-variable');
//var select_query = document.getElementById('select-query');
var canvas_dm = document.getElementById('canvas-dm');
var canvas_mm = document.getElementById('canvas-mm');
var canvas_hs = document.getElementById('canvas-hs');
var canvas_cross = document.getElementById('canvas-cross');
var file = document.getElementById('canvas').getAttribute("data-file");
//var project = canvas.getAttribute("data-project");
//var medium = canvas.getAttribute("data-medium");
//var windowsize = canvas.getAttribute("data-windowsize");
//var dicttimestamp = canvas.getAttribute("data-dicttimestamp");
var status_span = document.getElementById('status-msg');
//
//var nodeinfo_name = document.getElementById('nodeinfo-name');
//var nodeinfo_desc = document.getElementById('nodeinfo-desc');
//var edgeinfo_source_name = document.getElementById('edgeinfo-source-name');
//var edgeinfo_source_type = document.getElementById('edgeinfo-source-type');
//var edgeinfo_target_name = document.getElementById('edgeinfo-target-name');
//var edgeinfo_target_type = document.getElementById('edgeinfo-target-type');
//var edgeinfo_weight = document.getElementById('edgeinfo-weight');

//var edgeinfo_ddi_adr_di = document.getElementById('edgeinfo-ddi-adr-di');
//var btn_edgeinfo = document.getElementById('edgeinfo-btn');
var selected_nodes = document.getElementById('selected-nodes');
//var btn_mentions = document.getElementById('mentions-btn');
//var btn_comentions = document.getElementById('comentions-btn');
var btn_select = document.getElementById('btn-select');
//var btn_query = document.getElementById('query-btn');
var btn_simulation = document.getElementById('btn-simulation');
var btn_reset = document.getElementById('btn-reset');
var btn_orphans = document.getElementById('btn-orphans');
//var checkbox_node_drug = document.getElementById('checkbox-node-drug');
//var checkbox_node_symp = document.getElementById('checkbox-node-symp');
//var checkbox_node_herb = document.getElementById('checkbox-node-herb');
//var checkbox_edge_drug_drug = document.getElementById('checkbox-edge-drug-drug');
//var checkbox_edge_symp_symp = document.getElementById('checkbox-edge-symp-symp');
//var checkbox_edge_herb_herb = document.getElementById('checkbox-edge-herb-herb');
//var checkbox_edge_drug_symp = document.getElementById('checkbox-edge-drug-symp');
//var checkbox_edge_drug_herb = document.getElementById('checkbox-edge-drug-herb');
//var checkbox_edge_herb_symp = document.getElementById('checkbox-edge-herb-symp');

// Variables
var json;
var network_hs;
var network_mm;
var network_dm;
var network_cross;

// DataSet
var nodes = new vis.DataSet();
var edges = new vis.DataSet();


// DataViews
var nodes_view_selected = new vis.DataView(nodes, {
	filter: function(item) {
		return (item.selected == true);
	}
});
var nodes_hs_view_display = new vis.DataView(nodes, {
	filter: function(item) {
		return ((item.display == true) && (item.layer == 'HS'));
	}
});
var nodes_mm_view_display = new vis.DataView(nodes, {
	filter: function(item) {
		return ((item.display == true) && (item.layer == 'MM'));
	}
});
var nodes_dm_view_display = new vis.DataView(nodes, {
	filter: function(item) {
		return ((item.display == true) && (item.layer == 'DM'));
	}
});
var nodes_cross_view_display = new vis.DataView(nodes, {
	filter: function(item) {
		return true;
	}
});

//
var edges_hs_view_display = new vis.DataView(edges, {
	filter: function(item) {
		return ((item.display == true) && (item.layer == 'HS'));
	}
});
var edges_mm_view_display = new vis.DataView(edges, {
	filter: function(item) {
		return ((item.display == true) && (item.layer == 'MM'));
	}
});
var edges_dm_view_display = new vis.DataView(edges, {
	filter: function(item) {
		return ((item.display == true) && (item.layer == 'DM'));
	}
});
//
var edges_cross_view_display = new vis.DataView(edges, {
	filter: function(item) {
		return (item.type == 'cross');
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

function selectNodesInOtherLayers(network, nodeids) {
	nodeids.forEach(function(nodeid) {
		var neighborids = network_cross.getConnectedNodes(nodeid);
		network_cross.selectNodes(neighborids);
		//
		if (network != 'HS') {
			var tmp_hs = new vis.DataView(nodes_hs_view_display, {
				filter: function(item) {
					return (neighborids.includes(item.id));
				}
			});
			network_hs.selectNodes(tmp_hs.getIds());
		}
		//
		if (network != 'MM') {
			var tmp_mm = new vis.DataView(nodes_mm_view_display, {
				filter: function(item) {
					return (neighborids.includes(item.id));
				}
			});
			network_mm.selectNodes(tmp_mm.getIds());
		}
		//
		if (network != 'DM') {
			var tmp_dm = new vis.DataView(nodes_dm_view_display, {
				filter: function(item) {
					return (neighborids.includes(item.id));
				}
			});
			network_dm.selectNodes(tmp_dm.getIds());
		}
	});
}

function deselectNodesInOtherLayers(network) {
	network_hs.selectNodes([]);
	network_mm.selectNodes([]);
	network_dm.selectNodes([]);
	network_cross.selectNodes([]);
}

function returnNetworkObjectsFromString(s) {
	if (s == 'HS') {
		return {
			'network': network_hs,
			'nodes': nodes_hs_view_display,
			'edges': edges_hs_view_display,
			'select-find': select_hs_find,
			'select-variable': select_hs_variable
		};
	} else if (s == 'MM') {
		return {
			'network': network_mm,
			'nodes': nodes_mm_view_display,
			'edges': edges_mm_view_display,
			'select-find': select_mm_find,
			'select-variable': select_mm_variable
		};
	} else if (s == 'DM') {
		return {
			'network': network_dm,
			'nodes': nodes_dm_view_display,
			'edges': edges_dm_view_display,
			'select-find': select_dm_find,
			'select-variable': select_dm_variable
		};
	} else {
		throw new TypeError('String must be HS, MM or DM');
	}
}

function colorNetworkFromVariable(network_name) {
	objs = returnNetworkObjectsFromString(network_name);
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

var cmap_Reds = chroma.scale(chroma.brewer.Reds).nodata('white');

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


/*
// DataView Triggers
nodes_view_selected.on('add', function (event, params) {
	node = nodes.get(params.items[0]);
	// Create <li id="#"><span class="label label-default"></span></li>
	var li = document.createElement("li");
	li.id = 'node_selected_'+node.id
	li.className = "list-inline-item";
	var span = document.createElement("span");
	span.innerHTML = node.title;
	span.className = "badge badge-secondary";
	selected_nodes.appendChild(li);
	li.appendChild(span);

	// Enable Select button
	if (selected_nodes.children.length >= 1) {
		btn_select.className = "btn btn-sm btn-primary";
	}
});
nodes_view_selected.on('remove', function (event,params) {
	node = nodes.get(params.items[0]);
	a = document.getElementById('node_selected_'+node.id);
	a.parentElement.removeChild(a);
	
	// Disable Select button
	if (selected_nodes.children.length == 0){
		btn_select.className = "btn btn-sm btn-primary disabled";
	}
});
*/


// Network data
var network_dm_data;
var network_mm_data;
var network_hs_data;

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
			background: 'LightGray',
			border: 'Gray',
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
			color:'gray',
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
var network_cross_options = {
	autoResize:false,
	nodes: {
		fixed:true,
		shape:'dot',
		size:10,
		font: {
			size:12,
			color:'White',
		},
		borderWidth:0,
		color: {
			background:'LightGray',
			border:'Gray',
			highlight:'#00ff00',
		},
	},	
	edges: {
		labelHighlightBold:false,
		smooth: {
			enabled:false,
			type:'discrete',
			roundness:0.0,
		},
		color: {
			color:'gray',
			highlight:'#00ff00',
		},
		scaling: {
			min:1,
			max:30,
		},
	},
	interaction: {
		dragNodes:false,
		dragView:false,
		hover:false,
		multiselect:false,
		navigationButtons:false,
		multiselect:false,
		selectable:false,
		tooltipDelay:0,
		zoomView:false

	},
	layout: {
		improvedLayout:false,
	},
	physics: {
		//maxVelocity: 146,
		solver: 'forceAtlas2Based',
		//timestep: 0.35,
		stabilization: {
			enabled: false,
			iterations: 1,
			updateInterval:20
		},
	},
};

/*
function nodeDisplayGroup(btn, group) {
	var updates = [];
	display = (btn.checked==true) ? true : false;
	nodes.forEach(function(n) {
		if (n.group == group) {
			updates.push({id:n.id, display:display});
		}
	});
	nodes.update(updates);
}
*/
/*
function nodeDisplaySelected() {
	var updates = [];
	ids = nodes_view_selected.getIds();
	nodes.forEach(function(n) {
		if (!ids.includes(n.id)) {
			updates.push({id:n.id, display:false});
		}
	});
	nodes.update(updates);
}
*/
/*
function nodeDisplayOrphans(btn) {
	var updates = [];
	nodes.forEach(function(node) {
		//console.log("Node:",node);
		sid = node.id
		eids = network.getConnectedEdges(sid);
		//console.log("Node:",node, "Edges ids:",eids);
		orphan = true;
		for (var i=0; i<eids.length; i++) {
			eid = eids[i];
			edge = edges.get(eid);
			//console.log("edge:",eid, edge);
			if (edge.display) {
				orphan = false;
				break
			} else {
				tid = (edge.from == sid) ? edge.from : edge.to
				if (nodes.get(tid).display) {
					orphan = false;
					break
				}
			}
		}
		if (orphan) {
			updates.push({id:sid, display:false});
		}
	});
	nodes.update(updates);
}
*/
/*
function nodeEdgeDisplayReset() {
	var updates = [];
	nodes.forEach(function(n) {
		updates.push({id:n.id, display:true, label:undefined});
	},{ filter: function(n) {
		return ((n.display == false) || (n.label != undefined));
	}});
	nodes.update(updates);
	var updates = [];
	edges.forEach(function(e) {
		updates.push({id:e.id, display:true});
	}, {filter: function(e) {
		return (e.display == false);
	}});
	edges.update(updates);
}
*/
/*
function edgeDisplayGroup(btn, source, target) {
	var updates = []
	display = (btn.checked==true) ? true : false;
	edges.forEach(function (e) { 
		sid = e.from;
		tid = e.to;
		stype = nodes.get(sid).group;
		ttype = nodes.get(tid).group;
		if ( ((stype == source) && (ttype == target)) || ((stype == target) && (ttype == source)) ) {
			updates.push({id:e.id, display:display})
		}
	});
	edges.update(updates);
}
*/

function startStopSimulation(button) {
	isOn = network_dm.physics.options.enabled;
	if (isOn) {
		network_dm.setOptions( { physics: false } );
		network_mm.setOptions( { physics: false } );
		network_hs.setOptions( { physics: false } );
		button.innerHTML = 'Run';
		button.className = 'btn btn-sm btn-success';
	} else {
		network_dm.setOptions( { physics: true } );
		network_mm.setOptions( { physics: true } );
		network_hs.setOptions( { physics: true } );
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
			if (n['layer'] == 'HS') { 
				select_hs_find.appendChild(opt); 
			} else if (n['layer'] == 'MM') {
				select_mm_find.appendChild(opt);
			} else if (n['layer'] == 'DM') {
				select_dm_find.appendChild(opt);
			} else {
				throw new TypeError('node layer must be HS, MM or DM');
			}
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

		network_hs_data = {nodes: nodes_hs_view_display, edges: edges_hs_view_display};
		network_mm_data = {nodes: nodes_mm_view_display, edges: edges_mm_view_display};
		network_dm_data = {nodes: nodes_dm_view_display, edges: edges_dm_view_display};

		network_cross_data = {nodes: nodes_cross_view_display, edges: edges_cross_view_display};
		
		network_hs = new vis.Network(canvas_hs, network_hs_data, network_options);
		network_mm = new vis.Network(canvas_mm, network_mm_data, network_options);
		network_dm = new vis.Network(canvas_dm, network_dm_data, network_options);
		network_cross = new vis.Network(canvas_cross, network_cross_data, network_cross_options);

		// HS
		network_hs.on("stabilizationProgress", function(p) {
			progressbar_hs_stab.style.width = (p.iterations/p.total*100)+'%';
		});
		network_hs.on("stabilizationIterationsDone", function () {
			network_hs.setOptions( { physics: false } );
			progressbar_hs_stab.style.width = '100%';
			setTimeout(function () { progressbars_hs.style.display = 'none'; }, 500);			
			//network_hs.storePositions();
		});
		network_hs.on("selectNode", function(params) {
			node = nodes.get(params.nodes[0]);
			selectNodesInOtherLayers('HS', params.nodes);
		});
		network_hs.on("deselectNode", function(params) {
			deselectNodesInOtherLayers('HS');
		});
		
		// MM
		network_mm.on("stabilizationProgress", function(p) {
			progressbar_dm_stab.style.width = (p.iterations/p.total*100)+'%';
		});
		network_mm.on("stabilizationIterationsDone", function () {
			network_mm.setOptions( { physics: false } );
			progressbar_mm_stab.style.width = '100%';
			setTimeout(function () { progressbars_mm.style.display = 'none'; }, 500);			
			//network_mm.storePositions();
		});
		network_mm.on("selectNode", function(params) {
			node = nodes.get(params.nodes[0]);
			selectNodesInOtherLayers('MM', params.nodes);
		});
		network_mm.on("deselectNode", function(params) {
			deselectNodesInOtherLayers('MM');
		});
		
		// DM
		network_dm.on("stabilizationProgress", function(p) {
			btn_simulation.innerHTML = 'Stop';
			btn_simulation.className = 'btn btn-danger';
			progressbar_dm_stab.style.width = (p.iterations/p.total*100)+'%';
		});
		network_dm.on("stabilizationIterationsDone", function () {
			network_dm.setOptions( { physics: false } );
			progressbar_dm_stab.style.width = '100%';
			setTimeout(function () { progressbars_dm.style.display = 'none'; }, 500);			
			btn_simulation.innerHTML = 'Run';
			btn_simulation.className = 'btn btn-success';
			//network_dm.storePositions();
		});
		network_dm.on("selectNode", function(params) {
			node = nodes.get(params.nodes[0]);
			selectNodesInOtherLayers('DM', params.nodes);
		});
		network_dm.on("deselectNode", function(params) {
			deselectNodesInOtherLayers('DS');
		});
	});
	
	// Tooltips
	$(document).ready(function(){
		$('[data-toggle="tooltip"]').tooltip(); 
	});
	
});