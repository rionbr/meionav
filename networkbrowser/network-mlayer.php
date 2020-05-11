<?php
include_once("config.php");
$curpage = "Network";
$pagetitle = "MultiLayer Network : SpermNet";



?>
<!doctype html>
<html lang="en">
<head>
	<meta charset="utf-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge">
	<meta name="viewport" content="width=device-width, initial-scale=1">
	
	<!-- Bootstrap -->
	<link href="css/bootstrap.min.css" rel="stylesheet" type="text/css" />
	
	<!-- Vis-network.js -->
    <link href="css/vis-network.min.css" rel="stylesheet" type="text/css" />

	<!-- Oi Icons -->
	<link href="css/open-iconic-bootstrap.min.css" rel="stylesheet">

    <link href="css/spermnet.css" rel="stylesheet" type="text/css" />
    <style type="text/css">

    </style>
    
	<title><?php print($pagetitle);?></title>
</head>
<body>
	<header>
	<?php include "include/navigation.php"; ?>


	<main role="main">
		<div class="container">

			<nav aria-label="breadcrumb">
				<ol class="breadcrumb">
					<li class="breadcrumb-item"><a href="">Home</a></li>
					<li  class="breadcrumb-item active">MultiLayer Network</li>
				</ol>
			</nav>

			<div class="py-5 mb-4">
				<div class="row">
					<div class="col-md-9 col-ld-7">
						<h1>MultiLayer Network</h1>
						<p>
							The networks below represent individual layers of the multilayer network.
							To identify cross-layer edges, click to highlight a node in one of the layers. Their homolog genes will highlight in the other layers.
							Due to visualization restrictions, only intra-layer edges with experimental evidence are being shown.							
						</p>
					</div>
				</div>
			</div>
		</div>


		<!--
			Networks
		-->
		<div id="canvas" class="container" data-file="http://localhost/~rionbr/spermnet/networkbrowser/json/net_complete_core_mlayer_backbone.json">
		
			<!-- Left Canvas -->	
			<div class="row">
				<div class="col-9 p-2">

					<!-- DM Network Canvas -->
					<div class="position-relative mb-2">
						<h5 class="position-absolute m-2"><em>Drosophila melanogaster</em></h5>
						<div id="controls-dm" class="position-absolute fixed-bottom border rounded form-inline p-2 m-2 bg-light">
							<div class="input-group">
								<div class="input-group-prepend">
									<label for="select-dm-find" class="input-group-text">
										<span class="oi oi-magnifying-glass" aria-hidden="true" data-toggle="tooltip" title="Locates the node on the network."></span>
									</label>
								</div>
								<select id="select-dm-find" class="custom-select">

								</select>
								<div class="input-group-append">
									<button class="btn btn-primary" type="button" onClick="zoomInNode('DM');">Find</button>
								</div>
							</div>
							<div class="input-group mx-2">
								<div class="input-group-prepend">
									<label for="select-dm-variable" class="input-group-text">
										<span class="oi oi-beaker" aria-hidden="true" data-toggle="tooltip" title="Color the network according to specific variable."></span>
									</label>
								</div>
								<select id="select-dm-variable" class="custom-select">
									<option value="mean-fert-rate" data-type="continuous" data-cmap="cmap_Blues_r">DM Fertility Rate</option>
									<option value="new-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">New DM phenotype</option>
									<option value="known-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known DM phenotype</option>
									<option value="known-MM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known MM phenotype (homologs)</option>
									<option value="known-HS-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known HS phenotype (homologs)</option>
									<option value="logFPKM" data-type="continuous" data-cmap="cmap_Reds">log2(FPKM)</option>
									<option value="RNAi" data-type="discrete" data-cmap="cmap_RNAi">Validated RNAi</option>
									<option value="meiotic-entry" data-type="discrete" data-cmap="cmap_Boolean">Meiotic entry</option>
									<option value="meiotic-exit" data-type="discrete" data-cmap="cmap_Boolean">Meiotic exit</option>
									<option value="biotype" data-type="discrete" data-cmap="cmap_Biotype">Biotype</option>
									<option value="modules-infomap" data-type="discrete" data-cmap="cmap_modules_func(2)">Modules: Infomap</option>
									<option value="modules-DM-infomap" data-type="discrete" data-cmap="cmap_modules_func(28)">Modules: DM Infomap</option>
									<option value="modules-DM-louvain" data-type="discrete" data-cmap="cmap_modules_func(28)">Modules: DM Louvain</option>
								</select>
								<div class="input-group-append">
									<button class="btn btn-primary" type="button" onClick="colorNetworkFromVariable('DM');">View</button>
								</div>
							</div>
							<div class="input-group mx-2">
								<div class="input-group-prepend">
									<label for="input-dm-layout" class="input-group-text">
										<span class="oi oi-sun" aria-hidden="true" data-toggle="tooltip" title="Organize the network layout."></span>
									</label>
								</div>
								<div class="input-group-append">
									<button class="btn btn-light" type="button" id="btn-simulation" onclick='startStopSimulation(this)'>Run</button> 
								</div>
							</div>
						</div>
						<div id="canvas-dm" class="border rounded" style="height:700px"></div>
						<div id="progressbars-dm" class="position-absolute mx-auto my-auto fixed-top" style="width:600px; top:35%">
							<!-- Error Msg -->
							<div class="text-center mb-3">
								<p><span class="oi oi-dashboard"></span> <span id="status-msg">Loading scripts.</span></p>
							</div>
							Network force layout stabilization
							<div id="" class="progress">
								<div id="progressbar-dm-stab" class="progress-bar bg-success progress-bar-striped active" role="progressbar" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100">
									<span class="sr-only">0% Complete (success)</span>
								</div>
							</div>
							Loading network from json
							<div id="progress-file" class="progress">
								<div id="progressbar-file" class="progress-bar bg-warning progress-bar-striped active" role="progressbar" aria-valuenow="1" aria-valuemin="0" aria-valuemax="100" style="width:1%;">
									<span class="sr-only">0% Complete (success)</span>
								</div>
							</div>
						</div>
					</div>

					<!-- MM Network Canvas -->
					<div class="position-relative my-2">
						<h5 class="position-absolute m-2"><em>Mus musculus</em></h5>
						<div id="controls-dm" class="position-absolute fixed-bottom border rounded form-inline p-2 m-2 bg-light">
							<div class="input-group">
								<div class="input-group-prepend">
									<label for="select-mm-find" class="input-group-text">
										<span class="oi oi-magnifying-glass" aria-hidden="true" data-toggle="tooltip" title="Locates the node on the network."></span>
									</label>
								</div>
								<select id="select-mm-find" class="custom-select">

								</select>
								<div class="input-group-append">
									<button class="btn btn-primary" type="button" onClick="zoomInNode('MM');">Find</button>
								</div>
							</div>
							<div class="input-group mx-2">
								<div class="input-group-prepend">
									<label for="select-mm-variable" class="input-group-text">
										<span class="oi oi-beaker" aria-hidden="true" data-toggle="tooltip" title="Color the network according to specific variable."></span>
									</label>
								</div>
								<select id="select-mm-variable" class="custom-select">
									<option value="mean-fert-rate" data-type="continuous" data-cmap="cmap_Blues_r">DM Fertility Rate (homologs)</option>
									<option value="new-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">New DM phenotype (homologs)</option>
									<option value="known-MM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known MM phenotype</option>
									<option value="known-HS-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known HS phenotype (homologs)</option>
									<option value="known-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known DM phenotype (homologs)</option>
									<option value="logFPKM" data-type="continuous" data-cmap="cmap_Reds">log2(FPKM)</option>
									<option value="meiotic-entry" data-type="discrete" data-cmap="cmap_Boolean">Meiotic entry</option>
									<option value="meiotic-exit" data-type="discrete" data-cmap="cmap_Boolean">Meiotic exit</option>
									<option value="biotype" data-type="discrete" data-cmap="cmap_Biotype">Biotype</option>
									<option value="modules-infomap" data-type="discrete" data-cmap="cmap_modules_func(2)">Modules: Infomap</option>
									<option value="modules-MM-mammals-infomap" data-type="discrete" data-cmap="cmap_modules_func(36)">Modules: MM Mammals Infomap</option>
									<option value="modules-MM-mammals-louvain" data-type="discrete" data-cmap="cmap_modules_func(17)">Modules: MM Mammals Louvain</option>
								</select>
								<div class="input-group-append">
									<button class="btn btn-primary" type="button" onClick="colorNetworkFromVariable('MM');">View</button>
								</div>
							</div>
						</div>
						<div id="canvas-mm" class="canvas border rounded" style="height:700px"></div>
						<div id="progressbars-mm" class="position-absolute mx-auto my-auto fixed-top" style="width:600px; top:45%">
							Network force layout stabilization
							<div id="" class="progress">
								<div id="progressbar-mm-stab" class="progress-bar bg-success progress-bar-striped active" role="progressbar" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100">
									<span class="sr-only">0% Complete (success)</span>
								</div>
							</div>
						</div>
					</div>

					<!-- HS Network Canvas -->
					<div class="position-relative my-2">
						<h5 class="position-absolute m-2"><em>Homo sapiens</em></h5>
						<div id="controls-dm" class="position-absolute fixed-bottom border rounded form-inline p-2 m-2 bg-light">
							<div class="input-group">
								<div class="input-group-prepend">
									<label for="select-hs-find" class="input-group-text">
										<span class="oi oi-magnifying-glass" aria-hidden="true" data-toggle="tooltip" title="Locates the node on the network."></span>
									</label>
								</div>
								<select id="select-hs-find" class="custom-select">

								</select>
								<div class="input-group-append">
									<button class="btn btn-primary" type="button" onClick="zoomInNode('HS');">Find</button>
								</div>
							</div>
							<div class="input-group mx-2">
								<div class="input-group-prepend">
									<label for="select-hs-variable" class="input-group-text">
										<span class="oi oi-beaker" aria-hidden="true" data-toggle="tooltip" title="Color the network according to specific variable."></span>
									</label>
								</div>
								<select id="select-hs-variable" class="custom-select">
									<option value="mean-fert-rate" data-type="continuous" data-cmap="cmap_Blues_r">DM Fertility Rate (homologs)</option>
									<option value="new-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">New DM phenotype (homologs)</option>
									<option value="known-HS-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known HS phenotype</option>
									<option value="known-MM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known MM phenotype (homologs)</option>
									<option value="known-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known DM phenotype (homologs)</option>
									<option value="logFPKM" data-type="continuous" data-cmap="cmap_Reds">log2(FPKM)</option>
									<option value="meiotic-entry" data-type="discrete" data-cmap="cmap_Boolean">Meiotic entry</option>
									<option value="meiotic-exit" data-type="discrete" data-cmap="cmap_Boolean">Meiotic exit</option>
									<option value="biotype" data-type="discrete" data-cmap="cmap_Biotype">Biotype</option>
									<option value="modules-infomap" data-type="discrete" data-cmap="cmap_modules_func(2)">Modules: Infomap</option>
									<option value="modules-HS-mammals-infomap" data-type="discrete" data-cmap="cmap_modules_func(35)">Modules: HS Mammals Infomap</option>
									<option value="modules-HS-mammals-louvain" data-type="discrete" data-cmap="cmap_modules_func(14)">Modules: HS Mammals Louvain</option>
								</select>
								<div class="input-group-append">
									<button class="btn btn-primary" type="button" onClick="colorNetworkFromVariable('HS');">View</button>
								</div>
							</div>
						</div>
						<div id="canvas-hs" class="canvas border rounded" style="height:700px"></div>
						<div id="progressbars-hs" class="position-absolute mx-auto my-auto fixed-top" style="width:600px; top:45%">
							Network force layout stabilization
							<div id="" class="progress">
								<div id="progressbar-hs-stab" class="progress-bar bg-success progress-bar-striped active" role="progressbar" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100">
									<span class="sr-only">0% Complete (success)</span>
								</div>
							</div>
						</div>
					</div>
				</div>
	
				<!-- Right Canvas -->
				<div class="col-3 p-2">
					<!--
						Legend
					-->
					<div id="legend" class="">
						<div id="" class="border rounded p-2 mb-2">
							<div class="m-2">
								<h6>Fertility Rate</h6>
								<div class="gradient">
									<span class="grad-step" style="background-color:#08306b"></span><span class="grad-step" style="background-color:#08306b"></span><span class="grad-step" style="background-color:#08306b"></span><span class="grad-step" style="background-color:#08306b"></span><span class="grad-step" style="background-color:#08306b"></span><span class="grad-step" style="background-color:#083b7b"></span><span class="grad-step" style="background-color:#083b7b"></span><span class="grad-step" style="background-color:#083b7b"></span><span class="grad-step" style="background-color:#083b7b"></span><span class="grad-step" style="background-color:#08468c"></span><span class="grad-step" style="background-color:#08468c"></span><span class="grad-step" style="background-color:#08468c"></span><span class="grad-step" style="background-color:#08468c"></span><span class="grad-step" style="background-color:#08519c"></span><span class="grad-step" style="background-color:#08519c"></span><span class="grad-step" style="background-color:#08519c"></span><span class="grad-step" style="background-color:#08519c"></span><span class="grad-step" style="background-color:#105ca4"></span><span class="grad-step" style="background-color:#105ca4"></span><span class="grad-step" style="background-color:#105ca4"></span><span class="grad-step" style="background-color:#105ca4"></span><span class="grad-step" style="background-color:#1966ad"></span><span class="grad-step" style="background-color:#1966ad"></span><span class="grad-step" style="background-color:#1966ad"></span><span class="grad-step" style="background-color:#1966ad"></span><span class="grad-step" style="background-color:#2171b5"></span><span class="grad-step" style="background-color:#2171b5"></span><span class="grad-step" style="background-color:#2171b5"></span><span class="grad-step" style="background-color:#2171b5"></span><span class="grad-step" style="background-color:#2c7cbb"></span><span class="grad-step" style="background-color:#2c7cbb"></span><span class="grad-step" style="background-color:#2c7cbb"></span><span class="grad-step" style="background-color:#2c7cbb"></span><span class="grad-step" style="background-color:#3787c0"></span><span class="grad-step" style="background-color:#3787c0"></span><span class="grad-step" style="background-color:#3787c0"></span><span class="grad-step" style="background-color:#3787c0"></span><span class="grad-step" style="background-color:#4292c6"></span><span class="grad-step" style="background-color:#4292c6"></span><span class="grad-step" style="background-color:#4292c6"></span><span class="grad-step" style="background-color:#4292c6"></span><span class="grad-step" style="background-color:#509bcb"></span><span class="grad-step" style="background-color:#509bcb"></span><span class="grad-step" style="background-color:#509bcb"></span><span class="grad-step" style="background-color:#509bcb"></span><span class="grad-step" style="background-color:#5da5d1"></span><span class="grad-step" style="background-color:#5da5d1"></span><span class="grad-step" style="background-color:#5da5d1"></span><span class="grad-step" style="background-color:#5da5d1"></span><span class="grad-step" style="background-color:#6baed6"></span><span class="grad-step" style="background-color:#6baed6"></span><span class="grad-step" style="background-color:#6baed6"></span><span class="grad-step" style="background-color:#6baed6"></span><span class="grad-step" style="background-color:#7cb7da"></span><span class="grad-step" style="background-color:#7cb7da"></span><span class="grad-step" style="background-color:#7cb7da"></span><span class="grad-step" style="background-color:#7cb7da"></span><span class="grad-step" style="background-color:#8dc1dd"></span><span class="grad-step" style="background-color:#8dc1dd"></span><span class="grad-step" style="background-color:#8dc1dd"></span><span class="grad-step" style="background-color:#8dc1dd"></span><span class="grad-step" style="background-color:#9ecae1"></span><span class="grad-step" style="background-color:#9ecae1"></span><span class="grad-step" style="background-color:#9ecae1"></span><span class="grad-step" style="background-color:#9ecae1"></span><span class="grad-step" style="background-color:#abd0e6"></span><span class="grad-step" style="background-color:#abd0e6"></span><span class="grad-step" style="background-color:#abd0e6"></span><span class="grad-step" style="background-color:#abd0e6"></span><span class="grad-step" style="background-color:#b9d5ea"></span><span class="grad-step" style="background-color:#b9d5ea"></span><span class="grad-step" style="background-color:#b9d5ea"></span><span class="grad-step" style="background-color:#b9d5ea"></span><span class="grad-step" style="background-color:#c6dbef"></span><span class="grad-step" style="background-color:#c6dbef"></span><span class="grad-step" style="background-color:#c6dbef"></span><span class="grad-step" style="background-color:#c6dbef"></span><span class="grad-step" style="background-color:#cee0f2"></span><span class="grad-step" style="background-color:#cee0f2"></span><span class="grad-step" style="background-color:#cee0f2"></span><span class="grad-step" style="background-color:#cee0f2"></span><span class="grad-step" style="background-color:#d6e6f4"></span><span class="grad-step" style="background-color:#d6e6f4"></span><span class="grad-step" style="background-color:#d6e6f4"></span><span class="grad-step" style="background-color:#d6e6f4"></span><span class="grad-step" style="background-color:#deebf7"></span><span class="grad-step" style="background-color:#deebf7"></span><span class="grad-step" style="background-color:#deebf7"></span><span class="grad-step" style="background-color:#deebf7"></span><span class="grad-step" style="background-color:#e6f0fa"></span><span class="grad-step" style="background-color:#e6f0fa"></span><span class="grad-step" style="background-color:#e6f0fa"></span><span class="grad-step" style="background-color:#e6f0fa"></span><span class="grad-step" style="background-color:#eff6fc"></span><span class="grad-step" style="background-color:#eff6fc"></span><span class="grad-step" style="background-color:#eff6fc"></span><span class="grad-step" style="background-color:#eff6fc"></span><span class="grad-step" style="background-color:#f7fbff"></span><span class="grad-step" style="background-color:#f7fbff"></span><span class="grad-step" style="background-color:#f7fbff"></span><span class="grad-step" style="background-color:#f7fbff"></span>
									<span class="domain-min">0%</span><span class="domain-med">50%</span><span class="domain-max">100%</span>
								</div>
								<ul class="list-group list-group-flush">
									<li class="list-group-item p-1"><svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#ffffff" stroke="gray"></rect></svg> <small>no data</small></li>
								</ul>
							</div>
						</div>
						<div id="" class="border rounded p-2 mb-2">
							<div class="m-2">
								<h6>Phenotype Code</h6>
								<ul class="list-group list-group-flush">
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#d62728"></rect></svg>
										<small>Meiotic</small>
									</li>
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#ce6dbd"></rect></svg>
										<small>Post-meiotic</small>
									</li>
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#c7e9c0"></rect></svg>
										<small>Gametes</small>
									</li>
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#c7e9c0"></rect></svg>
										<small>Pre-meiotic</small>
									</li>
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#9edae5"></rect></svg>
										<small>General impairment of spermatogenesis</small>
									</li>
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#fdd0a2"></rect></svg>
										<small>Undetectable</small>
									</li>
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#dadaeb"></rect></svg>
										<small>Unspecified</small>
										</li>
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#bdbdbd"></rect></svg>
										<small>Non-germ cell autonomous</small>
									</li>
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#ffffff" stroke="gray"></rect></svg>
										<small>n/a</small>
									</li>
								</ul>
							</div>
						</div>
						<div id="" class="border rounded p-2 mb-2">
							<div class="m-2">
								<h6>Boolean variables</h6>
								<ul class="list-group list-group-flush">
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#2ca02c"></rect></svg>
										<small>Yes/True</small>
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#ffffff" stroke="gray"></rect></svg>
										<small>No/False/na</small>
									</li>
								</ul>
							</div>
						</div>
						<div id="" class="border rounded p-2 mb-2">
							<div class="m-2">
								<h6>log<sub>2</sub>(FPKM)</h6>
								<div class="gradient">
									<span class="grad-step" style="background-color:#fff5f0"></span><span class="grad-step" style="background-color:#fff5f0"></span><span class="grad-step" style="background-color:#fff5f0"></span><span class="grad-step" style="background-color:#fff5f0"></span><span class="grad-step" style="background-color:#ffeee6"></span><span class="grad-step" style="background-color:#ffeee6"></span><span class="grad-step" style="background-color:#ffeee6"></span><span class="grad-step" style="background-color:#ffeee6"></span><span class="grad-step" style="background-color:#fee7dc"></span><span class="grad-step" style="background-color:#fee7dc"></span><span class="grad-step" style="background-color:#fee7dc"></span><span class="grad-step" style="background-color:#fee7dc"></span><span class="grad-step" style="background-color:#fee0d2"></span><span class="grad-step" style="background-color:#fee0d2"></span><span class="grad-step" style="background-color:#fee0d2"></span><span class="grad-step" style="background-color:#fee0d2"></span><span class="grad-step" style="background-color:#fdd4c2"></span><span class="grad-step" style="background-color:#fdd4c2"></span><span class="grad-step" style="background-color:#fdd4c2"></span><span class="grad-step" style="background-color:#fdd4c2"></span><span class="grad-step" style="background-color:#fdc7b1"></span><span class="grad-step" style="background-color:#fdc7b1"></span><span class="grad-step" style="background-color:#fdc7b1"></span><span class="grad-step" style="background-color:#fdc7b1"></span><span class="grad-step" style="background-color:#fcbba1"></span><span class="grad-step" style="background-color:#fcbba1"></span><span class="grad-step" style="background-color:#fcbba1"></span><span class="grad-step" style="background-color:#fcbba1"></span><span class="grad-step" style="background-color:#fcad91"></span><span class="grad-step" style="background-color:#fcad91"></span><span class="grad-step" style="background-color:#fcad91"></span><span class="grad-step" style="background-color:#fcad91"></span><span class="grad-step" style="background-color:#fca082"></span><span class="grad-step" style="background-color:#fca082"></span><span class="grad-step" style="background-color:#fca082"></span><span class="grad-step" style="background-color:#fca082"></span><span class="grad-step" style="background-color:#fc9272"></span><span class="grad-step" style="background-color:#fc9272"></span><span class="grad-step" style="background-color:#fc9272"></span><span class="grad-step" style="background-color:#fc9272"></span><span class="grad-step" style="background-color:#fc8565"></span><span class="grad-step" style="background-color:#fc8565"></span><span class="grad-step" style="background-color:#fc8565"></span><span class="grad-step" style="background-color:#fc8565"></span><span class="grad-step" style="background-color:#fb7757"></span><span class="grad-step" style="background-color:#fb7757"></span><span class="grad-step" style="background-color:#fb7757"></span><span class="grad-step" style="background-color:#fb7757"></span><span class="grad-step" style="background-color:#fb6a4a"></span><span class="grad-step" style="background-color:#fb6a4a"></span><span class="grad-step" style="background-color:#fb6a4a"></span><span class="grad-step" style="background-color:#fb6a4a"></span><span class="grad-step" style="background-color:#f75a40"></span><span class="grad-step" style="background-color:#f75a40"></span><span class="grad-step" style="background-color:#f75a40"></span><span class="grad-step" style="background-color:#f75a40"></span><span class="grad-step" style="background-color:#f34b36"></span><span class="grad-step" style="background-color:#f34b36"></span><span class="grad-step" style="background-color:#f34b36"></span><span class="grad-step" style="background-color:#f34b36"></span><span class="grad-step" style="background-color:#ef3b2c"></span><span class="grad-step" style="background-color:#ef3b2c"></span><span class="grad-step" style="background-color:#ef3b2c"></span><span class="grad-step" style="background-color:#ef3b2c"></span><span class="grad-step" style="background-color:#e32f27"></span><span class="grad-step" style="background-color:#e32f27"></span><span class="grad-step" style="background-color:#e32f27"></span><span class="grad-step" style="background-color:#e32f27"></span><span class="grad-step" style="background-color:#d72422"></span><span class="grad-step" style="background-color:#d72422"></span><span class="grad-step" style="background-color:#d72422"></span><span class="grad-step" style="background-color:#d72422"></span><span class="grad-step" style="background-color:#cb181d"></span><span class="grad-step" style="background-color:#cb181d"></span><span class="grad-step" style="background-color:#cb181d"></span><span class="grad-step" style="background-color:#cb181d"></span><span class="grad-step" style="background-color:#be151a"></span><span class="grad-step" style="background-color:#be151a"></span><span class="grad-step" style="background-color:#be151a"></span><span class="grad-step" style="background-color:#be151a"></span><span class="grad-step" style="background-color:#b21218"></span><span class="grad-step" style="background-color:#b21218"></span><span class="grad-step" style="background-color:#b21218"></span><span class="grad-step" style="background-color:#b21218"></span><span class="grad-step" style="background-color:#a50f15"></span><span class="grad-step" style="background-color:#a50f15"></span><span class="grad-step" style="background-color:#a50f15"></span><span class="grad-step" style="background-color:#a50f15"></span><span class="grad-step" style="background-color:#900a12"></span><span class="grad-step" style="background-color:#900a12"></span><span class="grad-step" style="background-color:#900a12"></span><span class="grad-step" style="background-color:#900a12"></span><span class="grad-step" style="background-color:#7c0510"></span><span class="grad-step" style="background-color:#7c0510"></span><span class="grad-step" style="background-color:#7c0510"></span><span class="grad-step" style="background-color:#7c0510"></span><span class="grad-step" style="background-color:#67000d"></span><span class="grad-step" style="background-color:#67000d"></span><span class="grad-step" style="background-color:#67000d"></span><span class="grad-step" style="background-color:#67000d"></span><span class="grad-step" style="background-color:#67000d">								
									<span class="domain-min">0%</span><span class="domain-med">50%</span><span class="domain-max">100%</span>
								</div>							
							</div>
						</div>
						<div id="" class="border rounded p-2 mb-2">
							<div class="m-2">
								<h6>Biotype</h6>
								<ul class="list-group list-group-flush">
									<li class="list-group-item p-1">
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#d62728"></rect></svg>
										<small>Protein-coding</small>
										<svg class="bd-placeholder-img rounded mr-1" width="17" height="17" xmlns="http://www.w3.org/2000/svg" focusable="false" role="img"><rect width="100%" height="100%" fill="#ffffff" stroke="gray" stroke-linecap="round"></rect></svg>
										<small>Other</small>
									</li>
								</ul>
							</div>
						</div>
						<div id="" class="border rounded p-2 mb-2">
							<div class="m-2">
								<h6>Modules</h6>
								<div class="gradient border-bottom">
									<span class="grad-step" style="background-color:#a6cee3"></span><span class="grad-step" style="background-color:#a6cee3"></span><span class="grad-step" style="background-color:#a6cee3"></span><span class="grad-step" style="background-color:#a6cee3"></span><span class="grad-step" style="background-color:#a6cee3"></span><span class="grad-step" style="background-color:#589cc8"></span><span class="grad-step" style="background-color:#589cc8"></span><span class="grad-step" style="background-color:#589cc8"></span><span class="grad-step" style="background-color:#589cc8"></span><span class="grad-step" style="background-color:#589cc8"></span><span class="grad-step" style="background-color:#3688ad"></span><span class="grad-step" style="background-color:#3688ad"></span><span class="grad-step" style="background-color:#3688ad"></span><span class="grad-step" style="background-color:#3688ad"></span><span class="grad-step" style="background-color:#3688ad"></span><span class="grad-step" style="background-color:#8bc495"></span><span class="grad-step" style="background-color:#8bc495"></span><span class="grad-step" style="background-color:#8bc495"></span><span class="grad-step" style="background-color:#8bc495"></span><span class="grad-step" style="background-color:#8bc495"></span><span class="grad-step" style="background-color:#8acb6c"></span><span class="grad-step" style="background-color:#8acb6c"></span><span class="grad-step" style="background-color:#8acb6c"></span><span class="grad-step" style="background-color:#8acb6c"></span><span class="grad-step" style="background-color:#8acb6c"></span><span class="grad-step" style="background-color:#40a736"></span><span class="grad-step" style="background-color:#40a736"></span><span class="grad-step" style="background-color:#40a736"></span><span class="grad-step" style="background-color:#40a736"></span><span class="grad-step" style="background-color:#40a736"></span><span class="grad-step" style="background-color:#929d60"></span><span class="grad-step" style="background-color:#929d60"></span><span class="grad-step" style="background-color:#929d60"></span><span class="grad-step" style="background-color:#929d60"></span><span class="grad-step" style="background-color:#929d60"></span><span class="grad-step" style="background-color:#fa9392"></span><span class="grad-step" style="background-color:#fa9392"></span><span class="grad-step" style="background-color:#fa9392"></span><span class="grad-step" style="background-color:#fa9392"></span><span class="grad-step" style="background-color:#fa9392"></span><span class="grad-step" style="background-color:#ec494a"></span><span class="grad-step" style="background-color:#ec494a"></span><span class="grad-step" style="background-color:#ec494a"></span><span class="grad-step" style="background-color:#ec494a"></span><span class="grad-step" style="background-color:#ec494a"></span><span class="grad-step" style="background-color:#e83d2d"></span><span class="grad-step" style="background-color:#e83d2d"></span><span class="grad-step" style="background-color:#e83d2d"></span><span class="grad-step" style="background-color:#e83d2d"></span><span class="grad-step" style="background-color:#e83d2d"></span><span class="grad-step" style="background-color:#f89c5e"></span><span class="grad-step" style="background-color:#f89c5e"></span><span class="grad-step" style="background-color:#f89c5e"></span><span class="grad-step" style="background-color:#f89c5e"></span><span class="grad-step" style="background-color:#f89c5e"></span><span class="grad-step" style="background-color:#fea746"></span><span class="grad-step" style="background-color:#fea746"></span><span class="grad-step" style="background-color:#fea746"></span><span class="grad-step" style="background-color:#fea746"></span><span class="grad-step" style="background-color:#fea746"></span><span class="grad-step" style="background-color:#ff8206"></span><span class="grad-step" style="background-color:#ff8206"></span><span class="grad-step" style="background-color:#ff8206"></span><span class="grad-step" style="background-color:#ff8206"></span><span class="grad-step" style="background-color:#ff8206"></span><span class="grad-step" style="background-color:#e39a71"></span><span class="grad-step" style="background-color:#e39a71"></span><span class="grad-step" style="background-color:#e39a71"></span><span class="grad-step" style="background-color:#e39a71"></span><span class="grad-step" style="background-color:#e39a71"></span><span class="grad-step" style="background-color:#c0a6d0"></span><span class="grad-step" style="background-color:#c0a6d0"></span><span class="grad-step" style="background-color:#c0a6d0"></span><span class="grad-step" style="background-color:#c0a6d0"></span><span class="grad-step" style="background-color:#c0a6d0"></span><span class="grad-step" style="background-color:#8862ad"></span><span class="grad-step" style="background-color:#8862ad"></span><span class="grad-step" style="background-color:#8862ad"></span><span class="grad-step" style="background-color:#8862ad"></span><span class="grad-step" style="background-color:#8862ad"></span><span class="grad-step" style="background-color:#91709a"></span><span class="grad-step" style="background-color:#91709a"></span><span class="grad-step" style="background-color:#91709a"></span><span class="grad-step" style="background-color:#91709a"></span><span class="grad-step" style="background-color:#91709a"></span><span class="grad-step" style="background-color:#e7e099"></span><span class="grad-step" style="background-color:#e7e099"></span><span class="grad-step" style="background-color:#e7e099"></span><span class="grad-step" style="background-color:#e7e099"></span><span class="grad-step" style="background-color:#e7e099"></span><span class="grad-step" style="background-color:#deb969"></span><span class="grad-step" style="background-color:#deb969"></span><span class="grad-step" style="background-color:#deb969"></span><span class="grad-step" style="background-color:#deb969"></span><span class="grad-step" style="background-color:#deb969"></span><span class="grad-step" style="background-color:#b15928"></span><span class="grad-step" style="background-color:#b15928"></span><span class="grad-step" style="background-color:#b15928"></span><span class="grad-step" style="background-color:#b15928"></span><span class="grad-step" style="background-color:#b15928"></span><span class="grad-step" style="background-color:#b15928"></span>
									<span class="domain-min">0</span><span class="domain-med">...</span><span class="domain-max">n</span>
								</div>
								<p><small>Colors are dynamicaly assigned based on the total number of modules.</small></p>
							</div>
						</div>
					</div>
				</div>
			</div>

		<!-- This is an invisible canvas for the cross-layer nodes & edges -->
		<div id="canvas-cross" class="border rounded d-none" style="height:300px"></div>
				
	</main>

	<footer class="text-muted">
		<hr>
		<div class="container py-4">
			<p class="float-right"><a href="#">Back to top</a></p>
			<ul class="list-unstyled list-inline">
				<li class="list-inline-item"><a href="<?php print(DOCUMENT_ROOT);?>index">Home</a></li>
				<li class="list-inline-item"><a href="<?php print(DOCUMENT_ROOT);?>publications">Publications</a></li>
				<li class="list-inline-item"><a href="<?php print(DOCUMENT_ROOT);?>privacy">Privacy</a></li>
				<li class="list-inline-item"><a href="<?php print(DOCUMENT_ROOT);?>about">About</a></li>
			</ul>
			<p>
				<a href="www.igc.gulbenkian.pt" title="Instituto Gulbenkian de Ciência"> Instituto Gulbenkian de Ciência</a>
				|
				This application is still beta, <em>build: 0.1-beta</em>
				|
				<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/" target="_blank">
					<img alt="Creative Commons License" style="border-width:0" src="img/cc-80x15.png" />
				</a>
			</p>
		</div>
	</footer>

	<?php include "include/javascripts.php"; ?>
	<script type="text/javascript" src="js/network.js"></script>
</body>
</html>