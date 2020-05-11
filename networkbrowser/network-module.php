<?php
include_once("config.php");
$curpage = "Network Module";


$layer = (string)preg_replace('/[^A-Za-z]/', '', $_GET['layer']);
$algorithm = (string)preg_replace('/[^A-Za-z]/', '', $_GET['algorithm']);	
$module = (string)preg_replace('/[^0-9A-Za-z]/', '', $_GET['module']);	
#
switch($layer) {
	case 'HS':
		$species = 'Homo sapiens';
		break;
	case 'MM':
		$species = 'Mus musculus';
		break;
	case 'DM':
		$species = 'Drosophila melanogaster';
		break;
	default:
		$species = 'Error';
		break;
}
#
$pagetitle = "Network Module : " . $layer . " : " . $algorithm . " : " . $module . " : SpermNet";
#
$file =  DOMAIN.'json/net_threshold-0p5-' . $layer . '-' . $algorithm . '-' . $module . '.json';
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
					<li class="breadcrumb-item">Interactive Networks</li>
					<li class="breadcrumb-item"><?=$species?></li>
					<li  class="breadcrumb-item active">SVD module <?=$module?></li>
				</ol>
			</nav>

			<div class="py-5 mb-4">
				<div class="row">
					<div class="col-md-9 col-ld-7">
						<h1>Network SVD module <?=$module?></h1>
						<p>
							The network below shows a Singular Value Decomposition (SVD) module from the <em><?=$species?></em> network layer.
						</p>
					</div>
				</div>
			</div>

			<!-- Left Canvas -->	
			<div class="row">
				<div class="col-9 p-2">

					<!-- Network Module Canvas -->
					<div class="position-relative mb-2">
						<div class="position-absolute m-2"> 
							<h5 class="m-0 p-0"><em><?=$species?></em></h5>
							<h7 class="m-1 p-1">Algorithm: <?=$algorithm?> | Module: <?=$module?></h7>
						</div>
						<div id="controls" class="position-absolute fixed-bottom border rounded form-inline p-2 m-2 bg-light">
							<div class="input-group">
								<div class="input-group-prepend">
									<label for="select-find" class="input-group-text">
										<span class="oi oi-magnifying-glass" aria-hidden="true" data-toggle="tooltip" title="Locates the node on the network."></span>
									</label>
								</div>
								<select id="select-find" class="custom-select">

								</select>
								<div class="input-group-append">
									<button class="btn btn-primary" type="button" onClick="zoomInNode();">Find</button>
								</div>
							</div>
							<div class="input-group mx-2">
								<div class="input-group-prepend">
									<label for="select-variable" class="input-group-text">
										<span class="oi oi-beaker" aria-hidden="true" data-toggle="tooltip" title="Color the network according to specific variable."></span>
									</label>
								</div>
								<select id="select-variable" class="custom-select">
									<?php
									if ($layer == 'HS'):
									?>
									<option value="core" data-type="discrete" data-cmap="cmap_Boolean">Core gene</option>
									<option value="mammals" data-type="discrete" data-cmap="cmap_Boolean">Mammal gene</option>
									<option value="mean-fert-rate" data-type="continuous" data-cmap="cmap_Blues_r">DM Fertility Rate (homologs)</option>
									<option value="new-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">New DM phenotype (homologs)</option>
									<option value="known-HS-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known HS phenotype</option>
									<option value="known-MM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known MM phenotype (homologs)</option>
									<option value="known-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known DM phenotype (homologs)</option>
									<option value="logFPKM" data-type="continuous" data-cmap="cmap_Reds">log2(FPKM)</option>
									<option value="meiotic-entry" data-type="discrete" data-cmap="cmap_Boolean">Meiotic entry</option>
									<option value="meiotic-exit" data-type="discrete" data-cmap="cmap_Boolean">Meiotic exit</option>
									<option value="biotype" data-type="discrete" data-cmap="cmap_Biotype">Biotype</option>
									<?php
									elseif ($layer == 'MM'):
									?>
									<option value="core" data-type="discrete" data-cmap="cmap_Boolean">Core gene</option>
									<option value="mammals" data-type="discrete" data-cmap="cmap_Boolean">Mammal gene</option>
									<option value="mean-fert-rate" data-type="continuous" data-cmap="cmap_Blues_r">DM Fertility Rate (homologs)</option>
									<option value="new-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">New DM phenotype (homologs)</option>
									<option value="known-MM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known MM phenotype</option>
									<option value="known-HS-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known HS phenotype (homologs)</option>
									<option value="known-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known DM phenotype (homologs)</option>
									<option value="logFPKM" data-type="continuous" data-cmap="cmap_Reds">log2(FPKM)</option>
									<option value="meiotic-entry" data-type="discrete" data-cmap="cmap_Boolean">Meiotic entry</option>
									<option value="meiotic-exit" data-type="discrete" data-cmap="cmap_Boolean">Meiotic exit</option>
									<option value="biotype" data-type="discrete" data-cmap="cmap_Biotype">Biotype</option>
									<?php
									elseif ($layer == 'DM'):
									?>			
									<option value="core" data-type="discrete" data-cmap="cmap_Boolean">Core gene</option>
									<option value="mean-fert-rate" data-type="continuous" data-cmap="cmap_Blues_r">Fertility Rate</option>
									<option value="new-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">New DM phenotype</option>
									<option value="known-DM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known DM phenotype</option>
									<option value="known-MM-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known MM phenotype (homologs)</option>
									<option value="known-HS-phenotype" data-type="discrete" data-cmap="cmap_pheno_code">Known HS phenotype (homologs)</option>
									<option value="logFPKM" data-type="continuous" data-cmap="cmap_Reds">log2(FPKM)</option>
									<option value="RNAi" data-type="discrete" data-cmap="cmap_RNAi">Validated RNAi</option>
									<option value="meiotic-entry" data-type="discrete" data-cmap="cmap_Boolean">Meiotic entry</option>
									<option value="meiotic-exit" data-type="discrete" data-cmap="cmap_Boolean">Meiotic exit</option>
									<option value="biotype" data-type="discrete" data-cmap="cmap_Biotype">Biotype</option>
									<?php
									endif;
									?>
								</select>
								<div class="input-group-append">
									<button class="btn btn-primary" type="button" onClick="colorNetworkFromVariable();">View</button>
								</div>
							</div>
							<div class="input-group mx-2">
								<div class="input-group-prepend">
									<label for="input-layout" class="input-group-text">
										<span class="oi oi-sun" aria-hidden="true" data-toggle="tooltip" title="Organize the network layout."></span>
									</label>
								</div>
								<div class="input-group-append">
									<button class="btn btn-light" type="button" id="btn-simulation" onclick='startStopSimulation(this)'>Run</button> 
								</div>
							</div>
						</div>
						<div id="canvas" class="border rounded" style="height:900px" data-file="<?=$file?>" data-layer="<?=$layer?>" data-algorithm="<?=$algorithm?>" data-module="<?=$module?>"></div>
						<div id="progressbars" class="position-absolute mx-auto my-auto fixed-top" style="width:600px; top:35%">
							<!-- Error Msg -->
							<div class="text-center mb-3">
								<p><span class="oi oi-dashboard"></span> <span id="status-msg">Loading scripts.</span></p>
							</div>
							Network force layout stabilization
							<div id="" class="progress">
								<div id="progressbar-stab" class="progress-bar bg-success progress-bar-striped active" role="progressbar" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100">
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
				</div>
	
				<!-- Right Legend -->
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
									<span class="domain-min">0.0</span><span class="domain-med">5.0</span><span class="domain-max">10.0</span>
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
					</div>
				</div>
			</div>

	</main>

	<?php include "include/footer.php"; ?>
	
	<?php include "include/javascripts.php"; ?>
	<script type="text/javascript" src="js/network-module.js"></script>
</body>
</html>