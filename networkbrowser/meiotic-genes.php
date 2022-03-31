<?php
include_once("config.php");
$curpage = "meiotic-genes";
$pagetitle = "Functional analysis of meiotic genes : Meiotic Navigator";
$pagedescription = "Navigate the functional analysis of conserved meiotic genes.";

$filephp =  'json/page-meiotic-genes.json';
$filejs =  PATH . 'json/page-meiotic-genes.json';
?>
<!doctype html>
<html lang="en">
<head>
	<meta charset="utf-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge">
	<meta name="viewport" content="width=device-width, initial-scale=1">
	<meta name="description" content="<?php print($pagedescription);?>">
	
	<!-- Bootstrap -->
	<link href="<?=PATH?>css/bootstrap.min.css" rel="stylesheet" type="text/css" />

	<!-- Bootstrap Icons -->
	<link href="<?=PATH?>css/bootstrap-icons.css" rel="stylesheet">

    <link href="<?=PATH?>css/andronet.css" rel="stylesheet" type="text/css" />
    <style type="text/css">
	
		/* Points / Text / Errorbar / area */
		svg .names text {
			font-size: small;
		}
		svg text.label {
			font-size: smaller;
		}
		svg .items circle {
			stroke: black;
		}
		svg .items rect {
			stroke: black;
		}
		svg g.errorbar line.errorbar {
			stroke: black;
			stroke-width: 1px;
			fill: none;
		}
		svg g.selection #bar {
			fill: #ffc107;
			fill-opacity:0.2;
			stroke: black;
			stroke-opacity:0.2;
		}
		svg#chart .area {
			fill: #c7c7c7;
			fill-opacity: 0.5;
		}
		svg#chart .ui-rect {
			cursor: pointer;
		}
		svg#focus .area {
			fill: #7f7f7f;
			fill-opacity: 1.0;
		}
		
		/* Coding colors */
		.code-A { fill: #d62728; color: #d62728; }
		.code-B { fill: #ce6dbd; color: #ce6dbd; }
		.code-C { fill: #756bb1; color: #756bb1; }
		.code-D { fill: #c7e9c0; color: #c7e9c0; }
		.code-E { fill: #9edae5; color: #9edae5; }
		.code-F { fill: #fdd0a2; color: #fdd0a2; }
		.code-G { fill: #dadaeb; color: #dadaeb; }
		.code-H { fill: #bdbdbd; color: #bdbdbd; }
		.code-null { fill:white; color: #7f7f7f; stroke: #7f7f7f; stroke-width: 1 }
		.code-true { fill: #7f7f7f; color: #7f7f7f; }
	
		/* Grid */
		svg .grid line {
			fill: none;
			stroke: lightgrey;
			stroke-opacity: 0.6;
			shape-rendering: crispEdges;
		}
		svg .grid path {
			stroke-width:0;
		}
		
		/* Legend */
		svg .legend .cell rect {
			stroke: black;
		}
		
		/* Annotations */
		svg .annotation {
			font-size: small;
			fill: lightgray;
		}
		
    </style>

	<!-- Global site tag (gtag.js) - Google Analytics -->
	<script async src="https://www.googletagmanager.com/gtag/js?id=G-W81S0LJ1PT"></script>
	<script>
	  window.dataLayer = window.dataLayer || [];
	  function gtag(){dataLayer.push(arguments);}
	  gtag('js', new Date());
	  gtag('config', 'G-W81S0LJ1PT');
	</script>
    
	<title><?php print($pagetitle);?></title>
</head>
<body>
	<header>
		<?php include "include/navigation.php"; ?>
		<div class="container">
			<nav aria-label="breadcrumb">
				<ol class="breadcrumb">
					<li class="breadcrumb-item"><a href="<?=PATH?>index.php">Home</a></li>
					<li  class="breadcrumb-item active">Functional analysis of conserved meiotic genes</li>
				</ol>
			</nav>
		</div>
	</header>

	<main role="main">
		<div class="container">

			<div class="pt-5">
				<h1>Functional analysis of conserved meiotic genes </h1>
				<div class="row">
					<div class="col-12">
						<p>
							We used <em>in vivo</em> RNAi in insects (the fruit fly <em>Drosophila melanogaster</em>) to silence 920 evolutionarily-conserved meiotic genes specifically as male germ cells entered meiosis (GAL4/UAS system; driver: <em>bam-GAL4</em>).
							Fertility rate measures male reproductive fitness and corresponds to the average egg eclosion rate of four independent mating experiments of RNAi-silenced males with wild-type females. Gene colours depict the testicular phenotype (as assessed by phase-contrast microscopy, see phenotype legend).
						</p>
						<div class="row py-4">
							<div class="col-4">
								<div class="card py-2 mb-4 border-muted">
									<div class="card-body">
										<h5 class="card-title">Was my favourite Drosophila gene tested?</h5>
										Use the <span class="bg-warning">drop-down list</span> to select your gene of interest, and explore all associated data.
									</div>
								</div>
							</div>
							<div class="col-8">
								<div class="card py-2 mb-4 border-muted">
									<div class="card-body">
										<h5 class="card-title">I’m interested in a mouse or human gene. Was the corresponding Drosophila ortholog tested?</h5>
										Use CTRL-F (Command-F in Mac) to search the cross-species list for the corresponding Drosophila ortholog of your mammalian gene of interest. Then use the <span class="bg-warning">drop-down list</span> to select the identified insect ortholog, and explore all associated data.
									</div>
								</div>
							</div>
						</div>
					</div>
				</div>
			</div>

			<div class="row">
			<div class="input-group mb-3 offset-2 col-4">
				<div class="input-group-prepend">
					<label class="input-group-text bg-warning" for="selector">Fruit fly gene</label>
				</div>
				<select id="selector" name="selector" class="custom-select"><option value="">loading...</option></select>
			</div>
			</div>

			<div id="canvas" class="" data-file="<?=$filejs?>"></div>
			
			<div class="card my-4">
				<div class="card-body">
					<h4 class="card-title">Legend</h4>
					<h6 class="card-subtitle text-muted mb-2">Circle colors depict the observed testicular phenotype.</h6>
					<ul class="list-inline card-text">
						<li class="list-inline-item"><i class="bi bi-circle-fill code-D"></i> Pre-meiotic</li>
						<li class="list-inline-item"><i class="bi bi-circle-fill code-A"></i> Meiotic</li>
						<li class="list-inline-item"><i class="bi bi-circle-fill code-B"></i> Post-meiotic</li>
						<li class="list-inline-item"><i class="bi bi-circle-fill code-C"></i> Gametes</li>
						<li class="list-inline-item"><i class="bi bi-circle-fill code-F"></i> Undetectable</li>
						<li class="list-inline-item"><i class="bi bi-circle code-null"></i> No testicular analysis (fertility >75%)”</li>
					</ul>
					<h6 class="card-subtitle text-muted mb-2">Square colors depict RNAseq expression levels in wild-type testis (Vedelek et al., 2018).</h6>
					<div id="legend" class="mb-2"></div>
					<h6 class="card-subtitle text-muted mb-2">Gray squares denote previously reported spermatogenesis / male infertility phenotypes in humans ("HS pheno"), mice ("MM pheno") or fruit flies ("DM pheno").</h6>
				</div>
			</div>


			<div class="my-5">				
				<h3 class="my-3"><i class="bi bi-clipboard-data"></i> Gene Information</h3>
				<div class="row">
					<div class="col-sm-3 py-2 text-right font-weight-bold border-bottom">Gene:</div>
					<div id="gene-name" class="col-sm-9 py-2 lead border-bottom"><span class="text-muted">Select gene above.</span></div>

					<div class="col-sm-3 py-2 text-right font-weight-bold border-bottom">
						<i class="bi bi-info-circle" data-toggle="tooltip" data-placement="top" title="Stock identifiers refer to the BDSC (numerals only), VDRC (numerals preceded by 'v'), or in-house generated (indicated by 'ih') Drosophila lines."></i>
						Tested stock:
					</div>
					<div id="gene-stock" class="col-sm-9 py-2 border-bottom">-</div>
					
					<div class="col-sm-3 py-2 text-right font-weight-bold border-bottom">Summary (FlyBase):</div>
					<div id="gene-summary" class="col-sm-9 py-2 border-bottom">-</div>
					
					<div class="col-sm-3 py-2 text-right font-weight-bold border-bottom">
						<i class="bi bi-info-circle" data-toggle="tooltip" data-placement="top" title="Orthologs are defined based on eggNOG Orthologous Groups and on the expression of these genes in human spermatogenesis."></i>
						Human orthologs:
 					</div>
					<div id="gene-orthologs-human" class="col-sm-9 py-2 border-bottom">-</div>
					
					<div class="col-sm-3 py-2 text-right font-weight-bold border-bottom">
						<i class="bi bi-info-circle" data-toggle="tooltip" data-placement="top" title="Orthologs are defined based on eggNOG Orthologous Groups and on the expression of these genes in mouse spermatogenesis."></i>
						Mouse orthologs:
 					</div>
					<div id="gene-orthologs-mouse" class="col-sm-9 py-2 border-bottom">-</div>
					
					<div class="col-sm-3 py-2 text-right font-weight-bold border-bottom">See in FlyBase:</div>
					<div id="" class="col-sm-9 py-2 border-bottom"><a id="gene-flybase" class="btn btn-warning disabled" href="#" target="_blank"><i class="bi bi-link-45deg"></i>link</a></div>
				</div>
				
				
				<p>
				</p>
			</div>

			<div class="my-5">
				<h3 class="my-3"><i class="bi bi-search"></i> Testicular phenotype</h3>
				<div class="row my-3">
					<div class="col-md-3">		
						<figure class="figure">
							<figcaption class="figure-caption text-center mb-1" style="font-size:120%">Apical testis</figcaption>
							<a id="pheno-url-1" href="#" onClick="return false;">
								<img id="pheno-img-1" src="<?=PATH?>img/results/placeholder_240x240.jpg" width="240" alt="Apical testis" class="figure-img img-fluid rounded" />
							</a>
						</figure>
					</div>
					<div class="col-md-3">
						<figure class="figure">
							<figcaption class="figure-caption text-center mb-1" style="font-size:120%">Median testis</figcaption>
							<a id="pheno-url-2" href="#" onClick="return false;">
								<img id="pheno-img-2" src="<?=PATH?>img/results/placeholder_240x240.jpg" width="240" alt="Median testis" class="figure-img img-fluid rounded" />
							</a>
						</figure>
					</div>
					<div class="col-md-3">
						<figure class="figure">
							<figcaption class="figure-caption text-center mb-1" style="font-size:120%">Distal testis</figcaption>
							<a id="pheno-url-3" href="#" onClick="return false;">
								<img id="pheno-img-3" src="<?=PATH?>img/results/placeholder_240x240.jpg" width="240" alt="Distal testis" class="figure-img img-fluid rounded" />
							</a>
						</figure>
					</div>
					<div class="col-md-3">
						<figure class="figure">
							<figcaption class="figure-caption text-center mb-1" style="font-size:120%">Gametes</figcaption>
							<a id="pheno-url-4" href="#" onClick="return false;">
								<img id="pheno-img-4" src="<?=PATH?>img/results/placeholder_240x240.jpg" width="240" alt="Gametes" class="figure-img img-fluid rounded" />
							</a>
						</figure>
					</div>
				</div>
			</div>
			
			<div class="my-5">
				<h3 class="my-3"><i class="bi bi-clipboard-minus"></i> Control testis</h3>
				<div class="row mb-3">
					<div class="col-md-3">
						<figure class="figure">
							<figcaption class="figure-caption text-center mb-1" style="font-size:120%">Apical testis</figcaption>
							<a id="control-4" href="#" onClick="toggleModal('Control - Apical testis', '<?=PATH?>img/results/control/mC_1.jpg'); return false;">
								<img id="pheno-img-4" src="<?=PATH?>img/results/control/mC_1_thumb.jpg" width="240" alt="Apical testis" class="figure-img img-fluid rounded" />
							</a>
						</figure>
					</div>
					<div class="col-md-3">
						<figure class="figure">
							<figcaption class="figure-caption text-center mb-1" style="font-size:120%">Median testis</figcaption>
							<a id="control-4" href="#" onClick="toggleModal('Control - Median testis', '<?=PATH?>img/results/control/mC_2.jpg'); return false;">
								<img id="pheno-img-4" src="<?=PATH?>img/results/control/mC_2_thumb.jpg" width="240" alt="Median testis" class="figure-img img-fluid rounded" />
							</a>
						</figure>
					</div>
					<div class="col-md-3">
						<figure class="figure">
							<figcaption class="figure-caption text-center mb-1" style="font-size:120%">Distal testis</figcaption>
							<a id="control-4" href="#" onClick="toggleModal('Control - Distal testis', '<?=PATH?>img/results/control/mC_3.jpg'); return false;">
								<img id="pheno-img-4" src="<?=PATH?>img/results/control/mC_3_thumb.jpg" width="240" alt="Distal testis" class="figure-img img-fluid rounded" />
							</a>
						</figure>
					</div>
					<div class="col-md-3">
						<figure class="figure">
							<figcaption class="figure-caption text-center mb-1" style="font-size:120%">Gametes</figcaption>
							<a id="control-4" href="#" onClick="toggleModal('Control - Gametes', '<?=PATH?>img/results/control/mC_4.jpg'); return false;">
								<img id="pheno-img-4" src="<?=PATH?>img/results/control/mC_4_thumb.jpg" alt="Gametes" class="figure-img img-fluid rounded" />
							</a>
						</figure>
					</div>
				</div>
			</div>

<?php
# Load Data
$json = json_decode(file_get_contents($filephp), TRUE);
# number of Items
$n = count($json);

# Function that parses a certain proportion of results and results array of results
function iter_json($json, $i_start=0, $i_end=10) {
	#
	$arr_results = array();
	#
	for ($i = $i_start; $i < $i_end; $i++) {
		$arr = $json[$i];
		$gene_dm = $arr['gene'];
		$mean_fert_rate = $arr['mean(fert-rate)'];
		# MM orthologs
		$genes_mm = array();
		foreach ($arr["orthologs-MM"] as $oarr) {
			$oid = $oarr["id"];
			$ogene = $oarr["gene"];
			array_push($genes_mm, $ogene);
		}
		$genes_mm_html = implode(', ', $genes_mm);
		
		# HS orthologs
		$genes_hs = array();
		foreach ($arr["orthologs-HS"] as $oarr) {
			$oid = $oarr["id"];
			$ogene = $oarr["gene"];
			//
			array_push($genes_hs, $ogene);
		}
		$genes_hs_html = implode(', ', $genes_hs);
		#
		array_push($arr_results,
			array(
				'id' => $i,
				'mean-fert-rate' => $mean_fert_rate,
				'gene_dm' => $gene_dm,
				'genes_mm' => $genes_mm_html,
				'genes_hs' => $genes_hs_html
			));
	}
	return $arr_results;
}

?>
			<div class="my-5">
				<h3 class="my-3"><i class="bi bi-arrow-left-right"></i> Cross species gene correspondence <small>(ortholog list)</small></h3>
				<p>
				<div class="row my-3">
					<div class="col-md-4">		
						<table class="table-sm table-striped table-hover" style="font-size:80%">
							<thead>
								<tr>
									<th>Insect</th>
									<th>Mouse</th>
									<th>Human</th>
								</tr>
							</thead>
							<tbody>
								<?php
								// Items 0-289
								$list = iter_json($json, 0, 340);
								foreach($list as $key => $item) {
									$gene_dm = $item['gene_dm'];
									$mean_fert_rate = $item['mean-fert-rate'];
									$genes_mm = $item['genes_mm'];
									$genes_hs = $item['genes_hs'];
								?>
								<tr>
									<td><?=$gene_dm?></td>
									<td><?=$genes_mm?></td>
									<td><?=$genes_hs?></td>
								</tr>
								<?php
								} // end foreach
								?>
							</tbody>
						</table>
					</div>
					<div class="col-md-4">
						<table class="table-sm table-striped table-hover" style="font-size:80%">
							<thead>
								<tr>
									<th>Insect</th>
									<th>Mouse</th>
									<th>Human</th>
								</tr>
							</thead>
							<tbody>
								<?php
								// Items 289-600
								$list = iter_json($json, 340, 660);
								foreach($list as $key => $item) {
									$gene_dm = $item['gene_dm'];
									$mean_fert_rate = $item['mean-fert-rate'];
									$genes_mm = $item['genes_mm'];
									$genes_hs = $item['genes_hs'];
								?>
								<tr>
									<td><?=$gene_dm?></td>
									<td><?=$genes_mm?></td>
									<td><?=$genes_hs?></td>
								</tr>
								<?php
								} // end foreach
								?>
							</tbody>
						</table>
					</div>
					<div class="col-md-4">
						<table class="table-sm table-striped table-hover" style="font-size:80%">
							<thead>
								<tr>
									<th>Insect</th>
									<th>Mouse</th>
									<th>Human</th>
								</tr>
							</thead>
							<tbody>
								<?php
								// Items 600-919
								$list = iter_json($json, 660, 919);
								foreach($list as $key => $item) {
									$gene_dm = $item['gene_dm'];
									$mean_fert_rate = $item['mean-fert-rate'];
									$genes_mm = $item['genes_mm'];
									$genes_hs = $item['genes_hs'];
								?>
								<tr>
									<td><?=$gene_dm?></td>
									<td><?=$genes_mm?></td>
									<td><?=$genes_hs?></td>
								</tr>
								<?php
								} // end foreach
								?>							</tbody>
						</table>
					</div>

				</div>
			</div>

		</div> <!-- /end Container -->

	</main>

	<!-- Modal -->
	<div class="modal fade" id="modal" tabindex="-1" role="dialog" aria-labelledby="modal" aria-hidden="true">
		<div class="modal-dialog modal-dialog-centered modal-lg" role="document">
			<div class="modal-content">
				<div class="modal-header">
					<h5 class="modal-title" id="modal-title">Testicular phenotype</h5>
					<button type="button" class="close" data-dismiss="modal" aria-label="Close">
						<span aria-hidden="true">&times;</span>
					</button>
				</div>
				<div class="modal-body">
					<img id="modal-img" src="<?=PATH?>holder.js/800x800" class="figure-img img-fluid rounded" alt="Testicular phenotype" />
					<p class="text-right">Scale bar: 20 &mu;m.</p>
				</div>
				<div class="modal-footer">
					<button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
				</div>
			</div>
		</div>
	</div>
	<?php include "include/footer.php"; ?>
	<?php include "include/javascripts.php"; ?>
	<script type="text/javascript" src="<?=PATH?>js/page-meiotic-genes.js"></script>
	<script>
		/* Enable tooltips */
		$(function () {
			$('[data-toggle="tooltip"]').tooltip()
		})
	</script>
</body>
</html>