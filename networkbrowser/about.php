<?php
include_once("config.php");
$curpage = "About";
$pagetitle = "About : SpermNet";
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
		<div class="container">
			<nav aria-label="breadcrumb">
				<ol class="breadcrumb">
					<li class="breadcrumb-item"><a href="">Home</a></li>
					<li  class="breadcrumb-item active">About</li>
				</ol>
			</nav>
		</div>
	</header>

	<main role="main">
		<div class="container">

			<div class="py-5 mb-4">
				<div class="row">
					<div class="col-md-9 col-ld-7">
						<h1>About</h1>
					</div>
					<div class="col-8">
						<p>We are scientists working to solve complex issues in human health using biological and computational tools. <a href="/index">SpermNet<sup style="color:darkgreen; font-style:italic;">beta</sup></a> is a pipeline that leverages large scale computational data and wet-lab experimental & comparative biology to advane our understanding of male infertility.</p>
						<p>More text more text more text more text more text more text more text more text more text more text more text more text more text more text more text more text more text more text more text more text more text.</p>
						<p></p>
						<p>This app is a scientific endeavor hosted by <a href="http://igc.gulbenkian.pt/" target="_blank">Instituto Gulbenkian de Ciência</a>. It is currently developed by <a href="http://homes.soic.indiana.edu/rionbr/" target="_blank">Dr. Rion Brattig Correia</a>.</p>
					</div>
				</div>
			</div>


			<div class="row mb-4">
				<div class="col-2 offset-2 text-center">
					<figure class="figure mb-0">
						<a href="https://twitter.com/germcells" target="_blank" title="Paulo Navarro Costa">
							<img class="figure-img img-fluid rounded-circle" src="img/PauloNavarroCosta_160x160.jpg" alt="Paulo Navarro Costa">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="https://twitter.com/germcells" target="_blank" title="Paulo Navarro Costa">Dr. Paulo Navarro Costa, PI</a></figcaption>
				</div>
				<div class="col-2 offset-1 text-center">
					<figure class="figure mb-0">
						<a href="http://homes.soic.indiana.edu/rionbr/" target="_blank" title="Rion Brattig Correia">
							<img class="figure-img img-fluid rounded-circle" src="img/RionBrattigCorreia_160x160.jpg" alt="Rion Brattig Correia">
						</a>
					</figure>			
					<figcaption class="figure-caption"><a href="http://homes.soic.indiana.edu/rionbr/" target="_blank" title="Rion Brattig Correia">Dr. Rion Brattig Correia, Postdoc</a></figcaption>
				</div>
			</div>

			<div class="row mb-4">
				<div class="col-2 offset-2 text-center">
					<figure class="figure mb-0">
						<a href="https://gulbenkian.pt/ciencia/research/research-groups/plant-genomics/" target="_blank" title="Jörg Becker">
							<img class="figure-img img-fluid rounded-circle" src="img/JorgBecker_160x160.jpg" alt="Jörg Becker">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="http://www.informatics.indiana.edu/rocha" target="_blank" title="Jörg Becker">Dr. Jörg Becker, Co-PI</a></figcaption>
				</div>
				
				<div class="col-2 offset-1 text-center">
					<figure class="figure mb-0">
						<a href="http://www.informatics.indiana.edu/rocha" target="_blank" title="Luis M. Rocha">
							<img class="figure-img img-fluid rounded-circle" src="img/LuisMRocha_160x160.jpg" alt="Luis M. Rocha">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="http://www.informatics.indiana.edu/rocha" target="_blank" title="Luis M. Rocha">Prof. Luis M. Rocha, Co-PI</a></figcaption>
				</div>
			</div>

			<div class="row mb-4">
				<div class="col-3 text-center">
					<figure class="figure mb-0">
						<a href="https://gulbenkian.pt/ciencia/researcher/chandra-shekhar-misra/" title="Chandra Shekhar Misra">
							<img class="figure-img img-fluid rounded-circle" src="img/ChandraMisra_160x160.jpg" alt="Chandra Shekhar Misra">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="https://gulbenkian.pt/ciencia/researcher/chandra-shekhar-misra/" target="_blank" title="Jörg Becker">Msc. Chandra Shekhar Misra</a></figcaption>
				</div>
				<div class="col-3 text-center">
					<figure class="figure mb-0">
						<a href="#/" target="_blank" title="Joana Almeida">
							<img class="figure-img img-fluid rounded-circle" src="img/JoanaAlmeida_160x160.jpg" alt="Joana Almeida">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="#" title="Joana Almeida" target="_blank">Joana Almeida</a></figcaption>
				</div>
				<div class="col-3 text-center">
					<figure class="figure mb-0">
						<a href="#" target="_blank" title="Neide Silva">
							<img class="figure-img img-fluid rounded-circle" src="img/NeideSilva_160x160.jpg" alt="Neide Silva">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="#" title="Neide Silva" target="_blank">Neide Silva</a></figcaption>
				</div>
			</div>
		</div>

	</main>

	<?php include "include/footer.php"; ?>
	
	<?php include "include/javascripts.php"; ?>
	<script type="text/javascript" src="js/network-module.js"></script>
</body>
</html>