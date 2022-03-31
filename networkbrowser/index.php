<?php
include_once("config.php");
$curpage = "home";
$pagetitle = "Meiotic Navigator";
$pagedescription = "The Meiotic Navigator is a scientific research projet that leverages large scale computational data, wet-lab experiments, and comparative biology to advance our understanding of male infertility.";
?>
<!doctype html>
<html lang="en">
<?php include "include/head.php"; ?>
<body>
	<header>
		<?php include "include/navigation.php"; ?>
		<div class="container">
			<nav aria-label="breadcrumb">
				<ol class="breadcrumb">
					<li class="breadcrumb-item active">Home</li>
				</ol>
			</nav>
		</div>
	</header>

	<main role="main">
		<div class="container">
			<div class="py-5 mb-5">
				<h1 class="text-center">Welcome to the Meiotic Navigator</h1>
				<div class="row">
					<div class="col-md-10 col-ld-7">
						<img src="<?=PATH?>img/logo-meiotic-navigator.svg" alt="Meiotic Navigator" width="150" height="90" class="float-left">
						<p class="lead">The Meiotic Navigator is a comparative biology resource that leverages large scale computational data, new wet-lab experiments, and gene orthology analysis to advance our understanding of the evolutionarily-conserved basis of male meiosis.</p>
					</div>
				</div>
			</div>
		</div> <!-- /end Container -->

		<!-- Action link -->
		<div class="bg-info my-5 py-5">
			<div class="container text-center">
				<a class="btn btn-lg btn-outline-light my-5" href="<?=PATH?>meiotic-genes.php">Navigate the functional analysis of 920 conserved meiotic genes</a>
			</div>
		</div>
		
		<div class="container">
			<div class="pt-5 mb-5">
				<div class="row pt-4">
					<div class="offset-md-6">
						<p class="lead">The Meiotic Navigator was designed and built with <span class="oi oi-heart text-danger"></span> at the <a href="https://gulbenkian.pt/ciencia/" target="_blank">Instituto Gulbenkian de Ciência</a>, in Oeiras, Portugal.</p>
						<img src="<?=PATH?>img/gulbenkian.svg" alt="Instituto Gulbenkian de Ciência" width="260" height="90" class="float-right">
					</div>
				</div>
			</div>
		</div> <!-- /end Container -->
	</main>
	<?php include "include/footer.php"; ?>
	<?php include "include/javascripts.php"; ?>
</body>
</html>