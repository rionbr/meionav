<?php
include_once("config.php");
$curpage = "contacts";
$pagetitle = "Contacts : Meiotic Navigator";
$pagedescription = "The scientists behind the Meiotic Navigator.";
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
					<li class="breadcrumb-item"><a href="index.php">Home</a></li>
					<li  class="breadcrumb-item active">Contacts</li>
				</ol>
			</nav>
		</div>
	</header>

	<main role="main">
		<div class="container">

			<div class="py-5 mb-4">
				<div class="row">
					<div class="col-md-9 col-ld-7">
						<h1>Contacts</h1>
					</div>
					<div class="col-8">
						<p>
							We are an interdisciplinary research team working to solve complex issues in human health using biological and computational tools.
							The <a href="<?=PATH?>">Meiotic Navigator<sup style="color:darkgreen; font-style:italic;">beta</sup></a> is hosted by the <a href="https://gulbenkian.pt/ciencia" target="_blank">Instituto Gulbenkian de Ciência</a>.
							It is currently developed by <a href="https://rionbr.github.io/" target="_blank">Dr. Rion Brattig Correia</a> & <a href="https://twitter.com/germcells" target="_blank">Dr. Paulo Navarro-Costa</a>. Please contact us for any questions.</p>
					</div>
				</div>
			</div>


			<div class="row mb-4">
				<div class="col-4 text-center">
					<!-- Paulo Navarro-Costa -->
					<figure class="figure mb-0">
						<a href="https://twitter.com/germcells" target="_blank" title="Paulo Navarro-Costa">
							<img class="figure-img img-fluid rounded-circle" src="img/PauloNavarroCosta_160x160.jpg" alt="Paulo Navarro-Costa">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="https://twitter.com/germcells" target="_blank" title="Paulo Navarro-Costa">Dr. Paulo Navarro-Costa</a></figcaption>
				</div>
				<div class="col-4 text-center">
					<!-- Rion Brattig Correia -->
					<figure class="figure mb-0">
						<a href="https://rionbr.github.io/" target="_blank" title="Rion Brattig Correia">
							<img class="figure-img img-fluid rounded-circle" src="img/RionBrattigCorreia_160x160.jpg" alt="Rion Brattig Correia">
						</a>
					</figure>			
					<figcaption class="figure-caption"><a href="https://rionbr.github.io/" target="_blank" title="Rion Brattig Correia">Dr. Rion Brattig Correia</a></figcaption>
				</div>
			</div>

			<div class="row mb-4">
				<div class="col-4 text-center">
					<!-- Luis M. Rocha -->
					<figure class="figure mb-0">
						<a href="http://www.informatics.indiana.edu/rocha" target="_blank" title="Luis M. Rocha">
							<img class="figure-img img-fluid rounded-circle" src="img/LuisMRocha_160x160.jpg" alt="Luis M. Rocha">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="http://www.informatics.indiana.edu/rocha" target="_blank" title="Luis M. Rocha">Prof. Luis M. Rocha</a></figcaption>
				</div>
				<div class="col-4 text-center">
					<!-- Jörg Becker -->
					<figure class="figure mb-0">
						<a href="https://gulbenkian.pt/ciencia/research/research-groups/plant-genomics/" target="_blank" title="Jörg Becker">
							<img class="figure-img img-fluid rounded-circle" src="img/JorgBecker_160x160.jpg" alt="Jörg Becker">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="https://gulbenkian.pt/ciencia/research/research-groups/plant-genomics/" target="_blank" title="Jörg Becker">Dr. Jörg Becker</a></figcaption>
				</div>
			</div>

			<div class="row mb-4">
				<div class="col-4 text-center">
					<!-- Joana Almeida -->
					<figure class="figure mb-0">
						<a href="https://twitter.com/almeidajoana2" target="_blank" title="Joana Almeida">
							<img class="figure-img img-fluid rounded-circle" src="img/JoanaAlmeida_160x160.jpg" alt="Joana Almeida">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="https://twitter.com/almeidajoana2" target="_blank" title="Joana Almeida">MSc. Joana Almeida</a></figcaption>
				</div>
				<div class="col-4 text-center">
					<!-- Chandra Shekhar Misra -->
					<figure class="figure mb-0">
						<a href="https://www.linkedin.com/in/chandra-shekhar-misra-a9245217/" target="_blank" title="Chandra Shekhar Misra">
							<img class="figure-img img-fluid rounded-circle" src="img/ChandraMisra_160x160.jpg" alt="Chandra Shekhar Misra">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="https://www.linkedin.com/in/chandra-shekhar-misra-a9245217/" target="_blank" title="Jörg Becker">Dr. Chandra Shekhar Misra</a></figcaption>
				</div>
			</div>
			<div class="row mb-4">
				<div class="col-4 text-center">
					<!-- Gastón Guilgur -->
					<figure class="figure mb-0">
						<a href="#" target="_blank" title="Gastón Guilgur">
							<img class="figure-img img-fluid rounded-circle" src="img/GastonGuilgur_160x160.jpg" alt="Gastón Guilgur">
						</a>
					</figure>
					<figcaption class="figure-caption"><a href="https://www.linkedin.com/in/leonardo-gast%C3%B3n-guilgur-b0661960" title="Gastón Guilgur" target="_blank">Dr. Gastón Guilgur</a></figcaption>
				</div>
				<div class="col-4 text-center">
					<!-- Neide Silva -->
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
</body>
</html>