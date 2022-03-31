<?php
include_once("config.php");
$curpage = "publications";
$pagetitle = "Publications : Meiotic Navigator";
$pagedescription = "Scientific publications involving The Meiotic Navigator.";
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
					<li class="breadcrumb-item"><a href="/andronet/index.php">Home</a></li>
					<li  class="breadcrumb-item active">Publications</li>
				</ol>
			</nav>
		</div>
	</header>

	<main role="main">
		<div class="container">

			<div class="py-5 mb-4">
				<div class="row">
					<div class="col-md-9 col-ld-7">
						<h1>Publications</h1>
						<p>List of publications and data sources that were used to generate the <a href="<?=PATH?>">Meiotic Navigator<sup style="color:darkgreen; font-style:italic;">beta</sup></a>.</p>
					</div>
				</div>
			</div>

			<div class="card mb-5">
				<div class="card-header">Our publications</div>
				<div class="card-body">
					<table class="table table-hover">
						<thead>
							<tr>
								<th scope="col">Year</th>
								<th scope="col">Description</th>
								<th scope="col">Reference</th>
							</tr>
						</thead>
						<tbody>
							<tr class="">
								<td>2022</td>
								<td>MeioNav</td>
								<td>R.B. Correia, J.M. Almeida, M. Wyrwoll, I. Julca, D. Sobral, C.S. Misra, L.G. Guilgur, H. Schuppe, N. Silva, P. Prudêncio, A. Nóvoa, A.S. Leocádio, J. Bom, M. Mallo, S. Kliesch, M. Mutwil, L.M. Rocha, F. Tüttelmann, J.D. Becker, P. Navarro-Costa [2022]. "<a href="https://www.biorxiv.org/content/10.1101/2022.03.02.482557v1" target="_blank">An old transcriptional program in male germ cells uncovers new causes of human infertility</a>". bioRxiv: 10.1101/2022.03.02.482557v1</td>
							</tr>
							<tr class="">
								<td>2021</td>
								<td>LowEnDe</td>
								<td>R.B. Correia, P. Navarro-Costa, and L.M. Rocha [2020]. "<a href="https://rionbr.github.io/publications/2020-COMPLEX_NETWORKS_2020_paper_148.pdf" target="_blank">Extraction of overlapping modules in networks via spectral methods and information theory</a>". Complex Networks 2020. The 9th International Workshop on Complex Networks and Their Applications. Dec. 1-3, 2020, Madrid, Spain (Online).</td>
							</tr>
						</tbody>
					</table>
					
				</div>	
				<div class="card-footer">
					<span class="text-danger">If you are using the data provided here, please consider citing our work.</span>
				</div>
			</div>


			<div class="card">
				<div class="card-header">Publications from other groups used in our work</div>
				<div class="card-body">
					<table class="table table-hover">
						<thead>
							<tr>
								<th scope="col">Year</th>
								<th scope="col">Description</th>
								<th scope="col">Reference</th>
							</tr>
						</thead>
						<tbody>
							<tr>
								<td>2021</td>
								<td>FlyBase</td>
								<td>Aoife Larkin, Steven J Marygold, Giulia Antonazzo, Helen Attrill, Gilberto dos Santos, Phani V Garapati, Joshua L Goodman, L Sian Gramates, Gillian Millburn, Victor B Strelets, Christopher J Tabone, Jim Thurmond, FlyBase Consortium. [2021]. <a href="https://doi.org/10.1093/nar/gkaa1026" target="_blank">FlyBase: updates to the Drosophila melanogaster knowledge base</a>. <em>Nucleic Acids Res</em>, 49(<strong>D1</strong>), D899–D907</a>. doi: 10.1093/nar/gkaa1026</td>
							</tr>
							<tr>
								<td>2019</td>
								<td>EggNOG</td>
								<td>Jaime Huerta-Cepas, Damian Szklarczyk, Davide Heller, Ana Hernández-Plaza, Sofia K Forslund, Helen Cook, Daniel R Mende, Ivica Letunic, Thomas Rattei, Lars J Jensen, Christian von Mering, Peer Bork. [2019]. <a href="https://doi.org/10.1093/nar/gky1085" target="_blank">eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses</a>. <em>Nucleic Acids Res</em>, 47(<strong>D1</strong>), D309–D314. doi: 10.1093/nar/gky1085</td>
							</tr>
							<tr>
								<td>2005</td>
								<td>RNAi reagents</td>
								<td>Lizabeth A Perkins, Laura Holderbaum, Rong Tao, Yanhui Hu, Richelle Sopko, Kim McCall, Donghui Yang-Zhou, Ian Flockhart, Richard Binari, Hye-Seok Shim, Audrey Miller, Amy Housden, Marianna Foos, Sakara Randkelv, Colleen Kelley, Pema Namgyal, Christians Villalta, Lu-Ping Liu, Xia Jiang, Qiao Huan-Huan, Xia Wang, Asao Fujiyama, Atsushi Toyoda, Kathleen Ayers, Allison Blum, Benjamin Czech, Ralph Neumuller, Dong Yan, Amanda Cavallaro, Karen Hibbard, Don Hall, Lynn Cooley, Gregory J Hannon, Ruth Lehmann, Annette Parks, Stephanie E Mohr, Ryu Ueda, Shu Kondo, Jian-Quan Ni, and Norbert Perrimon. [2015]. <a href="https://doi.org/10.1534/genetics.115.180208" target="_blank">The Transgenic RNAi Project at Harvard Medical School: Resources and Validation</a>. <em>Genetics</em>, 201(<strong>3</strong>) 843-52. doi: 10.1534/genetics.115.180208</td>
							</tr>
							<tr>
								<td>2018</td>
								<td>RNA-Seq in insect testis</td>
								<td>Viktor Vedelek, László Bodai, Gábor Grézal, Bence Kovács, Imre M. Boros, Barbara Laurinyecz & Rita Sinka. [2018]. <a href="">Analysis of <em>Drosophila melanogaster</em> testis transcriptome</a>. <em>BMC Genomics</em>, 19, 697. doi: 10.1186/s12864-018-5085-z</td>
							</tr>
							<tr>
								<td>2020</td>
								<td>Human infertility genes</td>
								<td>Brendan J Houston, Antoni Riera-Escamilla, Margot J Wyrwoll, Albert Salas-Huetos, Miguel J Xavier, Liina Nagirnaja, Corinna Friedrich, Don F Conrad, Kenneth I Aston, Csilla Krausz, Frank Tüttelmann, Moira K O’Bryan, Joris A Veltman, Manon S Oud [2020]. <a href="https://doi.org/10.1093/humupd/dmab030">A systematic review of the validated monogenic causes of human male infertility: 2020 update and a discussion of emerging gene–disease relationships</a>. <em>Human Reproduction Update</em>, dmab030, doi: 10.1093/humupd/dmab030</td>
							</tr>
						</tbody>
					</table>
				</div>
			</div>

		</div>


	</main>
	<?php include "include/footer.php"; ?>
	<?php include "include/javascripts.php"; ?>
</body>
</html>