<?php
include_once("config.php");
$curpage = "privacy";
$pagetitle = "Privacy : Meiotic Navigator";
$pagedescription = "Privacy policy for The Meiotic Navigator.";
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
					<li  class="breadcrumb-item active">Privacy</li>
				</ol>
			</nav>
		</div>
	</header>

	<main role="main">
		<div class="container">
		
			<div class="py-5 mb-4">
				<div class="row">
					<div class="col-md-9 col-ld-7">
						<h1>Privacy</h1>
						<p>The Meiotic Navigator is subject to the privacy policy of the <a href="https://gulbenkian.pt/ciencia" target="_blank">Instituto Gulbenkian de Ciência</a> and the <a href="https://gulbenkian.pt/" target="_blank">Calouste Gulbenkian Foundation</a> website. The governing terms of use and privacy policy can be <a href="https://gulbenkian.pt/en/conditions-of-use/" target="_blank">found here</a>.</p>
					</div>
				</div>
			</div>
	</main>
	<?php include "include/footer.php"; ?>	
	<?php include "include/javascripts.php"; ?>
</body>
</html>