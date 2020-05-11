<?php
include_once("config.php");
$curpage = "Privacy";
$pagetitle = "Privacy : SpermNet";
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
						<p>text goes here.</p>
					</div>
				</div>
			</div>
	</main>

	<?php include "include/footer.php"; ?>
	
	<?php include "include/javascripts.php"; ?>
	<script type="text/javascript" src="js/network-module.js"></script>
</body>
</html>