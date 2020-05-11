<?php
include_once("config.php");
$curpage = "Publications";
$pagetitle = "Publications : SpermNet";
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
						<p>Here is the list of our publications as well as publications and data sources that were used to generate or analyze the work provided here.</p>
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
								<td>2020</td>
								<td>...</td>
								<td>...</td>
							</tr>
							<tr class="">
								<td>2020</td>
								<td>...</td>
								<td>...</td>
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
								<td>YYYY</td>
								<td>...</td>
								<td>...</td>
							</tr>
							<tr>
								<td>YYYY</td>
								<td>...</td>
								<td>...</td>
							</tr>
						</tbody>
					</table>
				</div>
			</div>

		</div>


	</main>

	<?php include "include/footer.php"; ?>
	
	<?php include "include/javascripts.php"; ?>
	<script type="text/javascript" src="js/network-module.js"></script>
</body>
</html>