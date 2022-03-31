	<nav class="navbar navbar-expand-lg navbar-dark bg-dark py-3 d-print-block">
		<div class="container">
			<a class="navbar-brand" href="<?=PATH?>" title="Meiotic Navigator">
				<img src="<?=PATH?>img/logo-meiotic-navigator.svg" alt="Meiotic Navigator" width="30" height="24" class="d-inline-block align-text-top">
				Meiotic Navigator<sup style="color:lightgreen; font-style:italic;">beta</sup>
			</a>
			<button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbar" aria-controls="navbar" aria-expanded="false" aria-label="Toggle navigation">
				<span class="navbar-toggler-icon"></span>
			</button>
			<div class="collapse navbar-collapse d-print-block" id="navbar">
				<ul class="nav navbar-nav mr-auto">
					<?php $a = ($curpage == "meiotic-genes") ? "active" : ""; ?>
					<li class="nav-item <?=$a?> "><a class="nav-link" href="<?=PATH?>meiotic-genes.php">Meiotic genes</a></li>
					<?php $a = ($curpage == "publications") ? "active" : ""; ?>
					<li class="nav-item <?=$a?> "><a class="nav-link" href="<?=PATH?>publications.php">Publications</a></li>
				</ul>
		  
				<ul class="navbar-nav">
					<?php $a = ($curpage == "contacts") ? "active" : ""; ?>
					<li class="nav-item <?=$a?>"><a class="nav-link" href="<?=PATH?>contacts.php">Contacts</a></li>
					<?php $a = ($curpage == "privacy") ? "active" : ""; ?>
					<li class="nav-item <?=$a?>"><a class="nav-link" href="<?=PATH?>privacy.php">Privacy</a></li>
				</ul>
			</div>
		</div>
	</nav>