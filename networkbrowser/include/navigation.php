	<nav class="navbar navbar-expand-lg navbar-dark bg-dark py-3">
		<div class="container">
			<a class="navbar-brand" href="/index" title="">SpermNet<sup style="color:lightgreen; font-style:italic;">beta</sup></a>
			<button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbar" aria-controls="navbar" aria-expanded="false" aria-label="Toggle navigation">
				<span class="navbar-toggler-icon"></span>
			</button>
			
			<div class="collapse navbar-collapse" id="navbar">
	
				<ul class="nav navbar-nav mr-auto">
					<li class="nav-item dropdown <?=$a?>">
						<a class="nav-link dropdown-toggle" href="" id="navbar-network-dropdown" data-toggle="dropdown">Interactive Networks</a>
						<ul class="dropdown-menu multi-level">
							<li><a class="dropdown-item" href="network-mlayer.php">Multi-layer network</a></li>
							<li class="dropdown-divider"></li>
							<li><a class="dropdown-item disabled " href="#">Single-layer networks</a></li>
							<li class="dropdown-submenu"><a class="dropdown-item" href="#">Homo sapiens</a>
								<ul class="dropdown-menu">
									<li class="dropdown-submenu"><a class="dropdown-item" href="#">SVD Modules</a>
										<ul class="dropdown-menu">
											<li><a class="dropdown-item" href="network-module.php?layer=HS&algorithm=SVD&module=1c">Mod. 1 - Uniquitination</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=HS&algorithm=SVD&module=2c">Mod. 2 - Translation + Splicing</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=HS&algorithm=SVD&module=3c">Mod. 3 - Splicing</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=HS&algorithm=SVD&module=4c">Mod. 4 - Cell Cycle</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=HS&algorithm=SVD&module=5c">Mod. 5 - Unnamed</a></li>
										</ul>
									</li>
								</ul>
							</li>
							<li class="dropdown-submenu"><a class="dropdown-item href="#">Mus musculus</a>
								<ul class="dropdown-menu">
									<li class="dropdown-submenu"><a class="dropdown-item" href="#">SVD Modules</a>
										<ul class="dropdown-menu">
											<li><a class="dropdown-item" href="network-module.php?layer=MM&algorithm=SVD&module=1c">Mod. 1 - Uniquitination</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=MM&algorithm=SVD&module=2c">Mod. 2 - Translation</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=MM&algorithm=SVD&module=3c">Mod. 3 - Splicing</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=MM&algorithm=SVD&module=4c">Mod. 4 - Cell Cycle</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=MM&algorithm=SVD&module=5c">Mod. 5 - Unnamed</a></li>
										</ul>
									</li>
								</ul>
							</li>
							<li class="dropdown-submenu"><a class="dropdown-item" href="#">Drosophila melanogaster</a>
								<ul class="dropdown-menu">
									<li class="dropdown-submenu"><a class="dropdown-item dropdown-toggle" href="#">SVD Modules</a>
										<ul class="dropdown-menu">
											<li><a class="dropdown-item" href="network-module.php?layer=DM&algorithm=SVD&module=1c">Mod. 1 - Translation</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=DM&algorithm=SVD&module=2c">Mod. 2 - Splicing</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=DM&algorithm=SVD&module=3c">Mod. 3 - Uniquitination + Cell cycle</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=DM&algorithm=SVD&module=4c">Mod. 4 - Post-meiotic development</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=DM&algorithm=SVD&module=5c">Mod. 5 - Unnamed</a></li>
											<li><a class="dropdown-item" href="network-module.php?layer=DM&algorithm=SVD&module=6c">Mod. 6 - Unnamed</a></li>
										</ul>
									</li>
								</ul>
							</li>
						</ul>
					</li>
					<li class="nav-item <?=$a?> "><a class="nav-link" href="publications.php">Publications</a></li>
				</ul>
		  
				<ul class="navbar-nav">
					<li class="nav-item <?=$a?>"><a class="nav-link" href="privacy.php">Privacy</a></li>
					<li class="nav-item <?=$a?>"><a class="nav-link" href="about.php">About</a></li>
				</ul>
		
			</div>
		</div>
	</nav>