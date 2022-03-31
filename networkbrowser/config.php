<?php

/*
 * ERROR DEBUGING
 */
error_reporting(E_ALL);
ini_set('display_errors', 1);


/*
 * CONSTANT VARIABLES
 */
if ($_SERVER['SERVER_NAME'] == 'localhost') {
	define("PATH", "");
} else{
	define("PATH", "/meionav/");
}

?>