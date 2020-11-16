<?php
    // set server variables
    $scriptDir = '/Users/nick/Documents/GitHub/pathway_searching/webserver_scripts/';
    $pyshebang = '#!/Users/nick/Documents/GitHub/motif_builder/py35env/bin/python';
    // remove any keys without values, set sessionDir
    $params = array_filter($_POST);
    // copy script to new file containing $pyshebang
    $script = $scriptDir . json_decode($_POST["pyscript"]);
    $arr = file($script); // load in pyscript to replace shebang
    $arr[0] = $pyshebang . "\n"; // replace shebang with $pyshebang
    file_put_contents($script, implode($arr)); // write back to file
    chmod($script, 0755); // add execution permission
    // execute script
    $script = $script . ' ' . escapeshellarg(json_encode($params));
    $output = shell_exec($script);
    echo $output;
?>
