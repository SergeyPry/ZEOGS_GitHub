<?php
if ($_GET['randomId'] != "SXzO91NsmlLk0rBdeIuuZI_EpMJFn2J8yh1Vi9x6DKqMkaP6Y6K8Wg94pMEcAHka") {
    echo "Access Denied";
    exit();
}

// display the HTML code:
echo stripslashes($_POST['wproPreviewHTML']);

?>  
