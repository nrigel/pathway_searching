var relative_dir; //relative_dir on the server, needed for xhttp call
var session_id;

var for_debug;

var key_flag = 0;

//var neo4j = require('neo4j-driver');
var target;


document.addEventListener('keydown', function(event) {

    if (event.keyCode === 85 && key_flag === 1) {
        key_flag = 2;
    } else if (event.keyCode === 77 && key_flag === 2) {
        alert("Advance Mode is On");
        document.getElementById("advance").style.display = "block";
    } else if (event.keyCode === 83) {
        key_flag = 1;
    } else {
        key_flag = 0;
    }
});


function query_database(flag){
    var oReq = new XMLHttpRequest();
    var base_url = window.location.href;

    var params = new FormData(document.getElementById("process_load"));
    
    var userfile = document.getElementById("fileinput").value

    //document.getElementById("infor").innerHTML = "Load your spectrum from the server. Please be patient.";
    if (document.getElementById("fileinput").value) {

        oReq.open("POST", base_url.concat('Scripts/pathway_search.py'), true);
        oReq.setRequestHeader("Authorization", null);
        
        oReq.onload = function(oEvent) {
            //result = JSON.parse(oReq.responseText.replace(/\bNaN\b/g, "null"));
            window.alert(oReq.responseText);
            if (oReq.status === 200) {
                //alert(oReq.responseText);
                //try{
                result = JSON.parse(oReq.responseText.replace(/\bNaN\b/g, "null"));
            }
        }
        oReq.send(params);
        window.alert(userfile)
        
   } else {
    window.alert('No file selected...')
   }
    
}



function check_change(flag) {
    if (flag === "motifs") {
        if (document.getElementById("motifs").checked) {
            document.getElementById("motif_0").checked = true;
            document.getElementById("motif_1").checked = true;
            document.getElementById("motif_2").checked = true;
        } else {
            document.getElementById("motif_0").checked = false;
            document.getElementById("motif_1").checked = false;
            document.getElementById("motif_2").checked = false;
        }
    }

    if (flag === "submotifs") {
        if (document.getElementById("submotifs").checked) {
            document.getElementById("submotif_1").checked = true;
            document.getElementById("submotif_2").checked = true;
            document.getElementById("submotif_3").checked = true;
            document.getElementById("submotif_4").checked = true;
        } else {
            document.getElementById("submotif_1").checked = false;
            document.getElementById("submotif_2").checked = false;
            document.getElementById("submotif_3").checked = false;
            document.getElementById("submotif_4").checked = false;
        }
    }
}
    
function applyfilters() {
    target = [];
}