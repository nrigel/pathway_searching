

//define base_url     
if (!window.location.origin) {
    window.location.origin = window.location.protocol + "//" + window.location.hostname + (window.location.port ? ':' + window.location.port : '');
}
var base_url = window.location.origin;


$(document).ready(function() {

    $('form#motif').submit(function(e) {
        var fd = new FormData(document.getElementById("motif"));
        var oReq = new XMLHttpRequest();
        var base_url = window.location.origin;


        oReq.open("POST", base_url.concat('/index.php/motif/run'), true);
        oReq.setRequestHeader("Authorization", null);
        oReq.onload = function(oEvent)
        {
            if (oReq.status === 200)
            {
                var result = JSON.parse(oReq.responseText);
                clear_list();
                for(var i=0;i<result.length;i++)
                {
                    add_item_to_list (result[i]); 
                }
                if(result.length===0)
                {
                    var ul = document.getElementById("result_list");
                    var div = document.createElement("li");
                    div.appendChild(document.createTextNode("No match"));
                    ul.appendChild(div);
                }
            }
        }
        oReq.send(fd);
        e.preventDefault();
    });

});

function clear_list()
{
    var ul = document.getElementById("result_list");
    while (ul.firstChild) {
        ul.removeChild(ul.firstChild);
    }
}

function add_item_to_list(x) {
    var ul = document.getElementById("result_list");
    var div = document.createElement("li");
   
    div.appendChild(document.createTextNode("Compound name: "));
    div.appendChild(document.createTextNode(x.name));
    div.appendChild(document.createTextNode(", matching rmsd:"));
    div.appendChild(document.createTextNode(x.rmsd));
    div.appendChild(document.createTextNode(", data source:"));

    var a = document.createElement('a');
    a.setAttribute('href', x.db_link);
    a.setAttribute('target',"_blank");
    a.appendChild(document.createTextNode(x.source));
    div.appendChild(a);
    div.appendChild(document.createElement("BR"));


    var newimg1 = document.createElement("img");
    newimg1.setAttribute("src", origin+"/example/"+x.msm_1st);
    var newimg2 = document.createElement("img");
    newimg2.setAttribute("src", origin+"/example/"+x.msm_2nd);
    
    div.appendChild(newimg1);
    div.appendChild(newimg2);
    

    ul.appendChild(div);
}

function load_example()
{
    document.getElementById("peaklist").value = "68.598 4.243\n63.141 3.575\n22.145 1.318";
}