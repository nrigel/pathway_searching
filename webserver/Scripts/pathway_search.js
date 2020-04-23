var relative_dir; //relative_dir on the server, needed for xhttp call
var session_id;


var pathway_data = {}; // pathway checkbox ids to nodes and edges
var pathway_checkboxes = {}; // doc ids for userfile to pathway checkbox ids
var cy_layout = "random"; // layout option for cytoscape graphing
var node_labels = "name"; // node attribute that is displayed on cytoscape graph
var edge_labels = "type"; // edge labels from Reactome

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

function check_change(flag) {
    if (flag === "motifs") {
       
        if (document.getElementById("motifs").checked) {
            console.log(document.getElementById("motifs").onchange);
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
    // pathway list check boxes
    if (flag.slice(0, 2) === "pl") {
        var checked = document.getElementById(flag).checked;
        if (flag in pathway_checkboxes) {
            var plist = pathway_checkboxes[flag];
            for (i = 0; i < plist.length; i++) { 
                if (checked) {
                    if (document.getElementById(plist[i]).checked === false) {
                        document.getElementById(plist[i]).checked = true;
                        check_change(plist[i]);
                    }
                } else {
                    if (document.getElementById(plist[i]).checked) {
                        document.getElementById(plist[i]).checked = false;
                        check_change(plist[i]);
                    }
                }
            }
        } else {
            if (checked) {
                visualize(flag, pathway_data[flag].nodes, pathway_data[flag].edges);
            } else {
                de_visualize(flag, pathway_data[flag].nodes);
            }            
        }
    }
}

function query_database(flag){
    var oReq = new XMLHttpRequest();
    var base_url = window.location.href;

    var params = new FormData(document.getElementById("process_load"));
    
    var userfile = document.getElementById("fileinput").value

    //document.getElementById("infor").innerHTML = "Load your spectrum from the server. Please be patient.";
    if (userfile) {

        oReq.open("POST", base_url.concat('Scripts/pathway_search.py'), true);
        oReq.setRequestHeader("Authorization", null);
        
        oReq.onload = function(oEvent) {
            //window.alert(oReq.responseText);
            if (oReq.status === 200) {
                //console.log(oReq.response);
                var result = JSON.parse(oReq.response);
                var pathway_ids = set_pathwaylist(userfile, result.pathway_list)
                for (i = 0; i < result.pathway_list.length; i++) { 
                    pathway_data[pathway_ids[i]] = result.pathway_data[result.pathway_list[i]];
                    var nodes = pathway_data[pathway_ids[i]].nodes;
                    var edges = pathway_data[pathway_ids[i]].edges;
                    for (j = 0; j < nodes.length; j++) {
                        nodes[j].data["js_pathway_id"] = pathway_ids[i];
                        // add specific pathway id into node ids
                        var node_id = nodes[j].data["id"].toString();
                        nodes[j].data["id"] = node_id.concat('_').concat(pathway_ids[i]);
                    }
                    for (j = 0; j < edges.length; j++) {
                        edges[j].data["js_pathway_id"] = pathway_ids[i];
                        // add specific pathway id into node ids
                        var source_id = edges[j].data['source'].toString();
                        var target_id = edges[j].data['target'].toString();
                        edges[j].data['source'] = source_id.concat('_').concat(pathway_ids[i]);
                        edges[j].data['target'] = target_id.concat('_').concat(pathway_ids[i]);
                    }
                }
        }   }
        oReq.send(params);
        window.alert("Searching for Pathways...");
        document.getElementById("fileinput").value = userfile
    } else {
        window.alert("No file selected...")
   }
    
}

function set_pathwaylist(userfile, pathway_list) {
    var ul = document.getElementById("pathway_list");
    var filename = userfile.split(/(\\|\/)/g).pop();
    if (document.getElementById("pl_".concat(filename)) == null) {
        var l = document.createElement("ol");
        l.setAttribute("id", "pl_".concat(filename));
        l.setAttribute("style", "list-style-type:none; padding-left: 15px;");

        var id = "plbox_".concat(filename)
        pathway_checkboxes[id] = [];
        var x = document.createElement("input");
        x.setAttribute("type", "checkbox");
        x.setAttribute("id", id);
        x.setAttribute("onclick", "check_change('".concat(id).concat("')"));
        var label = document.createElement("label");
        label.htmlFor = id;
        label.appendChild(document.createTextNode(filename));
        ul.appendChild(x);
        ul.appendChild(label);
        ul.appendChild(l);
        ul.appendChild(document.createElement("br"));
    }
    
    var pathway_ids = {};
    var l = document.getElementById("pl_".concat(filename));
    var i;
    for (i = 0; i < pathway_list.length; i++) { 
        var id = null;
        var c = 0;
        while (id == null) {
            var s = "plbox_".concat(filename).concat("_").concat(pathway_list[i])
            if (c) {
                s = s.concat("_").concat(c);
            }
            if (document.getElementById(s) == null) {
                var id = s;
            } else {
                c += 1;
            }
        }
        var x = document.createElement("input");
        x.setAttribute("type", "checkbox");
        x.setAttribute("id", id);
        x.setAttribute("onclick", "check_change('".concat(id).concat("')"));
        var label = document.createElement("label");
        label.htmlFor = id;
        label.appendChild(document.createTextNode(pathway_list[i]));
        l.appendChild(x);
        l.appendChild(label);
        l.appendChild(document.createElement("br"));
        pathway_ids[i] = id;
        pathway_checkboxes["plbox_".concat(filename)].push(id);
    }
    return pathway_ids
}

function clear_pathwaylist() {
    const ul = document.getElementById("pathway_list");
    ul.innerHTML = "";
}

function setup_cy(container_id) {
    var cy = window.cy = cytoscape({
                container: document.getElementById('cy'),

                style: [
                    {selector: 'node',
                        style: {'content': 'data('.concat(node_labels).concat(')')}},

                    {selector: 'edge',
                        style: {'label': 'data('.concat(edge_labels).concat(')'),
                            'curve-style': 'bezier',
                            'target-arrow-shape': 'triangle'
                        }}],
                elements: {nodes: [], edges: []}, layout: {name: cy_layout}});
}

function visualize(pathway, nodes, edges) {
    // this will overlap the nodes if the same node label already exists...
    cy.add(nodes);
    cy.add(edges);
    var layout = cy.layout({name: cy_layout});
    layout.run();
}

function de_visualize(pathway, nodes) {
    cy.remove(cy.elements('node[js_pathway_id = "'.concat(pathway).concat('"]')));
    cy.remove(cy.elements('edge[js_pathway_id = "'.concat(pathway).concat('"]')));
}