var relative_dir; //relative_dir on the server, needed for xhttp call
var session_id;

var structure_labels = {"metabolites": "Metabolites", 
                    "motif_0": "Motifs (0th Shell)", "motif_1": "Motifs (1st Shell)", "motif_2": "Motifs (2nd Shell)",
                    "submotif_1": "Sub-Motifs (1st Shell)", "submotif_2": "Sub-Motifs (2nd Shell)", "submotif_3": "Sub-Motifs (3rd Shell)", "submotif_4": "Sub-Motifs (4th Shell)"};
var clique_checkboxes = {};
var pathway_metabolites = new Set(); // pathway doc ids to list of metabolites
var pathway_data = {}; // pathway checkbox ids to nodes and edges
var metabolites_to_pathways = {}; // metabolite names to list of pathway ids
var pathway_checkboxes = {}; // doc ids for userfile to pathway checkbox ids
var pathwayid_pathwayname = {}; // all pathway ids to the displayname
var pathway_lists = {}; // pathway list id to child nodes

var cy_layout = "cose"; // layout option for cytoscape graphing
var node_labels = "displayName"; // node attribute that is displayed on cytoscape graph
var edge_labels = null;//"reaction"; // edge labels from Reactome
var connect_pathways = true; // connect shared metabolites when multiple pathways are drawn

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
    let state;
    if (flag === "motifs") {
        if (document.getElementById("motifs").checked) { state = true; } 
        else { state = false; }

        for (let id of ["motif_0", "motif_1", "motif_2"]) {
            document.getElementById(id).checked = state;
        }
    }

    if (flag === "submotifs") {
        if (document.getElementById("submotifs").checked) { state = true; } 
        else { state = false; }

        for (let id of ["submotif_1", "submotif_2", "submotif_3", "submotif_4"]) {
            document.getElementById(id).checked = state;
        }
    }
    // pathway list check boxes
    if (flag.slice(0, 2) === "pl") {
        var checked = document.getElementById(flag).checked;
        if (checked) {
            visualize(flag, pathway_data[flag].nodes, pathway_data[flag].edges);
        } else {
            de_visualize(flag);      
        }
    }
    // change color of pathways to indicate metabolite hit
    if (flag == "metabolite_filter") {
        let sel = document.getElementById(flag);
        let selections = [];
        for (let opt of sel.options) {
            if (opt.selected === true) {
                selections.push(opt.text);
            }
        }

        for (let node of cy.elements('node')) {
            node.json({ selected: false });
        }

        for (let pathway in pathway_data) {
            if (pathway == null) { continue; }
            document.getElementById(pathway.concat("_label")).style.color = "black";
        }        

        if (selections.includes("All Metabolites") == 0 && selections.length > 0) {
            // filter pathways
            for (let metabolite of selections) {
                for (let node of cy.elements('node[displayName = "'.concat(metabolite).concat('"]'))) {
                        node.json({ selected: true }); // select the node to highlight it in the cy plot
                    }
                for (let pathway of metabolites_to_pathways[metabolite]) {
                    document.getElementById(pathway.concat("_label")).style.color = "magenta";
                } 
            }
        }        
    }
}

function query_database(flag){

    if (flag == 'preliminary') {
        document.getElementById("search_opt").setAttribute("value", "true");
    } else {
        document.getElementById("search_opt").setAttribute("value", "false");
    }

    var oReq = new XMLHttpRequest();
    var base_url = window.location.href;

    var params = new FormData(document.getElementById("process_load"));
    var userfile = document.getElementById("fileinput").value

    if (userfile) {
        document.getElementById("loader").style.display = "block";
        oReq.open("POST", base_url.concat('Scripts/pathway_search.py'), true);
        oReq.setRequestHeader("Authorization", null);
        
        oReq.onload = function(oEvent) {
            if (oReq.status === 200) {
                //console.log(typeof oReq.response);
                try {
                    var result = JSON.parse(oReq.response);
                }
                catch(err) {
                    console.log(oReq.response);
                }
                
                let structure_options = result.structure_options;

                clique_checkboxes = {}
                
                let ptable = document.getElementById("clique_table");

                let clique_cutoff = parseInt(document.getElementById("clique_cutoff").value, 10);
                let pcount = document.getElementById("number_of_pathways");
                let cliques = {};
                
                for (let m of structure_options) {
                    cliques[m] = {};
                    for (i = 0; i < result['pathway_cliques'][m]['keys'].length; i++) { 
                        let clique = result['pathway_cliques'][m]['keys'][i].slice(1).join(', ');
                        cliques[m][clique] = result['pathway_cliques'][m]['values'][i];
                    }
                }

                for (let m of structure_options) {
                    let header = ['<tbody><tr>',
                            '<th style="width: 5%; word-wrap: break-all; border: none;">Select</th>',
                            '<th style="width: 76.5%; word-wrap: break-all; border: none;">'.concat(structure_labels[m]).concat(' Clique</th>'),
                            '<th style="width: 18.5%; word-wrap: break-all; border: none;">Pathway Count</th>',
                            '</tr>',
                            '</tbody>']
                    ptable.innerHTML = header.join('\n');

                    for (let [clique, count] of Object.entries(cliques[m])) {
                        let row = ptable.insertRow(-1);
                        let n = ptable.getElementsByTagName('td').length;
                        
                        let cell0 = row.insertCell(0);
                        let cbox = document.createElement("input");
                        cbox.setAttribute("type", "checkbox");
                        cbox.setAttribute("name", "cbox".concat(n.toString()));
                        cbox.setAttribute("value", m.concat(', ').concat(clique))
                        cbox.checked = true;
                        if (clique_cutoff < count) {
                            cbox.checked = false;
                        }
                        cell0.appendChild(cbox);
                        cell0.setAttribute("value", "cbox");
                        cell0.setAttribute("valign", "middle");
                        cell0.setAttribute("style", "border: none;");

                        let cell1 = row.insertCell(1);
                        cell1.appendChild(document.createTextNode(clique));
                        cell1.setAttribute("value", "clique");
                        cell1.setAttribute("valign", "middle");
                        cell1.setAttribute("style", "border: none;");
                        
                        let cell2 = row.insertCell(2);
                        cell2.appendChild(document.createTextNode(count.toString()));
                        cell2.setAttribute("value", "count");
                        cell2.setAttribute("style", "border: none;");
                        cell2.setAttribute("valign", "middle");
                    }
                }
                for (let col of ptable.getElementsByTagName('td')) {
                    if (col.value == "cbox") {
                        col.setAttribute("style", "width: 5%; word-wrap: break-all; border: none;");
                    }
                    if (col.value == "clique") {
                        col.setAttribute("style", "width: 76.5%; word-wrap: break-all; border: none;");
                    }
                    if (col.value == "count") {
                        col.setAttribute("style", "width: 18.5%; word-wrap: break-all; border: none;");
                    }
                }
            
                if (flag != 'preliminary') {
                    var pathway_ids = set_pathwaylist(userfile, result.pathway_list, structure_options);

                    for (let m of structure_options) {
                        if (!(m in result.pathway_list)) {
                            continue;
                        }

                        for (i = 0; i < result.pathway_list[m].length; i++) { 
                            if (pathway_ids[m][i] === null) {
                                continue;
                            }
                            pathwayid_pathwayname[pathway_ids[m][i]] = result.pathway_list[m][i];
                            pathway_data[pathway_ids[m][i]] = result.pathway_data[m][result.pathway_list[m][i]];
                            var nodes = pathway_data[pathway_ids[m][i]].nodes;
                            var edges = pathway_data[pathway_ids[m][i]].edges;

                            for (j = 0; j < nodes.length; j++) {
                                // add metabolite to metabolite pathway dict
                                if (metabolites_to_pathways[nodes[j].data["displayName"]] == null){
                                    metabolites_to_pathways[nodes[j].data["displayName"]] = [];
                                }
                                metabolites_to_pathways[nodes[j].data["displayName"]].push(pathway_ids[m][i]);
                                // add tag to node id so that we can keep searches separate
                                nodes[j].data["js_pathway_id"] = pathway_ids[m][i];
                                // add metabolite name to our selection list
                                pathway_metabolites.add(nodes[j].data["displayName"]);   
                                // add specific pathway id into node ids
                                var node_id = nodes[j].data["id"].toString();
                                nodes[j].data["id"] = node_id.concat('_').concat(pathway_ids[m][i]);
                                nodes[j].data["parent"] = pathway_ids[m][i];
                                nodes[j].data["class"] = "Metabolite";                            
                            }

                            for (j = 0; j < edges.length; j++) {
                                edges[j].data["js_pathway_id"] = pathway_ids[m][i];
                                // add specific pathway id into node ids
                                var source_id = edges[j].data['source'].toString();
                                var target_id = edges[j].data['target'].toString();
                                edges[j].data['source'] = source_id.concat('_').concat(pathway_ids[m][i]);
                                edges[j].data['target'] = target_id.concat('_').concat(pathway_ids[m][i]);
                            }
                        }
                    }

                    var met_menu = document.getElementById("metabolite_filter");
                    
                    for (i = met_menu.options.length - 1; i >= 1; i--) {
                        met_menu.remove(i);
                    }

                    var metabolitelist = Array.from(pathway_metabolites)
                    metabolitelist.sort(function (a, b) {return a.toLowerCase().localeCompare(b.toLowerCase());});
                    for (let item of metabolitelist) {
                        let option = document.createElement("option");
                        option.text = item; // metabolite name
                        met_menu.add(option);
                    }
                } 
            }
        
            document.getElementById("loader").style.display = "none";
        }
        oReq.send(params);
    } else {
        window.alert("No file selected...")
   } 
}

function set_pathwaylist(userfile, pathway_list, structure_options) {
    var ul = document.getElementById("pathway_list");

    var filename = userfile.split(/(\\|\/)/g).pop();
    var id = "pl_".concat(filename)
    if (document.getElementById(id) == null) {
        pathway_checkboxes[id] = {'Metabolites': [], 'Motifs': [], 'Sub-Motifs': []};
        
        let label = document.createElement("label");
        label.setAttribute("id", id);
        label.appendChild(document.createTextNode(filename));
        label.style = "font-weight:bold"
        ul.appendChild(label);
        ul.appendChild(document.createElement("br"));

        pathway_lists[id] = {};
        pathway_checkboxes[id] = {};
    } 

    let font_colors = {"met": "blue", "mot": "darkgreen", "sub": "red"}
    let cat_names = {"metabolites": "Metabolites", 
                    "motif_0": "Motifs (0th Shell)", "motif_1": "Motifs (1st Shell)", "motif_2": "Motifs (2nd Shell)",
                    "submotif_1": "Sub-Motifs (1st Shell)", "submotif_2": "Sub-Motifs (2nd Shell)", "submotif_3": "Sub-Motifs (3rd Shell)", "submotif_4": "Sub-Motifs (4th Shell)"}
    var pathway_ids = {};
    var i;
    for (let m of structure_options) {
        let m_id = "pl".concat(m).concat("_").concat(filename)

        if (!(m in pathway_lists[id])) {
            let label = document.createElement("label");
            label.appendChild(document.createTextNode(cat_names[m]));
            label.setAttribute("style", "padding-left: 15px; color: ".concat(font_colors[m.slice(0, 3)]).concat(";"));
            ul.appendChild(label);
            
            let m_div = document.createElement("div");
            m_div.setAttribute("id", m_id.concat("_div"));
            m_div.setAttribute("style", "padding-left: 25px;");
            ul.appendChild(m_div);
            ul.appendChild(document.createElement("br"));
            
            pathway_lists[id][m] = {};
            pathway_checkboxes[id][m] = [];
        }
        pathway_ids[m] = {};
        let l = document.getElementById(m_id.concat("_div"));
        for (i = 0; i < pathway_list[m].length; i++) { 
            let pid = m_id.concat("_").concat(pathway_list[m][i])
            if (!(pid in pathway_lists[id][m])) {  
                var x = document.createElement("input");
                x.setAttribute("type", "checkbox");
                x.setAttribute("id", pid);
                x.setAttribute("onclick", "check_change('".concat(pid).concat("')"));
                
                var label = document.createElement("label");
                label.htmlFor = pid;
                label.setAttribute("id", pid.concat("_label"));
                label.appendChild(document.createTextNode(pathway_list[m][i]));

                pathway_lists[id][m][pid] = [x, label];

                l.appendChild(x);
                l.appendChild(label);
                l.appendChild(document.createElement("br"));
                
                pathway_ids[m][i] = pid;
                pathway_checkboxes[id][m].push(pid);
            } else {
                pathway_ids[m][i] = null;
            }
        }
    }
    return pathway_ids
}

function clear_pathwaylist() {
    let ul = document.getElementById("pathway_list");
    //ul.innerHTML = "";

    while (ul.firstChild) {
        //console.log(ul.lastChild)
        ul.removeChild(ul.lastChild);
    }

    cy.remove(cy.elements('node'));
    cy.remove(cy.elements('edge'));

    // clear metabolite filter engine
    let filter = document.getElementById("metabolite_filter");
    $("#metabolite_filter").empty();
    let option = document.createElement("option");
    option.text = "All Metabolites";
    option.value = "metfilt_AllMetabolites";
    filter.add(option);
}

function set_nodesize(flag) {
    
    function sizeset(opt, fc1, fc2) {
        // opt either equals a color or the word "reset"
        if (!isNaN(fc1) && !isNaN(fc2)) {
            for (let node of cy.elements("node")) {
                
                if (opt = "set") {

                    if (parseFloat(fc1) > parseFloat(node.data()["size"])) {
                        node.data("userSize", parseFloat(fc1))
                    }

                    if (parseFloat(fc2) < parseFloat(node.data()["size"])) {
                        node.data("userSize", parseFloat(fc2))
                    }
                }

                if (opt = "reset") {
                    node.data("userSize", undefined);
                } 
            }
        }
    }

    if (flag == "reset") {
        sizeset("reset", 0, 0);
        document.getElementById("ns_min").value = "min";
        document.getElementById("ns_max").value = "max";
    } else {
        let fc1 = document.getElementById("ns_min").value;
        if (fc1 == "min") {
            fc1 = 0;
        }
        let fc2 = document.getElementById("ns_max").value;
        sizeset("set", fc1, fc2);
    }

}

function set_nodecolor(flag) {

    function colorset(opt, color, fc1, fc2) {
        // opt either equals a color or the word "reset"
        if (!isNaN(fc1) && !isNaN(fc2)) {
            for (let node of cy.elements("node")) {
                if (node.data()["Fold Change"]){
                    if (parseFloat(fc1) <= parseFloat(node.data()["Fold Change"]) && parseFloat(node.data()["Fold Change"]) <= parseFloat(fc2)) {
                        if (opt == "reset") {
                            if (node.data()["userColor"] == color) {
                                node.data("userColor", node.data()["color"]);
                            }
                        } else {
                            node.data("userColor", color);
                        }
                    }
                }
            }
        }
    }

    if (flag == "check") {
        // called because nodes were added to plot
        let l = document.getElementById("usercolorlist");
        for (let node of l.childNodes) {
            let id = node.id.split(" ");
            let fc1 = parseFloat(id[0]);
            let fc2 = parseFloat(id[1]);
            let color = id[2];
            colorset("set", color, fc1, fc2);
        }
    } else {

        let e = document.getElementById("fc_setcolor");
        let color = e.options[e.selectedIndex].value;
        let colorname = e.options[e.selectedIndex].text;
        let fc1 = document.getElementById("fc_1").value;
        let fc2 = document.getElementById("fc_2").value;
 
        colorset("set", color, fc1, fc2);
        
        let l = document.getElementById("usercolorlist");
        let node = document.createElement("LI");
        let text = fc1.concat(" - ").concat(fc2).concat(" (").concat(colorname).concat(") ");
        let id = fc1.concat(" ").concat(fc2).concat(" ").concat(color)
        node.setAttribute("id", id);
        let textnode = document.createTextNode(text);
        node.appendChild(textnode); 
        let button = document.createElement("button");
        button.innerHTML = "remove";
        button.onclick = function(){ 
            node.remove();
            colorset("reset", color, fc1, fc2);
        };
        node.appendChild(button);
        l.appendChild(node);
        document.getElementById("fc_1").value = "min";
        document.getElementById("fc_2").value = "max";
    }
}

function setup_cy(container_id) {
    // colors = [#E50002, #E02C00, #DC5900, #D88500, #D4AF00, #C8CF00, #99CB00, #6CC700, #41C300, #18BF00], #0074BF
    var cy = window.cy = cytoscape({
                container: document.getElementById('cy'),

                style: [
                    {selector: 'node',
                        style: {'content': 'data('.concat(node_labels).concat(')'),
                                'background-color': function(ele) {
                                    if (ele.data()['userColor']) {
                                        return ele.data()['userColor']
                                    } else {
                                        return ele.data()['color'];
                                    }},
                                'shape': function(ele) {return ele.data()['shape']},
                                'height': function(ele) {
                                    if (ele.data()['userSize']) {
                                        return ele.data()['userSize']
                                    } else {
                                        return ele.data()['size']}
                                    },
                                'width': function(ele) {
                                    if (ele.data()['userSize']) {
                                        return ele.data()['userSize']
                                    } else {
                                        return ele.data()['size']}
                                    },
                                'font-size': '14px',
                                "text-wrap": "wrap",
                                "text-max-width": 75
                            }},
                    
                    {selector: ':parent',
                        style: {'content': function(ele) {return ele.data()['pathwayName'];},
                                'font-size': '24px',
                                'text-valign': 'top',
                                'text-halign': 'center',
                                "text-max-width": 1000,
                                'background-opacity': 0.1,
                                'background-color': '#FF00FF'
                        }
                    },

                    {selector: ':selected',
                        style: {'background-color': 'yellow', 
                                'line-color': 'yellow', 
                                'target-arrow-color': 'black', 
                                'source-arrow-color': 'black',}},
                    
                    {selector: 'edge',
                        style: {'label': function(ele) {
                            if (edge_labels) {
                                return 'data('.concat(edge_labels).concat(')')
                            } else { return '' }},
                            'curve-style': 'bezier',
                            'target-arrow-shape': 'triangle',
                            'font-size': '18px'}}
                ],
                elements: {nodes: [], edges: []}, 
                layout: {name: cy_layout,
                        nodeDimensionsIncludeLabels: true}});

    cy.on('tapend', function() {
        for (let ele of cy.$(':selected')) {            
            if (ele.data()['type'] != 'Pathway') {
                ele.qtip({
                    content: function(){ 
                        if (ele.isNode()) {
                            if (ele.data()['type'] == 'Metabolite') {
                                let svg = ele.data()['SVG']; //.split('\n'); // most molecules have SVG stored
                                //console.log(svg)
                                var stats = ele.data()[node_labels].concat('<br/>').concat(svg);                 
                                for (let key of ['Metabolite Matches', 'Fold Change', 'P-Value', 'databaseName', 'name', 'displayName', 'CHEBI_ID', 'COLMARm', 'HMDB_ID', 'SMILES_2D', 'SMILES_3D', 'formula', 'labels', 'motif_0', 'motif_1', 'motif_2', 'submotif_1', 'submotif_2', 'submotif_3', 'submotif_4', 'spinsystems', 'nodes', 'url']) {
                                    if (ele.data()[key]) {       
                                        if (key == 'url') {
                                            var value = key.concat(': <a href="').concat(ele.data()[key]).concat('" target="_blank">').concat(ele.data()[key]).concat('</a>')
                                        } else {
                                            var value = key.concat(': ').concat(ele.data()[key].toString())
                                        }
                                    } else {
                                        var value = key.concat(': ')
                                    }
                                    var stats = stats.concat('<br/><br/>').concat(value)
                                }
                                return stats;

                            } else { return ''; }
                        } else if (ele.isEdge()) {
                            return ele.data()['reaction']['displayName']
                        }   
                    },
                    position: {
                        my: 'top center',
                        at: 'bottom center'
                    },
                    style: {
                        classes: 'qtip-bootstrap',
                        tip: {
                            width: 8,
                            height: 8
                        },
                        'max-width': '200px',
                        'max-height': '275px',
                    }
                });
            }          
        }
    })
}

function visualize(pathway, nodes, edges) {
    // this will overlap the nodes if the same node label already exists...
    cy.add([{'data': {'id': pathway, 'displayName': '', 'pathwayName': pathwayid_pathwayname[pathway], 'type': 'Pathway', 'size': 0, 'color': 'red', 'shape': 'ellipse', 'js_pathway_id': pathway}, 'selectable': false}]);
    cy.add(nodes);
    cy.add(edges);

    if (connect_pathways) {
        let disp_nodes = cy.elements('node');
        let pathway_connectors = [];
        var i;
        for (i = 0; i < disp_nodes.length; i++) {
            for (j = i+1; j < disp_nodes.length; j++) {
                if (disp_nodes[i].data()['displayName'] === disp_nodes[j].data()['displayName']) {
                    let edge = {'data': {}};
                    edge.data['source'] = disp_nodes[i].data()['id'];
                    edge.data['target'] = disp_nodes[j].data()['id'];
                    pathway_connectors.push(edge);
                }
            }
        }
        cy.add(pathway_connectors);
    }
    
    let options = {
      name: 'cose', //cose works best so far; cose-bilkent could be worth a try but requires more code

      // Called on `layoutready`
      ready: function(){},

      // Called on `layoutstop`
      stop: function(){},

      // Whether to animate while running the layout
      // true : Animate continuously as the layout is running
      // false : Just show the end result
      // 'end' : Animate with the end result, from the initial positions to the end positions
      animate: true,

      // Easing of the animation for animate:'end'
      animationEasing: undefined,

      // The duration of the animation for animate:'end'
      animationDuration: undefined,

      // A function that determines whether the node should be animated
      // All nodes animated by default on animate enabled
      // Non-animated nodes are positioned immediately when the layout starts
      animateFilter: function ( node, i ){ return true; },

      // The layout animates only after this many milliseconds for animate:true
      // (prevents flashing on fast runs)
      animationThreshold: 250,

      // Number of iterations between consecutive screen positions update
      refresh: 20,

      // Whether to fit the network view after when done
      fit: true,

      // Padding on fit
      padding: 20, //30

      // Constrain layout bounds; { x1, y1, x2, y2 } or { x1, y1, w, h }
      boundingBox: undefined,

      // Excludes the label when calculating node bounding boxes for the layout algorithm
      nodeDimensionsIncludeLabels: true,

      // Randomize the initial positions of the nodes (true) or use existing positions (false)
      randomize: false,

      // Extra spacing between components in non-compound graphs
      componentSpacing: 40, //40

      // Node repulsion (non overlapping) multiplier
      nodeRepulsion: function( node ){ return 2048; },

      // Node repulsion (overlapping) multiplier
      nodeOverlap: 1, //4

      // Ideal edge (non nested) length
      idealEdgeLength: function( edge ){ return 32; }, //32

      // Divisor to compute edge forces
      edgeElasticity: function( edge ){ return 32; }, //32

      // Nesting factor (multiplier) to compute ideal edge length for nested edges
      nestingFactor: 1.2,

      // Gravity force (constant)
      gravity: 1,

      // Maximum number of iterations to perform
      numIter: 1000,

      // Initial temperature (maximum node displacement)
      initialTemp: 1000,

      // Cooling factor (how the temperature is reduced between consecutive iterations
      coolingFactor: 0.99,

      // Lower temperature threshold (below this point the layout will end)
      minTemp: 1.0
    };

    var layout = cy.layout(options);
    layout.run();
    check_change("metabolite_filter"); // check for metabolite selections
    set_nodecolor("check"); // change node colors if any user color filters have been applied
}

function de_visualize(pathway) {
    cy.remove(cy.elements('node[js_pathway_id = "'.concat(pathway).concat('"]')));
    cy.remove(cy.elements('edge[js_pathway_id = "'.concat(pathway).concat('"]')));
}

function saveResults(flag) {
    // save results using options from save box
    if (flag == 'SVG') {
        var textToWrite = cy.svg({scale: 1, full: false});
        var textFileAsBlob = new Blob([textToWrite], { type: 'image/svg+xml;charset=utf-8' });
        var fileNameToSaveAs = document.getElementById("svg_name").value;
    } else if (flag == 'PNG') {
        var textToWrite = cy.png({'output': 'blob', scale: 1, full: false});
        var textFileAsBlob = new Blob([textToWrite], { type: 'image/png' });
        var fileNameToSaveAs = document.getElementById("png_name").value;
    } else if (flag == 'CSV') {
        var textToWrite = '';
        var textFileAsBlob = new Blob([textToWrite], { type: 'text/csv' });
        var fileNameToSaveAs = document.getElementById("csv_name").value;
    }

    var downloadLink = document.createElement("a");
    downloadLink.download = fileNameToSaveAs;
    downloadLink.innerHTML = "Download File";
    if (window.URL != null) {
        // Chrome allows the link to be clicked
        // without actually adding it to the DOM.
        downloadLink.href = window.webkitURL.createObjectURL(textFileAsBlob);
    } else {
        // Firefox requires the link to be added to the DOM
        // before it can be clicked.
        downloadLink.href = window.URL.createObjectURL(textFileAsBlob);
        downloadLink.onclick = destroyClickedElement;
        downloadLink.style.display = "none";
        document.body.appendChild(downloadLink);
    }

    downloadLink.click();
}
