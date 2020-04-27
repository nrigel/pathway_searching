var relative_dir; //relative_dir on the server, needed for xhttp call
var session_id;

var pathway_metabolites = new Set(); // pathway doc ids to list of metabolites
var pathway_data = {}; // pathway checkbox ids to nodes and edges
var metabolites_to_pathways = {}; // metabolite names to list of pathway ids
var pathway_checkboxes = {}; // doc ids for userfile to pathway checkbox ids
var pathwayid_pathwayname = {}; // all pathway ids to the displayname
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
                de_visualize(flag);
            }            
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

        for (let pathway in pathway_data) {
            if (pathway == null) {continue;}
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
                //console.log(oReq.response);
                var result = JSON.parse(oReq.response);
                var pathway_ids = set_pathwaylist(userfile, result.pathway_list);
                
                for (i = 0; i < result.pathway_list.length; i++) { 
                    if (pathway_ids[i] === null) {
                        continue;
                    }
                    pathwayid_pathwayname[pathway_ids[i]] = result.pathway_list[i];
                    pathway_data[pathway_ids[i]] = result.pathway_data[result.pathway_list[i]];
                    var nodes = pathway_data[pathway_ids[i]].nodes;
                    var edges = pathway_data[pathway_ids[i]].edges;

                    // calculate fold change color scale
                    let fold_changes = [];
                    for (j = 0; j < nodes.length; j++) { 
                        if (nodes[j].data['Fold Change']) {
                            fold_changes.push(nodes[j].data['Fold Change']);
                        }
                    }
                    let fc_min = Math.min.apply(Math, fold_changes);
                    let fc_max = Math.max.apply(Math, fold_changes);

                    if (fc_min !== fc_max) {
                        var fc_range = fc_max-fc_min;
                    } else {
                        var fc_range = fc_min;
                    }

                    var colors = ['#E50002', '#E02C00', '#DC5900', '#D88500', '#D4AF00', '#C8CF00', '#99CB00', '#6CC700', '#41C300', '#18BF00', '#0074BF'];
                    
                    for (j = 0; j < nodes.length; j++) {
                        // add metabolite to metabolite pathway dict
                        if (metabolites_to_pathways[nodes[j].data["displayName"]] == null){
                            metabolites_to_pathways[nodes[j].data["displayName"]] = [];
                        }
                        metabolites_to_pathways[nodes[j].data["displayName"]].push(pathway_ids[i]);
                        // add tag to node id so that we can keep searches separate
                        nodes[j].data["js_pathway_id"] = pathway_ids[i];
                        // add metabolite name to our selection list
                        pathway_metabolites.add(nodes[j].data["displayName"]);   
                        // add specific pathway id into node ids
                        var node_id = nodes[j].data["id"].toString();
                        nodes[j].data["id"] = node_id.concat('_').concat(pathway_ids[i]);
                        nodes[j].data["parent"] = pathway_ids[i];
                        nodes[j].data["class"] = "Metabolite";
                        if (nodes[j].data['Fold Change']) {
                            var fc = (nodes[j].data['Fold Change']-fc_min)/fc_range;
                            var color = null;
                            var diff = null;
                            for (k = 0; k < 10; k++) {
                                if (Math.abs(k-fc*10) < diff || diff === null) {
                                    diff = Math.abs(k-fc*10);
                                    color = colors[k];
                                } 
                            }
                            var size = 30*(nodes[j].data['P-Value']);
                        } else {
                            color = colors[10];
                            size = 30;
                        }
                        nodes[j].data['color'] = color;
                        nodes[j].data['size'] = size;
                        // add shape detail based on if in COLMAR
                        if (nodes[j].data["COLMAR"].length > 0) {
                            nodes[j].data['shape'] = "ellipse";
                        } else {
                            nodes[j].data['shape'] = "square";
                        }
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

                var met_menu = document.getElementById("metabolite_filter");
                
                for (i = met_menu.options.length - 1; i >= 1; i--) {
                    met_menu.remove(i);
                }

                var metabolitelist = Array.from(pathway_metabolites)
                metabolitelist.sort(function (a, b) {return a.toLowerCase().localeCompare(b.toLowerCase());});
                for (let item of metabolitelist) {
                    var option = document.createElement("option");
                    option.text = item; // metabolite name
                    //option.setAttribute("id", "plmetlist_".concat(filename))
                    met_menu.add(option);
                }
            } 
        document.getElementById("loader").style.display = "none";
        }
        oReq.send(params);
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
        label.setAttribute("id", id.concat("_label"));
        label.htmlFor = id;
        label.appendChild(document.createTextNode(filename));
        ul.appendChild(x);
        ul.appendChild(label);
        ul.appendChild(l);
    } 
    
    var pathway_ids = {};
    var l = document.getElementById("pl_".concat(filename));
    var i;
    for (i = 0; i < pathway_list.length; i++) { 
        var id = "plbox_".concat(filename).concat("_").concat(pathway_list[i])
        if (document.getElementById(id) == null) {
            var x = document.createElement("input");
            x.setAttribute("type", "checkbox");
            x.setAttribute("id", id);
            x.setAttribute("onclick", "check_change('".concat(id).concat("')"));
            var label = document.createElement("label");
            label.htmlFor = id;
            label.setAttribute("id", id.concat("_label"));
            label.appendChild(document.createTextNode(pathway_list[i]));
            l.appendChild(x);
            l.appendChild(label);
            l.appendChild(document.createElement("br"));
            pathway_ids[i] = id;
            pathway_checkboxes["plbox_".concat(filename)].push(id);
        } else {
            pathway_ids[i] = null;
        }
    }
    return pathway_ids
}

function clear_pathwaylist() {
    const ul = document.getElementById("pathway_list");
    ul.innerHTML = "";
    cy.remove(cy.elements('node'));
    cy.remove(cy.elements('edge'));
}

function setup_cy(container_id) {
    // colors = [#E50002, #E02C00, #DC5900, #D88500, #D4AF00, #C8CF00, #99CB00, #6CC700, #41C300, #18BF00], #0074BF
    var cy = window.cy = cytoscape({
                container: document.getElementById('cy'),

                style: [
                    {selector: 'node',
                        style: {'content': 'data('.concat(node_labels).concat(')'),
                                'background-color': function(ele) {return ele.data()['color'];},
                                'shape': function(ele) {return ele.data()['shape']},
                                'height': function(ele) {return ele.data()['size']},
                                'width': function(ele) {return ele.data()['size']},
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
                                'background-opacity': 0.2,
                                'background-color': '#FF00FF'
                        }
                    },

                    {selector: ':selected',
                        style: {'background-color': 'yellow', 'line-color': 'yellow', 
                        'target-arrow-color': 'black', 'source-arrow-color': 'black',}},
                    
                    {selector: 'edge',
                        style: {'label': function(ele) {
                            if (edge_labels) {
                                return 'data('.concat(edge_labels).concat(')')
                            } else {
                                return ''
                            }
                        },
                        
                            'curve-style': 'bezier',
                            'target-arrow-shape': 'triangle',
                            'font-size': '18px'}}
                ],
                elements: {nodes: [], edges: []}, 
                layout: {name: cy_layout,
                        nodeDimensionsIncludeLabels: true}});
}

function visualize(pathway, nodes, edges) {
    // this will overlap the nodes if the same node label already exists...
    cy.add([{'data': {'id': pathway, 'displayName': '', 'pathwayName': pathwayid_pathwayname[pathway], 'type': 'Pathway', 'size': 0, 'color': 'red', 'shape': 'ellipse', 'js_pathway_id': pathway}}]);
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
}

function de_visualize(pathway) {
    cy.remove(cy.elements('node[js_pathway_id = "'.concat(pathway).concat('"]')));
    cy.remove(cy.elements('edge[js_pathway_id = "'.concat(pathway).concat('"]')));
}