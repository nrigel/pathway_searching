var mainplot; //hsqc plot object
var mainplot1;
var mainplot3; //hsqc_tocsy plot object
var mainplot2; //tocsy plot object
var picked_peaks; //picked peaks from the server
var f1, f2, f3;

var mytrace = new Array(); //input peaks
var mytrace1 = new Array();
var mytrace2 = new Array();
var compounds = new Array(); //matched compounds
var current_compound; //index of current showing compound, in the details region
var current_data0 = new Array();
var scale = new Array(); //contour level control  , top hsqc panel
var scale1 = new Array(); //contour level control , hsqc
var scale2 = new Array(); //contour level control , tocsy
var scale3 = new Array(); //contour level control , hsqc_tocsy

var stage; //global var to define which stage we are in:  plotting or matching
var tooltiv; //tooltiv obj
var relative_dir; //relative_dir on the server, needed for xhttp call
var session_id;

var show_all_flag = 0; //flag of show all compounds or show no compounds

var compound_textToWrite;
var peak_textToWrite;

var for_debug;

var key_flag = 0;

//predefined shapes 
var triangle = "M0,-5.26429605180997L6.078685485212741,5.26429605180997 -6.078685485212741,5.26429605180997Z";
var diamond = "M0,-7.444838872816797L4.298279727294168,0 0,7.444838872816797 -4.298279727294168,0Z";
var cross = "M-5.366563145999495,-1.7888543819998317H-1.7888543819998317V-5.366563145999495H1.7888543819998317V-1.7888543819998317H5.366563145999495V1.7888543819998317H1.7888543819998317V5.366563145999495H-1.7888543819998317V1.7888543819998317H-5.366563145999495Z";
var square = "M-4,-4L4,-4 4,4 -4,4Z";
var circle = "M0,4.51351666838205A4.51351666838205,4.51351666838205 0 1,1 0,-4.51351666838205A4.51351666838205,4.51351666838205 0 1,1 0,4.51351666838205Z";

var shapes = [triangle, diamond, cross, square, circle];

//used for automatic re-referencing.
var dss = [
    { cs_x: 0.0001, cs_y: 0.0010, name: "dss" },
    { cs_x: 0.6258, cs_y: 17.683, name: "dss" },
    { cs_x: 1.7561, cs_y: 21.802, name: "dss" },
    { cs_x: 2.9074, cs_y: 57.068, name: "dss" }
];

var alanine = [
    { cs_x: 1.4661, cs_y: 18.868, name: "alanine" },
    { cs_x: 3.7669, cs_y: 53.216, name: "alanine" }
];

var leucine = [
    { cs_x: 0.9421, cs_y: 23.589, name: "leucine" },
    { cs_x: 0.9555, cs_y: 24.751, name: "leucine" },
    { cs_x: 1.7027, cs_y: 26.871, name: "leucine" },
    { cs_x: 1.7028, cs_y: 42.526, name: "leucine" },
    { cs_x: 3.7219, cs_y: 56.112, name: "leucine" }
];

var glucose = [
    { cs_x: 3.8849, cs_y: 63.454, name: "glucose" },
    { cs_x: 3.7129, cs_y: 63.454, name: "glucose" },
    { cs_x: 3.4533, cs_y: 78.690, name: "glucose" },
    { cs_x: 3.3973, cs_y: 72.309, name: "glucose" },
    { cs_x: 3.4741, cs_y: 78.473, name: "glucose" },
    { cs_x: 3.2334, cs_y: 76.854, name: "glucose" },
    { cs_x: 4.6359, cs_y: 98.642, name: "glucose" }
];

var lactic = [
    { cs_x: 1.3176, cs_y: 22.813, name: "lactic" },
    { cs_x: 4.1008, cs_y: 71.227, name: "lactic" }
];



var target;

var totalpeak = 20;



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
    target = [];


    if (document.getElementById("dss").checked) {
        target = target.concat(dss);
    }
    if (document.getElementById("alanine").checked) {
        target = target.concat(alanine);
    }
    if (document.getElementById("leucine").checked) {
        target = target.concat(leucine);
    }
    if (document.getElementById("glucose").checked) {
        target = target.concat(glucose);
    }
    if (document.getElementById("lactic").checked) {
        target = target.concat(lactic);
    }



    if (typeof mainplot !== "undefined") mainplot.dodss(target);

    totalpeak = target.length;

    updatetotalpeak();
}

function updatetotalpeak() {
    document.getElementById("totalpeak").value = totalpeak;

    var matchpeak = parseInt(document.getElementById("matchpeak").value);

    if (matchpeak > totalpeak)
        document.getElementById("matchpeak").value = totalpeak;

}




function doref() {
    var ref1 = parseFloat(document.getElementById("ref1").value);
    var ref2 = parseFloat(document.getElementById("ref2").value);

    if (ref1 >= 90 && ref2 >= 90) {
        alert('autoref calculation failed!');
        return;
    }



    apply_ref(ref1, ref2);


    mainplot.doref();
    if (typeof mainplot1 !== "undefined") { mainplot1.doref(); }
    if (typeof mainplot2 !== "undefined") { mainplot2.doref(); }
    if (typeof mainplot3 !== "undefined") { mainplot3.doref(); }


    $.ajax({ url: base_url.concat("/index.php/colmarm/saveref"), type: "post", data: { "ref1": ref1, "ref2": ref2, "relative_dir": relative_dir } });

    document.getElementById("ref1").value = "0.0";
    document.getElementById("ref2").value = "0.0";
}


function apply_ref(ref1, ref2) {

    for (var i = 0; i < picked_peaks.length; i++) {
        picked_peaks[i].cs_x += ref1;
        picked_peaks[i].cs_y += ref2;
    }

    for (var i = 0; i < f1.contour.length; i++) {
        for (var j = 0; j < f1.contour[i].data.length; j++) {
            for (var k = 0; k < f1.contour[i].data[j].length; k++) {
                f1.contour[i].data[j][k][0] += ref1;
                f1.contour[i].data[j][k][1] += ref2;
            }
            f1.contour[i].spans[j][0] += ref1;
            f1.contour[i].spans[j][1] += ref1;
            f1.contour[i].spans[j][2] += ref2;
            f1.contour[i].spans[j][3] += ref2;
        }
    }


    if (typeof f2 !== "undefined" && typeof f2.contour !== "undefined") {
        for (var i = 0; i < f2.contour.length; i++) {
            for (var j = 0; j < f2.contour[i].data.length; j++) {
                for (var k = 0; k < f2.contour[i].data[j].length; k++) {
                    f2.contour[i].data[j][k][0] += ref1;
                    f2.contour[i].data[j][k][1] += ref1;
                }
                f2.contour[i].spans[j][0] += ref1;
                f2.contour[i].spans[j][1] += ref1;
                f2.contour[i].spans[j][2] += ref1;
                f2.contour[i].spans[j][3] += ref1;
            }
        }
    }

    if (typeof f3 !== "undefined" && typeof f3.contour !== "undefined") {
        for (var i = 0; i < f3.contour.length; i++) {
            for (var j = 0; j < f3.contour[i].data.length; j++) {
                for (var k = 0; k < f3.contour[i].data[j].length; k++) {
                    f3.contour[i].data[j][k][0] += ref1;
                    f3.contour[i].data[j][k][1] += ref2;
                }
                f3.contour[i].spans[j][0] += ref1;
                f3.contour[i].spans[j][1] += ref1;
                f3.contour[i].spans[j][2] += ref2;
                f3.contour[i].spans[j][3] += ref2;
            }
        }
    }

}





function autoref() {
    if (target.length === 0)
        return;


    var t = "";

    for (var i = 0; i < target.length; i++) {
        t = t.concat(target[i].cs_x.toString(), " ", target[i].cs_y.toString(), "\r\n");
    }


    $.ajax({
        type: "POST",
        url: base_url.concat("/index.php/colmarm/autoref"),
        data: { "peaks1": $('#peaklist').val(), "peaks2": t, "relative_dir": relative_dir, "nmatch": $('#matchpeak').val(), "cut_c1": "7.5", "cut_h1": "0.5", "cut_c2": "0.3", "cut_h2": "0.03" },
        success: function(data) {
            console.log(data);
            result = JSON.parse(data);
            $('#ref2').val(result.c.toString())
            $('#ref1').val(result.h.toString())
            mainplot.ref(result.h, result.c);
        }
    });

}




//define base_url     
if (!window.location.origin) {
    window.location.origin = window.location.protocol + "//" + window.location.hostname + (window.location.port ? ':' + window.location.port : '');
}
var base_url = window.location.origin;

window.onbeforeunload = closingCode;



$(document).ready(function() {

    stage = 1;

    document.getElementById("button2").disabled = true;
    document.getElementById("button3").disabled = true;
    document.getElementById("plot").style.display = "none";

    document.getElementById("plot1").style.display = "none";
    document.getElementById("plot2").style.display = "none";
    document.getElementById("plot3").style.display = "none";


    document.getElementById('reporting').style.display = "none";
    document.getElementById('detail').style.display = "none";


    $(window).resize(function() {
        resize();
    });

    tooldiv = d3.select("body")
        .append("div")
        .attr("class", "tooltip2")
        .style("opacity", 0);

    $('input:radio[name=color]').change(function() { color_change(); });
    $('input:radio[name=peak]').change(function() { cross_change(); });


    $('form#process').submit(function(e) {
        var fd = new FormData(document.getElementById("process"));
        var oReq = new XMLHttpRequest();
        var base_url = window.location.origin;

        document.getElementById("button1").disabled = true;
        document.getElementById("button11").disabled = true;
        document.getElementById("fileinput").disabled = true;
        document.getElementById("infor").innerHTML = "Processing your spectrum on the server. Please be patient, it may take up to 2 minutes for large files.";

        oReq.open("POST", base_url.concat('/index.php/colmarm/getdata/1'), true);
        oReq.setRequestHeader("Authorization", null);
        oReq.onload = function(oEvent) {
            if (oReq.status === 200) {
                //alert(oReq.responseText);
                //try{
                var result = JSON.parse(oReq.responseText);
                document.getElementById("infor").innerHTML = "";
                if (typeof result.error !== "undefined") {
                    alert(result.error);
                    document.getElementById("button1").disabled = false;
                    document.getElementById("fileinput").disabled = false;
                    document.getElementById("infor").innerHTML = "";
                } else {
                    f1 = result.f1;
                    picked_peaks = result.picked_peaks;
                    document.getElementById("dir").value = result.dir;
                    relative_dir = result.dir;
                    session_id = result.session_id;
                    document.getElementById('id_info').innerHTML = "<p> You session ID is " + result.session_id + " and user chosen name is " + result.session_name + "</p>";
                    setscale(result.f1.r1);
                    draw(result.f1, result.picked_peaks);
                    setscale1(result.f1.r1);
                    draw1(result.f1);
                    document.getElementById("plot1").style.display = "block";
                    if (typeof result.f2 !== "undefined") {
                        f2 = result.f2;
                        setscale2(result.f2.r2);
                        draw2(result.f2);
                        document.getElementById("plot2").style.display = "block";
                    }
                    if (typeof result.f3 !== "undefined") {
                        f3 = result.f3;
                        setscale3(result.f3.r3);
                        draw3(result.f3);
                        document.getElementById("plot3").style.display = "block";
                    }
                    document.getElementById("plot").style.display = "block";
                    document.getElementById("button2").disabled = false;
                    document.getElementById("message").innerHTML = "You are now at peaks adjusting state. you can click on one peak to delete it";
                }
                /*}catch(e){
                    document.getElementById("infor").innerHTML=""; 
                    alert(e+" Please contact us."); //error in the above string(in this case,yes)!
                }*/



            }
        };
        oReq.send(fd);
        e.preventDefault();
    });



    $('form#process_load').submit(function(e) {
        var fd = new FormData(document.getElementById("process_load"));
        var oReq = new XMLHttpRequest();
        var base_url = window.location.origin;

        document.getElementById("button11").disabled = true;
        document.getElementById("button1").disabled = true;
        document.getElementById("infor").innerHTML = "Load your spectrum from the server. Please be patient.";

        oReq.open("POST", base_url.concat('/index.php/colmarm/loaddata1'), true);
        oReq.setRequestHeader("Authorization", null);
        oReq.onload = function(oEvent) {
            if (oReq.status === 200) {
                //alert(oReq.responseText);
                //try{
                result = JSON.parse(oReq.responseText.replace(/\bNaN\b/g, "null"));
                //console.log(oReq.responseText);
                document.getElementById("infor").innerHTML = "";
                if (typeof result.error !== "undefined") {
                    alert(result.error);
                    document.getElementById("infor").innerHTML = "";
                    document.getElementById("button11").disabled = false;
                    document.getElementById("button1").disabled = false;
                } else {
                    document.getElementById("dir").value = result.dir;
                    relative_dir = result.dir;
                    session_id = result.session_id;
                    document.getElementById('id_info').innerHTML = "<p> You session ID is " + result.session_id + " and user chosen name is " + result.session_name + "</p>";
                    picked_peaks = result.picked_peaks;
                    f1 = result.f1;
                    if (typeof result.f2 !== "undefined") { f2 = result.f2; }
                    if (typeof result.f3 !== "undefined") { f3 = result.f3; }

                    if (typeof result.ref1 !== "undefned" && typeof result.ref2 !== "undefned") {
                        apply_ref(result.ref1, result.ref2);
                    }

                    setscale(result.f1.r1);
                    draw(result.f1, result.picked_peaks);
                    setscale1(result.f1.r1);
                    draw1(result.f1);
                    document.getElementById("plot1").style.display = "block";
                    if (typeof result.f2 !== "undefined") {
                        setscale2(result.f2.r2);
                        draw2(result.f2);
                        document.getElementById("plot2").style.display = "block";
                    }
                    if (typeof result.f3 !== "undefined") {
                        setscale3(result.f3.r3);
                        draw3(result.f3);
                        document.getElementById("plot3").style.display = "block";
                    }
                    document.getElementById("plot").style.display = "block";
                    document.getElementById("button2").disabled = false;
                    document.getElementById("message").innerHTML = "You are now at peaks adjusting state. you can click on one peak to delete it";

                    if (typeof result.compound !== "undefined") {
                        stage = 2;
                        process_query(result.mytrace, result.compound);
                        show_one_compound(0);
                    }
                }

            }
        };
        oReq.send(fd);
        e.preventDefault();
    });

    /*
    $('form#process_load_query').submit(function(e) {

        document.getElementById("session_id").value = session_id;
        var fd = new FormData(document.getElementById("process_load_query"));
        var oReq = new XMLHttpRequest();
        var base_url = window.location.origin;


        oReq.open("POST", base_url.concat('/index.php/colmarm/loaddata2'), true);
        oReq.setRequestHeader("Authorization", null);
        oReq.onload = function(oEvent) {
            if (oReq.status === 200) {
                result = JSON.parse(oReq.responseText);
                if (typeof result.compound !== "undefined") {
                    stage = 2;
                    process_query(result.mytrace, result.compound);
                    show_one_compound(0);
                }
            }
        };
        oReq.send(fd);
        e.preventDefault();
    });
    */


    $('form#hsqc').submit(function(e) {
        e.preventDefault();

        var fd = new FormData(document.getElementById("hsqc"));
        var oReq = new XMLHttpRequest();

        document.getElementById("infor").innerHTML = "Wait… matching in progress…";

        var base_url = window.location.origin;
        oReq.open("POST", base_url.concat('/index.php/colmarm/query_json'), true);
        oReq.onload = function(oEvent) {
            stage = 2;
            if (oReq.status == 200) {
                for_debug = oReq.responseText;
                result = JSON.parse(oReq.responseText);
                document.getElementById("infor").innerHTML = "";
                if (typeof result.error !== "undefined") {
                    alert(result.error);
                } else {
                    process_query(result.mytrace, result.compound);
                    show_one_compound(0);
                }
            } else {
                document.getElementById("infor").innerHTML = "Error " + oReq.status + " occurred submit your job.<br \/>";
            }

        };


        mytrace = [];
        mytrace1 = [];
        mytrace2 = [];
        compounds = [];
        oReq.send(fd);
        e.preventDefault();
    });


});


function process_query(a, b) {

    mytrace = a;
    compounds = b;

    for (var i = compounds.length - 1; i >= 0; i--) {
        if (typeof compounds[i].ngood === "undefined" || compounds[i].ngood <= 0) {
            compounds.splice(i, 1);
        }
    }

    if (compounds.length === 0) {
        alert("no matching!");
        return;
    }

    var bts = document.getElementsByClassName("compound_button");
    for (var i = 0; i < bts.length; i++)
        bts[i].style.display = "inline";

    var x = document.getElementById("selected_compound");
    var tt = x.options.length;

    for (i = tt - 1; i >= 0; i--) {
        x.options[i] = null;
    }
    x.options.length = 0;

    for (var i = 0; i < compounds.length; i++) {
        var option = document.createElement("option");
        option.text = compounds[i]['name'];
        x.add(option);
    }


    var c1 = 0,
        c2 = 0;
    for (var i = 0; i < mytrace.length; i++) {
        if (mytrace[i].z.length >= 1)
            c1++;
        else
            c2++;
    }

    generate_compound_table();
    document.getElementById('peak_info').innerHTML = "<p> " + c1 + " HSQC peaks match compounds and " + c2 + " HSQC peaks don't match any compounds.</p>";
    mainplot.addpeaks(mytrace, compounds);
    mainplot1.picked_peaks(mytrace);

    document.getElementById('reporting').style.display = "block";
    document.getElementById('detail').style.display = "block";
    document.getElementById("button2").disabled = true;
    document.getElementById("button3").disabled = false;
    document.getElementById("doref").disabled = true;
    document.getElementById("message").innerHTML = "You are now at matching state.";
}

function set_quality(qua) {
    for (var i = 0; i < compounds.length; i++) {
        compounds[i].quality = qua[i];
    }
}

function set_limit(tt) {
    document.getElementById("contour-lower").value = tt;
    mainplot.update_contour(+tt);
    var c = scale[(+tt) - 1] / scale[0] * 2;
    document.getElementById("contour-value-lower").innerHTML = c.toFixed(2);

}

function process_uniqe(index) {
    var temp = "";
    var t = new Array();
    var peaks = compounds[index].peaks;
    for (var i = 0; i < peaks.length; i++) {
        if (typeof peaks[i].share !== "undefined" && peaks[i].share.length > 0) {
            for (var j = 0; j < peaks[i].share.length; j++)
                t.push(peaks[i].share[j]);
        }
    }

    t = t.removeDuplicates();

    compounds[index].shared = t;

    for (var i = 0; i < t.length; i++) {
        if (t[i] === index)
            continue;
        temp = temp + (t[i] + 1).toFixed(0) + " ";
    }

    return temp;
}

function generate_compound_table() {
    var thead = "ID Name Matching_ratio 13C_RMSD 1H_RMSD Uniqueness origin pubchemID Uniqueness_detail\n"
    var tdata = new Array();
    var temp = "";
    for (var i = 0; i < compounds.length; i++) {
        var compound = compounds[i];
        var t = new Array();
        t.push((i + 1).toFixed(0));
        t.push(compound['name']);       
        t.push(compound['matching'].toFixed(2));     
        t.push(compound['c_rmsd'].toFixed(2));  
        t.push(compound['p_rmsd'].toFixed(3));     
        t.push(compound['unique']);      
        t.push(compound['origin']);    
        t.push(compound['pubchemid']);
        //temp=temp+compound['amp_ratio'].toFixed(2);t.push(compound['amp_ratio'].toFixed(2));
        t.push(process_uniqe(i));
        tdata.push(t);
    }
    var table = createTable(thead.split(" "), tdata);
    var table_region = document.getElementById('table_region');
    while (table_region.firstChild) {
        table_region.removeChild(table_region.firstChild);
    }
    table_region.appendChild(createPara("Compounds information"));
    table_region.appendChild(table);

    for (var i = 1; i < table.rows.length; i++) {
        table.rows[i].cells[1].onclick = function() { tableclick(this); }
    }

    for (var i = 1; i < table.rows.length; i++) {
        var link = document.createElement('a');
        link.setAttribute('href', compounds[i-1]['link']);
        link.setAttribute('target',"_blank");
        var linkText = document.createTextNode(table.rows[i].cells[6].innerHTML);
        link.appendChild(linkText);
        table.rows[i].cells[6].innerHTML="";
        table.rows[i].cells[6].appendChild(link);
    }
    sorttable.makeSortable(table);
}


function tableclick(node) {
    var s = (node.innerText || node.textContent);
    s=s.trim();
    var index = compounds.findIndex(function(sample) { return sample.name === s });
    show_one_compound(index);
    document.getElementById("selected_compound").selectedIndex = index;
}


function compound_report() {
    
    var newline="\n";
    if(getOS()==='Windows')
    {
        newline="\r\n";
    }

    var temp;
    var delimiter;
    var radios = document.getElementsByName('delimiter');
    if (radios[0].checked) {
        delimiter = " ";
        temp = "ID Name Matching_ratio 13C_RMSD 1H_RMSD Uniqueness User_selected"+newline;
    } else {
        delimiter = ",";
        temp = "ID,Name,Matching_ratio,13C_RMSD,1H_RMSD,Uniqueness,User_selected"+newline;
    }


    for (var i = 0; i < compounds.length; i++) {
        var compound = compounds[i];
        temp = temp + (i + 1).toFixed(0) + delimiter;
        temp = temp + compound['name'] + delimiter;
        temp = temp + compound['matching'].toFixed(2) + delimiter;
        temp = temp + compound['c_rmsd'].toFixed(2) + delimiter;
        temp = temp + compound['p_rmsd'].toFixed(3) + delimiter;
        temp = temp + compound['unique'] + delimiter;
        temp = temp + compound['origin'] + delimiter;
        temp = temp + delimiter + trans(compound['quality']);
        temp = temp + newline;

    }
    compound_textToWrite = temp;



    //Save the peak report to a file named from editbox preak_report_name
    var textFileAsBlob = new Blob([compound_textToWrite], { type: 'text/plain' });
    var fileNameToSaveAs = document.getElementById("compound_report_name").value;

    var downloadLink = document.createElement("a");
    downloadLink.download = fileNameToSaveAs;
    downloadLink.innerHTML = "Download File";
    if (window.URL != null) {
        // Chrome allows the link to be clicked
        // without actually adding it to the DOM.
        downloadLink.href = window.URL.createObjectURL(textFileAsBlob);
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

function peak_report() {
    
    var newline="\n";
    if(getOS()==='Windows')
    {
        newline="\r\n";
    }

    var delimiter;
    var textToWrite;
    var radios = document.getElementsByName('delimiter');
    if (radios[0].checked) {
        delimiter = " ";
        textToWrite = "Proton Carbon Amplitude compound_name\n";
    } else {
        delimiter = ",";
        textToWrite = "Proton,Carbon,Amplitude,compound_name\n";
    }


    for (var i = 0; i < mytrace.length; i++) {
        var peak = mytrace[i];
        var temp = peak['x'].toFixed(3) + delimiter + peak['y'].toFixed(3) + delimiter + peak['a'].toFixed(3);
        for (var j = 0; j < peak['z'].length; j++) {
            temp = temp + delimiter + compounds[peak['z'][j]].name;
        }
        temp = temp + newline;
        textToWrite = textToWrite + temp;
    }


    //Save the peak report to a file named from editbox preak_report_name
    var textFileAsBlob = new Blob([textToWrite], { type: 'text/plain' });
    var fileNameToSaveAs = document.getElementById("peak_report_name").value;

    var downloadLink = document.createElement("a");
    downloadLink.download = fileNameToSaveAs;
    downloadLink.innerHTML = "Download File";
    if (window.URL != null) {
        // Chrome allows the link to be clicked
        // without actually adding it to the DOM.
        downloadLink.href = window.URL.createObjectURL(textFileAsBlob);
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





function setscale(r) {
    scale = r;
    var c = scale[1] / scale[0] * 2;
    document.getElementById("contour-value-lower").innerHTML = c.toFixed(2);
    document.getElementById("contour-lower").max = r.length - 2;
}


function setscale1(r) {
    scale1 = r;
    var c = scale1[1] / scale1[0] * 4;
    document.getElementById("contour1-value-lower").innerHTML = c.toFixed(2);
    document.getElementById("contour1-lower").max = r.length - 2;
}


function setscale2(r) {
    scale2 = r;
    var c = scale2[1] / scale2[0] * 4;
    document.getElementById("contour2-value-lower").innerHTML = c.toFixed(2);
    document.getElementById("contour2-lower").max = r.length - 2;
}


function setscale3(r) {
    scale3 = r;
    var c = scale3[1] / scale3[0] * 4;
    document.getElementById("contour3-value-lower").innerHTML = c.toFixed(2);
    document.getElementById("contour3-lower").max = r.length - 2;
}





function draw(temp, picked_peaks) {

    var wid = parseInt(document.getElementById('body').clientWidth) - 180;
    var wid2 = wid * 0.65 + 50;

    document.getElementById('visualisation').style.height = wid2.toString().concat('px');
    document.getElementById('visualisation').style.width = wid.toString().concat('px');


    /*
    for (var i = 0; i < picked_peaks.length; i++) {
        if (typeof picked_peaks[i].cs_x === "undefined") { picked_peaks[i].cs_x = picked_peaks[i].x; }
        if (typeof picked_peaks[i].cs_y === "undefined") { picked_peaks[i].cs_y = picked_peaks[i].y; }
    }
    */


    var input = new Object();

    input.xscale = temp.xscale;
    input.yscale = temp.yscale;
    input.contour = temp.contour;

    input.PointData = picked_peaks;
    input.WIDTH = wid;
    input.HEIGHT = wid2;
    input.MARGINS = { top: 20, right: 20, bottom: 50, left: 50 };
    input.drawto = "#visualisation";
    input.drawto_legend = "#legend";
    input.drawto_peak = "#peaklist";




    mainplot = new plotit(input);
    mainplot.draw();
    mainplot.draw_peaks();


    check_change(0);

    d3.select("#contour-lower").on("input", function() {
        mainplot.update_contour(+this.value);
        document.getElementById("contour-value-lower").innerHTML = (scale[+this.value - 1] / scale[0] * 2).toFixed(2);
    });
    var tt = document.getElementById("contour-lower").value;
    mainplot.update_contour(+tt);


};


function draw1(input) {

    //var wid=parseInt(document.getElementById('body').clientWidth)-10;
    var tt = document.getElementById("size").value;
    var wid = 300 + 100 * tt;
    var wid2 = wid * 0.65 + 50;

    document.getElementById('visualisation1').setAttribute("height", wid2.toString());
    document.getElementById('visualisation1').setAttribute("width", wid.toString());



    input.WIDTH = wid;
    input.HEIGHT = wid2;
    input.MARGINS = { top: 20, right: 20, bottom: 50, left: 50 };
    input.drawto = "#visualisation1";
    input.xtext = "Proton Chemical Shift (ppm)";
    input.ytext = "Carbon Chemical Shift (ppm)";
    input.clipid = "clip1";
    input.x_symbol_size = document.getElementById("cutoff_h").value;
    input.y_symbol_size = document.getElementById("cutoff_c").value;;
    input.color_flag = get_color_select();
    input.index = 1;
    input.type = "hsqc";

    mainplot1 = new plotit3(input);
    mainplot1.init_draw();
    mainplot1.draw_contour();
    mainplot1.init_drag();

    d3.select("#contour1-lower").on("input", function() {
        mainplot1.update_contour_level(+this.value);
        document.getElementById("contour1-value-lower").innerHTML = (scale1[+this.value - 1] / scale1[0] * 4).toFixed(2);
    });


    var tt = document.getElementById("contour1-lower").value;
    mainplot1.update_contour_level(+tt);
}

function get_color_select() {
    var n;
    var radios = document.getElementsByName('color');

    for (var i = 0, length = radios.length; i < length; i++) {
        if (radios[i].checked) {
            n = i;
            break;
        }
    }
    return n;
}

function get_cross_select() {
    var n;
    var radios = document.getElementsByName('peak');

    for (var i = 0, length = radios.length; i < length; i++) {
        if (radios[i].checked) {
            n = i;
            break;
        }
    }
    return n;
}


function color_change() {
    var t = get_color_select();
    if (typeof mainplot1 !== "undefined") mainplot1.contour_color(t);
    if (typeof mainplot2 !== "undefined") mainplot2.contour_color(t);
    if (typeof mainplot3 !== "undefined") mainplot3.contour_color(t);
}

function draw2(input) {

    //var wid=parseInt(document.getElementById('body').clientWidth)-10;
    var tt = document.getElementById("size").value;
    var wid = 300 + 100 * tt;
    var wid2 = wid * 0.65 + 50;

    document.getElementById('visualisation2').setAttribute("height", wid2.toString());
    document.getElementById('visualisation2').setAttribute("width", wid.toString());

    input.WIDTH = wid;
    input.HEIGHT = wid2;
    input.MARGINS = { top: 20, right: 20, bottom: 50, left: 50 };
    input.drawto = "#visualisation2";
    input.xtext = "Proton Chemical Shift (ppm)";
    input.ytext = "Proton Chemical Shift (ppm)";
    input.clipid = "clip2";
    input.x_symbol_size = document.getElementById("cutoff_h").value;;
    input.y_symbol_size = document.getElementById("cutoff_h").value;;
    input.color_flag = get_color_select();
    input.type = "tocsy";

    mainplot2 = new plotit3(input);
    mainplot2.init_draw();
    mainplot2.draw_contour();


    d3.select("#contour2-lower").on("input", function() {
        mainplot2.update_contour_level(+this.value);
        document.getElementById("contour3-value-lower").innerHTML = (scale2[+this.value - 1] / scale2[0] * 4).toFixed(2);
    });


    var tt = document.getElementById("contour2-lower").value;
    mainplot2.update_contour_level(+tt);

}



function draw3(input) {

    //var wid=parseInt(document.getElementById('body').clientWidth)-10;
    var tt = document.getElementById("size").value;
    var wid = 300 + 100 * tt;
    var wid2 = wid * 0.65 + 50;

    document.getElementById('visualisation3').setAttribute("height", wid2.toString());
    document.getElementById('visualisation3').setAttribute("width", wid.toString());


    input.WIDTH = wid;
    input.HEIGHT = wid2;
    input.MARGINS = { top: 20, right: 20, bottom: 50, left: 50 };
    input.drawto = "#visualisation3";
    input.xtext = "Proton Chemical Shift (ppm)";
    input.ytext = "Carbon Chemical Shift (ppm)";
    input.clipid = "clip3";
    input.x_symbol_size = document.getElementById("cutoff_h").value;;
    input.y_symbol_size = document.getElementById("cutoff_c").value;;
    input.color_flag = get_color_select();
    input.type = "hsqc-tocsy";

    mainplot3 = new plotit3(input);
    mainplot3.init_draw();
    mainplot3.draw_contour();

    d3.select("#contour3-lower").on("input", function() {
        mainplot3.update_contour_level(+this.value);
        document.getElementById("contour3-value-lower").innerHTML = (scale3[+this.value] / scale3[0] * 4).toFixed(2);
    });


    var tt = document.getElementById("contour3-lower").value;
    mainplot3.update_contour_level(+tt);

}



function resize() {

    var wid = parseInt(document.getElementById('body').clientWidth) - 180;
    var wid2 = wid * 0.65 + 50;

    document.getElementById('visualisation').style.height = wid2.toString().concat('px');
    document.getElementById('visualisation').style.width = wid.toString().concat('px');


    var input = {
        WIDTH: wid,
        HEIGHT: wid2
    };


    if ('undefined' !== typeof(mainplot)) {
        mainplot.update(input);
        //mainplot.draw();
        //mainplot.transition_data();
        if (stage >= 2) { mainplot.dolegend(); }
    }



}


function redo_step4() {
    stage = 1;
    document.getElementById("button2").disabled = false;
    document.getElementById("button3").disabled = true;
    document.getElementById("doref").disabled = false;

    document.getElementById('reporting').style.display = "none";
    document.getElementById('detail').style.display = "none";


    document.getElementById("message").innerHTML = "You are now at peaks adjusting state. you can click on one peak to delete it";
    mainplot.removepeaks();

    var bts = document.getElementsByClassName("compound_button");
    for (var i = 0; i < bts.length; i++)
        bts[i].style.display = "none";

    var x = document.getElementById("selected_compound");
    var tt = x.options.length;

    for (i = tt - 1; i >= 0; i--) {
        x.options[i] = null;
    }
    x.options.length = 0;

    compounds = [];
    mytrace = [];
    mytrace1 = [];
    mytrace2 = [];


    if (typeof mainplot1 !== "undefined") {
        mainplot1.addpeaks([]);
        mainplot1.update_range([
            [0, 10]
        ], [
            [0, 150]
        ]);
    }
    if (typeof mainplot2 !== "undefined") {
        mainplot2.addpeaks([]);
        mainplot2.update_range([
            [0, 10]
        ], [
            [0, 10]
        ]);
    }
    if (typeof mainplot3 !== "undefined") {
        mainplot3.addpeaks([]);
        mainplot3.update_range([
            [0, 10]
        ], [
            [0, 150]
        ]);
    }

    var infor = document.getElementById('detail_info');
    while (infor.firstChild) {
        infor.removeChild(infor.firstChild);
    }

    var tab = document.getElementById('table_region');
    while (tab.firstChild) {
        tab.removeChild(tab.firstChild);
    }

    var tt = document.getElementById("contour-lower").value;
    mainplot.update_contour(+tt);

    check_change(0);

}

function resetzoom() {
    mainplot.resetzoom();
}

function popzoom() {
    mainplot.popzoom();
}

function zoomout() {
    mainplot.zoomout();
}

function toggle_contour() {
    mainplot.toggle_contour();
}

function toggle_peak() {
    mainplot.toggle_peak();
}

function show_match() {
    show_all_flag = 1 - show_all_flag;
    mainplot.legendall(show_all_flag);
}


function data2range(input, t) {
    var out = [];
    input.sort(function(a, b) { return a - b; });
    for (var i = 0; i < input.length; i++) {
        out.push([input[i] - t, input[i] + t]);
    }

    for (var i = out.length - 1; i > 0; i--) {
        if (out[i][0] <= out[i - 1][1]) {
            out[i - 1][1] = out[i][1];
            out.splice(i, 1);
        }
    }
    return out;
}

function trans(input) {
    var out;
    switch (input) {
        case 1:
            out = "Good";
            break;
        case 2:
            out = "Fair";
            break;
        case 3:
            out = "Poor";
            break;
        case 4:
            out = "Unknown";
            break;
        default:
            out = "not processed yet";
    }
    return out;
}


function custom() {
    var contents = $("#custom_peak").val();
    var lines = contents.split(/[\r\n]+/g);
    current_data0 = [];

    var j = 0;
    for (var i = 0; i < lines.length; i++) {
        var temp = lines[i].replace(/^ +/g, "");
        temp = temp.replace(/\t+/g, " ");
        temp = temp.trim();
        var data = temp.split(/ +/);
        if (data.length == 2) {
            current_data0[j] = new Object();
            current_data0[j]['x'] = parseFloat(data[0]);
            current_data0[j]['y'] = parseFloat(data[1]);
            current_data0[j]['z'] = new Array();
            j++;
        }
    }

    var temp1 = [];
    var temp2 = [];
    for (var i = 0; i < current_data0.length; i++) {
        temp2.push(current_data0[i].y);
        temp1.push(current_data0[i].x);
    }

    var tt = document.getElementById("scale").value;
    var xr = data2range(temp1, 1.5 / temp1.length / tt);
    var yr = data2range(temp2, 15 / temp2.length / tt);




}

function show_one_compound_table(index) {
    var peaks = mytrace;


    var buff2 = "For compound " + compounds[index].name + ", matching ratio is " + compounds[index].matching.toFixed(2) + " ,carbon rmsd is " +
        compounds[index].c_rmsd.toFixed(2) + " proton rmsd is " + compounds[index].p_rmsd.toFixed(3) + " User_selected is " + trans(compounds[index].quality);



    var thead = ["Index", "database 1H", "database 13C", "exp 1H", "exp 13C"];
    var tdata = new Array();
    for (var i = 0; i < compounds[index].peaks.length; i++) {
        var t = new Array();
        t.push(i.toFixed(0));
        t.push(compounds[index].peaks[i].x.toFixed(2));
        t.push(compounds[index].peaks[i].y.toFixed(2));
        var id = compounds[index].peaks[i].id - 1;
        if (id >= 0) {
            t.push(peaks[id].x.toFixed(2));
            t.push(peaks[id].y.toFixed(2));
            //t.push(peaks[id].a.toExponential(2));
            //t.push(id.toFixed(0));
        } else {
            t.push("N.A.");
            t.push("N.A."); //t.push("N.A.");t.push("N.A.");
        }
        tdata.push(t);
    }



    var table = createTable(thead, tdata);
    var para2 = createPara(buff2);

    var infor = document.getElementById('detail_info');
    while (infor.firstChild) {
        infor.removeChild(infor.firstChild);
    }
    infor.appendChild(para2);
    infor.appendChild(table);
}


function show_one_compound(index) {
    if (compounds.length === 0) return;

    current_compound = index;

    show_one_compound_table(index);


    current_data0 = compounds[index].peaks;

    var temp1 = [];
    var temp2 = [];
    for (var i = 0; i < current_data0.length; i++) {
        temp2.push(current_data0[i].y);
        temp1.push(current_data0[i].x);
    }


    var tt = document.getElementById("scale").value;
    var xr = data2range(temp1, 1.5 / temp1.length / tt);
    var yr = data2range(temp2, 15 / temp2.length / tt);

    if (typeof mainplot1 !== "undefined") { mainplot1.update_range(xr, yr); }
    if (typeof mainplot2 !== "undefined") { mainplot2.update_range(xr, xr); }
    if (typeof mainplot3 !== "undefined") { mainplot3.update_range(xr, yr); }

    var t = get_cross_select();
    if (typeof mainplot1 !== "undefined") {
        mainplot1.update_compound_index(index);
        mainplot1.addpeaks(compounds[index].peaks);
    }


    if (typeof mainplot2 !== "undefined") {
        mainplot2.addpeaks(compounds[index].peaks);
        mainplot2.init_spinsystem(compounds[index].system);
        mainplot2.cross_change(t);
    }

    if (typeof mainplot3 !== "undefined") {
        mainplot3.addpeaks(compounds[index].peaks);
        mainplot3.init_spinsystem(compounds[index].system);
        mainplot3.cross_change(t);
    }


    mainplot.legendall(0);
    mainplot.toggle_peak(0);

    for (var i = 0; i < compounds[index].shared.length; i++)
        mainplot.legendClick(compounds[index].shared[i]);

    if (compounds[index].shared.length === 0) mainplot.legendClick(index);

}


function cross_change() {
    var t = get_cross_select();
    if (typeof mainplot2 !== "undefined") { mainplot2.cross_change(t); }
    if (typeof mainplot3 !== "undefined") { mainplot3.cross_change(t); }


}

function scale_change(tt) {

    if (current_data0.length === 0) return;


    var temp1 = [];
    var temp2 = [];
    for (var i = 0; i < current_data0.length; i++) {
        temp2.push(current_data0[i].y);
        temp1.push(current_data0[i].x);
    }



    var xr = data2range(temp1, 1.5 / temp1.length / tt);
    var yr = data2range(temp2, 15 / temp2.length / tt);

    if (typeof mainplot1 !== "undefined") { mainplot1.update_range(xr, yr); }
    if (typeof mainplot3 !== "undefined") { mainplot3.update_range(xr, yr); }
    if (typeof mainplot2 !== "undefined") { mainplot2.update_range(xr, xr); }

}

function size_change(tt) {
    var wid = 300 + 100 * tt;
    var wid2 = wid * 0.65 + 50;

    //document.getElementById('visualisation3').style.height = wid2.toString().concat('px');
    //document.getElementById('visualisation3').style.width = wid.toString().concat('px');


    document.getElementById('visualisation3').setAttribute("height", wid2.toString());
    document.getElementById('visualisation3').setAttribute("width", wid.toString());

    document.getElementById('visualisation2').setAttribute("height", wid2.toString());
    document.getElementById('visualisation2').setAttribute("width", wid.toString());

    document.getElementById('visualisation1').setAttribute("height", wid2.toString());
    document.getElementById('visualisation1').setAttribute("width", wid.toString());



    var input = {
        WIDTH: wid,
        HEIGHT: wid2
    };


    if (typeof mainplot1 !== "undefined") { mainplot1.update(input); }
    if (typeof mainplot3 !== "undefined") { mainplot3.update(input); }
    if (typeof mainplot2 !== "undefined") { mainplot2.update(input); }


}


//call funtion when user click "good, fair, poor or unknown".
function label(flag) {
    //flag==1: Good 
    //flag==2: Fair
    //flag==3: Poor
    //Flag==4: Unknown
    compounds[current_compound].quality = flag;
    compound_change(1);

}

//call back function when user select one from droplist
function show_compound() {
    var x = document.getElementById("selected_compound");
    //var value = x.options[x.selectedIndex].value;
    show_one_compound(x.selectedIndex);
}


function compound_change(flag) {
    var x = document.getElementById("selected_compound");

    if (x.selectedIndex === x.options.length - 1 && flag !== 0)
        return;
    if (x.selectedIndex === 0 && flag === 0)
        return;


    if (flag == 0) {
        x.selectedIndex--;
    } else {
        x.selectedIndex++;
    }

    show_one_compound(x.selectedIndex);
}


function saveTextAsFile() {
    var textToWrite = document.getElementById("peaklist").value;
    var textFileAsBlob = new Blob([textToWrite], { type: 'text/plain' });
    var fileNameToSaveAs = document.getElementById("inputFileNameToSaveAs").value;

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


function save_vector() {
    if (typeof mainplot1 !== "undefined")
    {
        saveAs(new Blob([document.getElementById("vis1").innerHTML]),"hsqc.svg");
    }
    if (typeof mainplot2 !== "undefined")
    {
        saveAs(new Blob([document.getElementById("vis2").innerHTML]),"tocsy.svg");
    }
    if (typeof mainplot3 !== "undefined")
    {
        saveAs(new Blob([document.getElementById("vis3").innerHTML]),"hsqc-tocsy.svg");
    }
}

function save() {
    if (typeof mainplot1 !== "undefined") save_png("visualisation1", "hsqc.png");
    if (typeof mainplot2 !== "undefined") save_png("visualisation2", "tocsy.png");
    if (typeof mainplot3 !== "undefined") save_png("visualisation3", "hsqc-tocsy.png");

}





function save_png(input, fname) {


    var iDiv1 = document.createElement('div');
    iDiv1.id = 'svgdataurl';
    document.body.appendChild(iDiv1);

    var iDiv2 = document.createElement('div');
    iDiv2.id = 'pngdataurl';
    document.body.appendChild(iDiv2);



    var wid = document.getElementById(input).getAttribute("width");
    var height = document.getElementById(input).getAttribute("height");

    document.getElementById('canvas').setAttribute("height", height);
    document.getElementById('canvas').setAttribute("width", wid);

    var html = d3.select("#" + input)
        .attr("version", 1.1)
        .attr("xmlns", "http://www.w3.org/2000/svg")
        .node().parentNode.innerHTML;

    console.log(html);
    var imgsrc = 'data:image/svg+xml;base64,' + btoa(html);
    var img = '<img src="' + imgsrc + '">';


    d3.select("#svgdataurl").html(img);


    var canvas = document.querySelector("canvas"),
        context = canvas.getContext("2d");

    var image = new Image;
    image.src = imgsrc;
    image.onload = function() {
        context.fillStyle = 'white';
        context.fillRect(0, 0, canvas.width, canvas.height);
        context.drawImage(image, 0, 0);


        var canvasdata = canvas.toDataURL("image/png");

        var pngimg = '<img src="' + canvasdata + '">';
        d3.select("#pngdataurl").html(pngimg);

        var a = document.createElement("a");
        a.download = fname;
        a.href = canvasdata;
        document.body.appendChild(a);
        a.click();

        document.body.removeChild(iDiv1);
        document.body.removeChild(iDiv2);

    };

};


//save user results to server when exiting
function closingCode() {

/*
    if (typeof relative_dir !== "undefined" && relative_dir.indexOf("uploads") >= 0 && typeof picked_peaks !== "undefined" && picked_peaks.length > 0) {
        var data = '"picked_peaks":'.concat(JSON.stringify(picked_peaks));
        var limit = document.getElementById("contour-lower").value;
        $.ajax({ url: base_url.concat("/index.php/colmarm/savepeaks"), type: "post", data: { "limit": limit, "data": data, "relative_dir": relative_dir } });

    }*/
    
    if (typeof relative_dir !== "undefined" && relative_dir.indexOf("uploads") >= 0 && typeof compounds !== "undefined" && compounds.length > 0) {
        var t = new Object();
        t.compound = compounds;
        t.mytrace = mytrace;
        var data = JSON.stringify(t);
        $.ajax({ url: base_url.concat("/index.php/colmarm/savequery"), type: "post", data: { "data": data, "relative_dir": relative_dir } });
    }
    return null;
}


function show_one_compound_update(no_use) {
    no_use = 0;

    var t = compounds[current_compound];
    t.ngood = 0;
    t.nbad = 0;
    cc = 0;
    pp = 0;
    for (var i = 0; i < t.peaks.length; i++) {
        if (t.peaks[i].id > 0) {
            t.ngood++;
            var p_error = t.peaks[i].exp_x - t.peaks[i].x;
            var c_error = t.peaks[i].exp_y - t.peaks[i].y;
            cc = cc + c_error * c_error;
            pp = pp + p_error * p_error;
        } else {
            t.nbad++;
        }
    }
    t.matching = t.ngood / (t.ngood + t.nbad);
    t.c_rmsd = Math.sqrt(cc / t.ngood);
    t.p_rmsd = Math.sqrt(pp / t.ngood);


    //update uniqueness, uniqueness detail. amp_ratio
    for (var i = 0; i < compounds.length; i++) {
        var peaks = compounds[i]['peaks'];
        var unique = 0;
        var non_unique = 0;
        var max_a = -1e20;
        var min_a = 1e20;

        for (var j = 0; j < peaks.length; j++) {
            if (peaks[j].id > 0) {
                if (mytrace[peaks[j].id - 1].z.length > 1) {
                    non_unique++;
                    compounds[i]['peaks'][j].share = mytrace[peaks[j].id - 1].z.slice();
                } else
                    unique++;
                if (mytrace[peaks[j].id - 1].a < min_a)
                    min_a = mytrace[peaks[j].id - 1].a;
                if (mytrace[peaks[j].id - 1].a > max_a)
                    max_a = mytrace[peaks[j].id - 1].a;
            }
        }
        compounds[i]['unique'] = unique.toFixed(0) + "/" + compounds[i].ngood.toFixed(0);
        compounds[i]['amp_ratio'] = max_a / min_a;
    }

    //need to redo one row of main table, (uniqueness and amp_ratio has been changed





    show_one_compound_table(current_compound);
    if (typeof mainplot2 !== "undefined") {
        mainplot2.addpeaks(compounds[current_compound].peaks);
        mainplot2.cross_change(get_cross_select());
    }
    if (typeof mainplot3 !== "undefined") {
        mainplot3.addpeaks(compounds[current_compound].peaks);
        mainplot3.cross_change(get_cross_select());
    }
}

