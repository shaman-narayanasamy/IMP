// format time
String.prototype.toHHMMSS = function () {
    var sec_num = parseInt(this, 10); // don't forget the second param
    var hours   = Math.floor(sec_num / 3600);
    var minutes = Math.floor((sec_num - (hours * 3600)) / 60);
    var seconds = sec_num - (hours * 3600) - (minutes * 60);

    if (hours   < 10) {hours   = "0"+hours;}
    if (minutes < 10) {minutes = "0"+minutes;}
    if (seconds < 10) {seconds = "0"+seconds;}
    var time    = hours+'h:'+minutes+'m:'+seconds +"s";
    return time;
}

// time chart rendering
function renderTimeChart(data, title, key, selector, total_runtime){
    var converted = [];
    for(item in data){
        converted.push({
            label: item,
            value: data[item][key],
        });
    }
    var pie = new d3pie(selector, {
        header: {
            title: {
                text: title
            },
            subtitle: {
                text: (total_runtime).toString().toHHMMSS() + " total."
            },
            location: "pie-center",
        },
        data: {
            //sortOrder: "value-asc",
            content: converted
        },
        size: {
        canvasHeight: 500,
        canvasWidth: 500,
        pieInnerRadius: "80%",
        pieOuterRadius: null
      },
        labels: {
            outer: {
                format: "label",
                hideWhenLessThanPercentage: 0,
                pieDistance: 30
        },
            inner: {
                format: "percentage",
                hideWhenLessThanPercentage: 4
            },
            lines: {
                enabled: true,
                style: "curved",
            }
        }
    });
}



// bar chart rendrering
function renderBarChart(data, title, key, selector, total_runtime){
    var converted = [];
    var total = 0.0;
    for(item in data){
        converted.push({
            label: item,
            value: data[item][key],
        });
        total += data[item][key];
    }
    converted.sort(function(a, b){
        return b.value - a.value;
    });

    function percent(value){
        var v = (value * 100 / total).toFixed(0);
        if(v == 0){
            return "<1%";
        }
        return v + "%";
    }
    var margin = {top: 30, right: 50, bottom: 10, left: 400},
    width = 960 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

    var x = d3.scale.linear()
    .range([0, width])

    var y = d3.scale.ordinal()
    .rangeRoundBands([0, height], .2);

    var xAxis = d3.svg.axis()
    .scale(x)
    .orient("top");

    var svg = d3.select(selector).append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    x.domain(d3.extent(converted, function(d) { return d.value; })).nice();
    y.domain(converted.map(function(d) { return d.label; }));

    var bar = svg.selectAll(".bar")
    .data(converted)
    .enter().append("g")
    .attr("transform", function(d){
        return "translate(" + x(Math.min(0, d.value)) + "," + y(d.label) + ")"});

    bar.append("rect")
    .attr("class", "bar")
    .attr("width", function(d) { return Math.abs(x(d.value) - x(0)); })
    .attr("height", y.rangeBand());

    bar.append("text")
    .attr("x", -6)
    .attr("y", y.rangeBand())
    .attr("dy", ".35em")
    .style("text-anchor", "end")
    .text(function(d) { return d.label; });

    bar.append("text")
    .attr("x", function(d){
        return width + 12;
    })
    .attr("dy", ".35em")
    .style("text-anchor", "end")
    .text(function(d) { return percent(d.value); });


    // axis
    svg.append("g")
    .attr("class", "x axis")
    .call(xAxis);

    svg.append("g")
    .attr("class", "y axis")
    .append("line")
    .attr("x1", x(0))
    .attr("x2", x(0))
    .attr("y2", height);


    function type(d) {
        d.value = +d.value;
        return d;
    }
}


function renderFiles(data, selector, total_runtime){
    $('#basic-info').html("IMP has run in <span class='label label-success'>" + (total_runtime).toString().toHHMMSS()+ "</span> and <span class='label label-success'>" + Object.keys(data).length + " files were generated</span>");
}

function renderConfiguration(data, selector){
    var child = null, val = null, params = null, p = null;
    // rendering raws
    child = $("<li class='list-group-item'><ul class='list-group'><span class='label label-info'>Data</span></ul></li>");
    child.append("<li class='list-group-item'><b>Metagenomics</b> : " + data.raws.Metagenomics + "</li>");
    child.append("<li class='list-group-item'><b>Metatranscriptomics</b> : " + data.raws.Metatranscriptomics + "</li>");
    $(selector).append(child);

    // rendering non-params
    child = $("<li class='list-group-item'><ul class='list-group'><span class='label label-info'>General parameters</span></ul></li>");
    for(var gconf in data){
        if(gconf != 'raws'){
            val = data[gconf];
            if(typeof val !== 'object'){
                child.append("<li class='list-group-item'><b>" + gconf + "</b> : " + val + "</li>");
            }
        }
    }
    $(selector).append(child);

    // rendering nested params
    for(var gconf in data){
        if(gconf != 'raws'){
            val = data[gconf];
            if(typeof val === 'object'){
                child = $("<li class='list-group-item'><ul class='list-group'><span class='label label-info'>" + gconf + "</span></ul></li>");
                for(p in val){
                    child.append("<li class='list-group-item'><b>" + p + "</b> : " + val[p] + "</li>");
                }
                $(selector).append(child);
            }

        }
    }
}

if(IMP_CONFIG){
    renderConfiguration(IMP_CONFIG, '#configuration');
}


// runtime statistics
if(typeof IMP_STATS !== 'undefined') {
    //renderTimeChart(IMP_STATS.rules, '% runtime per task ', 'mean-runtime', 'time-charts-mean', IMP_STATS.total_runtime);
    renderFiles(IMP_STATS.files, '#files-per-samples', IMP_STATS.total_runtime);
    renderBarChart(IMP_STATS.rules, '% runtime per task ', 'mean-runtime', '#bar-chart', IMP_STATS.total_runtime);
} else {
    $('#imp-stats').remove();
}


// BUTTONS
$('#but-raw').on('change', function(){
    var val = $('input[name="raw-data"]:checked').val();
    $("#ifr-raw-stat").attr('src', val);
});
$('#but-preprocess').on('change', function(){
    var val = $('input[name="preprocess-data"]:checked').val();
    $("#ifr-preprocess-stat").attr('src', val);
});

$('#but-annot').on('change', function(){
    var val = $('input[name="annot-data"]:checked').val();
    $("#ifr-kronaplot").attr('src', val);
});


$('#but-assembly').on('change', function(){
    var val = $('input[name="assembly-data"]:checked').val();
    if(val.endsWith('png')){
        $("#assembly-wrapper").replaceWith("<div id='assembly-wrapper' class='row'><img class='img-responsive' src='" + val + "'/></div>");
    } else {
        $("#assembly-wrapper").replaceWith("<div id='assembly-wrapper' class='row'><iframe style='position: absolute; width: 100%;height: 100%; border: none' src='" + val + "'></iframe></div>");
    }
});


$('#but-mapping').on('change', function(){
    var val = $('input[name="mapping-data"]:checked').val();
    $("#mapping-wrapper").attr('src', val);
});

$('#but-ratio').on('change', function(){
    var val = $('input[name="ratio-data"]:checked').val();
    $("#ratio-wrapper").attr('src', val);
});

$('#but-variant').on('change', function(){
    var val = $('input[name="variant-data"]:checked').val();
    $("#variant-wrapper").attr('src', val);
});



// set last tab in cookie
$('#tabs a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
  Cookies('tab', $(e.target).attr('href'));
})
// show lat iopenned tab or first one
if(Cookies('tab')){
    $('#tabs li > a[href=' + Cookies('tab') + ']').tab('show')
} else {
    $('#tabs li:first > a').tab('show');
}

// load kronaplot if people click on tab. Make a lot of errors if loaded from start.
$('#annottab').on('show.bs.tab', function (e) {
    if($("#MG-kronaplot").attr('src') === undefined){
        $("#ifr-kronaplot").attr('src', 'Analysis/results/MG.gene_kegg_krona.html');
    }
});
