var koopa = null;
// template
function carousel_tmpl(selector){
    return '<div id="${selector}" class="carousel slide" data-interval="false"><ol class="carousel-indicators"></ol><div class="carousel-inner" role="listbox"></div></div>';
}



function simple_carousel(selector){
  return '<div class="col-sm-3 col-md-3"></div><div class="col-sm-6 col-md-6"><div id="${selector}" class="carousel slide" data-interval="false"><ol class="carousel-indicators"></ol><div class="carousel-inner" role="listbox"></div></div><a class="left carousel-control" role="button" data-slide="prev"><span class="glyphicon glyphicon-chevron-left" aria-hidden="true"></span><span class="sr-only">Previous</span></a><a class="right carousel-control" role="button" data-slide="next"><span class="glyphicon glyphicon-chevron-right" aria-hidden="true"></span><span class="sr-only">Next</span></a></div><div class="col-sm-3 col-md-3"></div>';
}



// define functions
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
    var child = null, val = null, params = null;
    var general_config = data['general'];
    child = $("<li class='list-group-item'><ul class='list-group'><span class='label label-info'>General</span></ul></li>");
    for(var gconf in general_config){
        if(gconf == 'raws'){
            val = general_config[gconf]['Metagenomics'];
            child.append("<li class='list-group-item'><b>Metagenomics</b> : " + val + "</li>");
            val = general_config[gconf]['Metatranscriptomics'];
            child.append("<li class='list-group-item'><b>Metatranscriptomics</b> : " + val + "</li>");
        } else {
            val = general_config[gconf];
            child.append("<li class='list-group-item'><b>" + gconf + "</b> : " + val + "</li>");
        }
    }
    $(selector).append(child);
    for(item in data){
        if(item != 'general'){
            child = $("<li class='list-group-item'><ul class='list-group'><span class='label label-info'>" + item + "</span></ul></li>");
            params = data[item];
            for(p in params){
                val = params[p];
                child.append("<li class='list-group-item'><b>" + p + "</b> : " + val + "</li>");
            }
            $(selector).append(child);
        }
    }
}

function renderStats(selector, path, files){
    var fpath = null;
    var node = $("<ul class='media-list'></ul>");
    for(var idx in files){
        fpath = path + '/' + files[idx];
        node.append("<li class='media'><img class='img-responsive img-thumbnail' src='" + fpath + "'></img></li>");
    }
    $(selector).append(node);
}

function renderStatsCarousel(selector, path, files){
    var fpath = null, dt = null, sl = null;
    var carousel = $(carousel_tmpl(selector));
    var slides = $(carousel).find(".carousel-inner");
    for(var idx in files){
        fpath = path + '/' + files[idx];
        sl = $("<div class='item'><img class='img-responsive img-thumbnail' src='" + fpath + "'></img></div>");
        if(idx == 0){
            sl.addClass("active");
        }
        slides.append(sl);
    }
    $('#' + selector + "-wrapper").append(carousel);
}

function renderSimpleCarousel(selector, path, files){
  var fpath = null, dt = null, sl = null;
  var carousel = $(simple_carousel(selector));
  var slides = $(carousel).find(".carousel-inner");
  for(var idx in files){
      fpath = path + '/' + files[idx];
      sl = $("<div class='item'><img class='img-responsive img-thumbnail' src='" + fpath + "'></img></div>");
      if(idx == 0){
          sl.addClass("active");
      }
      slides.append(sl);
  }
  $('#' + selector + "-wrapper").append(carousel);
}



// Trigger functions


// configuration
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
    $('#imp-stats-alert').show();
//renderTimeChart(data['rules'], 'Max runtime', 'max-runtime', 'time-charts-max');
}


// ht_stats
renderStatsCarousel("carousel-mg-raw", "Preprocessing/stats/MG", [
    "cycle_composition_2.png",
    "cycle_quality_1.png",
    "cycle_quality_2.png",
    "cycle_quality_box_1.png",
    "cycle_quality_box_2.png",
    "lane_tile_quality_1.png",
    "lane_tile_quality_2.png",
    "quality_QQ.png",
    "reads_length.png",
    "reads_quality.png"
]);

renderStatsCarousel("carousel-mg-preprocessed", "Preprocessing/stats_after_preprocessing/MG", [
    "cycle_composition_2.png",
    "cycle_quality_1.png",
    "cycle_quality_2.png",
    "cycle_quality_box_1.png",
    "cycle_quality_box_2.png",
    "lane_tile_quality_1.png",
    "lane_tile_quality_2.png",
    "quality_QQ.png",
    "reads_length.png",
    "reads_quality.png"
]);

renderStatsCarousel("carousel-mt-raw", "Preprocessing/stats/MT", [
   "cycle_composition_2.png",
   "cycle_quality_1.png",
   "cycle_quality_2.png",
   "cycle_quality_box_1.png",
   "cycle_quality_box_2.png",
   "lane_tile_quality_1.png",
   "lane_tile_quality_2.png",
   "quality_QQ.png",
   "reads_length.png",
   "reads_quality.png"
]);

renderStatsCarousel("carousel-mt-preprocessed", "Preprocessing/stats_after_preprocessing/MT", [
    "cycle_composition_2.png",
    "cycle_quality_1.png",
    "cycle_quality_2.png",
    "cycle_quality_box_1.png",
    "cycle_quality_box_2.png",
    "lane_tile_quality_1.png",
    "lane_tile_quality_2.png",
    "quality_QQ.png",
    "reads_length.png",
    "reads_quality.png"
]);

renderSimpleCarousel('carousel-assembly', 'Analysis/results', [
    'IMP-vizbin_length.png',
    'IMP-vizbin_length_GC.png'
]);

renderSimpleCarousel('carousel-mapping', 'Analysis/results', [
    'IMP-reads_density.png',
    'IMP-rpkm_density.png',
    'IMP-coverage_density.png',
    'IMP-depth_density.png',
    'IMP-vizbin_length_MGcov.png',
    'IMP-vizbin_length_MTcov.png',
    'IMP-vizbin_length_MGdepth.png',
    'IMP-vizbin_length_MTdepth.png'
]);

renderSimpleCarousel('carousel-ration', 'Analysis/results', [
    'IMP-vizbin_length_depthRatio.png',
    'IMP-vizbin_length_rpkmRatio.png'
]);

renderSimpleCarousel('carousel-variant', 'Analysis/results', [
    'IMP-var_count.png',
    'IMP-var_density.png',
    'IMP-vizbin_length_MGvardens.png',
    'IMP-vizbin_length_MTvardens.png'
]);

// hide sections
//$("section.row").hide();
// $("section.row").prev("h2").click(function(){
//     $(this).find(".glyphicon").toggleClass("glyphicon-plus glyphicon-minus");
//     $(this).next("section").toggle();
// });

//carousel
$('.left.carousel-control').click(function(){
    $('.carousel').carousel('next');
});
$('.right.carousel-control').click(function(){
    $('.carousel').carousel('prev');
});

// show first tab
$('#tabs li:first > a').tab('show');

//load html tables
$(".table-result").each(function(idx, node) {
  var ident = node.id;
  //console.log(ident);
  $("#" + ident).load("MGMT/results/" + ident + ".html", function() {
    //console.log(ident);
  });
});
