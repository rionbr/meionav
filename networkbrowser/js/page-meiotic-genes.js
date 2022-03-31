var file = document.getElementById('canvas').getAttribute("data-file")
var data;
var debug;

/*
 * Toggle Modal
 */
function toggleModal(title, file) {
    console.log("togglemodal:",title, file);
    $('#modal-title').text(title);
    $('#modal-img').attr("src", file);
    $("#modal").modal();
    return false
}

function loadThumbAndSetClickActions(tag, title, image, thumb) {
    // Thumb
    $(tag).load(thumb, function(response, status, xhr) {
        if (status == "error") {
            $(this).attr("src", "img/results/not-available_240x240.jpg");
            // Remove action
            $(tag).off("click");
        } else {
            $(this).attr("src", thumb);
            // Add action
            $(tag).on("click", function() {
                toggleModal(title, image)
            });
        }
    })
}

/*
 * Update Gene Information Table
 */
function updateGeneInformation(id) {
    gene = data[id]['id_gene'];
    name = data[id]['gene'];
    stock = data[id]['stock'];
    ortho_hs = data[id]['orthologs-HS']
    ortho_hs_str = ortho_hs.map(d => d['gene']).join(', ')
    ortho_mm = data[id]['orthologs-MM']
    ortho_mm_str = ortho_mm.map(d => d['gene']).join(', ')
    // Name
    $("#gene-name").html(name);
    $("#gene-stock").html(stock);
    $("#gene-flybase").attr("class", "btn btn-warning").attr("href", `https://flybase.org/reports/${gene}`);
    // Orthologs
    $("#gene-orthologs-human").html(ortho_hs_str);
    $("#gene-orthologs-mouse").html(ortho_mm_str);
    // Phenotype
    path = `img/results/screened/${stock}/`
    var types = ['Apical testis', 'Median testis', 'Distal testis', 'Gametes']
    types.map((x, i) => {
        i = i + 1;
        urlId = `#pheno-url-${i}`;
        imgId = `#pheno-img-${i}`;
        //
        image = path +`${stock}_${i}.jpg`;
        thumb = path + `${stock}_${i}_thumb.jpg`;
        title = `${name} (${stock}) - ${x}`
        // Load thumb and set click action
        loadThumbAndSetClickActions(`#pheno-img-${i}`, title, image, thumb);
    })
    // FlyBase API Query
    requestFlyBaseSummary(gene);
}

/*
 * Async Request FlyBase Summary
 */
async function requestFlyBaseSummary(gene) {
    /* FlyBase API Request */
    d3.select("#gene-summary").text("Requesting summary to FlyBase.")
    url = `https://api.flybase.org/api/v1.0/gene/summaries/auto/${gene}`
    fbdata = await d3.json(url)
    summary = fbdata.resultset.result[0].summary
    d3.select("#gene-summary").text(summary)
}

/*
 * Load JSON Results
 */ 
async function loadData() {
    // load
    data = await d3.json(file)
    // sort
    data = data.slice().sort((a, b) => d3.ascending(a['mean(fert-rate)'], b['mean(fert-rate)']))
    // Map ID
    data.map((d,i) => d['id'] = i)
}

loadData().then(() => {
    plot_results();
})


var colors_fpkm = d3.scaleLog()
    .domain([1, 4000])
    .range(["white","red"]);


/*
 * Buils D3 Viz
 */
function plot_results() {
    width = 960;
    heightBand = 80;
    heightChart = 350;
    heightFocus = 80;
    margin = {top: 5, right: 20, bottom: 38, left: 90}
    curSelection = 155;
    
    // Filter that is being displayed
    numberOfPoints = 50
    /*
     * x, y scales
     */
    xBand = d3.scaleLinear()
        .domain([-.5, data.length])
        .range([margin.left, width - margin.right])
        .nice()

    xChart = d3.scaleLinear()
        .domain([-.5, data.length])
        .range([margin.left, width - margin.right])
        .nice()

    yChart = d3.scaleLinear()
        .domain([-0.03, 1.05])
        .range([heightChart - margin.bottom, margin.top])
        //.nice()

    xFocus = d3.scaleLinear()
        .domain([-.5, data.length])
        .range([margin.left, width - margin.right])
        //.nice()

    yFocus = d3.scaleLinear()
        .domain([0, 1])
        .range([heightFocus - margin.bottom, margin.top])
        //.nice()

    /*
     * x, y axes
     */
    xAxisChart = (g, x) => g.call(d3.axisBottom(x))
    yAxisChart = (g, y) => g.call(d3.axisLeft(y))
    y2AxisChart = (g, y) => g.call(d3.axisRight(y))
    xAxisFocus = (g) => g.call(d3.axisBottom(xFocus))
    yAxisFocus = (g) => g.call(d3.axisLeft(yFocus).ticks(3))

    // axes for Grid Lines
    xAxisChartGrid = (g, x) => g.call(d3.axisBottom(x)
        .ticks(5)
        .tickSize( - heightChart + margin.top + margin.bottom)
        .tickFormat(""))
    yAxisChartGrid = (g, x) => g.call(d3.axisLeft(x)
        .ticks(5)
        .tickSize( - width + margin.left + margin.right)
        .tickFormat(""))

    /*
     * Area
     */
    area = (x, y) => d3.area()
      .defined(d => !isNaN(d['mean(fert-rate)']))
      .x(d => x(d['id']))
      .y0(y(0))
      .y1(d => y(d['mean(fert-rate)']))

    /*
     * Mouse over
     */
    handleMouseOut = (d, i) => {
    	// Change radius size
        id = i['id']
        d3.select(`#circle-${id}`).node().r.baseVal.value = 6
    }
    handleMouseOver = (d, i) => {
    	id = i['id']
    	// Change radius size
    	d3.select(`#circle-${id}`).node().r.baseVal.value = 12
    }
    /*
     * Mouse Click
     */
    handleClick = (d) => {
        debug = d
        //id = d.target.id
        id = d.target.id.split("-")[2]

        updateGeneInformation(id);

        // Update Current Selection
        curSelection = id;
        // Move selection
        d3.select("#bar")
            .transition().duration(500)
            .attr("x", xChart(curSelection))

    }

    handleSelector = (e) => {
        debug = e;
        id = Number(d3.select(e.target).property("value"));
        updateGeneInformation(id)

        // Update Current Selection
        curSelection = id;
        // Move selection
        d3.select("#bar")
            .transition().duration(500)
            .attr("x", xChart(curSelection))

        // Move the brush (and window)
        if (id < 50) {
            minX = -1;
            maxX = 50;
        } else if (id > 870) {
            minX = 870;
            maxX = 920;
        } else {
            minX = id - 25;
            maxX = id + 25;
        }
        console.log(id,minX,maxX);
        gb.call(brush.move, [xFocus(minX), xFocus(maxX)]);
    }
    /*
     * Focus
     */
    svgFocus = d3.select("#canvas").append("svg")
      .attr("id", "focus")
      .attr("viewBox", [0, 0, width, heightFocus])
      .style("display", "block");

    // y-axis
    svgFocus.append("g")
      .attr("class", "y-axis focus")
      .attr("transform", `translate(${margin.left},0)`)
      .call(yAxisFocus);

    // x-axis
    svgFocus.append("g")
      .attr("class", "x-axis focus")
      .attr("transform", `translate(0,${heightFocus - margin.bottom})`)
      .call(xAxisFocus);

    // area
    svgFocus.append("g")
        .attr("class", "areas")
        .append("path")
        .attr("class", "area")
        .datum(data)
        .attr("d", area(xFocus, yFocus));
        
    // annotation
    svgFocus.append("text")
        .attr("class", "annotation")
        .attr("y", margin.top + 23)
        .attr("x", margin.left + 14)
        .text("Drag the gray box");

    /*
     * Chart
     */
    svgChart = d3.select('#canvas').append('svg')
      .attr("id", "chart")
      .attr("viewBox", [0, 0, width, heightChart])
      .style("display", "block");

    svgChart.append("clipPath")
        .attr("id", "clipChart")
      .append("rect")
        .attr("x", margin.left)
        .attr("y", margin.top)
        .attr("width", width - margin.left - margin.right)
        .attr("height", heightChart - margin.top - margin.bottom)

    // x-grid
    svgChart.append("g")
        .attr("class", "grid x-grid")
        .attr("transform", `translate(0, ${heightChart - margin.bottom})`)
        .call(xAxisChartGrid, xChart)

    // y-grid
    svgChart.append("g")
        .attr("class", "grid y-grid")
        .attr("transform", `translate(${margin.left},0)`)
        .call(yAxisChartGrid, yChart)

    // x-axis
    svgChart.append("g")
        .attr("class", "axis x-axis")
        .attr("transform", `translate(0,${heightChart - margin.bottom})`)
        .call(xAxisChart, xChart);

    // y-axis
    svgChart.append("g")
      .attr("class", "axis y-axis")
      .attr("transform", `translate(${margin.left},0)`)
      .call(yAxisChart, yChart);

    // x-label
    svgChart.append("text")
        .attr("class", "x-label")
        //.attr("transform", "")
        .attr("y", margin.top + heightChart)
        .attr("x", width / 2)
    	.attr("dy", "-.4em")
        .style("text-anchor", "middle")
        .text("Gene rank");

    // y-label
    svgChart.append("text")
        .attr("class", "y-label")
        .attr("transform", "rotate(-90)")
        .attr("y", 0 - margin.left)
        .attr("x", 0 - (heightChart / 2))
        .attr("dy", "9em")
        .style("text-anchor", "middle")
        .text("Fertility rate");

    // area
    svgChart.append("g")
        .attr("class", "areas")
        .attr("clip-path", "url(#clipChart)")
        // apth
        .append("path")
        .attr("class", "area")
        .datum(data)
        .attr("d", area(xChart, yChart));

    // selection-bar
    svgChart.append("g")
        .attr("class", "selection")
        .attr("clip-path", "url(#clipChart)")
        .append("rect")
        .attr("id", "bar")
        .attr("width", 18)
        .attr("height", heightChart - margin.top - margin.bottom)
        .attr("x", 0)  // start at zero
        .attr("y", margin.top)
        .attr("transform", `translate(-9,0)`)

    // circles
    svgChart.append("g")
        .attr("class", "items")
        .attr("clip-path", "url(#clipChart)")
        //
        .selectAll("circle")
        .data(data)
        .enter()
        //
        .append("circle")
        .attr("id", d => `circle-${d['id']}`)
        .attr("class", d => "item code-" + d['our-DM-code'])
        .attr("r", 6)
        .attr("cx", d => xChart(d['id']))
        .attr("cy", d => yChart(d['mean(fert-rate)']))

    // names
    svgChart.append("g")
        .attr("class", "names")
        .attr("clip-path", "url(#clipChart)")
        .selectAll("text")
        .data(data)
        .enter()
        //
        .append("text")
        .attr("class", "name")
        .attr("dominant-baseline", "middle")
        .attr("text-anchor", d => d['mean(fert-rate)'] < 0.75 ? "start" : "end") // if <0.75 ? left : right
        .attr("x", d => d['mean(fert-rate)'] < 0.75 ? xChart(d['id'] + 0.6) : xChart(d['id'] - 0.6))
        .attr("y", d => yChart(d['mean(fert-rate)']))
        .attr("transform", d => "rotate(-90," + xChart(d['id']) + "," + yChart(d['mean(fert-rate)']) + ")")
        .text(d => d['gene'])

	// interaction boxes
    svgChart.append("g")
        .attr("class", "ui-rects")
        .attr("clip-path", "url(#clipChart)")
        //
        .selectAll("rect")
        .data(data)
        .enter()
        //
        //
        .append("rect")
        .attr("id", d => `ui-rect-${d['id']}`)
        .attr("class", "ui-rect")
        .attr("width", 12)
        .attr("height", d => 10 * d['gene'].length)
        .attr("x", d => xChart(d['id']))
        .attr("y", d => yChart(d['mean(fert-rate)']))
        .attr("fill", "black")
        .attr("fill-opacity", 0.0)
        .attr("transform", d => d['mean(fert-rate)'] < 0.75 ? `translate(-6,-${(10 * d['gene'].length) - 6})` : `translate(-6,-${d['gene'].length + 3})`)
        //
        .on("mouseover", handleMouseOver)
        .on("mouseout", handleMouseOut)
        .on("click", handleClick)
	/*    
    // errorbars
    svgChart.append("g")
        .attr("class", "errorbars")
        .attr("clip-path", "url(#clipChart)")
        .selectAll("line")
        .data(data)
        .enter()
        //
        .append("line")
        .attr("class", "errorbar")
        .attr('x1', d => xChart(d['id']))
        .attr('x2', d => xChart(d['id']))
        .attr('y1', d => yChart(d['mean(fert-rate)'] - d['std(fert-rate)']))
        .attr('y2', d => yChart(d['mean(fert-rate)'] + d['std(fert-rate)']));
	*/

    /*
     * Band
     */
    svgBand = d3.select('#canvas').append('svg')
      .attr("id", "band")
      .attr("viewBox", [0, 0, width, heightBand])
      .style("display", "block");

    svgBand.append("clipPath")
        .attr("id", "clipBand")
      .append("rect")
        .attr("x", margin.left)
        .attr("width", width - margin.left - margin.right)
        .attr("height", heightBand)

    // FPKM rects
    svgBand.append("g")
        .attr("class", "items")
        .attr("clip-path", "url(#clipBand)")
        //
        .selectAll("rect")
        .data(data)
        .enter()
        //
        .append("rect")
        .attr("class", "item FPKM")
        .attr("width", 12)
        .attr("height", 12)
        .attr("x", d => xChart(d['id']))
        .attr("y", 16)
        .attr("fill", d => colors_fpkm(d['FPKM']))
        .attr("transform", `translate(-6,-6)`)

    // Ext-DM-Code
    svgBand.append("g")
        .attr("class", "items")
        .attr("clip-path", "url(#clipBand)")
        //
        .selectAll("rect")
        .data(data)
        .enter()
        .filter(function(d) { return d['ext-HS-code'] !== null; })
        //
        .append("rect")
        .attr("class", d => "item code-true") /* + d['ext-HS-code']) */
        .attr("width", 12)
        .attr("height", 12)
        .attr("x", d => xChart(d['id']))
        .attr("y", 22+2+12)
        .attr("transform", `translate(-6,-6)`)

    // Ext-MM-Code
    svgBand.append("g")
        .attr("class", "items")
        .attr("clip-path", "url(#clipBand)")
        //
        .selectAll("rect")
        .data(data)
        .enter()
        .filter(function(d) { return d['ext-MM-code'] !== null; })
        //
        .append("rect")
        .attr("class", d => "item code-true") /* + d['ext-MM-code']) */
        .attr("width", 12)
        .attr("height", 12)
        .attr("x", d => xChart(d['id']))
		.attr("y", 22+2+12+2+12)
        .attr("transform", `translate(-6,-6)`)

    // Ext-HS-Code
    svgBand.append("g")
        .attr("class", "items")
        .attr("clip-path", "url(#clipBand)")
        //
        .selectAll("rect")
        .data(data)
        .enter()
        .filter(function(d) { return d['ext-DM-code'] !== null; })
        //
        .append("rect")
        .attr("class", d => "item code-true") /* + d['ext-DM-code']) */
        .attr("width", 12)
        .attr("height", 12)
        .attr("x", d => xChart(d['id']))
        .attr("y", 22+2+12+2+12+2+12)
        .attr("transform", `translate(-6,-6)`)

    // Labels
    svgBand.append("text")
        .attr("class", "label fpkm")
        .attr("x", 82)
        .attr("y", 20)
        .attr("text-anchor", "end")
        .text("Log(FPKM+1)");

    svgBand.append("text")
        .attr("class", "label ext-hs-code")
        .attr("x", 80)
        .attr("y", 40)
        .attr("text-anchor", "end")
        .text("HS pheno.");

    svgBand.append("text")
        .attr("class", "label ext-mm-code")
        .attr("x", 80)
        .attr("y", 55)
        .attr("text-anchor", "end")
        .text("MM pheno.");

    svgBand.append("text")
        .attr("class", "label ext-dm-code")
        .attr("x", 80)
        .attr("y", 70)
        .attr("text-anchor", "end")
        .text("DM pheno.");

    /*
     * Populate the Select
     */
    selector = d3.select("#selector");
    selector.selectAll("option")
        .remove(); // first remove the loading

    selector.selectAll("option")
        .data(data.slice().sort((a, b) => d3.ascending(a['gene'], b['gene'])))
        .enter()
            .append("option")
            .attr("id", d => `option-${d["id"]}`)
            .attr("value", d => d["id"])
            .text(d => d["gene"]);

    d3.select("#option-155").attr("selected", true); // Initial state on CG9776
    selector.on("change", handleSelector);

	/*
	 * Focus Actions
	 */
    brush = d3.brushX()
      .extent([[margin.left, margin.top - 0.5], [width - margin.right, heightFocus - margin.bottom + 0.5]])
      .on("end", brushed);

    const defaultSelection = [xFocus(140), xFocus(190)];

    const gb = svgFocus.append("g")
      .call(brush)
      .call(brush.move, defaultSelection);

    updateGeneInformation(155);

    function brushed({selection}) {
        if (selection) {
            // Reset domain
            [minXFocus, maxXFocus] = selection
            xChart.domain(selection.map(xFocus.invert, xFocus));
            [minXChart, maxXChart] = xChart.domain()

            xBand.domain(selection.map(xFocus.invert, xFocus));

            // Update Area
            svgChart.selectAll(".area")
                .transition().duration(400)
                .attr("d", area(xChart, yChart));

            // Update Circles
            svgChart.selectAll(".item")
                .transition().duration(450)
                .attr("cx", d => xChart(d['id']))
                //.attr("cy", d => yChart(d['mean(fert-rate)']))

            // Update UI-Rects
            svgChart.selectAll(".ui-rect")
                .transition().duration(450)
                .attr("x", d => xChart(d['id']))

            // Update Text
            svgChart.selectAll(".name")
                .transition().duration(450)
                .attr("x", d => d['mean(fert-rate)'] < 0.75 ? xChart(d['id'] + 0.6) : xChart(d['id'] - 0.6))
                //.attr("y", d => yChart(d['mean(fert-rate)']))
                .attr("transform", d => "rotate(-90," + xChart(d['id']) + "," + yChart(d['mean(fert-rate)']) + ")")

            /*
            // Update errorbars
            svgChart.selectAll(".errorbar")
                .transition().duration(450)
                .attr("x", d => xChart(d['id']))
            */

            // Update selection-bar
            svgChart.selectAll("#bar")
                .transition().duration(450)
                .attr("x", d => xChart(curSelection))

            // Update x-axis
            svgChart.selectAll(".x-axis")
                .transition().duration(300)
                .call(xAxisChart, xChart)

            // Update y-axis
            svgChart.selectAll(".y-axis")
                .transition().duration(300)
                .call(yAxisChart, yChart)

            // Update x-grid
            svgChart.selectAll(".x-grid")
                .transition().duration(300)
                .call(xAxisChartGrid, xChart)

            // Update y-grid
            svgChart.selectAll(".y-grid")
                .transition().duration(300)
                .call(yAxisChartGrid, yChart)

            // Update Band
            svgBand.selectAll(".item")
                .transition().duration(450)
                .attr("x", d => xBand(d['id']))
                .attr("cx", d => xBand(d['id']))
    

        } else {
            gb.call(brush.move, defaultSelection);
        }
    }
	/* 
	 * Legend
	 */
	widthLegend = 300;
	heightLegend = 50;
	const svg = d3.select("#legend").append("svg")
	  .attr("width", widthLegend)
	  .attr("height", heightLegend)
	  .attr("viewBox", [0, 0, widthLegend, heightLegend])
	  /*.style("overflow", "visible")*/
	  .style("display", "block");

	svg.append("g")
	  .attr("class", "legend")
	  .attr("transform", "translate(1,10)");

	var legend = d3.legendColor()
	  .cells([1, 10, 100, 1000, 4000])
	  //.labels(['10e1', '10e2', '10e3', '40e3'])
	  .labelFormat(".0e")
	  .shapeWidth(30)
	  .orient("horizontal")
	  .shapePadding(10)
	  .scale(colors_fpkm);

	svg.select(".legend")
		.call(legend);
}