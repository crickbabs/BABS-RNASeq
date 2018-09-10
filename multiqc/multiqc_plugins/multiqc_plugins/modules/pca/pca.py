#!/usr/bin/env python

from __future__ import print_function
from collections import OrderedDict
import logging
import json
#from unicodedata import normalize

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# the multiqc logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):

	def __init__(self):

		super(MultiqcModule, self).__init__(
				name = "Principal Component Analysis",
				target = "",
				anchor = "pca",
				href = "",
				info = ""
				)

		self.data = list()
		
		#########################################################################
		for i, f in enumerate( self.find_log_files("pca") ):

			# the pca as a json
			j = json.loads(f['f'])

			# percentage of variance
			xtitle = "PC1 (" + str(j["percent_var"][0]) + " %)"
			ytitle = "PC2 (" + str(j["percent_var"][1]) + " %)"

			######################################################################
			# transform points into plotly traces

			points = list()

			for point in j["pca"]:

				d = dict()

				d["mode"] = "markers+text"
				d["type"] = "scatter"
				d["name"] = point["name"]
				d["x"] = [point["PC1"]]
				d["y"] = [point["PC2"]]
				d["text"] = [point["name"]]
				d["textposition"] = "bottom center"
				d["hovertext"] = point["hover"]

				###################################################################
				# marker features

				marker = dict()
				marker["size"] = 12

				if "color" in point:
					marker["color"] = point["color"]

				if "symbol" in point:
					marker["symbol"] = point["symbol"]

				if "line" in point:
					marker["line"] = {
							"color": point["line"],
							"width": 5
					}

				d["marker"] = marker
				###################################################################

				points.append(d)
			######################################################################

			######################################################################
			# PLOT

			description = """PC1 against PC2 plot. The experimental variables are
			represented by the color, the shape and the contour line of the
			markers."""

			mailto = """<a href="mailto:bioinformatics@crick.ac.uk">
			bioinformatics@crick.ac.uk</a>"""

			self.add_section(
					description = description,
					helptext = mailto,
					plot = self.pca_plot(i+1, points, xtitle, ytitle)
					)

		log.info( "Found {} reports".format(str(i+1)) )

	############################################################################
	def pca_plot(self, num, data, xtitle, ytitle):
		
		label = "pca_plot_" + str(num)
		
		html = """
		<div
			id='""" + label + """'
			style="width:100%"
			class="pyplot-plot"
		></div>

		<script type="text/javascript">

			var label = '""" + label + """';

			var data = """ + json.dumps(data)  + """;

			var layout = {
				xaxis: { title: '""" + xtitle + """' },
				yaxis: { title: '""" + ytitle + """' },
				showlegend: false,
				hovermode: "closest"
			};

			Plotly.newPlot(label, data, layout);

		</script>
		"""

		return html

	############################################################################
	def test_plot(self):
		
		html = """
		<div id="pca_test_plot" style="width:100%" class="pyplot-plot"></div>
		<script type="text/javascript">

			var trace1 = {
				x: [1, 2, 3, 4, 5],
				y: [1, 6, 3, 6, 1],
				mode: 'markers',
				type: 'scatter',
				name: 'Team A',
				text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
				marker: {
					symbol: 12,
					size: 12,
					line: {
						color: "#FF0000",
						width: 5
					},
				}
			};

			var trace2 = {
				x: [1.5, 2.5, 3.5, 4.5, 5.5],
				y: [4, 1, 7, 1, 4],
				mode: 'markers',
				type: 'scatter',
				name: 'Team B',
				text: ['B-a', 'B-b', 'B-c', 'B-d', 'B-e'],
				marker: { size: 12 }
			};

			var data = [ trace1, trace2 ];

			var layout = {
				xaxis: {
					range: [ 0.75, 5.25 ],
					yaxis: { range: [0, 8] },
					title: 'Data Labels Hover'
				}
			};

			Plotly.newPlot('pca_test_plot', data, layout);

		</script>
		"""

		return html

