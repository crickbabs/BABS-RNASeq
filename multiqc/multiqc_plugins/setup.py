#!/usr/bin/env python

from setuptools import setup, find_packages

# jinja templates
html_files = ["multiqc_plugins/templates/custom/includes.html"]

# the javascript libraries. The highcharts-more is not necessary at the moment,
# but we'll se in the future
js_files = [
	"multiqc_plugins/templates/custom/assets/js/highcharts-more.v5.0.6.js",
	"multiqc_plugins/templates/custom/assets/js/plotly-latest.min.js"
	]

setup(
		name = "multiqc_plugins",
		version = "0.1",
		author = "nourdine bah",
		author_email = "nourdine.bah@crick.ac.uk",
		description = "MultiQC plugins for special metrics",
		packages = find_packages(),
		include_package_data = True,
		install_requires = ["multiqc>=1.5"],
		data_files=[("html", html_files), ("js", js_files)],
		entry_points = {
			"multiqc.modules.v1": [
				"pca = multiqc_plugins.modules.pca:MultiqcModule"
				],
			"multiqc.templates.v1": [
				"custom = multiqc_plugins.templates.custom"
				]
			}
		)

