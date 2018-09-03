# coding: utf-8

# nourdine.bah@crick.ac.uk

import sys
import pandas as pd
from os.path import abspath
import pydot
from pydot import Edge
from pydot import Graph
from pydot import Node
from pydot import Dot
import argparse

###
def get_label_map(graph, labels=None):
    """Returns a dictionary that associates the node names to their label."""

    # the dictionary to fill
    d = {}

    # the all the labels for all nodes
    for n in graph.get_node_list():
        name = n.get_name()
        label = str(n.get_label()).replace('"', '')
        if name in d.keys():
            d[name].add(label)
        else:
            d[name] = set([label])

    # check if there is multiple labels for a node
    x = list(map( lambda x: len(x)!=1, [l for n, l in d.items()] ))
    assert True not in x, "Problem, there are multiples labels for a node"
    for n, s in d.items():
        d[n] = list(s)[0]

    # filter on the selected label
    if labels != None:
        d = dict( [(n, l) for n, l in d.items() if l in labels] )

    return d


###
def get_successors(graph):
    """Returns a dict of all successors of each node."""
    d = {}
    for e in graph.get_edge_list():
        src = e.get_source()
        dst = e.get_destination()
        if src in d.keys():
            d[src].add(dst)
        else:
            d[src] = set([dst])
    return d


###
def get_paths(node, successors, path, paths):
    """Returns all the paths starting from a particular node."""
    path.append(node)
    if node not in successors.keys():
        paths.add(tuple(path))
        path.pop()
    else:
        for successor in successors[node]:
            get_paths(successor, successors, path, paths)
        path.pop()


###
def get_map_from_df(df, key, value):
    return pd.Series(df[value].values, index=df[key]).fillna("").to_dict()


###############################################################################

description = """
Takes a nextflow output dot file and rebuilds it by just keeping the nodes
speficied in the design csv file. If you're not specifying --pdf or --png the
script will write the dot code on the standard output.
"""

if __name__ == "__main__":
	
	# arguments of the scripts
	ap = argparse.ArgumentParser(description=description)

	# mandatory arguments
	ap.add_argument("dot", help="nextflow output dot file path")
	ap.add_argument("csv", help="design csv path")

	# export arguments
	ap.add_argument("--png",
			help="path of the png that will be exported",
			required=False
			)
	ap.add_argument("--pdf",
			help="path of the pdf that will be exported",
			required=False
			)
	ap.add_argument("--svg",
			help="path of the svg that will be exported",
			required=False
			)

	# parsing
	args = ap.parse_args()

	# the design file
	design = pd.read_csv(args.csv)
	
	# the dot file
	(g,) = pydot.graph_from_dot_file(args.dot)
	
	# extract the selected dot node names by their labels
	labels = design["label"].tolist()
	name_to_label = get_label_map(g, labels)
	names = sorted(list( name_to_label.keys() ))
	
	# add the dot node names to the design
	tmp_df = pd.DataFrame.from_records(
	        list(tuple( name_to_label.items() )),
	        columns=["name", "label"]
	        )
	design = design.merge(tmp_df, left_on="label", right_on="label")
	
	# useful maps
	name_to_color = get_map_from_df(design, "name", "color")
	name_to_newlabel = get_map_from_df(design, "name", "new_label")
	name_to_software = get_map_from_df(design, "name", "software")
	
	# the successor of each node of the graph
	successors = get_successors(g)
	
	# associate all the paths starting from each node
	name_to_paths = {}
	for name in names:
	    paths = set()
	    get_paths(name, successors, [], paths)
	    name_to_paths[name] = paths
	
	# create a new graph and add the nodes to it
	graph = Dot()
	for name in names:
	
	    # the label of the node
	    if name_to_newlabel[name] != "":
	        label = str(name_to_newlabel[name]) + "\n" + name_to_software[name]
	    else:
	        label = name_to_software[name]
	
	    graph.add_node(
	            Node(
	                name,
	                fillcolor=name_to_color[name],
	                color=name_to_color[name],
	                style="filled",
	                fontname="cantarell bold",
	                fontcolor="#ffffff",
	                fontsize=22,
	                label=label
	                )
	            )
	
	# get all the edges connected two processes
	edges = set()
	for name, paths in name_to_paths.items():
	    for path in paths:
	        for descendant in path:
	            if (descendant in names) and (descendant != name):
	                edges.add( (name, descendant) )
	                break
	
	# add the edges to the graph
	for edge in edges:
	    graph.add_edge(
	            Edge(
	                edge[0],
	                edge[1],
	                penwidth=4,
	                color=name_to_color[ edge[0] ]
	                )
	            )
	
	#graph.set("size", '"10,7"')
	#graph.set("ratio", '"fill"')
	#graph.set("ratio", "1")
	#graph.set_simplify(True)
	
	# save the graph as a images or print the dot code
	if args.pdf != None:
		graph.write_pdf(args.pdf)
	if args.png != None:
		graph.write_png(args.png)
	if args.svg != None:
		graph.write_svg(args.svg)
	if args.pdf == None and args.png == None and args.svg == None:
		print(graph.to_string())

