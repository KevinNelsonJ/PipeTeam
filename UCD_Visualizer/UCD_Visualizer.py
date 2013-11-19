# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>


from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.cm import rainbow
from itertools import chain

import matplotlib.pyplot as plt
import numpy as np
import math
import sys


from collections import defaultdict

# a little trick to make a recursive default dict.  This structure autmatically
# intializes as many nested dictionaries as needed by the call.
rdefaultdict = lambda:defaultdict(rdefaultdict)

def convert(arry, args):
    # takes a list and casts the elements in that list according to the types given
    return [f(s) for f,s in zip(args,arry)]

class UCD_Parser:

    def __init__(self):
        
        self.entry_types = {
            'node':(int,float,float,float),
            'line':(int,int,str,int,int),
            'tri' :(int,int,str,int,int,int),
            'quad':(int,int,str,int,int,int,int),
            'hex' :(int,int,str,int,int,int,int,int,int,int,int)
        }
        
        self.state = 'initial'
        self.values = {}
        self.ctr = 0
        
        self.nodes = {}
        self.cells = rdefaultdict()

    def make_transition(self,index,line):
        
        if self.state == 'initial':
            self.state = 'comments'
            
        elif self.state == 'comments' and line.lstrip()[0] != '#':
            self.ctr = index
            self.state = 'header'
            
        elif self.state == 'header' and index == self.ctr+1:  
            self.ctr = index
            self.state = 'nodes'
            
        elif self.state == 'nodes' and index == self.ctr+self.values['numNodes']:
            self.ctr = index
            self.state = 'cells'
    
    def transition(self,index,line):
        self.index = index
        self.line = line
        
        while(True):
            oldstate = self.state
            self.make_transition(index,line)
            if self.state == oldstate: break
        
    def action(self):
        if self.state == 'initial':
            # Should never get here, throw an error
            raise RuntimeError('<Inital> is not a valid state to perform an action on')
            
        elif self.state == 'comments':
            # Just ignore these
            pass
        
        elif self.state == 'header':
            header = self.line.split()
            self.values['numNodes'] = int(header[0])
            self.values['numCells'] = int(header[1])
        
        elif self.state == 'nodes':
            entry = convert(self.line.split(), self.entry_types['node'])
            self.nodes[entry[0]] = entry[1:]
        
        elif self.state == 'cells':
            entry = self.line.split()
            entry = convert(entry, self.entry_types[entry[2]])
            
            self.cells[entry[2]][entry[0]] = [entry[1]] + entry[3:]


def riffleShuffle(pile):
    # given an indexable data structure, riffle shuffles that structure
    
    middle = len(pile)/2
    
    cut1 = pile[:(middle)]
    cut2 = pile[(middle):]
    
    riffle = []
    for i in range(middle):
        riffle.append(cut1[i])
        riffle.append(cut2[i])
        
    riffle.extend(cut1[middle:])
    riffle.extend(cut2[middle:])
                
    return riffle
    
def getBounds(parser):
    # finds the minum and maximum X Y and Z values for nodes in the UCD file
    
    x,y,z = zip(*[n for n in parser.nodes.values()])
    return [[min(x),min(y),min(z)],[max(x),max(y),max(z)]]
    
    
def plotElements(ax, elem, conn, color, text=False, lwidth=1):
    # actually does the work of drawing the wireframe of the element
    
    for i,e in parser.cells[elem].iteritems():
        verts = e[1:]
        for line in conn:
            # look up the nodes of the next line to draw in the shape feed that 
            # information to matplotlib
            ax.plot3D(*zip(parser.nodes[verts[line[0]]], parser.nodes[verts[line[1]]]),color=color, linewidth=lwidth)
        
        # draw the index number of the element in the middle of the element
        if text:
            txt = [sum(l)/len(verts) for l in zip(*[parser.nodes[n] for n in verts])]
            txt.append(i)
            ax.text(*txt,color=color) 


def showElements(ax, parser, hex_text, quad_text, line_text):
    ax.set_aspect('equal')
    
    # gather the connectivity information about UCD files.  This comes from the 
    # deal.ii documentation on what order verticies come in in a hex
    # Then, draw the element.
    hex_conn = [(0,1),(1,2),(2,3),(3,0),(0,4),(1,5),(2,6),(3,7),(4,5),(5,6),(6,7),(7,4)]
    plotElements(ax, 'hex', hex_conn, 'b', hex_text, lwidth=4)
    
    # do the same for quads
    quad_conn = [(0,1),(1,2),(2,3),(3,0)]
    plotElements(ax, 'quad', quad_conn, 'r', quad_text, lwidth=2)
    
    # and lines
    line_conn = [(0,1)]
    plotElements(ax, 'line', line_conn, 'g', line_text, lwidth=1)    
    
    return fig
        
        
def showBoundaries(ax, parser, dim, text):
    ax.set_aspect('equal')
    
    if dim==3:
    
        # calculate how many boundary terms we have
        # determines how many shades of colors we need and how the legend is drawn
        bounds = set([b[0] for b in parser.cells['quad'].values()])
        bounds_min = min(bounds)
        bounds_max = max(bounds)
        bounds_range = bounds_max-bounds_min
        
        # I riffle shuffle the colors to try to make sure that boundary terms
        # near each other will get contrasting colors
        colors = riffleShuffle(rainbow(np.linspace(0,1,bounds_range+1)))
    
        for i,e in parser.cells['quad'].iteritems():
            verts = e[1:]
            material = e[0]
            
            coords = [parser.nodes[verts[n]] for n in range (4)]
            tri = Poly3DCollection([coords])
            tri.set_color(colors[material-bounds_min])
            tri.set_edgecolor('k')
            ax.add_collection(tri)
        
        # create the legend manually, by createing and filling in the rectangles
        # and then appending the material number    
        guide = [plt.Rectangle((0, 0), 1, 1, fc=color) for color in colors]
        plt.legend(guide, bounds, loc='best')

    return fig
 
         
if __name__=='__main__':       
    with open(sys.argv[1], 'rb') as f:
        contents=f.readlines()
        
        parser = UCD_Parser()
        
        # parse the UCD file
        for index,line in enumerate(contents):
            parser.transition(index,line)
            parser.action()
       
       
        ranges = [list(x) for x in zip(*getBounds(parser))]
        
        fig = plt.figure()
        
        ax = fig.add_subplot(1,2,1,projection='3d')
        ax.set_xlim(ranges[0])
        ax.set_ylim(ranges[1])
        ax.set_zlim(ranges[2])
        ax.set_title('Element IDs')
        showElements(ax, parser, hex_text=True, quad_text=True, line_text=True)
        
        ax = fig.add_subplot(1,2,2,projection='3d')
        ax.set_xlim(ranges[0])
        ax.set_ylim(ranges[1])
        ax.set_zlim(ranges[2])
        ax.set_title('Boundary Labels')
        showBoundaries(ax, parser, 3, True)
        
        plt.show()
    



