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

# stop windows from showing up when we don't want them
plt.ioff()

# <codecell>

from collections import defaultdict

def convert(arry, args):
    return [f(s) for f,s in zip(args,arry)]

rdefaultdict = lambda:defaultdict(rdefaultdict)

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


def plotElements(ax, elem, conn, color, text=False, lwidth=1):
    for i,e in parser.cells[elem].iteritems():
        verts = e[1:]
        for line in conn:
            ax.plot3D(*zip(parser.nodes[verts[line[0]]], parser.nodes[verts[line[1]]]),color=color, linewidth=lwidth)
        
        if text:
            txt = [sum(l)/len(verts) for l in zip(*[parser.nodes[n] for n in verts])]
            txt.append(i)
            ax.text(*txt,color=color) 

def showElements(parser, hex_text, quad_text, line_text):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection='3d')
    ax.set_aspect('equal')
    
    hex_conn = [(0,1),(1,2),(2,3),(3,0),(0,4),(1,5),(2,6),(3,7),(4,5),(5,6),(6,7),(7,4)]
    plotElements(ax, 'hex', hex_conn, 'b', hex_text, lwidth=4)
    
    quad_conn = [(0,1),(1,2),(2,3),(3,0)]
    plotElements(ax, 'quad', quad_conn, 'r', quad_text, lwidth=2)
    
    line_conn = [(0,1)]
    plotElements(ax, 'line', line_conn, 'g', line_text, lwidth=1)    
    
    plt.show()
   
def riffleShuffle(pile):
    
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
        
def showBoundaries(parser, dim, text):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.set_aspect('equal')
    
    if dim==3:
    
        bounds = set([b[0] for b in parser.cells['quad'].values()])
        bounds_min = min(bounds)
        bounds_max = max(bounds)
        bounds_range = bounds_max-bounds_min
        
        colors = riffleShuffle(rainbow(np.linspace(0,1,bounds_range+1)))
    
        for i,e in parser.cells['quad'].iteritems():
            verts = e[1:]
            material = e[0]
            
            coords = [parser.nodes[verts[n]] for n in range (4)]
            tri = Poly3DCollection([coords])
            tri.set_color(colors[material-bounds_min])
            tri.set_edgecolor('k')
            ax.add_collection(tri)
            
        guide = [plt.Rectangle((0, 0), 1, 1, fc=color) for color in colors]
        plt.legend(guide, bounds, loc='best')

    plt.show()
         
        
with open('test_cube_two_materials.ucd', 'rb') as f:
    contents=f.readlines()
    
    parser = UCD_Parser()
    
    # parse the UCD file
    for index,line in enumerate(contents):
        parser.transition(index,line)
        parser.action()
   
    #showElements(parser, hex_text=True, quad_text=False, line_text=False)
    showBoundaries(parser, 3, True)
    



