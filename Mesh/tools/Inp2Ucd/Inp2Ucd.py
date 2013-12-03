from collections import defaultdict
import sys
import os
import textwrap

def printList(l):
  print '\n'.join([str(x) for x in l])

def convert(arry, args):
  # takes a list and casts the elements in that list according to the types given
  return [f(s) for f,s in zip(args,arry)]
    
def stripList(l):
  return [x.strip() for x in l]


class recursiveDict(dict):
  def __getitem__(self, key):
    try:
      return super(recursiveDict, self).__getitem__(key)
    except KeyError:
      super(recursiveDict, self).__setitem__(key, recursiveDict())
      return super(recursiveDict, self).__getitem__(key)

    
class InpFileParser:
  
  def p_node(self):
    # defines logic to perform if we're inside a 'node' command
    
    # split the node into individual cells, and strip cells of whitespace
    # cast those cells into the following set form [int, float, float, float]
    # the [int] is the node number, use as key in dictionary
    # the [float, float, float] are the xyz coordinates of node, store as value in dictionary 
    node = convert(stripList(self.current_line.split(',')), [int, float, float, float])
    self.parsed['nodes'][node[0]] = node[1:]
    
  def p_element(self):
    # defines logic to perform if we're inside a 'element' command

    # split the element into individual cells, and strip cells of whitespace
    # cast those cells into ints
    # the 1st entry is the element number, use as key in dictionary
    # all the rest of the entrys are nodes in element, store as value in dictionary
    element = map(int, stripList(self.current_line.split(',')))
    self.parsed['elements'][self.arguements['type']][element[0]] = element[1:]
  
  def i_elset(self):
    # initialization code for elset command
    self.parsed['elsets'][self.arguements['elset']] = []
    
  def p_elset(self):
    # defines the logic to perform if we're inside an 'elset' command
    
    # this code only works on non 'generate' elsets
    if 'generate' not in self.arguements:
      # that list will then have all the element numbers appended into it
      elements = [int(x) for x in stripList(self.current_line.split(','))]
      for e in elements:
        self.parsed['elsets'][self.arguements['elset']].append(e)
        
    # in the case where the elset is 'generated'
    else:
      # get the 3 numbers that represent the generate syntax 
      # which are (first element, last element, interval)
      generates = map(int, stripList(self.current_line.split(',')))
      
      # generate the list of elements from the above number.  
      # the goofy syntax at the end adds the last element because range() is inclusive-exclusive
      elements = range(*generates) + [generates[1]]
      
      # store the list of elements in the parsed data structure
      self.parsed['elsets'][self.arguements['elset']] = elements
       
  def i_surface(self):
    # initialization code for surface command
    self.parsed['surfaces'][self.arguements['name']] = []
        
  def p_surface(self):
    # defines the logic to perform if we're inside a 'surface' command
    if self.arguements['type'] == 'element':
      preparedLine = self.current_line.strip().lower() 
      arglist = stripList(preparedLine.split(','))
      elset = arglist[0]
      face = int(arglist[1][1:])
      
      self.parsed['surfaces'][self.arguements['name']].append((elset, face))

    else:
      print '>>Found surface that wasn\'t of type ELEMENT!'
      pass

  def p_unspecified(self):
    pass  
  
  def __init__(self):
  
    self.reserved_keywords = ['unspecified']
    
    self.keywords = [x[2:] for x in filter(lambda x: x.find('p_') == 0 and x not in [('p_' + x) for x in self.reserved_keywords], dir(InpFileParser))]
    
    self.current_lineNo = None
    self.current_line = None
    self.state = 'unspecified'
    self.arguements = None
    
    self.parsed = recursiveDict()
    
    
  def transition(self):
    # If this is a command line, then transition to parsing this command
    
    # try to find a command in the line
    newState = self.findKeyword(self.current_line, self.keywords)

    if newState is not None:
      # if a command is found then set the machine to parse the new state, and find any
      # relevant arguements
      self.state = newState
      self.arguements = self.findArguements(self.current_line)
      
      print '>>Transitioned line: %s to: <%s>' % (self.current_lineNo, self.state)
      print '>>Found Arguements: %s' % self.arguements
      
      try:
        getattr(self, 'i_' + self.state)()
      except AttributeError:
        print '>>No initialization method found'
      else:
        print '>>Initialized state'
        
      print ''
      return True
    return False
  
  def action(self):
    # don't try to parse lines only containing whitespace
    if self.current_line.strip() != '':
      getattr(self, 'p_' + self.state)()
  
  def findKeyword(self, line, keywords):
    # try to find any of the commands in the line
    preparedLine = line.lstrip().lower()
    for key in keywords:
      if preparedLine.find('*' + key) == 0:
        self.arguements = self.findArguements(line)
        return key
        
    # this is a command line, but not one we care about
    if preparedLine[0] == '*':
      return 'unspecified'
    return None
    
  def findArguements(self, line):
    # remove preceding whitespace, and change everything to lower case for processing
    preparedLine = line.lstrip().lower()
    
    # split out arguements by comma, discard the first one (which is the command)
    args = stripList(preparedLine.split(','))[1:]
    
    # split each arguement by '='. 
    # if arguement does not have a value, append 'None' (for a placeholder)
    # and put into a dictionary
    pairs = [stripList(x.split('=')) for x in args]
    for pair in pairs:
      if len(pair) == 1: pair.append(None)
    
    return dict(pairs)
  
  def parseLine(self, lineNo, line):
    # these will be used by the rest of the class for performing the correct action
    self.current_lineNo = lineNo
    self.current_line   = line
    
    # see if we're on a command line - if so, change state
    transitioned = self.transition()
    
    # if we are on a command line, do not try to perform an action on it
    # otherwise, parse the line according to the logic in the p_* functions
    if not transitioned:  
      self.action()
      

class UCDGenerator:
  def __init__(self):
    self.nodes = []
    self.elements = []
    self.boundaries = []
    
  def generateNodes(self, inpNodes):
  # expects a dictionary of the form {<nodenumber>:[x, y, z]}
    self.nodes.extend(sorted([[k] + v for k,v in inpNodes.iteritems()]))
  
  def generateHexes(self, inpElements):
    # expects a dictionary of the form {<nodenumber>:[1, 2, 3, 4, 5, 6, 7, 8]}  
    # expects nodes to be in the following order: (taken from ABAQUS documents)
    #         7-------6        7-------6
    #        /|       |       /       /|
    #       / |       |      /       / |
    #      /  |       |     /       /  |
    #     4   |       |    4-------5   |
    #     |   3-------2    |       |   2
    #     |  /       /     |       |  /
    #     | /       /      |       | /
    #     |/       /       |       |/
    #     0-------1        0-------1
    
    temp_elements = sorted([[k] + v for k,v in inpElements.iteritems()])
    self.elements.extend([([x[0], '1', 'hex', x[1], x[2], x[6], x[5], x[4], x[3], x[7], x[8]]) for x in temp_elements])
    
  def generateBoundaries(self, inpSurfaces, inpElsets, inpElements):
    face_defns = [[3,2,1,0],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]]
  
    ctr = 1
    for material,surface in enumerate(inpSurfaces):
      for elset,face in inpSurfaces[surface]:
        for element in inpElsets[elset]:
          self.boundaries.append([ctr, material+1, 'quad'] + [inpElements[element][x] for x in face_defns[face-1]])
          ctr += 1
      print '>>Chose material no <%s> for surface: %s' % (material+1, surface)
  
  def createUCD(self):
    numNodes = len(gen.nodes)
    numElements = len(gen.elements) + len(gen.boundaries)
    
    preamble = textwrap.dedent('''\
      # Inp2Ucd Converted Mesh 
      # Input file name: %s''' % sys.argv[1])
        
                  
    header = ("%s   %s   0   0   0" % (numNodes, numElements))
    
    nodes      = os.linesep.join(['  '.join(map(lambda x: str(x).ljust(14, ' '), node))    for node in self.nodes])
    elements   = os.linesep.join(['  '.join(map(lambda x: str(x).ljust(5, ' '), element))  for element in self.elements])
    boundaries = os.linesep.join(['  '.join(map(lambda x: str(x).ljust(5, ' '), boundary)) for boundary in self.boundaries])

    output = os.linesep.join([preamble, header, nodes, elements, boundaries])
    
    return output
    
if __name__ == '__main__':

  parser = InpFileParser()
  
  with open(sys.argv[1], 'rb') as f: 
    for lineNo,line in enumerate(f):
      parser.parseLine(lineNo+1, line)
     
  gen = UCDGenerator()
  
  gen.generateNodes(parser.parsed['nodes'])
  gen.generateHexes(parser.parsed['elements']['fc3d8'])
  gen.generateBoundaries(parser.parsed['surfaces'], parser.parsed['elsets'], parser.parsed['elements']['fc3d8']) 

  with open(sys.argv[2], 'wb') as f:
    f.write(gen.createUCD())

