import matplotlib.pyplot as plt

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

class UCDParser:

    def t_initial(self):
        return 'comments'

    def t_comments(self):
        if self.current_line.lstrip()[0] != '#':
            self.lineAnchor = self.current_lineNo
            return 'header'

    def t_header(self):
        if self.current_lineNo == self.lineAnchor+1:
            self.lineAnchor = self.current_lineNo
            return 'nodes'
            
    def t_nodes(self):
        if self.current_lineNo == self.lineAnchor + self.parsed['header']['numNodes']:
            self.lineAnchor = self.current_lineNo
            return 'elements'

    def t_elements(self):
        pass

    def p_initial(self):
        raise RuntimeError('>><%s> is not a valid state to perform actions in' % self.state)

    def p_comments(self):
        pass

    def p_header(self):
        header = map(int,stripList(self.current_line.split()))
        self.parsed['header']['numNodes']    = header[0]
        self.parsed['header']['numElements'] = header[1]

    def p_nodes(self):
        entry = convert(stripList(self.current_line.split()), (int,float,float,float))
        self.parsed['nodes'][entry[0]] = entry[1:]

    def p_elements(self):
        element_types = {
            'line':(int,int,str,int,int),
            'tri' :(int,int,str,int,int,int),
            'quad':(int,int,str,int,int,int,int),
            'hex' :(int,int,str,int,int,int,int,int,int,int,int)}
        
        entry = stripList(self.current_line.split())
        entry = convert(entry, element_types[entry[2]])

        self.parsed['elements'][entry[2]][entry[0]]['material'] = entry[1]
        self.parsed['elements'][entry[2]][entry[0]]['nodes'] = entry[3:]

    def __init__(self):
        self.state = 'initial'
        self.parsed = recursiveDict()

    def transition(self):
        try:
            newstate = getattr(self, 't_' + self.state)()
            if newstate != None:
                self.state = newstate
                print '>>Transitioned to <%s>' % self.state
        except:
            pass
            #print '>>Did not find transition <t_%s>' % self.state
        
    def action(self):
        try:
            getattr(self, 'p_' + self.state)()
        except:
            print '>>Did not find action <p_%s>' % self.state

    def parseLine(self, lineNo, line):
        self.current_lineNo = lineNo
        self.current_line   = line

        while(True):
            oldstate = self.state
            self.transition()
            if oldstate == self.state: break

        self.action()

if __name__ == '__main__':
    print 'Starting...'
    parser = UCDParser()

    with open('slide.ucd', 'rb') as f:
        for lineNo,line in enumerate(f):
            parser.parseLine(lineNo+1, line)

    colors = ['r', 'b', 'y', 'g']
    fig = plt.figure()
    for line in parser.parsed['elements']['line'].values():
        x,y,z = zip(*[parser.parsed['nodes'][node] for node in line['nodes']])
        plt.plot(x,y, colors[line['material']-1])

    guide = [plt.Rectangle((0, 0), 1, 1, fc=color) for color in colors]
    plt.legend(guide, ['1','2','3','4'], loc='best')
    plt.xlim([-1.2,1.2])
    plt.show()
        
    
