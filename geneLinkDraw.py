import math
from matplotlib import font_manager
from tkinter import font
from PIL import Image, ImageColor, ImageDraw, ImageFont

# Classes to represent genes, exons, links, ...

class Gene:
    """
    Class to represent a gene for drawing

    ...
    
    Attributes
    ----------
    id : str
        gene identifier
    species : str
        identifier of the species the gene belongs to
    len : int
        gene length in bp
    strand : str
        gene orientation, either '+' or '-'
    elements : dict
        gene elements like exons, has the structure {"type1": [(start0, end0), (start1, end1), ...], 
                                                     "type2": [...], 
                                                     ...}

    Methods
    -------
    addElement(elemtype: str, start: int, end: int):
        adds a gene element (e.g. an exon)
    """
    def __init__(self, geneID: str, speciesID: str, length: int, strand: str):
        """
        Constructor

        Parameters
        ----------
        geneID : str
            gene identifier
        geneID : str
            species identifier
        length : int
            gene length in bp
        strand : str
            gene orientation, either '+' or '-'
        """
        assert strand in ['+', '-'], "[ERROR] >>> Strand must be '+' or '-', not "+str(strand)
        self.id = geneID
        self.species = speciesID
        self.len = length
        self.strand = strand
        self.elements = dict()

    def __str__(self):
        return str({'id': self.id, 'species': self.species, 'len': self.len, 'strand': self.strand, 
                    'elements': self.elements})

    def addElement(self, elemtype: str, start: int, end: int):
        """
        Adds a gene element (e.g. an exon)

        Parameters
        ----------
        elemtype : str
            element type
        start : int
            relative element start position inside the gene, top strand, zero-based.
            I.e. start = 0 refers to the leftmost gene base when looking at the top strand of the chromosome,
                 start = gene.len-1 refers to the rightmost gene base
        end : int
            relative element end position inside the gene, top strand, zero-based.
        """
        if elemtype not in self.elements:
            self.elements[elemtype] = []

        self.elements[elemtype].append((start, end))



class Link:
    """
    Class to represent a link for drawing

    ...
    
    Attributes
    ----------
    genes : list
        list of gene identifiers for genes connected by the link
    pos : list
        list of positions for each respective gene where the link connects them
    strand : list
        optional list of strands on which the link connects the genes
    """
    def __init__(self, genes: list[Gene], pos: list[int], strands: list[str] = None):
        """
        Constructor

        Parameters
        ----------
        genes : list
            list of gene identifiers for genes connected by the link
        pos : list
            list of positions for each respective gene where the link connects them
        strand : list
            optional list of strands on which the link connects the genes
        """
        assert len(genes) == len(pos), "[ERROR] >>> genes and pos must be of equal length"
        assert len(genes) >= 2, "[ERROR] >>> A link must connect at least two genes"
        assert len(set(genes)) == len(genes), "[ERROR] >>> All genes in a link must be unique"
        self.genes = genes
        self.pos = pos
        if strands is not None:
            assert len(strands) == len(genes), "[ERROR] >>> genes and strands must be of equal length"
            assert all([s in ['+', '-'] for s in strands]), "[ERROR] >>> Strands must be '+' or '-'"
        
        self.strands = strands

    def __str__(self):
        return str({'genes': self.genes, 'pos': self.pos, 'strands': self.strands})



class GeneCoordinates:
    """
    Class to store drawing coordinates for a gene

    ...
    Attributes
    ----------
    color : str or tuple of RGB values
        color value (string or rgb) of drawn gene
    gene : Gene
        reference to the corresponding gene
    id : str
        id of the corresponding gene
    len : int
        gene length in pixels
    res : float
        drawing resolution in pixel / basepair
    x0 : int
        pixel coordinate of the first base
    y : int
        vertical pixel coordinate
    """
    def __init__(self, gene: Gene, x0: int, y: int, res: float, col = None):
        """
        Constructor

        Parameters
        ----------
        gene : Gene
            corresponding gene object
        x0 : int
            pixel coordinate of the first (leftmost) base
        xn : int
            pixel coordinate of the last (rightmost) base
        y : 0
            vertical pixel coordinate
        res : float
            drawing resolution in pixel / basepair
        col : str or tuple of RGB values
            optional color of drawn gene, defaults to blue for fwd strand genes and red otherwise
        """
        self.gene = gene
        self.id = gene.id
        assert x0 >= 0, '[ERROR] >>> pixel coordinate cannot be negative'
        assert y >= 0, '[ERROR] >>> pixel coordinate cannot be negative'
        assert res > 0, '[ERROR] >>> resolution must be positive and greater than zero'
        self.x0 = x0
        self.y = y
        self.updateResolution(res) # sets self.res, self.len and self.xn
        if col is None:
            self.color = 'darkblue' if self.gene.strand == "+" else 'darkorange'
        else:
            self.color = col

    def __str__(self):
        return str({'id': self.id, 'color': self.color, 'len': self.len, 'res': self.res, 'x0': self.x0, 'y': self.y})

    def updateResolution(self, res : float):
        """
        Set new drawing resolution, update pixel length accordingly

        Parameters
        ----------
        res : float
            new drawing resolution in pixel / basepair
        """
        self.res = res
        self.len = math.ceil(self.gene.len * self.res)
        self.xn = self.x0 + self.len - 1



# helper class containing array of handpicked color names for automatic coloring
class Palette:
    """
    Helper class containing an array of handpicked color names for automatic coloring

    ...
    Attributes
    ----------
    colors : list
        list of color names
    i : int
        current index

    Methods
    -------
    color():
        returns the current color
    inc():
        sets index to next color
    """
    def __init__(self, idx: int = 0):
        """
        Constructor

        Parameters
        ----------
        idx : int
            initial color index, defaults to 0, set to -1 if given value is too big
        """
        self.colors = ['darkblue', 'darkorange', 'blueviolet', 'burlywood', 'darkgreen', 'darkmagenta', 'cadetblue',
                       'chocolate', 'cyan', 'darkgoldenrod', 'darkcyan', 'darkkhaki', 'darkred', 'darkolivegreen', 
                       'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategrey', 'darkturquoise', 'deeppink',
                       'dodgerblue', 'gold', 'greenyellow', 'indigo', 'khaki', 'lightblue', 'maroon', 'plum', 'teal', 
                       'bisque', 'aquamarine']
        self.i = idx if idx < len(self.colors) else -1 # default to last color if idx too large

    def color(self):
        """
        Returns the current color name
        """
        return self.colors[self.i]

    def inc(self):
        """
        Sets the index to next color, starting over if color list has endet
        """
        self.i += 1
        if self.i >= len(self.colors):
            self.i = 0



# drawing functions

def draw(genes: list[Gene], links: list[Link], font = None,
         width: int = 1920, height: int = 1080, dpi: int = 100,
         outerMargin: int = None, genewidth: int = 5, linkwidth: int = 2, fontsize: int = 12,
         genecols: list = None, elementcols: dict = None, linkcol = None, show: bool = True):
    """
    Draws an image with genes as horizontal bars, possibly with gene elements like exons inside, 
      and links connecting them

    Returns a PIL.Image object

    Parameters:
        genes: list of Gene objects
        links: list of Link objects
        font: path to a truetype font (if omitted, default font will be used with fixed fontsize)
        width: image width in pixels
        height: image heigth in pixels
        dpi: image resolution in dpi
        outerMargin: distance from image edge to drawn content in pixels, automatically set if None
        genewidth: line width of drawn genes in pixels
        linkwidth: line witdh of drawn links in pixels
        fontsize: fontsize of text
        genecols: optional list of colors for each single gene, same length as genes
        elementcols: optional dict with element type names as keys and color values as values
        linkcol: optional color value (name or RGB tuple) or 2-tuple of color values to use for all links. If single
                 value, all links are colored this way, otherwise links between opposing strands are colored in the 
                 second color value. If None, links are colored automatically
        show: set to False to suppress image drawing
    """

    # some sanity checks
    assert len(genes) >= 2, "[ERROR] >>> Single genes cannot be linked"
    if outerMargin is not None:
        assert outerMargin >= 0, "[ERROR] >>> margin must be positive"
        assert outerMargin < min(width/2, height/2), "[ERROR] >>> margin too big"
    else:
        outerMargin = int(height*0.012)

    assert fontsize > 0, "[ERROR] >>> Fontsize must be bigger than 0"
    if font is not None:
        font = ImageFont.truetype(font, fontsize)

    if genecols is not None:
        assert len(genecols) == len(genes), "[ERROR] >>> genecols must have same length as genes"
    
    if linkcol is not None and type(linkcol).__name__ == 'tuple':
        assert len(linkcol) >= 2, "[ERROR] >>> linkcol must be a single color name or value or a 2-tuple of colors"

    # determine number of gene rows based on species
    speciesToRow = {} # {"species": row_num, ...}
    geneGrid = [] # [[gene_coord, gene_coord, ...], [...], ...]
    geneToCoord = {} # {"geneID": gene_coord, ...}
    for i in range(len(genes)):
        gene = genes[i]
        if genecols is not None:
            gc = GeneCoordinates(gene, 0, 0, 1, col = genecols[i]) # set coords and res later
        else:
            gc = GeneCoordinates(gene, 0, 0, 1) # set coords and res later

        if gene.species not in speciesToRow:
            speciesToRow[gene.species] = len(geneGrid)
            geneGrid.append([])

        geneGrid[speciesToRow[gene.species]].append(gc)
        geneToCoord[gene.id] = gc

    # special case: single species, draw all genes stacked
    if len(geneGrid) == 1:
        geneGrid = [[g] for g in geneGrid[0]]

    nrows = len(geneGrid)

    # assign y-coordinate (top) of each gene row
    rowheight = genewidth + int(genewidth/2) + fontsize + genewidth # last genewith as minimal space between rows
    rowpix = height - (2 * outerMargin)
    if nrows * rowheight > rowpix:
        height = (2 * outerMargin) + (nrows * rowheight)
        print("[WARNING] >>> Too many gene rows for image height, increasing image height to",
              height, "pixels. Adjust outerMargin, genewidth and fontsize if you want to keep the lower image height")

    vspace = int((height - (nrows*rowheight) - (2 * outerMargin)) / (nrows-1)) if nrows > 1 else 0
    #rowToY = {}
    y = outerMargin
    for row in geneGrid:
        #rowToY[r] = y
        for gc in row:
            gc.y = y

        y += (rowheight + vspace)
    
    # determine longest gene row and pixel/bp resolution
    maxRowLen = 0
    ngenes = 0
    for row in geneGrid:
        rowLen = sum([gc.gene.len for gc in row])
        if rowLen > maxRowLen:
            maxRowLen = rowLen
            ngenes = len(row)
        
    genesep = int(width * 0.006) # number of pixels between two genes in a row
    genepix = width - ((ngenes-1) * genesep) - (2 * outerMargin) # number of pixels for gene drawing
    pxpbp = genepix / maxRowLen # pixel per basepair (horizontal)

    # set for each gene the appropriate x0 and resolution
    for row in geneGrid:
        x0 = outerMargin+1
        for gc in row:
            gc.x0 = x0
            gc.updateResolution(pxpbp) # sets correct res, pixellen, xn
            x0 += (genesep + gc.len)
    
    # start drawing
    img = Image.new(mode = "RGB", size = (width, height), color = "white")
    drw = ImageDraw.Draw(img) # drawing context    
    palette = Palette(2) # default gene coloring are darkblue and darkorange, start at blueviolet for additional colors
    if elementcols is None:
        typeToPalette = {} # fill with default colors later
    else:
        typeToPalette = elementcols # use user-provided coloring

    # draw genes and elements
    for gid in geneToCoord:
        gc = geneToCoord[gid]
        drw.line(xy = ((gc.x0, gc.y), (gc.xn, gc.y)),
                    fill = gc.color,
                    width = genewidth)

        # write gene name below
        drw.text(xy = (gc.x0, gc.y+int(1.5*genewidth)), 
                    text = gc.gene.id,
                    font = font, fill = "black")

        for elemtype in gc.gene.elements:
            if elemtype not in typeToPalette:
                typeToPalette[elemtype] = palette.color()
                palette.inc()

            for elem in gc.gene.elements[elemtype]:
                a = gc.x0 + (gc.res * elem[0])
                b = gc.x0 + (gc.res * elem[1])
                drw.line(xy = ((a, gc.y), (b, gc.y)),
                         fill = typeToPalette[elemtype],
                         width = genewidth)

    # draw links
    if links is not None:
        if linkcol is None:
            sscol = palette.color()
            palette.inc()
            oscol = palette.color()
        else:
            if type(linkcol).__name__ == 'tuple' and len(linkcol) == 2:
                sscol = linkcol[0]
                oscol = linkcol[1]
            else: # assume 3-tuple of RGB values
                sscol = linkcol
                oscol = linkcol

        for link in links:
            for i in range(1, len(link.genes)):
                gc1 = geneToCoord[link.genes[i-1].id]
                gc2 = geneToCoord[link.genes[i].id]
                x1 = gc1.x0 + (gc1.res * link.pos[i-1])
                x2 = gc2.x0 + (gc2.res * link.pos[i])
                if link.strands is not None and link.strands[i-1] != link.strands[i]:
                    col = oscol
                else:
                    col = sscol

                drw.line(xy = ((x1, gc1.y), (x2, gc2.y)),
                        fill = col,
                        width = linkwidth)

    if show:
        img.show()

    return img



# user helper functions

def getAvailableFonts():
    """
    Returns and prints a list of font paths found on your system
    that you can use with the draw() function
    """
    
    system_fonts = font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
    print("\n".join(system_fonts))
    return system_fonts



def getColorSheet(font = None):
    """
    Draws a big image with all named colors available

    Parameters:
        font: path to a truetype font (if omitted, default font will be used with fixed fontsize)
    """

    fakegenes = []
    colors = []
    i = 0
    for name, code in ImageColor.colormap.items():
        fakegenes.append(Gene(geneID = str(i)+" - "+name+" // "+str(code),
                              speciesID = str(i // 4), # four colors per row
                              length = 100, strand = "+"))
        colors.append(name)
        i += 1

    # calculate image height
    genewidth = 10
    fontsize = 12
    outerMargin = 10
    rowheight = genewidth + int(genewidth/2) + fontsize + genewidth # last genewith as minimal space between rows
    height = (2 * outerMargin) + (((i // 4)+1) * rowheight)
        
    draw(fakegenes, None, genecols = colors, 
         genewidth = 10, outerMargin = outerMargin, fontsize = fontsize,
         height = height, font = font)