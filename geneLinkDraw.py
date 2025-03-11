from dataclasses import dataclass
import logging
import math
from matplotlib import font_manager
#import os
from pathlib import Path
#from tkinter import font
from PIL import Image, ImageColor, ImageDraw, ImageFont



class Palette:
    """
    Helper class containing an array of handpicked color names for automatic coloring

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
        """ Returns the current color name """
        return self.colors[self.i]
    
    def colorpp(self):
        """ Returns the current color name & increments the index, so that the next call will return the next color """
        c = self.colors[self.i]
        self.inc()
        return c

    def inc(self):
        """ Sets the index to next color, starting over if color list has endet """
        self.i += 1
        if self.i >= len(self.colors):
            self.i = 0



# Classes to represent genes, exons, links, ...

@dataclass
class Gene:
    """
    Class to represent a gene for drawing
    
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
    link_anchors: set[int]
        list of positions where links can connect to the gene, populated at Link object creation
    elements : dict
        gene elements like exons, has the structure {"type1": [(start0, end0), (start1, end1), ...], 
                                                     "type2": [...], 
                                                     ...}
    sites : dict
        gene sites like start and stop codons, has the structure {"type1": [pos1, pos2, ...], "type2": [...], ...}
    drawInfo : GeneDrawInfo 
        GeneDrawInfo object, contains drawing information like coordinates and colors, added later via draw() function

    Methods
    -------
    addElement(elemtype: str, start: int, end: int):
        adds a gene element (e.g. an exon)
    """

    id: str
    species: str
    len: int
    strand: str
    #elements: dict = field(default_factory=dict) # initialized empty dict, populate via addElement()
    #sites: dict = field(default_factory=dict) # initialized empty dict, populate via addSite()

    def __post_init__(self):
        assert self.strand in ['+', '-'], f"[ERROR] >>> Strand must be '+' or '-', not '{self.strand}'"
        assert self.len > 0, "[ERROR] >>> Gene length must be positive"
        self.link_anchors: set[int] = set()
        self.elements: dict[str, list[tuple[int, int]]] = {}
        self.sites: dict[str, list[int]] = {}
        self.drawInfo = None # added later in draw() function

    def __str__(self):
        return str({'id': self.id, 'species': self.species, 'len': self.len, 'strand': self.strand, 
                    'elements': self.elements, 'sites': self.sites, 'link_anchors': self.link_anchors})

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
            relative element end position inside the gene, top strand, zero-based, exclusive (i.e. NOT part of the element!).
        """
        if elemtype not in self.elements:
            self.elements[elemtype] = []

        self.elements[elemtype].append((start, end))

    def addSite(self, sitetype: str, positions: list[int]):
        """
        Adds one or more gene site (e.g. a start codon, occurrences of some sort, ...)

        Parameters
        ----------
        sitetype : str
            site type
        positions : int | list[int]
            (list of) position(s) where the site(s) occur(s) inside the gene, top strand, zero-based.
        """
        if sitetype not in self.sites:
            self.sites[sitetype] = []

        if type(positions) == int:
            self.sites[sitetype].append(positions)
        else:
            self.sites[sitetype].extend(positions)




@dataclass
class Link:
    """
    Class to represent a link for drawing. Can be in compressed mode to represent many links instead of individual ones.
    
    Attributes
    ----------
    genes : list
        list of Gene objects connected by the link
    pos : list | list[list]
        list of positions for each respective gene where the link connects them (uncompressed mode, default). If
        compressed mode is activated, this is a list of lists of positions, giving one or more positions for each
        respective gene
    strands : list
        optional list of strands on which the link connects the genes. Note that currently storing of differing strands
        per gene in compressed mode is not implemented.
    connect : bool
        True by default. If False, do not connect occurrences when drawing and only draw markers at occurrence sites
    compressed : bool
        False by default. If True, the object is storing many links among the same genes
    color : None or str or tuple of RGB values
        color of drawn link, set in draw(), otherwise is determined automatically
    """

    genes: list[Gene]
    pos: list[int]
    strands: list[str] = None
    connect: bool = True
    compressed: bool = False

    def __post_init__(self):
        assert len(self.genes) == len(self.pos), "[ERROR] >>> genes and pos must be of equal length"
        assert len(self.genes) >= 2, "[ERROR] >>> A link must connect at least two genes"
        if self.compressed:
            assert all([len(pl) > 0 for pl in self.pos]), \
                "[ERROR] >>> Need lists of positions for each gene in compressed mode, and each list must contain " \
                + f"at least one position. Got: {self.pos}"
            
        assert len(set([gene.id for gene in self.genes])) == len(self.genes), \
            "[ERROR] >>> All genes in a link must be unique"
        if self.strands is not None:
            assert len(self.strands) == len(self.genes), "[ERROR] >>> genes and strands must be of equal length"
            assert all([s in ['+', '-'] for s in self.strands]), "[ERROR] >>> Strands must be '+' or '-'"

        # populate link_anchors in genes
        for i, gene in enumerate(self.genes):
            if self.compressed:
                for p in self.pos[i]:
                    gene.link_anchors.add(p)
            else:
                gene.link_anchors.add(self.pos[i])

        self.color = None

    def __str__(self):
        return str({'genes': self.genes, 'pos': self.pos, 'strands': self.strands})



class GeneDrawInfo:
    """ 
    Class containing actual drawing information (coordinates and colors) for a gene, its label, elements, sites 
      and link anchors 
    """
    def __init__(self, gene: Gene, x0: int, y0: int, label: str, genewidth: int, 
                 font: ImageFont.FreeTypeFont, fontsize: int, res: float, 
                 elementColors: dict, siteColors: dict, geneColor=None):
        """
        Constructor

        Parameters
        ----------
        gene : Gene
            reference to the corresponding gene
        x0 : int
            leftmost pixel coordinate
        y0 : int
            topmost pixel coordinate
        label : str
            gene label to draw above the gene
        genewidth : int
            gene width (vertical) in pixels
        font : ImageFont.FreeTypeFont
            font object for drawing the label
        fontsize : int
            font size
        res : float
            drawing resolution in pixel / basepair
        elementColors : dict
            element type names as keys and color values as values
        siteColors : dict
            site type names as keys and color values as values
        geneColor : str or tuple of RGB values
            color value (string or rgb) of drawn gene, defaults to darkblue for '+' strand and darkorange for '-' strand
        
        """
        self.x0 = x0
        self.y0 = y0
        self.label = label
        self.genewidth = genewidth
        self.font = font
        self.fontsize = fontsize
        self.res = res
        self.elementColors = elementColors
        self.siteColors = siteColors
        if geneColor is None:
            self.geneColor = 'darkblue' if gene.strand == "+" else 'darkorange'
        else:
            self.geneColor = geneColor
        
        assert fontsize > 0, "[ERROR] >>> Fontsize must be bigger than 0"
        assert res > 0, "[ERROR] >>> Resolution must be positive and greater than zero"
        
        # set label coordinates
        try:
            textbbox_left, textbbox_top, textbbox_right, textbbox_bottom = font.getbbox(gene.id) # using top-left anchor
            textw = textbbox_right - textbbox_left
            texth = textbbox_bottom - textbbox_top
        except Exception as e:
            logging.warning("[geneLinkDraw.draw] >>> Could not determine text size via font.getbbox:" + str(e) \
                            + ". Trying deprecated font.getsize() method.")
            try:
                textw, texth = self.font.getsize(gene.id)
            except Exception as e:
                logging.error("[geneLinkDraw.draw] >>> Could not determine text size via font.getsize():" + str(e) \
                            + ". Using fixed text width of 100 pixels.")
                textw = 100
                texth = self.fontsize

        self.x0_label = self.x0
        self.x1_label = self.x0 + textw - 1
        self.y0_label = self.y0
        self.y1_label = self.y0 + texth - 1

        # set gene x-coordinates
        glen = math.ceil(gene.len * self.res)
        self.x0_gene = self.x0
        self.x1_gene = self.x0 + glen - 1
        self.y0_gene = self.y1_label + self.fontsize//2 # use half fontsize as space between label and gene
        self.y1_gene = self.y0_gene + self.genewidth - 1

        self.resetCoords(gene) # initialize element, site and link anchor coordinates

        # add this object to the gene
        gene.drawInfo = self

    
    def __str__(self):
        return f"""{{'x0': {self.x0}, 'y0': {self.y0}, 
                     'label (x0, y0, x1, y1)': ({self.x0_label}, {self.y0_label}, {self.x1_label}, {self.y1_label}),
                     'gene (x0, y0, x1, y1)': ({self.x0_gene}, {self.y0_gene}, {self.x1_gene}, {self.y1_gene}),
                     'elements': {self.elemCoords}, 'sites': {self.siteCoords}, 'linkAnchors': {self.linkAnchors},
                     'label': '{self.label}', 'genewidth': {self.genewidth}, 'fontsize': {self.fontsize}, 
                     'res': {self.res}, 'geneColor': '{self.geneColor}', 'elementColors': {self.elementColors}, 
                     'siteColors': {self.siteColors}}}"""

    
    def resetCoords(self, gene: Gene):
        """ Reset coordinates of gene elements, sites and anchor points when gene coordinates have changed """
        # set element coordinates
        self.elemCoords = {}
        for elemtype in gene.elements:
            self.elemCoords[elemtype] = []
            for elem in gene.elements[elemtype]:
                a = self.x0_gene + math.floor(self.res * elem[0])
                b = self.x0_gene + math.ceil(self.res * (elem[1]-1)) # end is exclusive, but b is inclusive
                self.elemCoords[elemtype].append((a, self.y0_gene, b, self.y1_gene))

        # set site coordinates
        self.siteCoords = {}
        for sitetype in gene.sites:
            self.siteCoords[sitetype] = []
            for pos in gene.sites[sitetype]:
                x = self.x0_gene + math.floor(self.res * pos)
                self.siteCoords[sitetype].append((x, self.y0_gene + self.genewidth//2))

        # set link anchor coordinates
        self.linkAnchors = {
            p: (self.x0_gene + math.floor(self.res * p), self.y0_gene + self.genewidth//2) for p in gene.link_anchors
        }


    def shiftToX(self, gene: Gene, x: int):
        """ Shift all coordinates to a new x-coordinate """
        xshift = x - self.x0
        self.x0 = x
        self.x0_label += xshift
        self.x1_label += xshift
        self.x0_gene += xshift
        self.x1_gene += xshift
        self.resetCoords(gene)


    def shiftToY(self, gene: Gene, y: int):
        """ Shift all coordinates to a new y-coordinate """
        yshift = y - self.y0
        self.y0 = y
        self.y0_label += yshift
        self.y1_label += yshift
        self.y0_gene += yshift
        self.y1_gene += yshift
        self.resetCoords(gene)



def optimizeGeneRow(genes: list[Gene], margin: int):
    """ If multiple genes are in the same row, optimize their positions. More precisely, gene labels can be longer than
        the gene itself, but in order to safe horizontal space, try to shift label heights a bit to still avoid label
        overlap but also not draw genes too far apart. """
    if len(genes) < 2:
        return
    
    assert all([g.drawInfo is not None for g in genes]), "[ERROR] >>> All genes must have a GeneDrawInfo object"
    assert all([g.drawInfo.y0 == genes[0].drawInfo.y0 for g in genes]), "[ERROR] >>> All genes must be in the same row"
    
    y0_orig = genes[0].drawInfo.y0 # save original y0 coordinate to re-align row later
        
    # first, move all genes close together, using the genome bar plus margin as reference
    x_left = min([g.drawInfo.x0 for g in genes])
    for g in genes:
        g.drawInfo.shiftToX(g, x_left)
        x_left = g.drawInfo.x1_gene + margin
        
    # then, try to shift labels up if they overlap
    def labelOverlap(g1, g2):
        return g1.drawInfo.x1_label >= g2.drawInfo.x0_label and g1.drawInfo.x0_label <= g2.drawInfo.x1_label \
                and g1.drawInfo.y0_label <= g2.drawInfo.y1_label and g1.drawInfo.y1_label >= g2.drawInfo.y0_label
    
    def shiftLabelUp(g1, g2) -> bool:
        if labelOverlap(g1, g2):
            logging.debug(f"[INFO] >>> Label overlap between {g1.id}: {g1.drawInfo}\n\nand {g2.id}: {g2.drawInfo}")
            # goal: g2.drawInfo.y1_label = g1.drawInfo.y0_label - 1 --> yshift = g1.drawInfo.y0_label - g2.drawInfo.y1_label
            yshift = g1.drawInfo.y0_label - g2.drawInfo.y1_label - 1
            g2.drawInfo.y0 += yshift
            g2.drawInfo.y0_label += yshift
            g2.drawInfo.y1_label += yshift
            logging.debug(f"[INFO] >>> After shift: {g1.id}: {g1.drawInfo}\n\nand {g2.id}: {g2.drawInfo}")
            assert not labelOverlap(g1, g2), "[ERROR] >>> Label overlap after shifting"
            return True
        return False

    for i in range(len(genes)-1):
        killswitch = 0
        while True and killswitch < 10:
            killswitch += 1 # in case something goes wrong, avoid infinite loop
            shifted = False 
            for rgene in genes[i+1:]:
                shifted = shifted or shiftLabelUp(genes[i], rgene)

            if not shifted:
                break

    # re-align row to original y0 coordinate
    min_y0 = min([g.drawInfo.y0 for g in genes])
    for g in genes:
        g.drawInfo.y0 = min_y0 # set all genes to the same y0 coordinate (top of row) to keep them aligned
        g.drawInfo.shiftToY(g, y0_orig) # shift all coordinates back to original y0 coordinate



# drawing function

def draw(genes: list[Gene], links: list[Link], fontpath: str = None,
         width: int = 1920, height: int = 1080, dpi: int = 100, forceDimensions: bool = False,
         outerMargin: int = None, genewidth: int = 5, linkwidth: int = 2, fontsize: int = 12,
         genecols: list = None, elementcols: dict = None, sitecols: dict = None, linkcols: list = None, 
         show: bool = True):
    """
    Draws an image with genes as horizontal bars, possibly with gene elements like exons inside, 
      and links connecting them

    Returns a PIL.Image object

    Parameters:
        genes: list of Gene objects. This determines the order of genes in the image.
        links: list of Link objects. Links must contain the same Gene objects as in genes list (same memory reference).
               Not all genes need to be linked, but all genes in links must be in genes list.
        fontpath: path (str) to a truetype font (if omitted, default font will be used with fixed fontsize)
        width: image width in pixels
        height: image heigth in pixels
        dpi: image resolution in dpi
        forceDimensions: set to True to force width and height, otherwise image will be resized to fit all genes with
                            approximately the desired width
        outerMargin: distance from image edge to drawn content in pixels, automatically set if None
        genewidth: line width of drawn genes in pixels
        linkwidth: line witdh of drawn links in pixels
        fontsize: fontsize of text
        genecols: optional list of colors for each single gene, same length as genes
        elementcols: optional dict with element type names as keys and color values as values
        sitecols: optional dict with site type names as keys and color values as values
        linkcols: optional list of color values (name or RGB tuple) for each link, same length as links
        show: set to False to suppress image drawing
    """

    # some sanity checks and preparations
    # -----------------------------------

    assert len(genes) >= 2, "[ERROR] >>> Single genes cannot be linked"
    _geneIDs = [gene.id for gene in genes]
    assert sorted(list(set(_geneIDs))) == sorted(_geneIDs), "[ERROR] >>> genes contain duplicate IDs"

    # assert that Gene objects are in links and genes are referencing the same objects in memory
    if links is not None:
        _gidToGene = {gene.id: gene for gene in genes}
        for link in links:
            for gene in link.genes:
                assert gene.id in _gidToGene, "[ERROR] >>> Link contains gene not in genes list"
                assert gene is _gidToGene[gene.id], "[ERROR] >>> Gene objects in links and genes must be the same"
    
    if outerMargin is not None:
        assert outerMargin >= 0, "[ERROR] >>> margin must be positive"
        assert outerMargin < min(width/2, height/2), "[ERROR] >>> margin too big"
    else:
        outerMargin = int(height*0.012)

    assert fontsize > 0, "[ERROR] >>> Fontsize must be bigger than 0"
    if fontpath is None:
        font = None # use default font later
    else:
        try:
            font = ImageFont.truetype(font, fontsize)
        except Exception as e:
            logging.error("[geneLinkDraw.draw] >>> Could not load font: " + str(e))
            logging.info("[geneLinkDraw.draw] >>> Using default font")
            font = None

    if font is None:
        try:
            # try to locate this module and use the included font in ./font/NugoSansLight-9YzoK.ttf
            modpath = Path(__file__).parent
            fontpath = modpath / "font" / "NugoSansLight-9YzoK.ttf"
            font = ImageFont.truetype(str(fontpath), fontsize)
        except Exception as e:
            logging.error("[geneLinkDraw.draw] >>> Could not load default font: " + str(e))
            logging.info("[geneLinkDraw.draw] >>> Using PIL default font") # fallback to default font
            font = ImageFont.load_default(fontsize)

    if genecols is not None:
        assert len(genecols) == len(genes), "[ERROR] >>> genecols must have same length as genes"
        for i in range(genecols):
            genes[i]._genecol = genecols[i] # use this hack to store the color in the gene object
    
    if linkcols is not None:
        assert len(linkcols) == len(links), "[ERROR] >>> linkcols must have same length as links"

    # default gene coloring are darkblue and darkorange, start at blueviolet for add. colors
    palette = Palette(2) if genecols is None else Palette() 

    elemtypes = set([elemtype for gene in genes for elemtype in gene.elements])
    if elementcols is None:
        elementcols = {elemtype: palette.colorpp() for elemtype in elemtypes}
    assert all([e in elemtypes for e in elementcols]), \
        "[ERROR] >>> Not all elementColors match an Gene.elements type"
    assert all(e in elementcols for e in elemtypes), \
        "[ERROR] >>> Not all Gene.elements types have a corresponding elementColor"
    
    sitetypes = set([sitetype for gene in genes for sitetype in gene.sites])
    if sitecols is None:
        sitecols = {sitetype: palette.colorpp() for sitetype in sitetypes}
    assert all([s in sitetypes for s in sitecols]), \
        "[ERROR] >>> Not all siteColors match an Gene.sites type"
    assert all(s in sitecols for s in sitetypes), \
        "[ERROR] >>> Not all Gene.sites types have a corresponding siteColor"
        

    # determine number of gene rows based on species and set gene coordinates
    # -----------------------------------------------------------------------

    generows = []
    speciesToRow = {} # {"species": row_num, ...}
    for gene in genes:
        if gene.species not in speciesToRow:
            speciesToRow[gene.species] = len(generows)
            generows.append([])

        generows[speciesToRow[gene.species]].append(gene)

    # special case: single species, draw all genes stacked
    if len(generows) == 1:
        generows = [[g] for g in generows[0]]

    # estimate pixel per basepair resolution to roughly fit all genes in the desired width
    maxgenelen = max([sum([g.len for g in rgenes]) for rgenes in generows])
    maxgenemargin = max([len(rgenes) for rgenes in generows])
    genepix = width - (2 * outerMargin) - ((maxgenemargin-1) * genewidth//2)
    res = genepix / maxgenelen # pixel per basepair (horizontal)

    killswitch = 0
    while True and killswitch < 10:
        killswitch += 1 # in case something goes wrong, avoid infinite loop
        x0 = outerMargin
        y0 = outerMargin
        for row in generows:
            for gene in row:
                genecol = gene._genecol if genecols is not None else None
                GeneDrawInfo(gene, x0, y0, gene.id, genewidth, font, fontsize, res, 
                             elementcols, sitecols, geneColor=genecol) # should add itself to gene

            optimizeGeneRow(row, genewidth) # optimize gene label positions
            y0 = max([g.drawInfo.y1_gene for g in row]) + fontsize # get next row start y coordinate

        # get current true width and height
        maxy = max([g.drawInfo.y1_gene for r in generows for g in r])
        trueheight = maxy + outerMargin
        truewidth = max([max(max(g.drawInfo.x1_gene, g.drawInfo.x1_label) for g in r) for r in generows]) + outerMargin

        if truewidth <= width and trueheight <= height:
            break
        elif not forceDimensions:
            break
        else:
            wres = res * width / truewidth if truewidth > width else res # new resolution to fit width
            hres = res * height / trueheight if trueheight > height else res # new resolution to fit height
            res = min(wres, hres) # set resolution to the smaller of the two and try again
            
    if not forceDimensions:
        logging.warning("[geneLinkDraw.draw] >>> Image dimensions adjusted to fit all genes")
        width = truewidth
        height = trueheight

    # draw image
    # ----------

    img = Image.new(mode = "RGB", size = (width, height), color = "white")
    drw = ImageDraw.Draw(img) # drawing context    
    
    # draw genes and elements
    siteradius = math.ceil(genewidth/3) #+ 1
    for row in generows:
        for gene in row:
            gdi: GeneDrawInfo = gene.drawInfo
            # gene bar
            drw.rectangle((gdi.x0_gene, gdi.y0_gene, gdi.x1_gene, gdi.y1_gene), fill=gdi.geneColor, outline=None, 
                          width=1)
            # gene elements
            for elemtype in gdi.elemCoords:
                for elem in gdi.elemCoords[elemtype]:
                    drw.rectangle(elem, fill=gdi.elementColors[elemtype], outline=None, width=1)
            # gene sites
            for sitetype in gdi.siteCoords:
                for site in gdi.siteCoords[sitetype]:
                    drw.ellipse((site[0]-siteradius, site[1]-siteradius, site[0]+siteradius, site[1]+siteradius), 
                                fill=gdi.siteColors[sitetype], outline=None, width=1)
            # gene label
            drw.text(xy=(gdi.x0_label, gdi.y0_label), text=gene.id, font=gdi.font, fill="black")

    # draw links
    if links is not None:
        lcol = palette.colorpp() if linkcols is None else None
        radius = math.ceil(linkwidth/2) + 1
        geneToRow = {gene.id: r for r, row in enumerate(generows) for gene in row}
        for li, link in enumerate(links):
            lcol = linkcols[li] if linkcols is not None else lcol
            anchorsByRow = {}
            for i in range(len(link.genes)):
                gene = link.genes[i]
                gid = gene.id
                gdi: GeneDrawInfo = gene.drawInfo
                r = geneToRow[gid]
                if r not in anchorsByRow:
                    anchorsByRow[r] = []
                if link.compressed:
                    logging.debug(f"[INFO] >>> Compressed link: {link.genes[i]} -> {link.pos[i]}")
                    logging.debug(f"[INFO] >>> Link anchors: {gdi}")
                    anchorsByRow[r].extend([gdi.linkAnchors[p] for p in link.pos[i]])
                else:
                    anchorsByRow[r].append(gdi.linkAnchors[link.pos[i]])

                for a in anchorsByRow[r]:
                    drw.ellipse((a[0]-radius, a[1]-radius, a[0]+radius, a[1]+radius), fill=lcol, outline=lcol, width=1)

            if link.connect:
                lrows = sorted(anchorsByRow.keys())
                for i in range(len(lrows)-1):
                    for a1 in anchorsByRow[lrows[i]]:
                        for a2 in anchorsByRow[lrows[i+1]]:
                            drw.line((a1, a2), fill=lcol, width=linkwidth)

    if show:
        img.show()

    return img



# user helper functions

def getAvailableFonts():
    """ Returns and prints a list of font paths found on your system that you can use with the draw() function """
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
        fakegenes.append(Gene(id=f"{i} - {name} // {code}",
                              species = str(i // 4), # four colors per row
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