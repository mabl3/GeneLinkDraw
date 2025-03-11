import logging
import unittest
import geneLinkDraw as gld

logging.basicConfig(level=logging.DEBUG)

class TestGeneLinkDraw(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_geneLinkDraw(self):
        genes = [
            gld.Gene('someverylonggenenname_'*10, 'speciesA', 2000, '+'),
            gld.Gene('gene2', 'speciesA', 10, '-'),
            gld.Gene('gene3'*20, 'speciesA', 500, '+'),
            gld.Gene('gene4', 'speciesA', 1000, '+'),
            gld.Gene('gene5', 'speciesB', 1500, '-'),
            gld.Gene('gene6', 'speciesC', 3000, '+'),
            gld.Gene('gene5a', 'speciesB', 1000, '-'),
        ]
        genes[0].addElement("element1", 100, 200)
        genes[0].addElement("element2", 500, 1000)
        genes[1].addElement("element1", 5, 10)
        genes[1].addElement("element2", 10, 20)
        genes[2].addElement("element2", 300, 400)
        genes[3].addElement("element1", 1000, 1000)
        genes[4].addElement("element2", 150, 250)
        genes[5].addElement("element2", 1050, 2050)
        genes[5].addElement("element1", 2500, 2501)
        genes[0].addSite("site1", 100)
        genes[0].addSite("site2", [500, 1000])
        genes[1].addSite("site1", 5)
        genes[1].addSite("site2", 10)
        genes[2].addSite("site2", [300, 400])
        genes[3].addSite("site1", [1000])
        genes[4].addSite("site2", [150, 250])

        links = [
            gld.Link(genes=[genes[0], genes[1], genes[6], genes[5]], 
                     pos=[[123, 234, 345], [5], [900], [2700]], compressed=True),
            gld.Link(genes=[genes[2], genes[5]],
                     pos=[300, 940], compressed=False),
            gld.Link(genes=[genes[2], genes[5]],
                     pos=[350, 990], compressed=False, connect=False),
        ]
        img = gld.draw(genes, links, genewidth=10, forceDimensions=False)
        img.save('test_geneLinkDraw.png')
        self.assertTrue(True)