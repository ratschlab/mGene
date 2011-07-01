"""
#######################################################################################
#                                                                                     #
#    sdimapper.py is a command-line tool to map sdi coordinates to existing gff       #
#    annotation coordinates                                                           #
#                                                                                     #
#    Copyright (C) 2010 Sebastian J. Schultheiss <sebi@umich.edu>                     #
#                                                                                     #
#    This program is free software; you can redistribute it and/or modify             #
#    it under the terms of the GNU General Public License as published by             #
#    the Free Software Foundation; either version 3 of the License, or                #
#    (at your option) any later version.                                              #
#                                                                                     #
#    This program is distributed in the hope that it will be useful,                  #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                   #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                     # 
#    GNU General Public License for more details.                                     #
#                                                                                     #
#    You should have received a copy of the GNU General Public License                # 
#    along with this program; if not, see http://www.gnu.org/licenses                 #
#    or write to the Free Software Foundation, Inc., 51 Franklin Street,              #
#    Fifth Floor, Boston, MA 02110-1301  USA                                          #
#                                                                                     #
#######################################################################################
#                                                                                     #
#  Please add a notice of any modifications here:                                     #
#                                                                                     #
#                                                                                     #
#######################################################################################
"""

__version__ = "0.1.0"
__license__  = "GNU General Public License"
__author__ = """Sebastian J. Schultheiss <sebi@umich.edu>"""
__usage__ = "%prog " + __version__ + """
  %prog maps coordinates from sdi files to a gff annotation, 
  updating it to reflect the new coordinates. 
  Usage: 
        %prog -d directory/to/chromosomes/ [options]"""
        
from Inputs import Gff3File
import sys
from optparse import OptionParser
from pickle import load









class SDIMapper(object):
    """Creates a mapper object that loads the sdi coordinates file 
    so you only have to wait once for that to complete -- send 
    multiple lookups afterwards without delay."""
    
    def __init__(self, sdipckfile = None):
        """Initialize the coordinates file"""    
        self.chromosomes = {}
        if sdipckfile is not None:
            self.loadSDIpickle(sdipckfile)
  
    def loadSDIpickle(self, sdipckfile):
        """Unpickle a coordinates file and load the index into memory"""
        try: 
            pckfilehandle = open(sdipckfile)
        except:
            print "Could not find or access specified SDI pickle file."
            raise
        try:
            chroms = load(pckfilehandle)
        except:
            print "SDI pickle file seems to be damaged."
            raise
        for chrom in chroms:
            self.chromosomes[chrom.chromno] = chrom
            
    def map(self, chromno, position):
        """Returns the chromosome name and corresponding position of the
        original sequence file at any specified position in the 
        virtual chromosome."""
        return self.chromosomes[int(chromno)].map(int(position))

    def mapMultiple(self, poslist):
        """Uses a list of tuples (chromosomenumber, position)
        and maps each back to its original file and position"""
        results = []
        for chromno, position in poslist:
            results.append(self.map(chromno, position))
        return results

def optionparse(parser):
    """Parse all command line options and set defaults"""
    parser.add_option("-s", "--sdi", dest = "sdifile", type = "string", 
                      help = "Path to SDI file")
    parser.add_option("-g", "--gff", dest = "gfffile", type = "string", 
                      help = "GFF file with original coordinates")
    parser.add_option("-o", "--output", dest = "output", type = "string", 
                      help = "Output file to write with new coordinates [default %default]")
    parser.add_option("-v", "--verbose", dest = "verbose", action = "store_true",
                      help = "Detailed program output on [default: %default]")
    parser.set_defaults(output = "new.gff",
                        verbose = False)

def main(argv = None):     
    """    
    #####################
    # the main function #
    #####################"""
    if argv is None:         
        argv = sys.argv    
    parser = OptionParser(version = "%prog " + __version__, usage = __usage__)
    optionparse(parser)
    (options, args) = parser.parse_args()
    if options.verbose:
        print "Reading SDI file", options.sdifile
    chroms = {}
    sdi = open(options.sdifile)
    for line in sdi:
        array = line.split('\t')
        chrom = array[0]
        if chrom not in chroms:
            chroms[chrom] = {}
        offset = int(array[2])
        if offset is not 0:
            chroms[chrom][int(array[1])] = offset
    if options.verbose:
        print "Reading GFF file", options.gfffile
    gff = Gff3File()
    max = gff.read(options.gfffile)
    if options.verbose:
        print "Adjusting coordinates"
    mchrom = {}    
    for chrom, offsets in chroms.iteritems():
        mapc = [0] * (max + 2)
        globaloffset = 0
        #waxingoffset = 0
        #waningoffset = 0
        for i in xrange(1, max + 1):
            if i in offsets:
                globaloffset += offsets[i]
            #if diff > 0: 
            #    globaloffset += diff
            #if diff < 0:
            #    waxingoffset.append(i - diff)
            mapc[i] = globaloffset
        mchrom[chrom] = mapc
    if options.verbose:
        print "Writing output to", options.output
    for line in gff.gff:
        chrom = line.getSeqId()
        start = line.getStart()
        end = line.getEnd()
        newstart = start + mchrom[chrom][start]
        newend = end + mchrom[chrom][end]
        line.setStartEnd(newstart, newend)
    gffout = open(options.output, "w")
    gffout.write(gff.getFormattedGff())
    gffout.close()
        
if __name__ == "__main__":
    main()