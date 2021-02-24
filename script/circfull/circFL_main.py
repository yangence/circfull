#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
author: Zelin Liu
email: zlliu@bjmu.edu.cn
license: GPL3
detail: Analyze circFL-seq data
"""
from docopt import docopt
import sys,warnings
from circfull.version import __version__
warnings.filterwarnings("ignore")
#__version__=0.01

helpInfo = '''
Usage: circfull <command> [options]
Command:
    RG             Reference guide detection
    DNSC           De novo self-correction to get consensus sequence
    cRG            Use RG algorithem with sequences corrected by DNSC
    strand         Identify the strand origion of circFL-seq reads
    sRG            RG algorithem with stranded circFL-seq reads.
    mRG            Merge circFL results from RG and cRG/sRG.
    anno           Annotate full-length circRNA.
    geneExp        Quantify gene expression to evaluate linear RNA residual
    var            Calling variants in circRNA.
'''

def main():
    command_log = 'circfull parameters: ' + ' '.join(sys.argv)
    if len(sys.argv) == 1:
        sys.exit(helpInfo)
    elif sys.argv[1] == '--version' or sys.argv[1] == '-v':
        sys.exit(__version__)
    elif sys.argv[1] == 'RG':
        import circfull.RG as RG
        RG.RG(docopt(RG.__doc__, version=__version__))
    elif sys.argv[1] == 'DNSC':
        import circfull.DNSC as DNSC
        DNSC.DNSC(docopt(DNSC.__doc__, version=__version__))
    elif sys.argv[1] == 'cRG':
        import circfull.cRG as cRG
        cRG.cRG(docopt(cRG.__doc__, version=__version__))
    elif sys.argv[1] =='mRG':
        import circfull.mRG as mRG
        mRG.mRG(docopt(mRG.__doc__, version=__version__))
    elif sys.argv[1] =='sRG':
        import circfull.sRG as sRG
        sRG.sRG(docopt(sRG.__doc__, version=__version__))
    elif sys.argv[1] == 'strand':
        import circfull.strand as strand
        strand.strand(docopt(strand.__doc__, version=__version__))
    elif sys.argv[1] == 'geneExp':
        import circfull.geneExp as geneExp
        geneExp.geneExp(docopt(geneExp.__doc__, version=__version__))
    elif sys.argv[1] == 'anno':
        import circfull.annoFL as annoFL
        annoFL.annoFL(docopt(annoFL.__doc__, version=__version__))
    elif sys.argv[1] == 'var':
        import circfull.var as var
        var.var(docopt(var.__doc__, version=__version__))
    else:
        sys.exit(helpInfo)

if __name__ == '__main__':
    main()