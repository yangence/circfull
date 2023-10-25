'''
Usage: circfull FL2BED -c circ [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -c circ                     circFL full-length file.
    -o output                   Output file [default: FL2bed.bed].
'''
import pandas as pd
from .RG_circFL_output import FL2bed
def FL2BED(options):
    circ_file=options['-c']
    output_file=options['-o']
    circ_df=pd.read_csv(circ_file,sep='\t')
    if circ_df.columns[0]!='circID':
    	circ_df=pd.read_csv(circ_file,sep='\t',header=None)
    	FL_colname=['circID','isoID','chr','start','end','len','exonNum','exon_start','exon_end','motif','leftSeq','rightSeq','exon_leftSeq','exon_rightSeq','strand','geneName','readCount','readID']
    	circ_df.columns=FL_colname[0:circ_df.shape[1]]
    #circ_df=filterCircFL(circ_df)
    circ_df_bed=FL2bed(circ_df)
    circ_df_bed.to_csv(output_file,header=None,index=None,sep='\t')
