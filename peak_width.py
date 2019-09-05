import argparse
import sys
import logging
import pandas as pd

def adjust_width(peak, width):
    interval = int((int(peak[2]) + int(peak[1]))/2)
    peak[2] = interval + int(width/2)
    peak[1] = interval - int(width/2)
    return peak
    

def main():
    #################### Setting arguments ########################
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-i","--inputBed", help="An input bed file to adjust its intervals")
    arg_parser.add_argument("-w","--width",default=500,help="A width of intervals")
    
    args = arg_parser.parse_args()

    if args.inputBed == None:
        logger.error("'-i' : An input file was missed.")
        exit()
    ###############################################################
    
    file_name = args.inputBed
    peak_width = int(args.width)
    
    peaks = pd.read_csv(file_name, sep='\t', header=None)
    adjusted_peaks = peaks.apply(lambda x: adjust_width(x, peak_width), axis=1)
    
    adjusted_peaks.to_csv(sys.stdout, index=False, header=False, sep='\t')

    
if __name__ == '__main__':
    logger = logging.getLogger("ConvLog")
    logger.setLevel(logging.DEBUG)               # The logger object only output logs which have
                                                # upper level than INFO.
    log_format = logging.Formatter('%(asctime)s:%(message)s')

    stream_handler = logging.StreamHandler()    # Log output setting for the command line.
    stream_handler.setFormatter(log_format)     # The format of stream log will follow this format.
    logger.addHandler(stream_handler)

    #file_handler = logging.FileHandler()        # Log output setting for the file.
    #logger.addHandler(file_handler)

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("USER INTERRUPT. \n")
        sys.exit()