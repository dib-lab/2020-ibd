#! /usr/bin/env python
"""
Given a catlas, extract a dominator to cDBG node map for level 1 dominators. 
Write to a csv file.
"""
import os
import sys
import argparse
from collections import defaultdict
import pandas as pd
from spacegraphcats.search.catlas import CAtlas

def main():
    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument('cdbg_prefix', help='cdbg prefix')
    p.add_argument('catlas_prefix', help='catlas prefix')
    args = p.parse_args()
    
    output = "cdbg_to_pieces.csv"
    output = os.path.join(args.catlas_prefix, output)
    
    catlas = CAtlas(args.cdbg_prefix, args.catlas_prefix)
    layer1_to_cdbg_df = pd.DataFrame([{"dominator":k, "cdbg_node": v} for k,v in catlas.layer1_to_cdbg.items()])
    layer1_to_cdbg_df = layer1_to_cdbg_df.explode('cdbg_node') 
    # write to csv
    layer1_to_cdbg_df.to_csv(output, index = False)

    return 0
    
if __name__ == "__main__":
    sys.exit(main())
