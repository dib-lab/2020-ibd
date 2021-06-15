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
    cdbg_to_pieces = defaultdict(set)
    for node_id in catlas:
        level = catlas.levels[node_id]
        if level == 1:
            pieces = catlas.layer1_to_cdbg.get(node_id)
            for cdbg_node in pieces:
                cdbg_to_pieces[cdbg_node] = set(pieces)

    cdbg_to_pieces_df = pd.DataFrame([{"dominator":k, "cdbg_node": v} for k,v in cdbg_to_pieces.items()])
    # make dict lines one per line
    cdbg_to_pieces_df = cdbg_to_pieces_df.explode('cdbg_node') 
    # write to csv
    cdbg_to_pieces_df.to_csv(output, index = False)

    return 0
    
if __name__ == "__main__":
    sys.exit(main())
