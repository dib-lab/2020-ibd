#! /usr/bin/env python
import csv
import sourmash
import argparse

from sourmash.lca.command_index import load_taxonomy_assignments
from collections import defaultdict


def pop_to_rank(lin, rank):
    "Remove lineage tuples from given lineage `lin` until `rank` is reached."
    lin = list(lin)

    txl = sourmash.lca.taxlist(include_strain=False)
    before_rank = []
    for txl_rank in txl:
        if txl_rank != rank:
            before_rank.append(txl_rank)
        else:
            break

    # are we already above rank?
    if lin and lin[-1].rank in before_rank:
        return tuple(lin)

    while lin and lin[-1].rank != rank:
        lin.pop()

    return tuple(lin)


def get_ident(ident):
    "Hack and slash identifiers."
    ident = ident.split()[0]
    ident = ident.split('.')[0]
    return ident


def summarize_gather_at(rank, tax_assign, gather_results):
    # collect!

    sum_uniq_weighted = defaultdict(float)
    for row in gather_results:
        match_ident = row['name']
        match_ident = get_ident(match_ident)
        lineage = tax_assign[match_ident]
        lineage = pop_to_rank(lineage, rank)
        assert lineage[-1].rank == rank, lineage[-1]

        f_uniq_weighted = row['f_unique_weighted']
        f_uniq_weighted = float(f_uniq_weighted)
        sum_uniq_weighted[lineage] += f_uniq_weighted

    items = list(sum_uniq_weighted.items())
    items.sort(key = lambda x: -x[1])

    for k, v in items:
        print(rank, f'{v:.3f}', sourmash.lca.display_lineage(k))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('gather_csv')
    p.add_argument('lineages_csv')
    args = p.parse_args()

    gather_results = []

    with open(args.gather_csv, 'rt') as fp:
        r = csv.DictReader(fp)
        for row in r:
            gather_results.append(row)
    print(f'loaded {len(gather_results)} gather results.')

    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                          start_column=3)
    print(f'loaded {len(tax_assign)} tax assignments.')

    n_missed = 0
    for row in gather_results:
        match_ident = row['name']
        match_ident = get_ident(match_ident)
        if match_ident not in tax_assign:
            n_missed += 1

    print(f'of {len(gather_results)}, missed {n_missed} lineage assignments.')

    assert n_missed == 0
    print(n_missed)

    for rank in sourmash.lca.taxlist(include_strain=False):
        summarize_gather_at(rank, tax_assign, gather_results)


if __name__ == '__main__':
    main()
