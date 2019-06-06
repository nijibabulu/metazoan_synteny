#! /usr/bin/env python3

from __future__ import print_function
import scprep
import magic
import click
import matplotlib.pyplot as plt

@click.command()
@click.argument("COUNT_TABLE")
@click.argument("OUT_FILE")
def main(count_table, out_file):
    print("started main")
    data = scprep.io.load_csv(count_table, cell_axis='column', delimiter='\t')
    print("loaded csv")

    # normalize with our method
    gtot = data.apply(sum, 0)
    ctot = data.apply(sum, 1)

    data_filt = data.loc[ctot >= 200, gtot >= 0]

    totu = data_filt.apply(lambda c: max(1,sum(c)), 1)
    data_norm = data_filt.div(totu, axis=0) * 1000

    print(data_norm.apply(sum, 0).head())
    print("normalized.")

    magic_op = magic.MAGIC()
    fig, ax = plt.subplots()
    magic_op.fit_transform(data_norm, plot_optimal_t=True, ax=ax)
    plt.savefig(out_file)


if __name__ == '__main__':
    main()
