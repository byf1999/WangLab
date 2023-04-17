import numpy as np
import matplotlib.pyplot as plt


def check_seq_depth(seq: list, sep=int(1e5)):
    seq = np.array(seq)
    seq = seq[seq != 0]
    prop = seq / sum(seq)
    Xs = list(range(1, sum(seq), int(sep)))
    Ys = [np.mean(np.count_nonzero(np.random.multinomial(n, prop, size=10), axis=1)) for n in Xs]

    ax = plt.subplot()
    ax.scatter(Xs, Ys)
    ax.set_xlabel('# of insertions',size=12)
    ax.set_ylabel('# of unique insertions',size=12)
    plt.show()
    return ax
