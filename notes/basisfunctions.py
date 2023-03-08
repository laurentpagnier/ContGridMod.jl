#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


def f1(x, y):
    return np.copy(x)


def f2(x, y):
    return np.copy(y)


def f3(x, y):
    return np.copy(1 - x - y)


def main():
    n = 10_001
    x, y = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
    excl = np.where(x > (1 - y))
    p1 = f1(x, y)
    p2 = f2(x, y)
    p3 = f3(x, y)
    p1[excl] = np.nan
    p2[excl] = np.nan
    p3[excl] = np.nan

    fig, axs = plt.subplots(
        1, 4, width_ratios=[1, 1, 1, 0.05], figsize=(9, 2.75))

    axs[0].pcolormesh(x, y, p1, cmap="cividis")
    axs[1].pcolormesh(x, y, p2, cmap="cividis")
    plot = axs[2].pcolormesh(x, y, p3, cmap="cividis")
    axs[0].text(0.95, 0.95, "$u_1$", ha="right", va="top")
    axs[1].text(0.95, 0.95, "$u_2$", ha="right", va="top")
    axs[2].text(0.95, 0.95, "$u_3$", ha="right", va="top")
    plt.colorbar(plot, cax=axs[3], label="$u_i(x,y)$")

    for ax in axs[:3]:
        ax.set_aspect('equal')
        ax.set_xticks([0, 1], minor=False)
        ax.set_yticks([0, 1], minor=False)
        ax.scatter([1/6, 1/6, 2/3], [1/6, 2/3, 1/6], c="k", marker='x')

    fig.tight_layout()
    fig.savefig("basisfunctions.png")
    # plt.show()


if __name__ == "__main__":
    main()
