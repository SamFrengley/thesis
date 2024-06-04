import HMS.ZNr as ZNr
import matplotlib.pyplot as plt


def cusp_sketch(Nr):
    for N, r in Nr:
        Z = ZNr.Ztil(N, r)
        cc = Z.sketch_cusps(
            tex=True, compile_tex=True, open_pdf=True, disp=False
        )
    return cc


def main():
    Nr = [(N, 1) for N in range(15,19)]
    cusp_sketch(Nr)
    plt.close("all")


if __name__ == "__main__":
    main()