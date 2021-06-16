import numpy as np
from matplotlib import pyplot as plt


def cost_breakdown(data):
    """
    Plots the cost breakdown of a value chain and saves an svg in the output folder.

    Args:
        data (obj) : Output of the src.read.get_data function

    """

    # INIT
    N = data["Y"]
    ind = np.arange(N)  # the x locations for the groups
    fig = plt.figure(figsize=(12, 12))
    pms = {  # params
        "sp1": {
            "width": .5,
            "colors": ['#e7e1ef', '#c994c7', '#dd1c77', 'grey'],
            "scale": 1e-6
        },
        "sp2": {
            "colors": ["red", "purple"],
            "scale": 1e-6
        },
        "sp3": {
            "width": 1/7,
            "colors": ["blue", "yellow", "green"]
        },
        "sp4": {
            "width": 1/7,
            "colors": ["#a8322d", "#edf8b1"]
        }
    }

    # SUBPLOT 1 - EMISSIONS
    w, c, s = pms["sp1"]["width"], pms["sp1"]["colors"], pms["sp1"]["scale"]
    sp1 = fig.add_subplot(411)
    sp1b = sp1.twinx()
    sp1.xaxis.grid(True)
    sp1.set_ylabel("Emissions [Mt$\textrm{CO}_2$/yr]")
    sp1.set_xticks(ind)
    sp1.set_xlim([0, N])
    sp1.set_xticklabels(["" for i in ind])
    sp1.set_ylim([0, data["emissions"]["wte"] * 6/5 * s])

    for y, capt in enumerate(data["emissions"]["captured"]):
        _ = sp1.hlines(
            capt/1e6, y, y+1, linestyle='-', color="red", zorder=0)
        if y == 0:
            _.set_label("Emissions Captured")
    sp1.bar(ind+1/2, (data["emissions"]["wte"] - data["emissions"]["stored"] + data["emissions"]["network"] + data["emissions"]
                      ["condition"]["capture"] + data["emissions"]["condition"]["storage"])*s, w, color=c[2], label="Storage Emissions", zorder=1)
    sp1.bar(ind+1/2, (data["emissions"]["wte"] - data["emissions"]["stored"] + data["emissions"]["network"] +
                      data["emissions"]["condition"]["capture"])*s, w, color=c[1], label="Network Emissions", zorder=1)
    sp1.bar(ind+1/2, (data["emissions"]["wte"] - data["emissions"]["stored"] + data["emissions"]
                      ["condition"]["capture"])*s, w, color=c[0], label="Capture Emissions", zorder=1)
    sp1.bar(ind+1/2, (data["emissions"]["wte"] - data["emissions"]
                      ["stored"])*s, w, color=c[-1], label="WtE Emissions", zorder=1)

    sp1b.set_yticks([data["emissions"]["wte"] / 5 * s * i for i in range(6)])
    sp1b.set_yticklabels([f"{20 * i}%" for i in range(6)])
    sp1b.set_ylabel("% of total WtE emissions")
    sp1b.set_ylim([0, data["emissions"]["wte"] * 6/5 * s])
    sp1.legend()

    # SUBPLOT 2 - COSTS
    c = pms["sp2"]["colors"]
    sp2 = fig.add_subplot(412)

    if "LCAC" in data["costs"] and "LCSC" in data["costs"]:
        LCAC, LCSC = data["costs"]["LCAC"], data["costs"]["LCSC"]
    else:
        LCAC = data["costs"]["yearly"] * \
            data["r"] / data["emissions"]["avoided"]
        LCSC = data["costs"]["yearly"] * \
            data["r"] / data["emissions"]["stored"]
    if type(LCAC) != list:
        LCAC, LCSC = [LCAC], [LCSC]
    for y, cost in enumerate(LCAC):
        _ = sp2.hlines(cost, y, y+1, linestyle='--', color="purple")
        if y == 0:
            _.set_label("Carbon avoided cost")
    for y, cost in enumerate(LCSC):
        _ = sp2.hlines(cost, y, y+1, color="red")
        if y == 0:
            _.set_label("Carbon stored cost")
    sp2.xaxis.grid(True)
    sp2.set_ylabel("Cost [EUR/tCO2]")
    sp2.set_xticks(ind)
    sp2.set_xlim([0, N])
    sp2.set_xticklabels(["" for i in np.arange(N)])
    sp2.legend()

    # SUBPLOT 3 - COSTS
    w, c = pms["sp3"]["width"], pms["sp3"]["colors"]
    sp3 = fig.add_subplot(413)

    def lvl(cost):
        return cost * data["r"] / data["emissions"]["stored"]

    capt = {
        "C": lvl(data["costs"]["capture"]["capex"]),  # capital
        "O": lvl(data["costs"]["capture"]["opex"]["total"]),  # operational
        "M": lvl(data["costs"]["capture"]["mainex"])  # maintenance
    }

    trans = {
        "C": lvl(data["costs"]["network"]["capex"]),  # capital
        "O": lvl(data["costs"]["network"]["opex"]["total"]),  # operational
        "M": lvl(data["costs"]["network"]["mainex"])  # maintenance
    }

    stor = {
        "C": lvl(data["costs"]["storage"]["capex"]),  # capital
        "O": lvl(data["costs"]["storage"]["opex"]["total"]),  # operational
        "M": lvl(data["costs"]["storage"]["mainex"])  # maintenance
    }

    # Capture
    sp3.bar(ind + 3/14, capt["C"] + capt["O"] +
            capt["M"], w, color=c[2], label="Maintenance")
    sp3.bar(ind + 3/14, capt["C"] + capt["O"],
            w, color=c[1], label="Operation")
    sp3.bar(ind + 3/14, capt["C"], w, color=c[0], label="Investment")
    # Transport
    sp3.bar(ind + 1/2,  trans["C"] + trans["O"] +
            trans["M"], w, color=c[2])
    sp3.bar(ind + 1/2, trans["C"] + trans["O"], w, color=c[1])
    sp3.bar(ind + 1/2, trans["C"], w, color=c[0])
    # Storage
    sp3.bar(ind + 11/14, stor["C"] + stor["O"] +
            stor["M"], w, color=c[2])
    sp3.bar(ind + 11/14, stor["C"] + stor["O"], w, color=c[1])
    sp3.bar(ind + 11/14, stor["C"], w, color=c[0])

    sp3.set_ylabel('Cost [EUR/tCO2]')
    sp3.xaxis.grid(True)
    sp3.set_xticks(ind)
    sp3.set_xticklabels(["" for i in ind])
    sp3.set_xlim([0, N])
    sp3.legend()

    # SUBPLOT 4 - CONDITIONING
    w, c = pms["sp4"]["width"], pms["sp4"]["colors"]
    sp4 = fig.add_subplot(414)
    sp4.bar(ind + 3/14, capt["O"], w,
            color=c[0], label="Conditioning")
    sp4.bar(ind + 3/14, capt["O"] - lvl(data["costs"]["capture"]
                                        ["opex"]["condition"]), w, color=c[1], label="Other")
    sp4.bar(ind + 1/2, trans["O"], w,
            color=c[0])
    sp4.bar(ind + 1/2, trans["O"] - lvl(data["costs"]["network"]
                                        ["opex"]["condition"]), w, color=c[1])
    sp4.bar(ind + 11/14, stor["O"], w,
            color=c[0])
    sp4.bar(ind + 11/14, stor["O"] - lvl(data["costs"]["storage"]
                                         ["opex"]["condition"]),  w, color=c[1])

    sp4.set_xlim([0, N])
    sp4.set_ylabel('Cost [EUR/tCO2]')
    sp4.set_xlabel('Years')
    sp4.set_xticks(ind+1/2)
    sp4.set_xlim([0, N])
    sp4.set_xticklabels([f"{2025 + i}" if not (i % 5 and i !=
                                               0 and i != N-1) else "" for i in np.arange(0, N)])
    sp4.legend()

    plt.savefig("output/cost-breakdown.svg")

    plt.show()


def R1_pareto(*datalists: list):
    """
    Plots the R1-cost pareto for multiple value chains and saves svg in output folder.

    Args: (one or more)
        datalist (list of [str, list(src.read.get_data)])
    """
    plt.figure(figsize=(12, 7))
    for dtl in datalists:
        points = [[float(samp["costs"]["LCOC"]), float(samp["risk"]["ECNS"]), round(
            samp["emissions"]["emitted"])] for samp in dtl[1]]
        points = sorted(points, key=lambda x: x[1])
        x = np.array([max(p[0], 1) for p in points])
        y = np.array([max(p[1], 0) for p in points])
        z = np.array([p[2] for p in points])
        r1 = 1 - y/np.max(y)
        plt.plot(x, r1, linewidth=1, marker="o", linestyle="dashed",
                 markersize=8, label=dtl[0])
    plt.xlabel("LCOC (EUR/tCO2)")
    plt.ylabel("$R_1$")
    plt.legend()
    plt.savefig("output/r1-pareto.svg")
    plt.show()
