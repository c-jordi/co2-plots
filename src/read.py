import scipy.io
import numpy as np
from matplotlib import pyplot as plt

from .constants import *


def get_capture_costs(struct, r):
    """Get capture components cost.
    """
    costs = {
        "capex": 0.,
        "mainex": 0.,
        "opex": {
            "total": 0.,
            "capture": 0.,
            "condition": 0.
        }
    }
    components = struct["component"][()]["CCS_CO2"][()]
    for i, comp in enumerate(components):
        if comp != [] and "investment_costY" in comp.dtype.names:
            costs["capex"] += comp["investment_costY"].item()
            costs["mainex"] += comp["maintenance_costY"].item()
            try:
                costs["opex"]["capture"] += comp["input"][()][1] * \
                    ELECTRICITY_PRICE/r
            except:
                costs["opex"]["capture"] += comp["input"]*ELECTRICITY_PRICE/r

    nbr_nodes = len(struct["network"][()]["nodeNames"].item())
    a = np.array([1 if comp != [] else 0 for comp in components])
    Import = struct["carrier"][()]["electricity"][()]["import"].item()
    Qtransp = struct["carrier"][()]["electricity"][()]["Q_transp"].item()

    index = np.zeros(nbr_nodes)
    index[:len(a)] = np.array(a)

    try:
        costs["opex"]["total"] = np.sum(Import[:, index == 1],
                                        axis=1)*ELECTRICITY_PRICE/r
        costs["opex"]["condition"] = np.sum(
            Qtransp[:, index == 1], axis=1)*ELECTRICITY_PRICE/r
    except:
        costs["opex"]["total"] = np.sum(Import[index == 1])*ELECTRICITY_PRICE/r
        costs["opex"]["condition"] = np.sum(
            Qtransp[index == 1])*ELECTRICITY_PRICE/r

    return costs


def get_storage_costs(struct, r):
    """Get storage components cost.
    """
    costs = {
        "capex": 0.,
        "mainex": 0.,
        "opex": {
            "total": 0.,
            "capture": 0.,
            "condition": 0.
        }
    }
    components = struct["component"][()]["CO2sinkOffshore"][()]
    for i, comp in enumerate(components):
        if comp != [] and "investment_costY" in comp.dtype.names:
            costs["capex"] += comp["investment_costY"].item()
            costs["mainex"] += comp["maintenance_costY"].item()
            try:
                costs["opex"]["storage"] += comp["input"][()][1] * \
                    ELECTRICITY_PRICE/r
            except:
                costs["opex"]["capture"] += comp["input"]*ELECTRICITY_PRICE/r

    nbr_nodes = len(struct["network"][()]["nodeNames"].item())
    a = np.array([1 if comp != [] else 0 for comp in components])
    Import = struct["carrier"][()]["electricity"][()]["import"].item()
    Qtransp = struct["carrier"][()]["electricity"][()]["Q_transp"].item()

    index = np.zeros(nbr_nodes)
    index[:len(a)] = np.array(a)

    try:
        costs["opex"]["total"] = np.sum(Import[:, index == 1],
                                        axis=1) * ELECTRICITY_PRICE/r
        costs["opex"]["condition"] = np.sum(
            Qtransp[:, index == 1], axis=1) * ELECTRICITY_PRICE/r
    except:
        costs["opex"]["total"] = np.sum(Import[index == 1])*ELECTRICITY_PRICE/r
        costs["opex"]["condition"] = np.sum(
            Qtransp[index == 1])*ELECTRICITY_PRICE/r
    return costs


def get_network_costs(struct, r):
    """Get network cost.
    """
    costs = {
        "opex": {
            "total": 0.,
            "fuel": 0.,
            "elec": 0.,
            "condition": 0.
        },
        "capex": 0.,
        "mainex": 0.,
        "detail": {}
    }
    networks = struct["network"][()]
    nws = [i for i in networks.dtype.names if i != 'nodeNames']
    for nw in nws:
        costs["detail"][nw] = {
            "capex": networks[nw][()]["investment_costY"].item(),
            "opex": networks[nw][()]["operation_costY"].item()
        }
        costs["capex"] += costs["detail"][nw]["capex"]
        costs["opex"]["fuel"] += costs["detail"][nw]["opex"]

    nbr_nodes = len(struct["network"][()]["nodeNames"].item())
    Import = struct["carrier"][()]["electricity"][()]["import"].item()
    Qtransp = struct["carrier"][()]["electricity"][()]["Q_transp"].item()

    # Get ids of nodes that are neither used for capture or storage
    captureIds = np.array(
        [1 if comp != [] else 0 for comp in struct["component"][()]["CCS_CO2"][()]])
    storageIds = np.array(
        [1 if comp != [] else 0 for comp in struct["component"][()]["CO2sinkOffshore"][()]])
    icapture, istorage = np.zeros(nbr_nodes), np.zeros(nbr_nodes)
    icapture[:len(captureIds)] = np.array(captureIds)
    istorage[:len(storageIds)] = np.array(storageIds)
    index = np.ones(nbr_nodes) - np.logical_or(icapture, istorage)

    try:
        costs["opex"]["elec"] = np.sum(
            Import[:, index == 1], axis=1) * ELECTRICITY_PRICE/r
        costs["opex"]["condition"] = np.sum(
            Qtransp[:, index == 1], axis=1) * ELECTRICITY_PRICE/r
    except:
        costs["opex"]["elec"] = np.sum(Import[index == 1])*ELECTRICITY_PRICE/r
        costs["opex"]["condition"] = np.sum(
            Qtransp[index == 1])*ELECTRICITY_PRICE/r

    costs["opex"]["total"] += costs["opex"]["fuel"] + costs["opex"]["elec"]
    return costs


def get_network_emissions(networks):
    """Get network emissions.
    """
    emissions = 0
    nws = [i for i in networks.dtype.names if i != 'nodeNames']
    for nw in nws:
        emissions += networks[nw][()]["emissionsY"].item() + \
            networks[nw][()]["emissions_indY"].item()
    return emissions


def get_plant_details(struct):
    """Get plant details.
    """
    Y = struct["horizon"][()]["number"].item()
    nbr_plants = len(struct["component"][()]["CCS_CO2"][()])
    connected = np.zeros([nbr_plants, Y])
    forgotten_emissions = np.zeros([nbr_plants, Y])
    for i, plant in enumerate(struct["component"][()]["CCS_CO2"][()]):
        connected[i, :] = plant["emissionscaptured"][()] > 1
        forgotten_emissions[i, :] = np.logical_not(
            connected[i, :]).astype(int) * plant["input"][()][0]

    return {
        "connected_to_co2_network": np.sum(connected, axis=0),
        "forgotten_emissions": np.sum(forgotten_emissions, axis=0)
    }


def get_condition_data(struct, r):
    """Get conditioning cost and emission.
    """

    Y = struct["horizon"][()]["number"].item()
    nbr_nodes = len(struct["network"][()]["nodeNames"].item())

    a = np.array(
        [1 if comp != [] else 0 for comp in struct["component"][()]["CCS_CO2"][()]])
    b = np.array([1 if comp != [] else 0 for comp in struct["component"]
                  [()]["CO2sinkOffshore"][()]])

    capt = np.zeros(nbr_nodes)
    capt[:len(a)] = a
    stor = np.zeros(nbr_nodes).astype(np.int)
    stor[:len(b)] = b
    capt = capt.astype(np.bool)
    stor = stor.astype(np.bool)
    other = np.ones(nbr_nodes) - np.logical_or(capt, stor).astype(np.int)

    Import = struct["carrier"][()]["electricity"][()
                                                  ]["import"].item() * ELECTRICITY_PRICE
    Emissions = struct["carrier"][()]["electricity"][()]["emissionsY"].item()
    if Y == 1:
        Import = np.reshape(Import, [1, Import.shape[0]])
        Emissions = np.reshape(Emissions, [1, Emissions.shape[0]])
    return {"opex": {
        "all": np.sum(Import[:, :], axis=1)/r,
        "capture": np.sum(Import[:, capt == 1], axis=1)/r,
        "storage": np.sum(Import[:, stor == 1], axis=1)/r,
        "others": np.sum(Import[:, other == 1], axis=1)/r,

    }
    }, {
        "all": np.sum(Emissions[:, :], axis=1),
        "capture": np.sum(Emissions[:, capt == 1], axis=1),
        "storage": np.sum(Emissions[:, stor == 1], axis=1),
        "others": np.sum(Emissions[:, other == 1], axis=1)
    }


def get_data(path):
    struct = scipy.io.loadmat(path, squeeze_me=True)["object"]

    data = {}

    # Params
    data["Y"] = struct["horizon"][()]["number"].item()
    data["discRate"] = struct["analysis"][()]["eco"][()]["discRate"].item()
    data["r"] = (1 + data["discRate"])**np.arange(0, data["Y"])
    data["modes"] = ", ".join(
        [i for i in struct["network"][()].dtype.names if i != 'nodeNames'])

    # Costs
    data["costs"] = {
        "total": struct["total_costs"].item(),
        "LCOC": struct["LCOC"],
        "capex": struct["investment_costs"].tolist(),
        "capexY": struct["investment_costsY"].tolist(),
        "opex": struct["operation_costs"].tolist(),
        "opexY": struct["operation_costsY"].tolist(),
        "mainex": struct["maintenance_costs"].tolist(),
        "mainexY": struct["maintenance_costsY"].tolist(),
    }

    data["costs"]["yearly"] = (
        data["costs"]["capexY"] + data["costs"]["opexY"] + data["costs"]["mainexY"])

    try:
        data["costs"]["LCSC"] = struct["LCSC"]
        data["costs"]["LCSC_nom"] = struct["LCSC_nom"]
        data["costs"]["LCAC_nom"] = struct["LCAC_nom"]
        data["costs"]["LCAC"] = struct["LCAC"]
    except:
        pass

    data["costs"]["capture"] = get_capture_costs(struct, data["r"])
    data["costs"]["storage"] = get_storage_costs(struct, data["r"])
    data["costs"]["network"] = get_network_costs(struct, data["r"])

    # Emissions
    data["emissions"] = {
        "wte": struct["industrialEmissions"].item(),
        "stored": np.sum(struct["sink_inputY"]),
        "emitted": struct["emissionsY"].item(),
        "captured": struct["emissionsCaptured"].item(),
        "network": get_network_emissions(struct["network"][()]),
        "scenario": struct["analysis"][()]["CO2"][()]["maxEmissionsY"].tolist()
    }
    data["emissions"]["avoided"] = data["emissions"]["wte"] - \
        data["emissions"]["emitted"]

    # Combined
    data["costs"]["condition"], data["emissions"]["condition"] = get_condition_data(
        struct, data["r"])

    data["network"] = {
        "plants": get_plant_details(struct)
    }

    # Risk
    data["risk"] = {}

    try:
        data["risk"]["ECNS"] = struct["risk"][()]['ECNS'][()]
        data["risk"]["ECS"] = struct["risk"][()]["ECS"][()]
        data["risk"]["NECNS"] = struct["risk"][()]["NECNS"][()]
        data["risk"]["LCSC"] = struct["risk"][()]["LCSC"][()]
        data["risk"]["RFEE"] = struct["risk"][()]["RFEE"][()]
        data["risk"]["p_nofail"] = struct["risk"][()]["p_nofail"][()]
    except ValueError:
        pass

    try:
        data["risk"]["f_emissions"] = struct["f_emissions"][()]
        data["risk"]["f_emissionsdY"] = struct["f_emissionsdY"][()]
    except ValueError:
        pass

    # Optional
    data["struct"] = struct

    return data


def check_validity(path):
    """Checks if the output file is valid and is not the result of unfeasible problem
    """
    struct = scipy.io.loadmat(path, squeeze_me=True)["object"]
    diagnostic = struct["diagnostic"][()]
    info = diagnostic["info"]
    problem = diagnostic["problem"]
    if problem == 1:
        print(f"> {info}: " + path)
        return False
    return True
