import requests
import json
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, rcParams
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import os
import csv
import numpy as np

class OrgSolvent():
    """Organic solvent used in liquid-liquid extraction in a flow system"""

    def __init__(self, smiles='', name='orgSlvt'):
        self.smiles = smiles
        self.name = name

    def getSmiles(self):
        return self.smiles

    def getName(self):
        return self.name
        # if not self.name:
        #     return self.smiles
        # else:
        #     return self.name


class OrgPhase():
    """Organic phase with two solvents (OrgSolvent) used in liquid-liquid extraction in a flow system"""

    def __init__(self, solvent1=None, solvent2=None, x2=0.0):
        self.one = solvent1
        self.two = solvent2
        self.frac = x2

    def getOne(self):
        return self.one

    def getTwo(self):
        return self.two

    def getFrac(self):
        return self.frac

    def getName(self):
        if self.two:
            return self.one.getName() + '-' + self.two.getName() + '({})'.format(self.frac)
        else:
            return self.one.getName()


class AquPhase():
    """Aqueous phase used in liquid-liquid extraction in a flow system"""

    def __init__(self, pH=7.0, ionicStrength=None):
        self.ionic = ionicStrength
        self.ph = pH

    def getpH(self):
        return self.ph

    def getIonicStrength(self):
        return self.ionic


def solventRepr(org):
    """
    :param org: an OrgPhase object, can be a simple solvent or binary solvents
    :return: A dict for the solvent part for the task
    """
    if (not org.getTwo()) or org.getFrac() == 0.0:
        solvent = {
          "smiles": org.getOne().getSmiles(),
          "solvents": None,
          "name": org.getOne().getName(),
          "mixedSolvent": False
        }
    else:
        solvent = {
          "smiles": None,
          "solvents": [{
            "volumeFraction": 1.0 - float(org.getFrac()),
            "smiles": org.getOne().getSmiles(),
            "name": org.getOne().getName()
          },
          {
            "volumeFraction": org.getFrac(),
            "smiles": org.getTwo().getSmiles(),
            "name": org.getTwo().getName()
          }],
          "name": org.getName(),
          "mixedSolvent": True
        }
    return solvent


def opt_cond(tgtArry, conds, idx=0, maxi=True):
    """
    :param tgtArry: the array for the target property
    :param conds: a list of the conditions (list or array)
    :param idx: the idx for the target molecule, default=0 for the SeparationDesigner
    :param maxi: whether to choose the maximum in tgtArray, default=True
    :return: the list of optimal conditions [opt cond A (np.array), opt cond B], [A1, B1] conditions for optimum 1
    """
    tgt = tgtArry[idx] # data for the target molecule
    ts = tgt.shape
    n_err = 0
    for c, cond in enumerate(conds):
        if cond.shape[0] != ts[c]:
            print('Condition {} provided mismatches the corresponding dimension of the tgtArray'.format(c))
            n_err += 1
            break
    if n_err > 0:
        return None
    else:
        if maxi is True:
            value = np.nanmax(tgt)
        else:
            value = np.nanmin(tgt)
        loc = np.where(tgt == value) # a tuple of np.array, the first index for all the hits, the second index, the third...
        # loc: a tuple of np.array, len(loc) = len(conds), (index among the first cond, index among the second, ...)
        try:
            opt_conds = []
            for c, cond in enumerate(conds):
                opt_conds.append(cond[loc[c]].tolist())
            return opt_conds
        except Exception as e:
            print('Condition {} has en exception {}'.format(c, e))


def opt_cond_requested(f_org_list=None, p_org_list=None, p_aqu_list=None,
                       tgt_idx=0, conds=[], max_phase='', min_phase='', tgt_prop=''):
    """
    :return: a list of list of optimal conditions for the opt goals
    For both yield and purity, [yield related, purity related]
    For either yield and purity, only yield or only purity, [max_phase, mix_phase], min_phase results are optional
    For each phase, [condA, condB, ...], each cond is an np.array, [condA1, condB1] for optimum 1
    """

    conds_max_forg = opt_cond(f_org_list, idx=tgt_idx, conds=conds, maxi=True)
    conds_min_forg = opt_cond(f_org_list, idx=tgt_idx, conds=conds, maxi=False)
    conds_max_porg = opt_cond(p_org_list, idx=tgt_idx, conds=conds, maxi=True)
    conds_min_porg = opt_cond(p_org_list, idx=tgt_idx, conds=conds, maxi=False)
    conds_max_paqu = opt_cond(p_aqu_list, idx=tgt_idx, conds=conds, maxi=True)
    conds_min_paqu = opt_cond(p_aqu_list, idx=tgt_idx, conds=conds, maxi=False)

    def opt_cond_yield(conds_max_forg, conds_min_forg, max_phase, min_phase):
        result = []
        if max_phase == 'aqu':
            result.append(conds_min_forg)
            if min_phase == 'aqu':
                result.append(conds_max_forg)
            elif min_phase == 'org':
                result.append(conds_min_forg)
            else:
                pass
        elif max_phase == 'org':
            result.append(conds_max_forg)
            if min_phase == 'aqu':
                result.append(conds_max_forg)
            elif min_phase == 'org':
                result.append(conds_min_forg)
            else:
                pass
        else:
            result.append([])
            if min_phase == 'aqu':
                result.append(conds_max_forg)
            elif min_phase == 'org':
                result.append(conds_min_forg)
            else:
                pass
        return result

    def opt_cond_purity(conds_max_porg, conds_min_porg, conds_max_paqu, conds_min_paqu, max_phase, min_phase):
        result = []
        if max_phase == 'aqu':
            result.append(conds_max_paqu)
            if min_phase == 'aqu':
                result.append(conds_min_paqu)
            elif min_phase == 'org':
                result.append(conds_min_porg)
            else:
                pass
        elif max_phase == 'org':
            result.append(conds_max_porg)
            if min_phase == 'aqu':
                result.append(conds_min_paqu)
            elif min_phase == 'org':
                result.append(conds_min_porg)
            else:
                pass
        else:
            result.append([])
            if min_phase == 'aqu':
                result.append(conds_min_paqu)
            elif min_phase == 'org':
                result.append(conds_min_porg)
            else:
                pass
        return result

    if tgt_prop == 'y':
        return opt_cond_yield(conds_max_forg, conds_min_forg, max_phase, min_phase)
    elif tgt_prop == 'p':
        return opt_cond_purity(conds_max_porg, conds_min_porg, conds_max_paqu, conds_min_paqu, max_phase, min_phase)
    elif tgt_prop == 'b':
        result = []
        result.append(opt_cond_yield(conds_max_forg, conds_min_forg, max_phase, min_phase))
        result.append(opt_cond_purity(conds_max_porg, conds_min_porg, conds_max_paqu, conds_min_paqu, max_phase, min_phase))
        return result
    else:
        err_msg = "Invalid target property {} {}, options are yield 'y', purity 'p', or both 'b'".format(tgt_prop, type(tgt_prop))
        return err_msg


def sep_2D_plot(x, y, x_name, y_name, molecules, semilogx=False, savePic=False):
    """2D plot to compare all the molecules
    :param x: data for x-axis
    :param y: data for y-axis, y[i,j], where i is for a molecule, j is for a data point corresponding to a x value
    :return: if savePic is True, then .png image, otherwise a 'matplotlib.figure.Figure' object
    """
    fig = plt.figure()
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    lineList = []
    labelList = []
    for i in range(len(molecules)):
        if semilogx:
            line, = plt.semilogx(x, y[i], 'o-')
        else:
            line, = plt.plot(x, y[i], 'o-')
        lineList.append(line)
        lineLabel = 'Comp {0}'.format(molecules[i])
        labelList.append(lineLabel)
        lgd = fig.legend(lineList, labelList, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=len(molecules) / 2)
    if savePic:
        file_name = y_name.replace(' ', '_')
        fig.savefig('{}.png'.format(file_name), bbox_extra_artists=(lgd,), bbox_inches='tight')
    return fig


def SPARC_logD(molecule_smi, org, aqu, temp=25.0, pHmin=0.0, pHincrement=0.1):
    """
    Calculate logD of a molecule in a org phase w.r.t. water from pHmin to pH = 14 with a pHincrement as specified
    """

    #  URL for the SPARC local server
    url = 'http://localhost:8080/sparc-integration/rest/calc/logd'

    task = {
        "type": "LOGD",
        "solvent": solventRepr(org),
        "temperature": temp,
        "pH_minimum": pHmin,
        "pH_increment": pHincrement,
        "ionic_strength": aqu.getIonicStrength(),
        "smiles": molecule_smi
    }
    # print(task)

    myResponse = requests.post(url, json=task)
    # print 'Status', myResponse.status_code
    # print(myResponse.content)
    # print myResponse.headers['content-type']

    # For successful API call, response code will be 200 (OK)
    if (myResponse.ok):
        # Loading the response data into a dict variable
        # json.loads takes in only binary or string variables so using content to fetch binary content
        # Loads: takes a Json file and converts into python data structure (dict or list, depending on JSON)
        jData = json.loads(myResponse.content)
        logdValuePair = jData['plotCoordinates'] # (pH, logD)

        return logdValuePair

    else:
        # If response code is not ok (200), print the resulting http error code with description
        myResponse.raise_for_status()


def logD_pH(molecules, org, aqu, temp=25.0, pHmin=0.0, pHincrement=0.1, savePic=False, saveCSV=False):
    """
    :param molecules: a list of SMILES for molecules
    :param savePic, saveCSV: Boolean variable, default False
    :return: logD at various pH for all the molecules & plot pic  & csv file (optional)
    """
    # Initialize np.array for collecting data
    n_molecules = len(molecules)
    n_pHPts = int((14.0 - pHmin)/pHincrement) + 1
    logD_list = np.zeros((n_molecules, n_pHPts))

    fig = plt.figure()
    plt.xlabel('pH')
    yLabel = 'logD ({0}/water)'.format(org.getName())
    plt.ylabel(yLabel)
    lineList = []
    labelList = []

    for i in range(n_molecules):
        valuePair = SPARC_logD(molecules[i], org, aqu, temp, pHmin, pHincrement)

        pHlist = []
        valueList = []
        for j in range(min(len(valuePair), n_pHPts)):
            pHlist.append(valuePair[j][0])
            valueList.append(valuePair[j][1])
            logD_list[i,j] = valuePair[j][1]

        line, = plt.plot(pHlist,valueList,'o-')
        lineList.append(line)
        lineLabel = 'Comp {0}'.format(molecules[i].getName())
        labelList.append(lineLabel)
        lgd = fig.legend(lineList, labelList, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=len(molecules)/2)

    if savePic:
        fig.savefig(os.path.dirname(__file__) + 'logD_{}.png'.format(org.getName()), bbox_extra_artists=(lgd,), bbox_inches='tight')

    if saveCSV:
        csvFile = open('logD_{}.csv'.format(org.getName()), "wb")
        csvOut = csv.writer(csvFile)
        csvOut.writerow(('Molecule', 'Solvent', 'pH', 'logD'))
        for i in range(n_molecules):
            for j in range(n_pHPts):
                csvOut.writerow((molecules[i].getName(), org.getName(), pHlist[j], logD_list[i,j]))
        csvFile.close()

    return pHlist, logD_list


def sep_at_different_pH(org, aqu, molecules=[], c_org=[], c_aqu=[], flowRatio_oa=None, temp=25.0, pHmin=0.0,
                        pHincrement=0.1, logDs_cal=None, saveCSV=False, plotPic=False, savePic=False):
    """
    Given solvents, flow rate ratio, concentrations in both phases, calculate the extraction efficiency
    and purity of each compound at different pH values
    :param molecules: a list of SMILES for molecules
    :param c_org, c_aqu: a list of concentrations of each molecule in two phases
    :param savePic, saveCSV, plotPic: Boolean variable, default False
    """

    if flowRatio_oa is None:
        flowRatio_oa = 1.0
        print('Missing flow ratio, set v_aqu = v_org')

    # Initialize np.array for collecting data
    n_molecules = len(molecules)
    n_pHPts = int((14.0 - pHmin)/pHincrement) + 1
    f_aqu_list = np.zeros((n_molecules, n_pHPts))
    f_org_list = np.zeros((n_molecules, n_pHPts))
    n_aqu_list = np.zeros((n_molecules, n_pHPts))
    n_org_list = np.zeros((n_molecules, n_pHPts))
    c_aqu_list = np.zeros((n_molecules, n_pHPts))
    c_org_list = np.zeros((n_molecules, n_pHPts))
    p_aqu_list = np.zeros((n_molecules, n_pHPts))
    p_org_list = np.zeros((n_molecules, n_pHPts))

    # When logD_list is not calculated
    if logDs_cal is None:
        logD_list = np.zeros((n_molecules, n_pHPts))
        for i in range(len(molecules)):
            valuePair = SPARC_logD(molecules[i], org, aqu, temp, pHmin, pHincrement)
            pHlist = []
            for j in range(min(len(valuePair), n_pHPts)):
                pHlist.append(valuePair[j][0])
                logD_list[i, j] = valuePair[j][1]
            pHlist = np.array(pHlist)
    else: # when logD_list is provided
        logD_list = logDs_cal
        pHlist = np.linspace(pHmin, 14.0, n_pHPts)

    # calculated extraction efficiency at various pH
    v_aqu = 1.0
    v_org = flowRatio_oa * v_aqu
    for i in range(len(molecules)):
        # total moles of molecule i
        try:
            tot_mol = c_aqu[i] * v_aqu + c_org[i] * v_org
        except Exception as e:
            tot_mol = 0
            print(e)
            print('Missing initial concentration? set tot_mol = 0 for mol {}: {}'.format(i, molecules[i]))

        for j in range(n_pHPts):
            f_aqu = v_aqu / (v_aqu + v_org * (10 ** logD_list[i,j]))
            f_org = 1 - f_aqu
            f_aqu_list[i,j] = f_aqu
            f_org_list[i,j] = f_org
            n_aqu_list[i,j] = tot_mol * f_aqu
            n_org_list[i,j] = tot_mol * f_org
            c_aqu_list[i,j] = tot_mol * f_aqu / v_aqu
            c_org_list[i,j] = tot_mol * f_org / v_org

    # Calculate purity of compound i in each phase at different pH
    for j in range(n_pHPts):
        for i in range(n_molecules):
            p_aqu_list[i,j] = n_aqu_list[i,j]/np.sum(n_aqu_list[:,j])
            p_org_list[i,j] = n_org_list[i,j]/np.sum(n_org_list[:,j])

    # Save csv file
    if saveCSV:
        csvFile = open("sep_{}.csv".format(org.getName()), "wb")
        csvOut = csv.writer(csvFile)
        csvOut.writerow(('Molecule', 'Solvent', 'pH', 'logD', 'f_org', 'f_aqu', 'p_org', 'p_aqu','c_org','c_aqu'))
        for i in range(n_molecules):
            for j in range(n_pHPts):
                csvOut.writerow((molecules[i].getName(), org.getName(), pHlist[j], logD_list[i, j], f_org_list[i, j],
                                 f_aqu_list[i, j], p_org_list[i, j], p_aqu_list[i, j], c_org_list[i,j], c_aqu_list[i,j]))
        csvFile.close()

    if plotPic:
        # logD w.r.t. pH
        y_name = 'logD {}-water'.format(org.getName())
        fd = sep_2D_plot(pHlist, logD_list, 'pH', y_name, molecules, semilogx=False, savePic=savePic)

        # Purity w.r.t. pH
        y_name = 'Purity in org phase {}'.format(org.getName())
        po = sep_2D_plot(pHlist, p_org_list, 'pH', y_name, molecules, semilogx=False, savePic=savePic)

        y_name = 'Purity in aqu phase'
        pa = sep_2D_plot(pHlist, p_aqu_list, 'pH', y_name, molecules, semilogx=False,savePic=savePic)

        # Extraction efficiency w.r.t. pH
        y_name = 'Extraction efficiency in org phase {}'.format(org.getName())
        fo = sep_2D_plot(pHlist, f_org_list, 'pH', y_name, molecules, semilogx=False, savePic=savePic)

        y_name = 'Extraction efficiency in aqu phase'
        fa = sep_2D_plot(pHlist, f_aqu_list, 'pH', y_name, molecules, semilogx=False, savePic=savePic)

        return fd, fo, fa, po, pa, pHlist, logD_list, f_org_list, c_org_list, c_aqu_list, p_org_list, p_aqu_list
    else:
        return pHlist, logD_list, f_org_list, f_aqu_list, c_org_list, c_aqu_list, p_org_list, p_aqu_list


def sep_pH_flowratio(org, aqu, molecules=[], c_orgs=[], c_aqus=[], temp=25.0, pHincrement=0.1,
                    min_fr=0.1, max_fr=10, num_fr=10, plotPlot=True, showPlot=False, savePlot=True):
    """Separation results for a set of molecules with a org/aqu flowratio_in (min_fr, max_fr) in pH (0, 14)
    """
    # flow rate ratio (org/aqu) in (min_fr, max_fr)

    flow_oa = np.linspace(np.log10(min_fr), np.log10(max_fr), num_fr)
    flow_oa_ratio = np.zeros(num_fr)
    for i in range(num_fr):
        flow_oa_ratio[i] = 10 ** flow_oa[i]

    # Initialize np.array for collecting data
    pHmin = 0.0
    n_molecules = len(molecules)
    n_pHPts = int((14.0 - pHmin)/pHincrement) + 1
    f_aqu_list = np.zeros((n_molecules, num_fr, n_pHPts))
    f_org_list = np.zeros((n_molecules, num_fr, n_pHPts))
    c_aqu_list = np.zeros((n_molecules, num_fr, n_pHPts))
    c_org_list = np.zeros((n_molecules, num_fr, n_pHPts))
    p_aqu_list = np.zeros((n_molecules, num_fr, n_pHPts))
    p_org_list = np.zeros((n_molecules, num_fr, n_pHPts))

    logDs_cal = None
    for i in range(num_fr):
        r = flow_oa_ratio[i]
        pHlist, logD_list, f_org, f_aqu, c_org, c_aqu, p_org, p_aqu = \
            sep_at_different_pH(org, aqu, molecules, c_orgs, c_aqus, flowRatio_oa=r, temp=temp, pHmin=pHmin,
                                pHincrement=pHincrement, logDs_cal=logDs_cal, saveCSV=False,
                                plotPic=False, savePic=False)

        logDs_cal = logD_list
        f_aqu_list[:, i, :] = f_aqu
        f_org_list[:, i, :] = f_org
        c_aqu_list[:, i, :] = c_aqu
        c_org_list[:, i, :] = c_org
        p_aqu_list[:, i, :] = p_aqu
        p_org_list[:, i, :] = p_org

    # Plot results
    if plotPlot:
        n_col = 2
        n_row = int(round(n_molecules/float(n_col)))

        def surface_plt(data, dataName='', showP=False, saveP=False):
            fig = plt.figure(figsize=plt.figaspect(n_row/float(n_col)))
            rcParams.update({'axes.labelsize': 'small'})
            rcParams.update({'xtick.labelsize': 'small'})
            rcParams.update({'ytick.labelsize': 'small'})
            x = pHlist
            y = flow_oa
            x, y = np.meshgrid(x, y)
            for m in range(n_molecules):
                ax = fig.add_subplot(n_row, n_col, m+1, projection='3d')
                z = data[m]
                # surf = ax.plot_surface(x, y, z, rstride=pHincrement, cstride=fr_step, cmap=cm.coolwarm,
                #                        label='{}'.format(molecules[m]))
                surf = ax.plot_surface(x, y, z, cmap=cm.jet, linewidth=0)
                ax.set_xlabel('pH')
                ax.set_ylabel('log10(flowRatio)')
                ax.set_zlabel('{}'.format(dataName))
                ax.set_zlim(-0.01, 1.01)
                # ax.set_title('{}'.format(molecules[m]))
                ax.set_title('mol {}'.format(m+1))
                ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
                fig.colorbar(surf, shrink=0.5)
            if showP:
                plt.show()
            if saveP:
                fig.savefig('{}'.format(dataName), bbox_inches='tight')
            return fig

        fa = surface_plt(f_aqu_list, 'f_aqu', showP=showPlot, saveP=savePlot)
        fo = surface_plt(f_org_list, 'f_org', showP=showPlot, saveP=savePlot)
        pa = surface_plt(p_aqu_list, 'p_aqu', showP=showPlot, saveP=savePlot)
        po = surface_plt(p_org_list, 'p_org', showP=showPlot, saveP=savePlot)

        return fo, fa, po, pa, pHlist, logD_list, flow_oa_ratio, f_org_list, f_aqu_list,\
               c_org_list, c_aqu_list, p_org_list, p_aqu_list
    else:
        return pHlist, logD_list, flow_oa_ratio, f_org_list, f_aqu_list, c_org_list, c_aqu_list, p_org_list, p_aqu_list

def sep_at_different_flowratio(org, aqu, molecules=[], c_orgs=[], c_aqus=[], temp=25.0,
                               min_fr=0.1, max_fr=10, num_fr=10, plotPic=True, savePlot=False):
    pH_list, logD_list, flow_oa_ratio, f_org_list, f_aqu_list, c_org_list, c_aqu_list, p_org_list, p_aqu_list = \
        sep_pH_flowratio(org, aqu, molecules, c_orgs, c_aqus, temp, pHincrement=0.1,
                         min_fr=min_fr, max_fr=max_fr, num_fr=num_fr, plotPlot=False, showPlot=False, savePlot=False)
    # (n_molecules, num_fr, n_pHPts)
    pH = aqu.getpH()
    pHlist = list(pH_list)
    if pH in pHlist:
        idx = pHlist.index(pH)
        logD = logD_list[:, idx]
        f_org = f_org_list[:, :, idx]
        f_aqu = f_aqu_list[:, :, idx]
        # c_org = c_org_list[:, :, idx]
        # c_aqu = c_aqu_list[:, :, idx]
        p_org = p_org_list[:, :, idx]
        p_aqu = p_aqu_list[:, :, idx]
    else:
        pHround = round(pH, 1)
        idx = pHlist.index(pHround)
        # when pH is smaller than the closest pH in the list
        if pHround > pH:
            highIndex = idx
            lowIndex = idx - 1
        # when pH is larger than the closest pH in the list
        else:
            highIndex = idx + 1
            lowIndex = idx

        def find_in_list(vArray, hi=highIndex, li=lowIndex, m=pH, mList=pHlist):
            vh = vArray[:, :,  hi]
            vl = vArray[:, :, li]
            return vh - (vh - vl) / (mList[hi] - mList[li]) * (mList[hi] - m)

        logD = find_in_list(logD_list)
        f_org = find_in_list(f_org_list)
        f_aqu = find_in_list(f_aqu_list)
        # c_org = find_in_list(c_org_list)
        # c_aqu = find_in_list(c_aqu_list)
        p_org = find_in_list(p_org_list)
        p_aqu = find_in_list(p_aqu_list)

    if plotPic:
        # # logD w.r.t. flow rate ratio
        # y_name = 'logD {}-water'.format(org.getName())
        # fd = sep_2D_plot(flow_oa_ratio, logD, 'flow rate ratio', y_name, molecules, semilogx=True, savePic=savePlot)

        # Purity w.r.t. pH
        y_name = 'Purity in org phase {}'.format(org.getName())
        po = sep_2D_plot(flow_oa_ratio, p_org, 'flow rate ratio', y_name, molecules, semilogx=True, savePic=savePlot)

        y_name = 'Purity in aqu phase'
        pa = sep_2D_plot(flow_oa_ratio, p_aqu, 'flow rate ratio', y_name, molecules, semilogx=True, savePic=savePlot)

        # Extraction efficiency w.r.t. pH
        y_name = 'Extraction efficiency in org phase {}'.format(org.getName())
        fo = sep_2D_plot(flow_oa_ratio, f_org, 'flow rate ratio', y_name, molecules, semilogx=True, savePic=savePlot)

        y_name = 'Extraction efficiency in aqu phase'
        fa = sep_2D_plot(flow_oa_ratio, f_aqu, 'flow rate ratio', y_name, molecules, semilogx=True, savePic=savePlot)

        return fo, fa, po, pa, flow_oa_ratio, logD, f_org, f_aqu, p_org, p_aqu
    else:
        return flow_oa_ratio, logD, f_org, f_aqu, p_org, p_aqu


class SeparationDesigner():
    """
    Flow Separation Design using SPARC
    """
    def __init__(self, logDcal = SPARC_logD):

        # logD calculator used
        self.cal = logDcal

        # Separation target
        self.tgtMol = '' # Target molecules
        self.tgtProp = '' # Target property, e.g. high yield 'y', high purity 'p' or both 'b'
        self.maxPhase = '' # In which phase the target molecule should be concentrated
        self.minPhase = '' # In which phase the target molecule should be minimized

        # A list of other molecules in the mixture
        self.molecules = []

        # Separation conditions
        self.aqu = AquPhase() # AquPhase
        self.org = OrgPhase() # OrgPhase
        self.flowRatio = 1.0 # OrgPhase/AquPhase flow rate ratio
        self.temp = 25.0

        # Initial conditions
        self.aquInitConc = [] # initial concentration of each molecules in the aqueous phase
        self.orgInitConc = [] # initial concentration of each molecules in the organic phase

        # Optimization
        self.optVar = ''
        self.pHStep = 0.1
        self.minFlowRatio = 0.1
        self.maxFlowRatio = 10
        self.num_flowratio = 10

        # Score
        self.scoreFunc = None

    def load_designer(self, userInput):
        """load the designer with necessary information from user input"""
        self.tgtMol = str(userInput['tgt_smiles'])
        other_mol = str(userInput['other_smiles']).split(',')
        self.molecules = [mol for mol in other_mol]
        self.tgtProp = str(userInput['opt_tgt_choice'])
        self.optVar = str(userInput['opt_option'])
        self.maxPhase = str(userInput['max_phase'])
        if userInput['min_phase']: self.minPhase = str(userInput['min_phase'])

        # Organic phase defined by the user
        slvt1 = OrgSolvent(userInput['org1'])
        if userInput['org2']:
            slvt2 = OrgSolvent(str(userInput['org2']).split(',')[0])
            self.org = OrgPhase(solvent1=slvt1, solvent2=slvt2, x2=float(str(userInput['org2']).split(',')[1]))
        else:
            self.org = OrgPhase(solvent1=slvt1)

        # if pH or ionic strength is given by the user
        if userInput['pH'] and userInput['ionic_strength']:
            self.aqu = AquPhase(pH=userInput['pH'], ionicStrength=userInput['ionic_strength'])
        elif userInput['pH']:
            self.aqu = AquPhase(pH=userInput['pH'])
        elif userInput['ionic_strength']:
            self.aqu = AquPhase(ionicStrength=userInput['ionic_strength'])
        else:
            self.aqu = AquPhase()

        if userInput['flow_ratio']: self.flowRatio = userInput['flow_ratio'] # if user provides flow ratio
        if userInput['temp']: self.temp = userInput['temp']
        if userInput['pH_step']: self.pHStep = userInput['pH_step']
        if userInput['flow_ratio_min']: self.minFlowRatio = userInput['flow_ratio_min']
        if userInput['flow_ratio_max']: self.maxFlowRatio = userInput['flow_ratio_max']
        if userInput['num_flowratio']: self.num_flowratio = userInput['num_flowratio']
        aquC = userInput['aqu_conc'].split(',')
        self.aquInitConc = [float(x) for x in aquC]
        orgC = userInput['org_conc'].split(',')
        self.orgInitConc = [float(x) for x in orgC]
        # if userInput['score_func']: self.scoreFunc = userInput['score_func']

    #### To do: def org_solvent_selection(), about how to go through the list of solvents available

    def sep_performance(self, at_pH=False):
        """
        Separation performance with all the separation conditions available
        """
        if not self.tgtMol:
            print('No target molecules input!')
            pass
        if not self.molecules:
            print('No other molecules present, no separation needed!')
            pass
        if not (self.org.getOne() and self.aquInitConc and self.orgInitConc):
            print('The separation conditions are not complete')
            if not self.org.getOne(): print('Missing organic phase')
            if not self.aquInitConc: print('Missing initial concentration in aqueous phase')
            if not self.orgInitConc: print('Missing initial concentration in organic phase')
            pass
        else:
            mol_list = [self.tgtMol]
            for mol in self.molecules:
                mol_list.append(mol)
            fd, fo, fa, po, pa, pHlist, logD_list, f_org_list, c_org_list, c_aqu_list, p_org_list, p_aqu_list = \
                sep_at_different_pH(org=self.org, aqu=self.aqu, molecules=mol_list, c_org=self.orgInitConc,
                                    c_aqu=self.aquInitConc, flowRatio_oa=self.flowRatio, temp=25.0, pHmin=0.0,
                                    pHincrement=0.1, logDs_cal=None, saveCSV=False, plotPic=True, savePic=False)
            if at_pH:
                pH = self.aqu.getpH()
                if pH in pHlist:
                    idx = pHlist.index(pH)
                    logD = logD_list[:, idx]
                    f_org = f_org_list[:, idx]
                    c_org = c_org_list[:, idx]
                    c_aqu = c_aqu_list[:, idx]
                    p_org = p_org_list[:, idx]
                    p_aqu = p_aqu_list[:, idx]
                else:
                    pHround = round(pH, 1)
                    idx = pHlist.index(pHround)
                    # when pH is smaller than the closest pH in the list
                    if pHround > pH:
                        highIndex = idx
                        lowIndex = idx - 1
                    # when pH is larger than the closest pH in the list
                    else:
                        highIndex = idx + 1
                        lowIndex = idx

                    def find_in_list(vArray, hi=highIndex, li=lowIndex, m=pH, mList=pHlist):
                        vh = vArray[:, hi]
                        vl = vArray[:, li]
                        return vh - (vh - vl) / (mList[hi] - mList[li]) * (mList[hi] - m)

                    logD = find_in_list(logD_list)
                    f_org = find_in_list(f_org_list)
                    c_org = find_in_list(c_org_list)
                    c_aqu = find_in_list(c_aqu_list)
                    p_org = find_in_list(p_org_list)
                    p_aqu = find_in_list(p_aqu_list)
                print('At pH {}'.format(pH))
                print('molecule\t' + 'logD\t' + 'f_org\t' + 'f_aqu\t' + 'c_org\t' + 'c_aqu\t' + 'p_org\t' + 'p_aqu')
                for m, mol in enumerate(mol_list):
                    print('{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}'.format(
                        mol, logD[m], f_org[m], 1-f_org[m], c_org[m], c_aqu[m], p_org[m], p_aqu[m]))
            else:
                return fd, fo, fa, po, pa, pHlist, logD_list, f_org_list, c_org_list, c_aqu_list, p_org_list, p_aqu_list
                # Figure objects: logD, fraction_org, fraction_aqu, purity_org, purity_aqu


    def opt_sep_pH_flow(self):
        """pH and flow rate ratio optimization given the solvents
        :return: a list of list of optimal conditions for the opt goals, and figures
        """
        # global f_org_list, p_org_list, p_aqu_list, conds, fo, fa, po, pa
        if not self.tgtMol:
            print('No target molecules input!')
            pass
        if not self.molecules:
            print('No other molecules present, no separation needed!')
            pass
        if not self.org.getOne():
            print('No organic phase specified!')
            pass

        mol_list = [self.tgtMol]
        for mol in self.molecules:
            mol_list.append(mol)

        if self.optVar == 'b':
            fo, fa, po, pa, pHlist, logD_list, flow_oa_ratio, \
            f_org_list, f_aqu_list, c_org_list, c_aqu_list, p_org_list, p_aqu_list = \
                sep_pH_flowratio(org=self.org, aqu=self.aqu, molecules=mol_list, c_orgs=self.orgInitConc,
                                 c_aqus=self.aquInitConc, temp=self.temp, pHincrement=self.pHStep,
                                 min_fr=self.minFlowRatio, max_fr=self.maxFlowRatio, num_fr=self.num_flowratio,
                                 plotPlot=True, showPlot=False, savePlot=False)
            conds = [flow_oa_ratio, pHlist]

        elif self.optVar == 'h':
            fd, fo, fa, po, pa, pHlist, logD_list,f_org_list, c_org_list, c_aqu_list, p_org_list, p_aqu_list = \
                sep_at_different_pH(org=self.org, aqu=self.aqu, molecules=mol_list, c_org=self.orgInitConc,
                                    c_aqu=self.aquInitConc, flowRatio_oa=self.flowRatio, temp=self.temp, pHmin=0.0,
                                    pHincrement=self.pHStep, logDs_cal=None, saveCSV=False, plotPic=True, savePic=False)
            conds = [pHlist]

        else: # self.optVar is 'f'
            fo, fa, po, pa, flow_oa_ratio, logD_list, f_org_list, f_aqu_list, p_org_list, p_aqu_list = \
                sep_at_different_flowratio(org=self.org, aqu=self.aqu, molecules=mol_list, c_orgs=self.orgInitConc,
                                           c_aqus=self.aquInitConc, temp=self.temp, min_fr=self.minFlowRatio,
                                           max_fr=self.maxFlowRatio, num_fr=self.num_flowratio,
                                           plotPic=True, savePlot=False)
            conds = [flow_oa_ratio]

        if self.tgtProp == 'y':
            result = opt_cond_requested(f_org_list=f_org_list, p_org_list=p_org_list, p_aqu_list=p_aqu_list,
                                        tgt_idx=0, conds=conds, max_phase=self.maxPhase, min_phase=self.minPhase,
                                        tgt_prop=self.tgtProp)
            return result, fo, fa
        elif self.tgtProp == 'p':
            if self.aquInitConc and self.orgInitConc:
                result = opt_cond_requested(f_org_list=f_org_list, p_org_list=p_org_list, p_aqu_list=p_aqu_list,
                                            tgt_idx=0, conds=conds, max_phase=self.maxPhase, min_phase=self.minPhase,
                                            tgt_prop=self.tgtProp)
                return result, po, pa
            else:
                err_msg = 'Please provide initial concentration for all the species'
                return err_msg
        else:
            if self.aquInitConc and self.orgInitConc:
                result = opt_cond_requested(f_org_list=f_org_list, p_org_list=p_org_list, p_aqu_list=p_aqu_list,
                                            tgt_idx=0, conds=conds, max_phase=self.maxPhase, min_phase=self.minPhase,
                                            tgt_prop=self.tgtProp)
                return result, fo, fa, po, pa
            else:
                result = opt_cond_requested(f_org_list=f_org_list, p_org_list=p_org_list, p_aqu_list=p_aqu_list,
                                            tgt_idx=0, conds=conds, max_phase=self.maxPhase, min_phase=self.minPhase,
                                            tgt_prop='y')
                print('Only fraction yield optimized due to missing initial concentrations')
                return result, fo, fa

    def result_to_html_table(self, result):
        """
        Convert numpy.array results into a html table
        Not really working now in the .html part
        """
        html = '<table>'
        if self.tgtProp == 'y' or self.tgtProp == 'p':
            for m in result: # maximization and minimization
                for c in m: # flow rate ratio and pH
                    html += '<tr>'
                    for n in c: # actual conditions
                        html += "<td>{}</td>".format(n)
                    html += '</tr>'
        else:
            for p in result: # yield and purity
                for m in p:
                    for c in m:  # flow rate ratio and pH
                        html += '<tr>'
                        for n in c:  # actual conditions
                            html += "<td>{}</td>".format(n)
                        html += '</tr>'
        html += '</table>'
        return html