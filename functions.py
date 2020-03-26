import os, shutil            #Operating System
#from pyswmm import Simulation, raingages, Nodes, SystemStats, Links
import pandas as pd
from swmmio.utils.modify_model import replace_inp_section
from swmmio import create_dataframeINP
import swmmio
from swmmio.graphics import swmm_graphics as sg
import math
import numpy as np
import datetime as DT
from functions import *


def prop2plot(dataFrame, Property, Ranges = [.2,.4,.6,.8,1.0] , FS=1, 
              ColorPalette='red_gradient', create_legend = False):
    # Funzione per plottare l'analisi del modello
    # si tratta di una funzione che viene richamata da plotModel
    # per versioni di swmmio successive alla 3.6 il parametro create_legend potrebbe dare errori, tenerlo False
    # ColorPalette possibili sono 'red_gradient', 'grey_scale', 'rainbow', 'blue_scale', oppure una palette personalizzata come una lista di 5 elementi
    # FS è il fattore di scala
    
    df = dataFrame
    pr = dataFrame[Property]
    delta = max(pr) - min(pr)
    
    a = min(pr) + delta * Ranges[0]
    b = min(pr) + delta * Ranges[1]
    c = min(pr) + delta * Ranges[2]
    d = min(pr) + delta * Ranges[3]
    e = min(pr) + delta * Ranges[4]
    
    if ColorPalette == 'red_gradient':
        df['draw_color'] = pr.apply(sg.gradient_grey_red, xmin = min(pr), xmax = max(pr)) # gradient grey-> red
    else:
        if ColorPalette == 'grey_scale':
            ColorPalette = ['#C8C8C8','#A0A0A0','#787878','#505050', '#282828','#000'] # grey scale
        elif ColorPalette == 'rainbow':
            ColorPalette = ['#C8C8C8','#00F','#0FF','#0F0', '#FF0', '#F00'] # rainbow
        elif ColorPalette == 'blue_scale':
            ColorPalette = ['#C8C8FF','#A0A0FF','#7878FF','#5050FF', '#2828F','#00F'] # rainbow
        df['draw_color'] = ColorPalette[0] #grey default
        df.loc[pr >= a, 'draw_color'] = ColorPalette[1] #blue
        df.loc[pr >= b, 'draw_color'] = ColorPalette[2] #cyano
        df.loc[pr >= c, 'draw_color'] = ColorPalette[3] #verde
        df.loc[pr >= d, 'draw_color'] = ColorPalette[4] #giiallo
        df.loc[pr >= d, 'draw_color'] = ColorPalette[5] #red
    
    #df['draw_size'] = pr * FS
    df['draw_size'] = 1 * FS
    df.loc[pr >= a, 'draw_size'] = 2 * FS
    df.loc[pr >= b, 'draw_size'] = 3 * FS
    df.loc[pr >= c, 'draw_size'] = 4 * FS
    df.loc[pr >= d, 'draw_size'] = 5 * FS
    df.loc[pr >= e, 'draw_size'] = 5 * FS
    
    if create_legend: 
        txt_legend = [min(pr),a,b,c,d,e]
        return(txt_legend)

## Plot dell'analisi della rete
def plotModel(inpfilepath, title, annotation, 
              Property_nodes, Property_conds, 
              outfile, ColorPalette='red_gradient', 
              legend=False):
    
    # plotModel è una funzione per la stampa del grafico delle proprietà di una rete
    # inpfilepath: str file .inp del modello
    # title: str del titolo
    # annotation: str delle note in fondo a dx
    # Property_nodes: str proprietà dei nodi da plottare
    # Property_conds: str proprietà delle condotte da plottare, se non si vuole usare ''
    # outfile: str file di output (.png)
    # ColorPalette: possibili sono 'red_gradient', 'grey_scale', 'rainbow', 'blue_scale', oppure una palette personalizzata come una lista di 5 elementi ['#C8C8C8','#A0A0A0','#787878','#505050', '#282828','#000']
    # per versioni di swmmio successive alla 3.6 il parametro create_legend potrebbe dare errori, tenerlo False

    if os.path.isfile(outfile):
        print('Grafico non plottato: ' + str(outfile) + ' esiste già')
        return

    mymodel = swmmio.Model(inpfilepath)
    nodes = mymodel.nodes()
    conds = mymodel.conduits()
    #mymodel = swmmio.Model(inpfilepath)
    
    inpfilebasename = inpfilepath[:-4] # Nome file senza estensione
    
    if legend: create_legend = True
    else:      create_legend = False
    
    if Property_nodes != '': 
        txt_legend_nodes = prop2plot(dataFrame=nodes, 
                                     Property=Property_nodes, 
                                     FS = 5, 
                                     ColorPalette=ColorPalette,
                                     create_legend = create_legend)
        #minNode = round(min(nodes[Property_nodes]), 1)
        #maxNode = round(max(nodes[Property_nodes]), 1)
        #print(nodes[Property_nodes])
    else:
        nodes['draw_size'] = 1
        #minNode = ''
        #maxNode = ''
        
    if Property_conds != '':
        txt_legend_conds = prop2plot(dataFrame=conds, 
                                     Property=Property_conds, 
                                     FS=5, 
                                     ColorPalette=ColorPalette,
                                     create_legend = create_legend)
        #minCond = round(min(conds[Property_conds]), 1)
        #maxCond = round(max(conds[Property_conds]), 1)
        #print(conds[Property_conds])
    else:
        conds['draw_size'] = 1
        #minCond = ''
        #maxCond = ''

    if legend: 
        #print(txt_legend_conds)
        annotation = 'Conds: grey -> ' +   str(round(txt_legend_conds[0],2)) + '\n' + \
                     'Conds: blue -> ' +   str(round(txt_legend_conds[1],2)) + '\n' + \
                     'Conds: cyan -> ' +   str(round(txt_legend_conds[2],2)) + '\n' + \
                     'Conds: green -> ' +  str(round(txt_legend_conds[3],2)) + '\n' + \
                     'Conds: yellow -> ' + str(round(txt_legend_conds[4],2)) + '\n' + \
                     'Conds: red -> '  +   str(round(txt_legend_conds[5],2))
    
    #annotation = 'Percentuale di riempimento\nnei nodi e nelle condotte'
    
    drw = sg.draw_model(nodes=nodes,
                        conduits=conds,
                        background_color = '#FFFFFF00', ## trasparente bianco swmmio modificato
                        title = title,
                        annotation=annotation, 
                        #legend = legend,
                        file_path = outfile)


def plotGeomNetwork(inpfile):
    # Plot della Geometria della rete
    inpfilebasename = inpfile[:-4]
    outfile = inpfilebasename + '_geom.png'
    plotModel(inpfilepath = inpfile, 
              title = 'Diametri della Rete\n' + inpfilebasename, 
              annotation = '',
              outfile = outfile,
              Property_nodes = '',
              Property_conds = 'Geom1',
              #ColorPalette='blue_scale',
              ColorPalette='rainbow',
              legend = True,
             )
    print(outfile)


def projectNet (inpfile, newinpfile):
    # metodo bruto per il dimensionamento
    # Prendo in analisi lo scenario meno svantaggioso,
    # in questo caso tp 5 giorni per un TR di 10 anni
    # controllo quali condotte hanno un riempimento elevato
    # queste le vado ad ingrandiere.
    # poi via via verso gli scrosci 

    # esegue verifica della rete
    swmm5solver( inpfile )

    # crea i df delle condotte
    model   = swmmio.Model(inpfile)
    #nodes   = model.nodes()
    conds   = model.conduits()
    
    # Diminuisco le dimensioni delle condutture troppo vuote
    #print(conds.loc[conds['MaxDPerc'] > 0.8, 'Geom1'] )
    CondLimitMin = 0.15
    Dmin = 0.10
    print('numero Condotte quasi vuote: ' + str(len(conds.loc[conds['MaxDPerc'] < CondLimitMin, 'Geom1'])) )
    conds.loc[conds['MaxDPerc'] < CondLimitMin, 'Geom1'] = conds.loc[conds['MaxDPerc'] < CondLimitMin, 'Geom1']  - 0.1 
    conds.loc[(conds['MaxDPerc'] < CondLimitMin)&(conds['Geom1'] < Dmin), 'Geom1'] = Dmin # imposto un diametro minimo
    #print(conds.loc[conds['MaxDPerc'] < CondLimitMin, 'Geom1'] )
    # se non ci sono cose  blocco tutto
    
    # Aumento le dimensioni delle condutture troppo piene
    CondLimitMax = .99
    Dmax = 1.8
    print('numero Condotte quasi piene: ' + str(len(conds.loc[conds['MaxDPerc'] > CondLimitMax, 'Geom1'])) )
    conds.loc[conds['MaxDPerc'] > CondLimitMax, 'Geom1'] = conds.loc[conds['MaxDPerc'] > CondLimitMax, 'Geom1']  + 0.1
    conds.loc[(conds['MaxDPerc'] > CondLimitMin)&(conds['Geom1'] > Dmax ), 'Geom1'] = Dmax # imposto un diametro massimo

    #print(all(conds['MaxDPerc'] <= 0.8))
    #print(conds.loc[conds['MaxDPerc'] > 0.8, 'Geom1'] )
    # se non ci sono cose  blocco tutto
    if all(conds['MaxDPerc'] <= CondLimitMax):
        print()
        print('Nuova Rete ' + inpfile +' corretamente dimensionata')
        return

    # creo un df da pre sostituire XSECTIONS nell'.inp file
    df = pd.DataFrame()
    df['Shape'] = conds.Shape
    df['Geom1'] = conds.Geom1
    df['Geom2'] = conds.Geom2
    df['Geom3'] = conds.Geom3
    df['Geom4'] = conds.Geom4
    df['Barrels'] = conds.Barrels

    # Aggiorno il file di progetto con le nuove dimensioni
    shutil.copyfile(inpfile, newinpfile)
    replace_inp_section(newinpfile, '[XSECTIONS]', df)

def theta(G):
    # Restituisce l'angolo dell'arco di circonferenza bagnato dato il grado di riempimento G
    return(2 * math.acos(1-2*G))

def ceil2list(x, lst):
    # Approssima al numero superiore nella lista, oppure prende il massimo della lista
    i = 0
    diff = lst[i] - x 
    while (diff < 0 and i < len(lst)-1):
        i += 1
        diff = lst[i] - x
    return(lst[i])

def round2list(x, lst):
    # Approssima al numero più vicino nella lista, oppure prende il massimo della lista
    i = 0
    diff = lst[i] - x 
    while (diff < 0 and i < len(lst)-1):
        i += 1
        diff = lst[i] - x
    N = x - lst[i-1]
    D = lst[i] - lst[i-1]
    if N/D < 0.5:
        i -= 1
    return(lst[i])

def LinkTimeFrame(node, timeend):
    # crea un dataframe 
    df = pd.DataFrame(columns=['Date', 'Time', 'Flow', 
                               'Velocity', 'Depth', 'Capacity'])
    for time in range(0,timeend):
        delta = DT.timedelta(hours=0, minutes = time, seconds=0)
        total_seconds = delta.total_seconds()
        hours = int(total_seconds // 3600)
        minutes = int(total_seconds // 60) - hours * 60
        seconds = int(total_seconds % 60)

        str_time = '{hours:02.0f}:{minutes:02.0f}:{seconds:02.0f}'.format(hours = hours,
                                                   minutes = minutes, 
                                                   seconds = seconds)

        entry = rpt.returnDataAtDTime(node, dtime= str_time)
        """
        ------------------------------------------------------------------
                                     Flow  Velocity     Depth  Capacity/
          Date        Time            LPS     m/sec    meters   Setting
        ------------------------------------------------------------------
        """
        df.loc[str_time]=(entry)

    return(df)

def PropCond(Conduit, Nodes, Conds, Dividers):
    # Ritorna un dataframe con le proprietà di una condotta
    # Conduit: str, indice della condotta
    # Nodes: df dei nodi
    # Conds: df delle condotte
    # Dividers: df dei dividers
    Cond = Conds[Conds.index == Conduit]
    
    #print(Cond)
    
    Length = Cond.Length[0]
    
    
    InNode = Cond.InletNode[0]
    try:
        InElevation = Nodes.InvertElev[Nodes.index == InNode ][0]
    except IndexError:
        InElevation = float(Dividers.Elevation[InNode])
        
    OutNode = Cond.OutletNode[0]
    try:
        OutElevation = Nodes.InvertElev[Nodes.index == OutNode ][0]
    except IndexError:
        OutElevation = float(Dividers.Elevation[OutNode])
    
    Pendenza = (InElevation - OutElevation) / Length
    
    Diametro = Cond.Geom1[0]
    
    #Rh05 = Diametro / 4 
    
    Ks = Cond.ManningN[0] ** -1
    
    #Velocity05 = Ks * Rh05 ** (2/3) * Pendenza ** .5
    
    return({'InNode'      : InNode,
            'OutNode'     : OutNode,
            'Length'      : Length,
            'InElevation' : InElevation,
            'OutElevation': OutElevation,
            'Pendenza'    : Pendenza,
            'Diametro'    : Diametro,
            'Ks'          : Ks,
      #      'Rh05'        : Rh05,
      #      'Velocity05'  : Velocity05,
           })

def D_opt(Q,Ks,i,Teta):
    # funzione di dimensionamento ottimale di una condotta
    b = 10 ** (-9/8) # coefficente per la conversione delle unità di misura
    Dopt = b * ( 2 ** (13/3) * Qmax / Ks / (i ** .5) *   \
            (1 - math.sin(Teta) / Teta) ** (-2/3) * \
            (Teta - math.sin(Teta)) ** -1 ) ** (3/8)
    return(Dopt)
    
    
def Cp(d,H):
    # Funzione di Costo Cp
    # d e H in metri
    # d diametro condotte [m]
    # H profondità media condotta [m]
    # ritorna in $ / m
    
    d = d / 0.3048
    H = H / 0.3048
    
    if d > 3.:
        CpUS = 30 *d + 4.9 * H - 105.9 
    else: 
        if H < 10:
            CpUS = 10.98 * d + 0.8 * H - 5.98
        else:
            CpUS = 5.94 * d + 1.166 * H + 0.504 *H* d- 9.64
    
    # $/ft = $/.3048m
    
    Cp = CpUS * 0.3048 # $/m
    return(Cp)

def Cm(h):
    # Funzione di costo Cm
    # h profondità pozzetti in m
    
    h = h / .3048
    Cm = 250 + h ** 2
    
    return(Cm)

def rStar(n):
    # calcolo di r con il metodo iterativo di Newton-Rapshon
    r = 1
    rr = 0.6
    toll = 10**-6
    while abs(r - rr) > toll:
        r = rr
        rr = r - ((1-n)*(np.exp(r) -1)-r)/((1-n)*np.exp(r)-1)
    return r
        
    

