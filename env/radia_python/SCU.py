
#############################################################################
# RADIA Python Test Example: Parallel (MPI based) Computation of Magnetic Field
# created by SuperConducting Undulator (APS-type)
# v 0.01
#############################################################################

from __future__ import absolute_import, division, print_function #Py 2.*/3.* compatibility
import radia as rad
from radia_ids import *
from uti_plot import *
from math import *
import os
import sys
import time

#*********************************Coils
def Coil():

    #Main ~racetrack part of the Coil
    zcCoil = 0
    strPart1 = rad.ObjRecCur([0.5*wNotch+0.25*(wPole-wNotch),0.5*gap+gapOffset+1.5*ttCoil+hYoke,zcCoil], [0.5*(wPole-wNotch),ttCoil,tlCoil], [jCoil,0,0])
    rMin = 0.5*hYoke
    arcPart1 = rad.ObjArcCur([0.5*wPole,0.5*gap+gapOffset+ttCoil+rMin,zcCoil], [rMin,rMin+ttCoil], [1.5*pi,2.5*pi], tlCoil, nsCoil, -jCoil)
    strPart2 = rad.ObjRecCur([0,0.5*gap+gapOffset+0.5*ttCoil,zcCoil], [wPole,ttCoil,tlCoil], [-jCoil,0,0])
    arcPart2 = rad.ObjArcCur([-0.5*wPole,0.5*gap+gapOffset+ttCoil+rMin,zcCoil], [rMin,rMin+ttCoil], [0.5*pi,1.5*pi], tlCoil, nsCoil, -jCoil)
    strPart3 = rad.ObjRecCur([-(0.5*wNotch+0.25*(wPole-wNotch)),0.5*gap+gapOffset+1.5*ttCoil+hYoke,zcCoil], [0.5*(wPole-wNotch),ttCoil,tlCoil], [jCoil,0,0])

    #trArcPart1 = rad.ObjArcCur([0.5*wPole,0.5*gap+ttCoil+rMin,zcCoil], [rMin,rMin+ttCoil], [1.5*pi,2.5*pi], tlCoil, nsCoil, -jCoil)
    #rtrPart = rad.ObjCnt([strPart1,arcPart1,strPart2,arcPart2,strPart3])

    #Transition part of the coil (from one loop to another), using ObjMltExtPgnCur
    iCoil = jCoil*ttCoil*tlCoil
    y1Tr = 0.5*gap+gapOffset+1.5*ttCoil+hYoke-0.5*ttCoil
    z1Tr = 0.5*tlCoil
    y2Tr = y1Tr+ttCoil
    z2Tr = z1Tr
    y3Tr = y2Tr
    z3Tr = z2Tr-tlCoil
    y4Tr = y3Tr-ttCoil
    z4Tr = z3Tr
    
    trPartR1 = [[0.5*wNotch,0,zcCoil+0.5*tlCoil+roCoil], [0,1,0], atCoil/ntCoil]
    trPartAll = [[trPartR1]]*ntCoil #[[trPartR1],[trPartR1],[trPartR1],[trPartR1]]
    trPartAll.append([[-slCoil*cos(atCoil),0,slCoil*sin(atCoil)]])

    trPartR2 = [[-0.5*wNotch,0,zcCoil+per-0.5*tlCoil-roCoil], [0,1,0], -atCoil/ntCoil]

    #print(trPartAll)
    trPartAll = trPartAll + [[trPartR2]]*ntCoil
    #print(trPartAll)

    trPart = rad.ObjMltExtPgnCur(0.5*wNotch, 'x', [[y1Tr,z1Tr],[y2Tr,z2Tr],[y3Tr,z3Tr],[y4Tr,z4Tr]], trPartAll, -iCoil, 'Frame->Lab')

    #trPart = rad.ObjArcCur([0.5*wNotch,0.5*gap+1.5*ttCoil+hYoke,zcCoil+0.5*tlCoil+roCoil], [roCoil,roCoil+tlCoil], [pi,atCoil+pi], ttCoil, ntCoil, -jCoil, 'man', 'y')

    #coil = rad.ObjCnt([strPart1,arcPart1,strPart2,arcPart2,strPart3,trPart])
    coilAux = rad.ObjCnt([strPart1,arcPart1,strPart2,arcPart2,strPart3])
    #coilAux = rad.ObjDpl(coil)
    #coil = rad.ObjAddToCnt(coil, [trPart])

    #coil = rad.ObjCnt([raceTrkPart, trPart])

    halfNumPer = int(0.5*nPer)
    coil = rad.ObjCnt([])

    for i in range(halfNumPer):
        trslOnePerR = rad.TrfTrsl([0,0,per*(i+0.25)])
        trslOnePerL = rad.TrfTrsl([0,0,-per*(i+0.75)])

        coilAuxR = rad.ObjDpl(coilAux)
        #coilAuxR = rad.ObjAddToCnt(coilAuxR, [trPart])

        coilAuxL = rad.ObjDpl(coilAux)
        coilAuxL = rad.ObjAddToCnt(coilAuxL, [rad.ObjDpl(trPart)])

        if(i < (halfNumPer - 1)): 
            coilAuxR = rad.ObjAddToCnt(coilAuxR, [rad.ObjDpl(trPart)])
            #coilAuxL = rad.ObjAddToCnt(coilAuxL, [trPart])

        coilAuxR = rad.TrfOrnt(coilAuxR, trslOnePerR)
        coilAuxL = rad.TrfOrnt(coilAuxL, trslOnePerL)

        coil = rad.ObjAddToCnt(coil, [coilAuxR, coilAuxL])
        #coil = rad.ObjAddToCnt(coil, [rad.TrfOrnt(coilAuxD,trslOnePerR),rad.TrfOrnt(coilAux,trslOnePerL)])
        #print('one more coil added')

    coil2 = rad.ObjDpl(coil)
    coil2 = rad.TrfOrnt(coil2, rad.TrfCmbL(rad.TrfTrsl([0,0,0.5*per]), rad.TrfInv()))

    rad.ObjDrwAtr(coil, cCoil1)
    rad.ObjDrwAtr(coil2, cCoil2)

    #coil = rad.ObjCnt([coil,coil2])

    #First and Second Terminating Coils
    coilTL1 = rad.ObjRaceTrk([0,0.5*gap+gapOffset+ttCoil+rMin,-per*(halfNumPer+0.25)], [rMin,rMin+ttCoilT1], [wPole,0], tlCoil, nsCoil, jCoilT1)
    rad.ObjDrwAtr(coilTL1, cCoilTL1)
    coilTL2 = rad.ObjRaceTrk([0,0.5*gap+gapOffset+ttCoil+rMin,-per*(halfNumPer+0.75)], [rMin,rMin+ttCoilT2], [wPole,0], tlCoil, nsCoil, -jCoilT2)
    rad.ObjDrwAtr(coilTL2, cCoilTL2)
    coilTR1 = rad.ObjRaceTrk([0,0.5*gap+gapOffset+ttCoil+rMin,per*(halfNumPer+0.25)], [rMin,rMin+ttCoilT1], [wPole,0], tlCoil, nsCoil, -jCoilT1)
    rad.ObjDrwAtr(coilTR1, cCoilTR1)
    coilTR2 = rad.ObjRaceTrk([0,0.5*gap+gapOffset+ttCoil+rMin,per*(halfNumPer+0.75)], [rMin,rMin+ttCoilT2], [wPole,0], tlCoil, nsCoil, jCoilT2)
    rad.ObjDrwAtr(coilTR2, cCoilTR2)
    
    coil = rad.ObjCnt([coil, coilTL1, coilTL2, coil2, coilTR1, coilTR2])
    coil = rad.TrfZerPara(coil, [0,0,0], [0,1,0])

    return coil

#*********************************Yoke Slice
def YokeSlice(_zc, _tlPole, _ttPole):

    #ycAll = 0.5*gap+_ttPole+0.5*hYoke
    ycAll = 0.5*gap+gapOffset+ttCoil+0.5*hYoke

    rSidePole = _ttPole+0.5*hYoke
    sidePolePart = rad.ObjArcPgnMag([0.5*wPole,ycAll], 'z', 
        [[0,_zc-0.5*_tlPole],[rSidePole,_zc-0.5*_tlPole],[rSidePole,_zc+0.5*_tlPole],[0,_zc+0.5*_tlPole]], [1.5*pi,2.5*pi], nsYoke)
   
    dxStr = 0.5*wPole - 0.5*dpYoke - rpYoke
    xcStr = 0.5*(0.5*wPole + (0.5*dpYoke + rpYoke))
    strPolePart = rad.ObjRecMag([xcStr, ycAll, _zc], [dxStr, 2*_ttPole+hYoke, _tlPole])

    dxStrUp = 0.5*wPole - 0.5*wNotch - dxStr
    xcStrUp = 0.5*(0.5*wNotch + (0.5*wNotch + dxStrUp))
    dyStrUp = rSidePole - rpYoke
    ycStrUp = ycAll + rpYoke + 0.5*dyStrUp
    strUpPolePart = rad.ObjRecMag([xcStrUp, ycStrUp, _zc], [dxStrUp, dyStrUp, _tlPole])
   
    dxStrLow = 0.5*wPole - dxStr
    xcStrLow = 0.5*dxStrLow
    ycStrLow = ycAll - rpYoke - 0.5*dyStrUp
    strLowPolePart = rad.ObjRecMag([xcStrLow, ycStrLow, _zc], [dxStrLow, dyStrUp, _tlPole])

    dxStrCen = 0.5*wPole - dxStr - 2*rpYoke
    xcStrCen = 0.5*dxStrCen
    dyStrCen = 2*rpYoke
    strCenPolePart = rad.ObjRecMag([xcStrCen, ycAll, _zc], [dxStrCen, dyStrCen, _tlPole])

    dxStrCenUp = 0.5*wPole - dxStr - dxStrUp
    xcStrCenUp = 0.5*dxStrCenUp
    dyStrCenUp = rSidePole - 0.5*dyStrCen - _ttPole
    ycStrCenUp = ycAll + 0.5*dyStrCen + 0.5*dyStrCenUp
    strCenUpPolePart = rad.ObjRecMag([xcStrCenUp, ycStrCenUp, _zc], [dxStrCenUp, dyStrCenUp, _tlPole])

    #Making "Pipes" cylindrical (filling-in corners by extruded triangles)
    cornExtTrgs = rad.ObjCnt([])
    xcp = 0.5*dpYoke; ycp = ycAll
    alp0 = 0.5*pi/npcYoke
    xCorn = xcp + rpYoke; yCorn = ycp + rpYoke
    x1 = xcp + rpYoke; y1 = ycp
    alp = alp0
    for i in range(npcYoke):
        x2 = xcp + rpYoke*cos(alp); y2 = ycp + rpYoke*sin(alp)
        cornExtTrgs = rad.ObjAddToCnt(cornExtTrgs, [rad.ObjThckPgn(_zc, _tlPole, [[x1,y1],[xCorn,yCorn],[x2,y2]], 'z')])
        x1 = x2; y1 = y2; alp += alp0

    xCorn = xcp - rpYoke
    x1 = xcp; y1 = ycp + rpYoke
    for i in range(npcYoke):
        x2 = xcp + rpYoke*cos(alp); y2 = ycp + rpYoke*sin(alp)
        cornExtTrgs = rad.ObjAddToCnt(cornExtTrgs, [rad.ObjThckPgn(_zc, _tlPole, [[x1,y1],[xCorn,yCorn],[x2,y2]], 'z')])
        x1 = x2; y1 = y2; alp += alp0

    yCorn = ycp - rpYoke
    x1 = xcp - rpYoke; y1 = ycp
    for i in range(npcYoke):
        x2 = xcp + rpYoke*cos(alp); y2 = ycp + rpYoke*sin(alp)
        cornExtTrgs = rad.ObjAddToCnt(cornExtTrgs, [rad.ObjThckPgn(_zc, _tlPole, [[x1,y1],[xCorn,yCorn],[x2,y2]], 'z')])
        x1 = x2; y1 = y2; alp += alp0

    xCorn = xcp + rpYoke
    x1 = xcp; y1 = ycp - rpYoke
    for i in range(npcYoke):
        x2 = xcp + rpYoke*cos(alp); y2 = ycp + rpYoke*sin(alp)
        cornExtTrgs = rad.ObjAddToCnt(cornExtTrgs, [rad.ObjThckPgn(_zc, _tlPole, [[x1,y1],[xCorn,yCorn],[x2,y2]], 'z')])
        x1 = x2; y1 = y2; alp += alp0

    slice = rad.ObjCnt([sidePolePart, strPolePart, strUpPolePart, strLowPolePart, strCenPolePart, strCenUpPolePart, cornExtTrgs])
    return slice

#*********************************Yoke
def Yoke():

    #Case of Rotation by pi around [0,1,0] vector
    sliceHalfPole = YokeSlice(0.25*tlPole, 0.5*tlPole, ttPole)
    sliceHalfPole = rad.ObjDivMag(sliceHalfPole, [sbPole[0],sbPole[1],(int(round(0.5*sbPole[2][0])))], 'Frame->LabTot')
    
    yoke = rad.ObjCnt([sliceHalfPole])

    zc = 0.25*per
    tlBwPoles = 0.5*per - tlPole

    for i in range(nPer + 2):
    #for i in range(1):

        sliceBwPoles = YokeSlice(zc, tlBwPoles, 0)
        sliceBwPoles = rad.ObjDivMag(sliceBwPoles, sbYokeBwPoles, 'Frame->LabTot')
        zc += 0.25*per

        slicePole = YokeSlice(zc, tlPole, ttPole)
        slicePole = rad.ObjDivMag(slicePole, sbPole, 'Frame->LabTot')

        zc += 0.25*per
        #yoke = rad.ObjAddToCnt(yoke, [sliceBwPoles])
        yoke = rad.ObjAddToCnt(yoke, [sliceBwPoles, slicePole])

    rad.MatApl(yoke, mtYoke)

    symYZ = rad.TrfPlSym([0,0,0], [1,0,0])
    yoke1 = rad.TrfOrnt(rad.ObjDpl(yoke), symYZ)
    yoke = rad.ObjCnt([yoke, yoke1])

    rot = rad.TrfRot([0,0,0], [0,1,0], pi)
    yoke = rad.TrfMlt(yoke, rot, 2)

    #yoke = rad.TrfZerPerp(yoke, [0,0,0], [1,0,0]) #This one is not allowed because of coils
    yoke = rad.TrfZerPara(yoke, [0,0,0], [0,1,0])

    return yoke

#*********************************Undulator
def SCU():

    coil = Coil()
    yoke = Yoke()

    und = rad.ObjCnt([coil, yoke])

    return und

#*********************************All Calculations
if __name__=="__main__":

    #Initialize MPI
    rank = rad.UtiMPI('on')

    #General Undulator Parameters
    per = 14. #Undulator Period [mm]
    nPer = 4 #30 #20 #Number of Periods (may need to even?)

    gap = 5.0 #4.0 #Magnetic Gap [mm]

    hYoke = 33 #30 #Height of Return Yoke [mm]
    
    wPole = 66.66 #50 #Width of Pole [mm]
    
    #wNotch = 10 #Width of Notch in Yoke for winding
    #hPole = #4.65 #Height of Notch in Yoke for winding

    nsYoke = 12 #16 #24 #Number of segments in arc parts of the Yoke
    
    air = 0.1 #Air gap between Coil and Yoke [mm]

    gapOffset = 0.

    tlCoil = 4.4 #Main Coil Thickness in Longitudinal Direction [mm]
    ttCoil = 3.9 #Main Coil Thickness in Transverse Plane [mm]

    tlPole = 0.5*per - tlCoil - 2*air  #Main Pole Thickness in Longitudinal Direction [mm]
    ttPole = ttCoil + gapOffset #Main Pole Height in Transverse Plane [mm]

    dpYoke = 30 #Horizontal distance bw cooling "Pipe" centers in Yoke
    rpYoke = 17.27/2. #Radius of cooling "Pipe" in Yoke
    
    npcYoke = 5 #Number of extruded triangles filling-in corners to make "Pipe"

    ttCoilT1 = 3.6 #First Terminating Coil Thickness in Transverse Plane [mm]
    ttCoilT2 = 2.1 #Second Terminating Coil Thickness in Transverse Plane [mm]

    nsCoil = 24 #Main Coil Number of Segments (approximating the coil)

    roCoil = 0.5 #Offset of Rotation Center of Coil Transition Region [mm]

    rtCoil = roCoil + 0.5*tlCoil #Average Radius of Coil Turning at Transition
    #print('rtCoil=', rtCoil)

    #slCoil = 16. #Length of Straight Coil part at Transition
    #cos_at = (per - 2*rtCoil)/(slCoil - 2*rtCoil)
    #atCoil = acos((per - 2*rtCoil)/(slCoil - 2*rtCoil)) #*(180/pi) #Angle of Coils Turning at Transition
    #wNotch = (slCoil + 2*rtCoil)*sin(atCoil) #Width of Notch in Yoke for winding
    
    atCoil = (45/180)*pi #Angle of Coils Turning at Transition
    slCoil = (per - 2*rtCoil*(1 - cos(atCoil)))/sin(atCoil) #Length of Straight Coil part at Transition
    wNotch = 2*rtCoil*sin(atCoil) + slCoil*cos(atCoil)

    #print('atCoil=', atCoil*(180/pi), 'deg.', ' wNotch=', wNotch, ' slCoil=', slCoil)

    ntCoil = 5 #Number of Segments in Transition Turning Part of Coil

    jCoil = -1300 #Coil Current Density [A/mm^2]
    jCoilT1 = jCoil
    jCoilT2 = jCoil

    cCoil1 = [0.2,0.1,0.8] #Coil Color (for one part of the coil)
    cCoil2 = [0.8,0.2,0.1] #Coil Color (for one part of the coil)
    cCoilTL1 = [0.8,0.2+0.2,0.1]
    cCoilTL2 = [0.2,0.1+0.6,0.8]
    cCoilTR1 = [0.2,0.1+0.3,0.8]
    cCoilTR2 = [0.8,0.2+0.6,0.1] 

    #Iron Material data
    arHM_T = [
        [0.,0.],[0.000056,0.099944],[0.000101,0.199899],[0.000143,0.299857],[0.000186,0.399814],[0.000203,0.499797],[0.000275,0.599725],
        [0.000324,0.699676],[0.000377,0.799623],[0.000433,0.899567],[0.000496,0.999504],[0.000571,1.09943],[0.000667,1.19933],[0.000735,1.24927],
        [0.000823,1.29918],[0.000964,1.34904],[0.001166,1.39883],[0.001479,1.44852],[0.001875,1.49813],[0.002348,1.54765],[0.002963,1.59704],
        [0.003446,1.62655],[0.003929,1.64607],[0.005151,1.69485],[0.007,1.743],[0.010001,1.79],[0.013215,1.83679],[0.017273,1.88273],[0.021,1.904],
        [0.025305,1.9247],[0.029319,1.94568],[0.033334,1.96667],[0.040017,1.98498],[0.0467,2.0033],[0.059999,2.04],[0.07,2.08],[0.08,2.095],
        [0.09,2.11],[0.119999,2.13],[0.149999,2.15],[0.198999,2.151],[0.247999,2.152],[0.297499,2.1525],[0.346999,2.153],[0.396499,2.1535],
        [0.445999,2.154],[0.554999,2.145],[0.644999,2.155],[0.744499,2.1555],[0.843999,2.156]]
    mtYoke = rad.MatSatIsoTab(arHM_T)
    #print('Pol. mat. index:', mtYoke)

    sbPole = [[6,4],[5,10],[4,1]] #Pole Subdivision Params
    #sbPole = [[6,5],[6,10],[4,1]] #Pole Subdivision Params
    sbYokeBwPoles = [[3,1],[2,1],[2,1]] #Pole Subdivision Params

    #sbPole = [[7,3],[7,5],[3,1]] #Pole Subdivision Params
    #sbYokeBwPoles = [[7,3],[7,5],[2,1]] #Pole Subdivision Params

    #Instantiating Undulator Model
    #und = Coil()
    #und = Yoke()
    und = SCU()

    #Displaying Undulator Model in 3D
    if(rank <= 0):
        from uti_radia_vtk import ObjDrwPyVista
        ObjDrwPyVista(und)

    #id_solve(und, _prc=0.0001, _max_it=5000, _print=2)
    id_solve(und, _print=(2 if(rank <= 0) else 0))
    if(rank <= 0): sys.stdout.flush()

    #Checking Current Density
    #x0 = -20; y0 = 38; z0 = 0
    #x0 = 20; y0 = 38; z0 = 0
    #x0 = 40; y0 = 21; z0 = 0
    #x0 = 0; y0 = 5; z0 = 0
    #x0 = -40; y0 = 21; z0 = 0
    #x0 = 0; y0 = 38; z0 = 0.5*per
    #x0 = -16.236753/2 + 2; y0 = 38; z0 = 0
    #x0 = 0; y0 = 5.5; z0 = -135
    #x0 = 0; y0 = 5.5; z0 = -142
    #x0 = 0; y0 = 5.5; z0 = -150
    #Jtst = rad.Fld(und, 'J', [x0,y0,z0])
    #print(Jtst)

    By0 = rad.Fld(und, 'by', [0,0,0])
    if(rank <= 0):
        print('Peak magnetic field: By0=', By0)
        sys.stdout.flush()
    
    #Calculating Magnetic Field
    arBh, arBv, ns, sRange, Keff, E1 = mag_fld_vs_long_pos(_id=und, _ns_per=50, _per=per, _nper=nPer, _long_ax='z', _rank=rank, _print=2, _half_nper_ext=8)

    #zRange = (nPer + 10)*per #Range for Magnetic Field calculation along Axis
    #nz = 50*(nPer + 10) #Number of Points for Magnetic Field calculation along Axis
    #zStep = zRange/(nz - 1)
    #ByOnAx = rad.Fld(und, 'by', [[0,0,-0.5*zRange+iz*zStep] for iz in range(nz)])

    if(rank <= 0):
        #print('arBh=', arBh)
        #print('arBv=', arBv)
        #sys.stdout.flush()

        #Plot the calculated Magnetic Field
        uti_plot1d(arBv, [-0.5e-03*sRange,0.5e-03*sRange,ns], ['Longitudinal Position [m]', 'By [T]', 'Vertical Magnetic Field']) #, ['mm', 'T'])
        #uti_plot1d(ByOnAx, [-0.5*zRange,0.5*zRange,nz], ['Longitudinal Position [m]', 'By [T]', 'Vertical Magnetic Field']) #, ['mm', 'T'])

        uti_plot_show()

    rad.UtiMPI('off')

