from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


au2ev = 27.211386245988
au2cm = 219474.6313705
au2K = 315777.663
amu2au = 1822.888486209
pi = math.pi


elemMass = {
    "H": 1.00782503223, "D": 2.01410177812, "He": 4.002602,
    "Li": 6.938, "N": 14.00307400443, "O": 15.99491461957,
    "F": 18.99840316273, "S": 32.06, "Cl": 34.968852682,
    "Ar": 39.9623831237,
}


@dataclass(frozen=True)
class para1D:
    atoms: tuple[str, str, str]
    energyUnit: int
    Zrange: tuple[float, float]
    nZ: int
    Zc: float
    delta: float
    Ec: float
    EcEv: float
    Ecol: np.ndarray
    EcolEv: np.ndarray
    Clr: float
    ZlrRange: float
    mu: float
    dZ: float


@dataclass(frozen=True)
class prop1D:
    Zgrid: np.ndarray
    dZ: float
    ZminAll: float
    ZmaxAll: float
    psi0: np.ndarray
    ae0: np.ndarray
    coeffVec: np.ndarray
    phase0: complex
    kinFact: float
    Hplus: float
    Hminus: float
    timeStep: float
    Zf: float
    iZf: int
    leftRange: float
    leftC: float
    fluxFloor: float


@dataclass(frozen=True)
class scanRes:
    deltaX: float
    Cabs: float
    ZabsStart: float
    timeTot: float
    nStep: int
    PofE: np.ndarray
    Pmax: float
    probFlux: np.ndarray
    probCurrent: np.ndarray
    accept: bool


# --------------- helpers: parsing, grid, packet ---------------

def getArg() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="1D absorber test (NIPtest logic)")
    ap.add_argument("--abcInf", type=Path, default=Path("ABC.inf"))
    ap.add_argument("--deltaX", type=str, default=None)
    ap.add_argument("--Cscan", type=str, default=None)
    ap.add_argument("--timeStep", type=float, default=400.0)
    ap.add_argument("--leftRange", type=float, default=20.0)
    ap.add_argument("--leftC", type=float, default=0.05)
    ap.add_argument("--probeSigma", type=float, default=6.0)
    ap.add_argument("--fluxFloor", type=float, default=1.0e-12)
    ap.add_argument("--coeffTol", type=float, default=1.0e-12)
    ap.add_argument("--maxOrder", type=int, default=512)
    ap.add_argument("--smoke", action="store_true")
    return ap.parse_args()

def parseNum(t: str) -> float:
    return float(t.replace("D","E").replace("d","e"))

def parseVec(t: str) -> np.ndarray:
    return np.asarray([parseNum(s.strip()) for s in t.split(",") if s.strip()], dtype=float)

def readNamelist(fp: Path) -> dict[str,str]:
    d: dict[str,str] = {}
    for raw in fp.read_text(encoding="utf-8").splitlines():
        ln = raw.split("!",1)[0].strip()
        if not ln or ln.startswith("&") or ln == "/": continue
        if "=" not in ln: continue
        k,v = ln.split("=",1); d[k.strip()] = v.strip().rstrip(",")
    return d

def parseAtoms(t: str) -> tuple[str,str,str]:
    a = [s.strip().strip("'").strip('"') for s in t.split(",")]
    if len(a)!=3: raise ValueError(f"bad atom line: {t}")
    return a[0],a[1],a[2]

def energyTrans(u: int) -> float:
    m = {1:1.0/au2cm, 2:1.0/au2ev, 3:1.0/au2K, 4:1.0}
    if u not in m: raise ValueError(f"unsupported energyUnit={u}")
    return m[u]

def getMassABC(atoms: tuple[str,str,str]) -> float:
    mA=elemMass[atoms[0]]*amu2au; mB=elemMass[atoms[1]]*amu2au; mC=elemMass[atoms[2]]*amu2au
    return mA*(mB+mC)/(mA+mB+mC)

def getEcol(Er: np.ndarray, dE: float, tr: float) -> tuple[np.ndarray,np.ndarray]:
    n=int((Er[1]-Er[0])/dE)+1; raw=Er[0]+np.arange(n,dtype=float)*dE; raw=np.minimum(raw,Er[1])
    return raw*tr, raw

def initPara(fp: Path) -> para1D:
    d=readNamelist(fp); atoms=parseAtoms(d["Atoms"])
    eu=int(parseNum(d["energyUnit"])); tr=energyTrans(eu)
    Zr=tuple(parseVec(d["IALR%Z_range"])); nZ=int(parseNum(d["IALR%nZ_IALR"]))
    dZ=(Zr[1]-Zr[0])/(nZ+1); Ec=parseNum(d["initWP%Ec"])*tr
    Ecol,EcolEv=getEcol(parseVec(d["E_range"]),parseNum(d["dE"]),tr)
    return para1D(atoms=atoms,energyUnit=eu,Zrange=(float(Zr[0]),float(Zr[1])),nZ=nZ,
        Zc=parseNum(d["initWP%Zc"]),delta=parseNum(d["initWP%delta"]),
        Ec=Ec,EcEv=Ec*au2ev,Ecol=Ecol,EcolEv=EcolEv,
        Clr=parseNum(d["Vabs%Clr"]),ZlrRange=parseNum(d["Vabs%Zlr_range"]),
        mu=getMassABC(atoms),dZ=dZ)

def sinDVRGrid(n: int, Zmin: float, Zmax: float) -> tuple[np.ndarray,float]:
    dZ=(Zmax-Zmin)/(n+1); return Zmin+np.arange(1,n+1,dtype=float)*dZ, dZ

def travelGaussian(E0: float, delta: float, Z0: float, mass: float,
                   Z: np.ndarray, dZ: float) -> np.ndarray:
    """Right-moving traveling Gaussian: exp(+ik0 Z) * envelope, like NIPtest INIT."""
    k0 = math.sqrt(2.0*E0*mass)
    fact = (1.0/(pi*delta*delta))**0.25
    env = np.exp(-((Z-Z0)**2)/(2.0*delta*delta))
    return (fact * env * np.exp(1.0j*k0*Z) * math.sqrt(dZ)).astype(np.complex128)


# --------------- AE0: incident spectral amplitude (NIPtest INIT) ---------------

def buildAe0(psi0: np.ndarray, Zgrid: np.ndarray, dZ: float,
             Ecol: np.ndarray, mu: float) -> np.ndarray:
    """AE0(E) = sum_i sqrt(2mu/pi)/sqrt(k) * conj(exp(ikZ_i)) * psi0_DVR(i) * sqrt(dZ) / 2.
    For right-moving packet, project onto exp(-ikZ) (incoming convention)."""
    nE = len(Ecol)
    ae0 = np.zeros(nE, dtype=np.complex128)
    kVec = np.sqrt(np.maximum(2.0*mu*Ecol, 0.0))
    fac0 = math.sqrt(2.0*mu/pi)
    sqDZ = math.sqrt(dZ)
    for iE in range(nE):
        if kVec[iE] <= 0.0: continue
        ae0[iE] = (fac0/math.sqrt(kVec[iE])) * np.sum(
            np.exp(-1.0j*kVec[iE]*Zgrid) * psi0 * sqDZ) * 0.5
    return ae0


# --------------- absorber (potent.f90 Vabs_EXP) ---------------

def VabsExp(Cabs: float, x0: float, xmax: float, x: np.ndarray) -> np.ndarray:
    return np.exp(-Cabs * ((x-x0)/(xmax-x0))**2)

def setVabs(Zgrid: np.ndarray, rangeAbs: float, Cabs: float) -> tuple[np.ndarray,float]:
    dZ=Zgrid[1]-Zgrid[0]; rangeAll=Zgrid[-1]+dZ; rNA=rangeAll-rangeAbs
    mask=np.ones_like(Zgrid,dtype=float)
    if rangeAbs>0 and Cabs>0:
        idx=Zgrid>=rNA
        if np.any(idx): mask[idx]=VabsExp(Cabs,rNA,rangeAll,Zgrid[idx])
    return mask, float(rNA)

def setLeftVabs(Zgrid: np.ndarray, rangeAbs: float, Cabs: float) -> tuple[np.ndarray,float]:
    dZ=Zgrid[1]-Zgrid[0]; rAll=Zgrid[0]-dZ; rNA=rAll+rangeAbs
    mask=np.ones_like(Zgrid,dtype=float)
    if rangeAbs>0 and Cabs>0:
        idx=Zgrid<=rNA
        if np.any(idx): mask[idx]=np.exp(-Cabs*((rNA-Zgrid[idx])/(rNA-rAll))**2)
    return mask, float(rNA)

def buildVabsMask(Zgrid, rangeAbs, Cabs, leftRange, leftC):
    rM,zS=setVabs(Zgrid,rangeAbs,Cabs); lM,lE=setLeftVabs(Zgrid,leftRange,leftC)
    return rM*lM, zS, lE


# --------------- probe, grid, Chebyshev ---------------

def chooseProbeLeft(Zgrid: np.ndarray, Zc: float, delta: float,
                    leftRange: float, probeSigma: float) -> tuple[float,int]:
    """Place probe to the LEFT of the packet (NIPtest POSFLUX geometry)."""
    dZ=Zgrid[1]-Zgrid[0]; ZminAll=Zgrid[0]-dZ
    leftEnd=ZminAll+leftRange
    Zf=max(Zc - probeSigma*delta, leftEnd + 4.0*dZ)
    if Zf >= Zc - 2.0*delta:
        raise ValueError("probe too close to packet; increase probeSigma or extend grid left")
    iZf=int(np.argmin(np.abs(Zgrid-Zf)))
    iZf=max(1,min(len(Zgrid)-2,iZf))
    return float(Zgrid[iZf]), iZf

def autoGrid(para: para1D, deltaXMax: float, probeSigma: float) -> tuple[np.ndarray,float,int]:
    dZ=para.dZ
    Zneed=deltaXMax+para.Zc+probeSigma*para.delta+4.0*dZ
    Zmax=max(para.Zrange[1],Zneed)
    nZ=max(para.nZ,int(math.ceil((Zmax-para.Zrange[0])/dZ))-1)
    Zgrid,dZn=sinDVRGrid(nZ,para.Zrange[0],Zmax)
    return Zgrid,dZn,nZ

def getCoeff(zVal: float, tol: float, maxOrd: int) -> np.ndarray:
    order=32
    while order<=maxOrd:
        nQ=max(2048,8*(order+1)); th=pi*(np.arange(nQ,dtype=float)+0.5)/nQ
        val=np.exp(-1.0j*zVal*np.cos(th))
        basis=np.cos(np.outer(np.arange(order+1,dtype=float),th))
        c=(2.0/nQ)*(basis@val); c[0]*=0.5
        idx=np.where(np.abs(c)>tol)[0]
        if idx.size==0: return c[:1]
        last=int(idx[-1])
        if last<=order-8 or order==maxOrd: return c[:last+1]
        order=min(2*order,maxOrd)
    raise RuntimeError("failed to build coefficient vector")

def Hpsi(psi, kF, Hp, Hm, out):
    out[1:-1]=(2.0*psi[1:-1]-psi[:-2]-psi[2:])*kF
    out[0]=(2.0*psi[0]-psi[1])*kF; out[-1]=(2.0*psi[-1]-psi[-2])*kF
    out-=Hp*psi; out/=Hm

def propCheby(psi, cV, ph0, kF, Hp, Hm):
    r=cV[0]*psi
    if len(cV)==1: return ph0*r
    p0=psi.copy(); p1=np.empty_like(psi); Hpsi(p0,kF,Hp,Hm,p1); r+=cV[1]*p1
    if len(cV)==2: return ph0*r
    p2=np.empty_like(psi)
    for c in cV[2:]:
        Hpsi(p1,kF,Hp,Hm,p2); p2*=2.0; p2-=p0; r+=c*p2; p0,p1,p2=p1,p2,p0
    return ph0*r


# --------------- buildProp, getTimeTot ---------------

def buildProp(para: para1D, arg, deltaXVec: np.ndarray) -> prop1D:
    Zgrid,dZ,_=autoGrid(para,float(np.max(deltaXVec)),arg.probeSigma)
    psi0=travelGaussian(para.Ec,para.delta,para.Zc,para.mu,Zgrid,dZ)
    ae0=buildAe0(psi0,Zgrid,dZ,para.Ecol,para.mu)
    Zf,iZf=chooseProbeLeft(Zgrid,para.Zc,para.delta,arg.leftRange,arg.probeSigma)
    Hmax=2.0/(para.mu*dZ*dZ); Hp=0.5*Hmax; Hm=0.5001*Hmax
    kF=1.0/(2.0*para.mu*dZ*dZ)
    cV=getCoeff(Hm*arg.timeStep,arg.coeffTol,arg.maxOrder)
    ph0=np.exp(-1.0j*Hp*arg.timeStep)
    return prop1D(Zgrid=Zgrid,dZ=dZ,ZminAll=float(Zgrid[0]-dZ),ZmaxAll=float(Zgrid[-1]+dZ),
        psi0=psi0,ae0=ae0,coeffVec=cV,phase0=ph0,kinFact=kF,Hplus=Hp,Hminus=Hm,
        timeStep=arg.timeStep,Zf=Zf,iZf=iZf,
        leftRange=arg.leftRange,leftC=arg.leftC,fluxFloor=arg.fluxFloor)

def getTimeTot(Zturn,Zc,Zstop,delta,k0,mu,dt):
    vg=k0/mu; path=(Zturn-Zc)+(Zturn-Zstop)+12.0*delta
    tTot=max(dt,1.1*path/max(vg,1e-12)); nS=int(math.ceil(tTot/dt))
    return nS*dt, nS


# --------------- NIPtest core: S2S_EX2, S2S_EXM, runNIPCase ---------------

def probS2S2(PHIE0, ae0, dt, Ecol, mu, Zf):
    """S2S_EX2: S(E)=PHIE0*dt/AE0 * sqrt(k/(2pi*mu)) * exp(i*k*Zf)."""
    k=np.sqrt(np.maximum(2.0*mu*Ecol,0.0))
    a0=np.where(np.abs(ae0)>1e-30,ae0,1e-30)
    S=PHIE0*dt/a0*np.sqrt(k/(2.0*pi*mu))*np.exp(1.0j*k*Zf)
    return np.abs(S)**2

def probCurrentStyle(PHIE0, PHIE0p, ae0, dt, mu):
    """S2S_EXM: P=-Im(conj(PHIE0)*PHIE0p)*dt^2/(mu*|AE0|^2)."""
    a2=np.maximum(np.abs(ae0)**2,1e-30)
    return -np.imag(np.conj(PHIE0)*PHIE0p)*dt*dt/(mu*a2)

def runNIPCase(para: para1D, prop: prop1D, deltaX: float, Cabs: float) -> scanRes:
    mask,zAS,lE=buildVabsMask(prop.Zgrid,deltaX,Cabs,prop.leftRange,prop.leftC)
    if zAS<=para.Zc+2.0*para.delta:
        raise ValueError(f"Δx={deltaX:.3f} too large: absorber overlaps packet")

    k0=math.sqrt(2.0*para.mu*para.Ec)
    tTot,nStep=getTimeTot(prop.ZmaxAll,para.Zc,prop.Zf,para.delta,k0,para.mu,prop.timeStep)

    nE=len(para.Ecol)
    PHIE0=np.zeros(nE,dtype=np.complex128)
    PHIE0p=np.zeros(nE,dtype=np.complex128)
    DR2=math.sqrt(prop.dZ)
    psi=prop.psi0.copy(); iZf=prop.iZf; ttot=0.0

    for _ in range(nStep):
        psi=propCheby(psi,prop.coeffVec,prop.phase0,prop.kinFact,prop.Hplus,prop.Hminus)
        psi*=mask
        ttot+=prop.timeStep
        C1=psi[iZf]/DR2
        CD1=(psi[iZf+1]-psi[iZf-1])/(2.0*prop.dZ*DR2)
        ef=np.exp(1.0j*para.Ecol*ttot)
        PHIE0+=C1*ef; PHIE0p+=CD1*ef

    dt=prop.timeStep
    pF=probS2S2(PHIE0,prop.ae0,dt,para.Ecol,para.mu,prop.Zf)
    pC=probCurrentStyle(PHIE0,PHIE0p,prop.ae0,dt,para.mu)

    ae0sq=np.abs(prop.ae0)**2; ae0mx=float(np.max(ae0sq))
    sup=ae0sq>ae0mx*prop.fluxFloor
    if np.any(sup):
        Pm=float(np.max(pF[sup])); acc=bool(np.all(pF[sup]<1e-3))
    else:
        Pm=float(np.max(pF)); acc=False

    return scanRes(deltaX=deltaX,Cabs=Cabs,ZabsStart=zAS,timeTot=tTot,nStep=nStep,
        PofE=pF,Pmax=Pm,probFlux=pF,probCurrent=pC,accept=acc)


# --------------- scan helpers ---------------

def getDeltaXScan(arg, para):
    if arg.smoke: return np.asarray([80.0,100.0,120.0],dtype=float)
    if arg.deltaX: return parseVec(arg.deltaX)
    v=np.arange(40.0,171.0,20.0,dtype=float)
    if para.ZlrRange not in v: v=np.unique(np.sort(np.append(v,para.ZlrRange)))
    return v

def getCscan(arg, para):
    if arg.smoke: return np.asarray([0.005,0.01,0.02],dtype=float)
    if arg.Cscan: return parseVec(arg.Cscan)
    v=np.asarray([0.001,0.002,0.005,0.01,0.02,0.05,0.1],dtype=float)
    if para.Clr not in v: v=np.unique(np.sort(np.append(v,para.Clr)))
    return v

def compressVal(vl):
    if not vl: return "none"
    vl=sorted(vl); return ", ".join(f"{x:.6g}" for x in vl)


# --------------- plotting ---------------

def plotScan(xV, yV, xL, outF, xLog=False):
    fig,ax=plt.subplots(figsize=(8.0,5.2),constrained_layout=True)
    ax.plot(xV,yV,marker="o",markersize=6.5,linewidth=2.0,color="blue",
            markerfacecolor="white",markeredgewidth=1.4)
    ax.axhline(1e-3,color="tab:red",linestyle="--",linewidth=1.4,label=r"$10^{-3}$")
    ax.set_xlabel(xL,fontsize=13); ax.set_ylabel(r"$P_{\max}$",fontsize=13)
    ax.set_yscale("log")
    if xLog: ax.set_xscale("log")
    ax.grid(False)
    ax.tick_params(axis="both",which="major",labelsize=11,width=1.2,length=6)
    ax.tick_params(axis="both",which="minor",width=1.0,length=3)
    for s in ax.spines.values(): s.set_linewidth(1.2)
    ax.legend(frameon=False,fontsize=11)
    fig.savefig(outF,dpi=500,facecolor="white",bbox_inches="tight"); plt.close(fig)

def plotBestPE(best, EcolEv, outF):
    fig,ax=plt.subplots(figsize=(8.0,5.2),constrained_layout=True)
    ax.plot(EcolEv,best.PofE,linewidth=2.0,color="blue")
    ax.axhline(1e-3,color="tab:red",linestyle="--",linewidth=1.4)
    ax.set_xlabel("Collision energy (eV)",fontsize=13)
    ax.set_ylabel(r"$P(E)$",fontsize=13); ax.set_yscale("log"); ax.grid(False)
    ax.tick_params(axis="both",which="major",labelsize=11,width=1.2,length=6)
    ax.tick_params(axis="both",which="minor",width=1.0,length=3)
    for s in ax.spines.values(): s.set_linewidth(1.2)
    ax.set_title(rf"Best NIP test: $\Delta x={best.deltaX:.3g}$, $C={best.Cabs:.3g}$",fontsize=13)
    fig.savefig(outF,dpi=500,facecolor="white",bbox_inches="tight"); plt.close(fig)


# --------------- summary ---------------

def writeSummary(outF, para, prop, dR, cR, pickDX, best):
    pD=[r.deltaX for r in dR if r.accept]; pC=[r.Cabs for r in cR if r.accept]
    L=["Direct 1D absorber summary (NIPtest logic)","="*44,"",
       f"System                 : {para.atoms[0]} + {para.atoms[1]}{para.atoms[2]}",
       f"Reduced mass mu (a.u.) : {para.mu:.9f}",
       f"Collision grid (eV)    : {para.EcolEv[0]:.6f} ... {para.EcolEv[-1]:.6f}",
       f"Grid points            : {len(prop.Zgrid)}",
       f"dZ (bohr)              : {prop.dZ:.9f}",
       f"Zf / POSFLUX (bohr)    : {prop.Zf:.9f}  (LEFT of packet)",
       f"Initial norm check     : {np.sum(np.abs(prop.psi0)**2):.9f}",
       f"max |AE0|^2            : {float(np.max(np.abs(prop.ae0)**2)):.6e}",
       f"Current C from ABC.inf : {para.Clr:.6f}",
       f"Current Δx from ABC.inf: {para.ZlrRange:.6f}",
       f"Fixed left guard       : Δx = {prop.leftRange:.6f}, C = {prop.leftC:.6f}","",
       "This script follows NIPtest logic (FLUX0/FLUX2 -> S2S_EX2):",
       "  1. right-moving traveling Gaussian packet (exp(+ik0 Z) * envelope)",
       "  2. probe/POSFLUX placed to the LEFT of the packet center",
       "  3. propagate with Chebyshev; apply Vabs_EXP absorber every step",
       "  4. at each step accumulate PHIE0(E) and PHIE0p(E) at the probe",
       "  5. S(E) = PHIE0*dt/AE0 * sqrt(k/(2pi*mu)) * exp(i*k*Zf)  [S2S_EX2]",
       "  6. P(E) = |S(E)|^2  (reflection probability)",
       "  7. cross-check: P_current from S2S_EXM","",
       "Δx scan","-------"]
    for r in dR:
        L.append(f"Δx = {r.deltaX:8.3f}  P_max = {r.Pmax:11.4e}  accept = {r.accept!s:5s}  steps = {r.nStep:5d}")
    L+=["",f"Accepted Δx values     : {compressVal(pD)}",
        f"Chosen Δx for C scan   : {pickDX:.6f}","","C scan","------"]
    for r in cR:
        L.append(f"C  = {r.Cabs:8.4f}  P_max = {r.Pmax:11.4e}  accept = {r.accept!s:5s}  steps = {r.nStep:5d}")
    L+=["",f"Accepted C values      : {compressVal(pC)}","",
        "Recommended pair","----------------",
        f"Δx = {best.deltaX:.6f}",f"C  = {best.Cabs:.6f}",
        f"P_max = {best.Pmax:.6e}",f"Accepted = {best.accept}"]
    outF.write_text("\n".join(L)+"\n",encoding="utf-8")


# --------------- main ---------------

def main():
    arg=getArg(); para=initPara(arg.abcInf)
    dXV=getDeltaXScan(arg,para); cV=getCscan(arg,para)
    prop=buildProp(para,arg,dXV)

    print(f"Loaded {arg.abcInf}")
    print(f"System: {para.atoms[0]} + {para.atoms[1]}{para.atoms[2]}")
    print(f"Reduced mass mu = {para.mu:.9f} a.u.")
    print(f"Ecol: {para.EcolEv[0]:.6f} ... {para.EcolEv[-1]:.6f} eV ({len(para.EcolEv)} pts)")
    print(f"Norm check = {np.sum(np.abs(prop.psi0)**2):.9f}")
    print(f"POSFLUX = {prop.Zf:.9f} bohr (LEFT of packet)")
    print(f"max |AE0|^2 = {float(np.max(np.abs(prop.ae0)**2)):.6e}")
    print(f"Chebyshev order {len(prop.coeffVec)-1} per step")
    print("NIPtest logic: FLUX0/FLUX2 -> S2S_EX2\n")

    dR=[]
    print("Scanning Δx at fixed C ...")
    for dx in dXV:
        r=runNIPCase(para,prop,float(dx),para.Clr); dR.append(r)
        print(f"Δx={r.deltaX:7.3f} C={r.Cabs:7.4f} P_max={r.Pmax:11.4e} accept={r.accept!s:5s}")

    pD=[r.deltaX for r in dR if r.accept]
    pickDX=float(min(pD)) if pD else float(min(dR,key=lambda x:x.Pmax).deltaX)

    cR=[]
    print(f"\nChosen Δx={pickDX:.6f}; scanning C ...")
    for c in cV:
        r=runNIPCase(para,prop,pickDX,float(c)); cR.append(r)
        print(f"Δx={r.deltaX:7.3f} C={r.Cabs:7.4f} P_max={r.Pmax:11.4e} accept={r.accept!s:5s}")

    pCR=[r for r in cR if r.accept]
    best=min(pCR,key=lambda x:(x.Pmax,abs(x.Cabs-para.Clr))) if pCR else min(cR,key=lambda x:x.Pmax)

    wd=arg.abcInf.resolve().parent
    plotScan(np.asarray([r.deltaX for r in dR]),np.asarray([max(r.Pmax,1e-16) for r in dR]),
             r"$\Delta x$ (bohr)",wd/"Vabs_deltaXScan.png",False)
    plotScan(np.asarray([r.Cabs for r in cR]),np.asarray([max(r.Pmax,1e-16) for r in cR]),
             "C",wd/"Vabs_Cscan.png",True)
    plotBestPE(best,para.EcolEv,wd/"Vabs_bestPE.png")
    writeSummary(wd/"Vabs_summary.log",para,prop,dR,cR,pickDX,best)

    print(f"\nDone.")
    print(f"Accepted Δx: {compressVal([r.deltaX for r in dR if r.accept])}")
    print(f"Accepted C : {compressVal([r.Cabs for r in cR if r.accept])}")
    print(f"Recommended: Δx={best.deltaX:.6f}, C={best.Cabs:.6f}, P_max={best.Pmax:.4e}, accept={best.accept}")

if __name__=="__main__":
    main()
