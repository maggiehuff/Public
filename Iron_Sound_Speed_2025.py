import numpy as np
import matplotlib.pyplot as plt

#this code is meant to take the data from https://doi.org/10.1103/PhysRevB.109.184311
#tables I and II and calculate sound speed and then Gruneisen parameter
#goal is to be able to recreate table III

class experiment():
    def __init__(self,Us_qz,Us_fe_avg,Us_fe,F_value,P,rho):
        #tell the program the quartz shock velocity (from table I) 
        self.Us_qz=Us_qz
        #tell the program the average iron shock velocity (from table I)
        self.Us_fe_avg=Us_fe_avg
        #tell the program the nonsteady waves corrected iron shock velocity (from table I)
        self.Us_fe=Us_fe
        #tell the program the F value for nonsteady waves (from table I)
        self.F_value=F_value
        #tell the program the pressure, density, and particle velocity for the iron (from table II)
        self.P=P
        self.rho=rho
        # Initial Densities in kg/m^3
        self.rho0 = {"qz": 2649, "fe": 7874}
        self.S1=1.46
        self.c0_1=4.54
        self.S2=1.26
        self.c0_2=5.82


    def calc_qz_snd_spd(self,rho_qz):
        #put soundspeedknown here
        #will have to use fit to cs-rho for qz high pressure or knudson desjarlais cannot use LANL Qz stuff
        #fit to rho-cs for LANL qz in relevant region: from Compare_LANLQz_MGLR_soundspeed.py
        return (4.393*rho_qz)-12.404
    def calc_qz_reshock_snd_spd(self,rho_qz):
        #put soundspeedknown here
        #will have to use fit to cs-rho for qz high pressure or knudson desjarlais cannot use LANL Qz stuff
        #fit to rho-cs for LANL qz in relevant region: from Compare_LANLQz_MGLR_soundspeed.py
        return (2.82*rho_qz)-3.45
    def calc_qz_reshock_pres(self,upusher):
        ureshock=0.62*upusher-3.84
        preshock=self.rho0['fe'] * self.Us_fe * ureshock * 10 ** (-3)
        rhoreshock=0.1795*upusher+4.3355
        return ureshock,preshock,rhoreshock
    
    def qzreshock(self):
        '''RESHOCK MODEL'''
        # Knudson/Desjarlais
        a = 5.477  # km/s
        b = 1.242
        c = 2.453
        d = 0.4336  # 1/km/s

        uparray=self.Us_Up_fit('qz')[0]
        u=self.qz_IM()[0]
        v1 = (self.Us_qz - u) / (self.rho0['qz'] * self.Us_qz)  #shocked specific volume at interface
        v0 = 1.0 / self.rho0['qz'] # 1/initial aluminum density, initial specific volume
        h = .00000001  # spacing

        v = np.linspace(float(v1), v0, round(float((v0 - v1) / h)))  # specific volume array
        gamma=0.64

        us = a + (b * uparray) - (c* uparray * (np.exp(-1 * d * uparray)))


        pressurearrayhug = us*uparray*self.rho0['qz']

        P2ray = self.Us_qz * self.rho0['qz'] * uparray
        i = np.argmin(np.abs(pressurearrayhug[2:] - P2ray[2:]))

        intersecindex = i
        upq = uparray[intersecindex]
        qpres = pressurearrayhug[intersecindex]
        varrayhug = (us - uparray) / (self.rho0['qz'] * us)
        qvol = varrayhug[intersecindex]

        # quartz density at interface
        qdensity = self.rho0['qz'] * self.Us_qz / (self.Us_qz- u)

        # Get the reshock line
        # The fraction in front of qvol specifies how far out to go from the qz volume!
        # if you change the ratio, and still don't get enough range, then add more to the max up in the uparray
        # difflist = varrayhug - ((4 / 5) * qvol)
        difflist = varrayhug - ((4/5) * qvol)

        minvindex = list(abs(difflist)).index(min(abs(difflist)))
        intrvarray = varrayhug[intersecindex:minvindex]

        intrpHugarray = pressurearrayhug[intersecindex:minvindex]

        # I changed the name ps from p2array and the name ups from upreshockarray so they would work with the if statements
        ps = (intrpHugarray - (gamma / intrvarray) * .5 * (intrpHugarray - qpres) * ((1 / self.rho0['qz']) - intrvarray)) / (
                1 - (gamma / intrvarray) * .5 * (qvol - intrvarray))
        ups = upq - (((ps - qpres) * (qvol - intrvarray)) ** .5)

        vzero = intrvarray[0]
        vone = intrvarray[len(intrvarray) - 1]
        hreshock = (vzero - vone) / len(ps)  # spacing
        vreshock = np.linspace(vone, vzero, len(ps))  # specific volume array

        return ps, v, ups,qdensity,hreshock, qpres, upq, qvol,h,vreshock,intrvarray
    
    def calc_sample_snd_spd(self):
        Usfit=self.Us_Up_fit('fe')[1]
        uparray=self.Us_Up_fit('qz')[0]
        Ppush=self.qz_IM()[2]
        upusher=self.qz_IM()[0]
        rhopush=self.qz_IM()[1]
        Pwitness=self.qz_IM()[5]
        uwitness=self.qz_IM()[3]
        rhowit=self.qz_IM()[4]
        upreshockedqzarray=self.qzreshock()[2]
        uqzreshock=self.calc_qz_reshock_pres(self.Us_qz)[0]
        reshockpres=self.calc_qz_reshock_pres(self.Us_qz)[1]
        reshockrho=self.calc_qz_reshock_pres(self.Us_qz)[2]
        #raylieghsam=self.rho0['fe']*uparray*self.Us_fe
        #rayleighsamnew = self.rho0['fe'] * self.Us_fe* (-upreshockedqzarray + upusher)
        #raylieghsamnewinterp=np.interp(upreshockedqzarray, (-upreshockedqzarray+ upusher), rayleighsamnew)
        #Preshockedqz=self.qzreshock()[0]
        #rhoreshockedqzarray=1/self.qzreshock()[9]
        #rhoreshockedqz=rhoreshockedqzarray[np.argmin(np.abs(Preshockedqz - raylieghsamnewinterp))]*10**(-3)
        #creshockedqz1=self.calc_qz_reshock_snd_spd(rhoreshockedqz)
        creshockedqz1=self.calc_qz_reshock_snd_spd(reshockrho)
        cwitness1=self.calc_qz_snd_spd(rhowit)
        cpush1=self.calc_qz_snd_spd(rhopush)
        Mpush=(reshockpres-Ppush)/(rhopush*cpush1*(upusher-uqzreshock))
        #Mreshock=(reshockpres-Ppush)/(rhoreshockedqz*creshockedqz1*(upusher-uqzreshock))
        Mreshock=(reshockpres-Ppush)/(reshockrho*creshockedqz1*(upusher-uqzreshock))
        #print("uqzreshock test",uqzreshock)
        #print("upusher test", upusher)
        #print("rhoreshockedqz test",rhoreshockedqz)
        #print("new reshocked rho",reshockrho)
        #print("creshockedqz1 test",creshockedqz1)
        #print("reshockpres test",reshockpres)
        #print("Ppush test",Ppush)
        #print("Mpush",Mpush)
        #print("Mreshock",Mreshock)
        Mwitness = (Pwitness) / (rhowit  * cwitness1 * uwitness)
        OpaqueF = self.F_value
        intermediatereshock=(1+Mreshock)/(1+Mpush)
        reshockTC=intermediatereshock
        #print("reshockTC",reshockTC)
        Mwitness = (Pwitness) / (rhowit  * cwitness1 * uwitness)
        Msample = 1 - (((1 - Mwitness) / OpaqueF) * reshockTC)
        #print("Msample",Msample)
        #want the usample used in the sound speed calculation to be from fit
        usample_avg=uparray[np.argwhere(self.Us_fe_avg<=Usfit)[0][0]]
        csample1 = (self.P) / (self.rho * usample_avg * Msample)
        return csample1,usample_avg, self.P, self.rho, Msample, Mwitness, Mreshock, Mpush, OpaqueF,reshockTC
    
    def calc_gruneisen_parameter(self):
        gammashot = (2/self.rho)*((self.calc_sample_snd_spd()[0] ** 2)*(self.rho**2) - self.dPdRho()*(self.rho**2)) / (( self.P ) -  self.dPdRho() * (self.rho**2) * (-(1 / self.rho) + (1 / (self.rho0['fe'] * 10 ** (-3)))))
        return gammashot
    
    def Us_Up_fit(self,mat):
        uparray=np.linspace(0.00, 49.9, 10000)
        #for quartz witness and pusher, use a linear fit to the LANL sjostrom and crockett data
        if mat=="qz":
            m=1.265
            b=4.828
            Usfit=(m*uparray)+b
            Pfit=self.rho0[mat]*Usfit*uparray
            part1=np.nan
            part2=np.nan
            breakpoint=np.nan
        #for iron, use the bilinear fit reported in the main text
        if mat=="fe":
            #identified breakpoint in particle velocity in main text
            breakpoint=np.argwhere(uparray>=6.0)[0][0]
            part1=(uparray[:breakpoint]*self.S1)+self.c0_1
            part2=(uparray[breakpoint:]*self.S2)+self.c0_2
            Usfit=np.concatenate((part1,part2))
            Pfit=self.rho0[mat]*Usfit*uparray
        return uparray,Usfit,Pfit,part1,part2,breakpoint
    
    def qz_IM(self):
        uparray=self.Us_Up_fit('qz')[0]
        Usfit=self.Us_Up_fit('qz')[1]
        upusher=uparray[np.argwhere(self.Us_qz<=Usfit)][0][0]
        rhopush = ((self.rho0['qz']*10**(-3)) / (1 - (upusher / self.Us_qz)))
        Ppush=(self.rho0['qz']*10**(-3))*self.Us_qz*upusher
        #make the witness density and particle velocity and pressure the same as the pusher
        uwitness=upusher
        rhowitness=rhopush
        Pwitness=Ppush

        return upusher,rhopush,Ppush,uwitness,rhowitness,Pwitness

    
    def dPdRho(self):
        bilinearvoluparray1=(self.Us_Up_fit('fe')[3] - self.Us_Up_fit('fe')[0][:self.Us_Up_fit('fe')[5]]) / (self.rho0['fe'] * 10 ** (-3) * self.Us_Up_fit('fe')[3])
        bilinearvoluparray2 = (self.Us_Up_fit('fe')[4] - self.Us_Up_fit('fe')[0][self.Us_Up_fit('fe')[5]:]) / (self.rho0['fe'] * 10 ** (-3) * self.Us_Up_fit('fe')[4])
        # bilinear dp/drho
        liquidgbilinear1 = np.polynomial.polynomial.Polynomial([0, -(self.c0_1 ** 2) * ((self.rho0['fe']*10**(-3)) ** 2), (self.c0_1 ** 2) * ((self.rho0['fe']*10**(-3)))])
        liquidhbilinear1 = np.polynomial.polynomial.Polynomial([(self.S1 ** 2) * ((self.rho0['fe']*10**(-3)) ** 2), 2 * self.S1 * (self.rho0['fe']*10**(-3)) - 2 * (self.S1 ** 2) * (self.rho0['fe']*10**(-3)),1 - 2 * self.S1 + (self.S1) ** 2])
        # do the quotient rule to take the derivative of this terrible expression (ratio of polynomials)
        liquidgprimebilinear1 = liquidgbilinear1.deriv(1)
        liquidhprimebilinear1 = liquidhbilinear1.deriv(1)
        liquidhsquaredbilinear1 = liquidhbilinear1 ** 2
        bilinear1dpdrhonum = liquidgprimebilinear1 * liquidhbilinear1 - liquidgbilinear1 * liquidhprimebilinear1
        bilinear1dpdrhodenom = liquidhsquaredbilinear1
        bilinear1dpdrho = bilinear1dpdrhonum(1 / bilinearvoluparray1) / bilinear1dpdrhodenom(1 / bilinearvoluparray1)

        # bilinear dp/drho
        liquidgbilinear2 = np.polynomial.polynomial.Polynomial([0, -(self.c0_2 ** 2) * ((self.rho0['fe']*10**(-3)) ** 2), (self.c0_2 ** 2) * (self.rho0['fe']*10**(-3))])
        liquidhbilinear2 = np.polynomial.polynomial.Polynomial([(self.S2 ** 2) * ((self.rho0['fe']*10**(-3)) ** 2), 2 * self.S2  * (self.rho0['fe']*10**(-3)) - 2 * (self.S2 ** 2) * (self.rho0['fe']*10**(-3)),1 - 2 * self.S2 + (self.S2) ** 2])
        # do the quotient rule to take the derivative of this terrible expression (ratio of polynomials)
        liquidgprimebilinear2 = liquidgbilinear2.deriv(1)
        liquidhprimebilinear2 = liquidhbilinear2.deriv(1)
        liquidhsquaredbilinear2 = liquidhbilinear2 ** 2
        bilinear2dpdrhonum = liquidgprimebilinear2 * liquidhbilinear2 - liquidgbilinear2 * liquidhprimebilinear2
        bilinear2dpdrhodenom = liquidhsquaredbilinear2
        bilinear2dpdrho = bilinear2dpdrhonum(1 / bilinearvoluparray2) / bilinear2dpdrhodenom(1 / bilinearvoluparray2)
        #identified breakpoint in particle velocity in main text
        if self.calc_sample_snd_spd()[1] <= 6.0:
            argslope = int(np.argwhere(1 / bilinearvoluparray1 >= self.rho)[0][0])
            slopeshot = bilinear1dpdrho[int(argslope)]
        if self.calc_sample_snd_spd()[1]>6.0:
            argslope = int(np.argwhere(1 / bilinearvoluparray2 > self.rho)[0][0])
            slopeshot = bilinear2dpdrho[int(argslope)]
        
        return slopeshot



shot_31381=experiment(34.65,24.5,26.20,1.1133,2870,20.0)
print("sample sound speed 31381",shot_31381.calc_sample_snd_spd())
print("sample gruneisen 31381",shot_31381.calc_gruneisen_parameter())


plt.show()