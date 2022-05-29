# -*- coding: utf-8 -*-
"""
Created on Thu May 26 23:28:29 2022

@author: USER
"""

from __future__ import unicode_literals # obsluga polskich znaków diaktrtycznych
import sys
from PyQt5.QtWidgets import QDialog, QApplication
from program import * # import kodu pythona ze schematem GUI
from math import sin, cos, sqrt, atan, atan2, degrees, radians
import math
import numpy as np
import PyQt5.QtGui


class MyForm(QDialog): # QDialog jako klasa nadrzedna

    def __init__(self):

        super().__init__()
        self.ui =  Ui_program() # Nazwa klasy z pliku przekonwertowanego z UI
        self.ui.setupUi(self)
        self.setWindowIcon(QtGui.QIcon('icon.jpg'))
        
        self.ui.GRS80.clicked.connect(self.grs80)
        self.ui.WGS84.clicked.connect(self.wgs84)
        
        self.ui.s5.clicked.connect(self.s5)
        self.ui.s6.clicked.connect(self.s6)
        self.ui.s7.clicked.connect(self.s7)
        self.ui.s8.clicked.connect(self.s8)
        self.ui.sA.clicked.connect(self.sA)
        
        self.ui.Oblicz00.clicked.connect(self.u2000)
        self.ui.Oblicz92.clicked.connect(self.u1992)
        self.ui.ObliczXYZ.clicked.connect(self.filamh2XYZ)
        self.show()
        
    def wgs84(self):
        self.a = 6378137.0 # semimajor_axis
        self.b = 6356752.31424518 # semiminor_axis
        
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) 
        self.ecc2 = (2 * self.flat - self.flat ** 2) 

    def grs80(self):
        self.a = 6378137.0
        self.b = 6356752.31414036

        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) 
        self.ecc2 = (2 * self.flat - self.flat ** 2) 
       
    def s5(self):
        self.ss = 1 
        self.s = 5
         
    def s6(self):
        self.ss = 1  
        self.s = 6
    
    def s7(self):
        self.ss = 1  
        self.s = 7
         
    def s8(self):
        self.ss = 1  
        self.s = 8
        
    def sA(self):
        self.ss = 0
        
    def u2000(self):
        '''
        funkcja dokonuje transformacji współrzędnych do układu 2000.
        '''

        if len(self.ui.stFI.text())!=0:
            fi_s = int(self.ui.stFI.text())
        if len(self.ui.mFI.text())!=0:
            fi_m = int(self.ui.mFI.text())
        if len(self.ui.secFI.text())!=0:
            fi_sec = int(self.ui.secFI.text())
        if len(self.ui.stLAM.text())!=0:
            lam_s = int(self.ui.stLAM.text())
        if len(self.ui.mLAM.text())!=0:
            lam_m = int(self.ui.mLAM.text())
        if len(self.ui.secLAM.text())!=0:
            lam_sec = int(self.ui.secLAM.text())
            
        fi = fi_s + fi_m/60 + fi_sec/3600 
        fi = math.radians(fi)
        lam = lam_s + lam_m/60 + lam_sec/3600 
        
        e_2 = self.ecc2 / (1 - self.ecc2)
        m_0 = 0.999923
        N = self.a / (math.sqrt(1 - self.ecc2 * np.sin(fi) ** 2))
        t = np.tan(fi)
        n2 = e_2 * np.cos(lam) ** 2  
        if self.ss == 1:         
            if self.s == 5:
                lam_0 = 15    
            elif self.s == 6:
                lam_0 = 18
            elif self.s == 7:
                lam_0 = 21
            elif self.s == 8:
                lam_0 = 24
        elif self.ss == 0:
            if lam > 13.5 and lam < 16.5:
                self.s = 5
                lam_0 = 15
            elif lam > 16.5 and lam < 19.5:
                self.s = 6
                lam_0 = 18
            elif lam > 19.5 and lam < 22.5:
                self.s = 7
                lam_0 = 21
            elif lam > 22.5 and lam < 25.5:
                self.s = 8
                lam_0 = 24


        
        lam = math.radians(lam)
        lam_0 = math.radians(lam_0)
        l = lam - lam_0
        A_0 = 1 - (self.ecc2 / 4) - (3 * (self.ecc2 ** 2)) / 64 - (5 * (self.ecc2 ** 3)) / 256
        A_2 = 3 / 8 * (self.ecc2 + ((self.ecc2 ** 2) / 4) + ((15 * self.ecc2 ** 3) / 128))
        A_4 = 15 / 256 * (self.ecc2 ** 2 + (3 * (self.ecc2 ** 3)) / 4)
        A_6 = (35 * (self.ecc2 ** 3)) / 3072
        sigma = self.a * ((A_0 * fi) - (A_2 * np.sin(2 * fi)) + (A_4 * np.sin(4 * fi)) - (A_6 * np.sin(6 * fi)))
        x = sigma + ((l ** 2) / 2) * (N * np.sin(fi) * np.cos(fi)) * (
                    1 + ((l ** 2) / 12) * ((np.cos(fi)) ** 2) * (5 - t ** 2 + 9 * n2 + (4 * n2 ** 2)) + ((l ** 4) / 360) * (
                        (np.cos(fi)) ** 4) * (61 - (58 * (t ** 2)) + (t ** 4) + (270 * n2) - (330 * n2 * (t ** 2))))
        y = l * (N * np.cos(fi)) * (1 + ((((l ** 2) / 6) * (np.cos(fi)) ** 2) * (1 - (t ** 2) + n2)) + (
                    ((l ** 4) / (120)) * (np.cos(fi) ** 4)) * (
                                                5 - (18 * (t ** 2)) + (t ** 4) + (14 * n2) - (58 * n2 * (t ** 2))))
        x00 = round(x * m_0, 3)
        y00 = round(y * m_0 + (self.s * 1000000) + 500000, 3)
        
        self.ui.wynikX.setText('x00 =' + str(x00))
        self.ui.wynikY.setText('y00 =' + str(y00))
        self.ui.wynikZ.setText('')
        
    def u1992(self):

        if len(self.ui.stFI.text())!=0:
            fi_s = int(self.ui.stFI.text())
        if len(self.ui.mFI.text())!=0:
            fi_m = int(self.ui.mFI.text())
        if len(self.ui.secFI.text())!=0:
            fi_sec = int(self.ui.secFI.text())
        if len(self.ui.stLAM.text())!=0:
            lam_s = int(self.ui.stLAM.text())
        if len(self.ui.mLAM.text())!=0:
            lam_m = int(self.ui.mLAM.text())
        if len(self.ui.secLAM.text())!=0:
            lam_sec = int(self.ui.secLAM.text())
            
        fi = fi_s + fi_m/60 + fi_sec/3600 
        fi = math.radians(fi)
        lam = lam_s + lam_m/60 + lam_sec/3600  
        lam = math.radians(lam)
        
        m_0 = 0.9993
        N = self.a/(math.sqrt(1-self.ecc2 * np.sin(fi)**2))
        t = np.tan(fi)
        n2 = self.ecc2 * np.cos(lam)**2
        lam_0 = math.radians(19)
        l = lam - lam_0
        
        A_0 = 1 - (self.ecc2/4) - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
        A_2 = 3/8 * (self.ecc2 + ((self.ecc2**2)/4) + ((15*self.ecc2**3)/128))
        A_4 = 15/256 * (self.ecc2**2 + (3*(self.ecc2**3))/4)
        A_6 = (35*(self.ecc2**3))/3072
        
        sigma = self.a* ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
        
        x = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x92 = round(x * m_0 - 5300000, 3)
        y92 = round(y * m_0 + 500000, 3)   
        
        self.ui.wynikX.setText('x92 =' + str(x92))
        self.ui.wynikY.setText('y92 =' + str(y92))
        self.ui.wynikZ.setText('')
                
    def filamh2XYZ(self):
        
        if len(self.ui.stFI.text())!=0:
            fi_s = int(self.ui.stFI.text())
        if len(self.ui.mFI.text())!=0:
            fi_m = int(self.ui.mFI.text())
        if len(self.ui.secFI.text())!=0:
            fi_sec = int(self.ui.secFI.text())
        if len(self.ui.stLAM.text())!=0:
            lam_s = int(self.ui.stLAM.text())
        if len(self.ui.mLAM.text())!=0:
            lam_m = int(self.ui.mLAM.text())
        if len(self.ui.secLAM.text())!=0:
            lam_sec = int(self.ui.secLAM.text())
            
        h = int(self.ui.h.text())
            
        fi = fi_s + fi_m/60 + fi_sec/3600 
        fi = math.radians(fi)
        lam = lam_s + lam_m/60 + lam_sec/3600  
        lam = math.radians(lam)
        
        N = self.a/math.sqrt(1-self.ecc2*math.sin(fi)**2)
        X = round((N + h) * math.cos(fi) * math.cos(lam), 3)
        Y = round((N + h) * math.cos(fi) * math.sin(lam), 3)
        Z = round((N*(1-self.ecc2) + h) * math.sin(fi), 3)
            
        self.ui.wynikX.setText('X =' + str(X))
        self.ui.wynikY.setText('Y =' + str(Y))
        self.ui.wynikZ.setText('Z =' + str(Z))
        
        
        
if __name__=="__main__":
    app = QApplication(sys.argv)
    w = MyForm()
    w.show()
    sys.exit(app.exec_())        
            