# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'program.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_program(object):
    def setupUi(self, program):
        program.setObjectName("program")
        program.resize(489, 597)
        self.podajFI = QtWidgets.QLabel(program)
        self.podajFI.setGeometry(QtCore.QRect(70, 30, 121, 41))
        self.podajFI.setTextFormat(QtCore.Qt.RichText)
        self.podajFI.setObjectName("podajFI")
        self.podajLAM = QtWidgets.QLabel(program)
        self.podajLAM.setGeometry(QtCore.QRect(60, 90, 131, 41))
        self.podajLAM.setTextFormat(QtCore.Qt.RichText)
        self.podajLAM.setObjectName("podajLAM")
        self.label = QtWidgets.QLabel(program)
        self.label.setGeometry(QtCore.QRect(90, 250, 101, 20))
        self.label.setObjectName("label")
        self.s5 = QtWidgets.QRadioButton(program)
        self.s5.setGeometry(QtCore.QRect(320, 260, 82, 17))
        self.s5.setAcceptDrops(False)
        self.s5.setToolTipDuration(-1)
        self.s5.setAutoFillBackground(False)
        self.s5.setCheckable(True)
        self.s5.setChecked(False)
        self.s5.setAutoRepeat(False)
        self.s5.setAutoExclusive(True)
        self.s5.setAutoRepeatDelay(300)
        self.s5.setObjectName("s5")
        self.buttonGroup = QtWidgets.QButtonGroup(program)
        self.buttonGroup.setObjectName("buttonGroup")
        self.buttonGroup.addButton(self.s5)
        self.s6 = QtWidgets.QRadioButton(program)
        self.s6.setGeometry(QtCore.QRect(320, 290, 82, 17))
        self.s6.setAcceptDrops(False)
        self.s6.setToolTipDuration(-1)
        self.s6.setAutoFillBackground(False)
        self.s6.setChecked(False)
        self.s6.setAutoRepeat(False)
        self.s6.setAutoExclusive(True)
        self.s6.setAutoRepeatDelay(300)
        self.s6.setObjectName("s6")
        self.buttonGroup.addButton(self.s6)
        self.label_2 = QtWidgets.QLabel(program)
        self.label_2.setGeometry(QtCore.QRect(280, 230, 121, 20))
        self.label_2.setObjectName("label_2")
        self.s7 = QtWidgets.QRadioButton(program)
        self.s7.setGeometry(QtCore.QRect(320, 320, 82, 17))
        self.s7.setAcceptDrops(False)
        self.s7.setToolTipDuration(-1)
        self.s7.setAutoFillBackground(False)
        self.s7.setCheckable(True)
        self.s7.setChecked(False)
        self.s7.setAutoRepeat(False)
        self.s7.setAutoExclusive(True)
        self.s7.setAutoRepeatDelay(300)
        self.s7.setObjectName("s7")
        self.buttonGroup.addButton(self.s7)
        self.s8 = QtWidgets.QRadioButton(program)
        self.s8.setGeometry(QtCore.QRect(320, 350, 82, 17))
        self.s8.setAcceptDrops(False)
        self.s8.setToolTipDuration(-1)
        self.s8.setAutoFillBackground(False)
        self.s8.setChecked(False)
        self.s8.setAutoRepeat(False)
        self.s8.setAutoExclusive(True)
        self.s8.setAutoRepeatDelay(300)
        self.s8.setObjectName("s8")
        self.buttonGroup.addButton(self.s8)
        self.wynikX = QtWidgets.QLabel(program)
        self.wynikX.setGeometry(QtCore.QRect(60, 530, 111, 31))
        self.wynikX.setFrameShape(QtWidgets.QFrame.Box)
        self.wynikX.setFrameShadow(QtWidgets.QFrame.Plain)
        self.wynikX.setMidLineWidth(0)
        self.wynikX.setText("")
        self.wynikX.setWordWrap(False)
        self.wynikX.setObjectName("wynikX")
        self.label_4 = QtWidgets.QLabel(program)
        self.label_4.setGeometry(QtCore.QRect(60, 493, 47, 20))
        self.label_4.setObjectName("label_4")
        self.Oblicz00 = QtWidgets.QPushButton(program)
        self.Oblicz00.setGeometry(QtCore.QRect(110, 450, 75, 23))
        self.Oblicz00.setObjectName("Oblicz00")
        self.GRS80 = QtWidgets.QCheckBox(program)
        self.GRS80.setGeometry(QtCore.QRect(100, 290, 70, 17))
        self.GRS80.setCheckable(True)
        self.GRS80.setAutoExclusive(True)
        self.GRS80.setTristate(False)
        self.GRS80.setObjectName("GRS80")
        self.buttonGroup_2 = QtWidgets.QButtonGroup(program)
        self.buttonGroup_2.setObjectName("buttonGroup_2")
        self.buttonGroup_2.addButton(self.GRS80)
        self.WGS84 = QtWidgets.QCheckBox(program)
        self.WGS84.setGeometry(QtCore.QRect(100, 320, 70, 17))
        self.WGS84.setCheckable(True)
        self.WGS84.setAutoExclusive(True)
        self.WGS84.setTristate(False)
        self.WGS84.setObjectName("WGS84")
        self.buttonGroup_2.addButton(self.WGS84)
        self.wynikY = QtWidgets.QLabel(program)
        self.wynikY.setGeometry(QtCore.QRect(190, 530, 111, 31))
        self.wynikY.setFrameShape(QtWidgets.QFrame.Box)
        self.wynikY.setFrameShadow(QtWidgets.QFrame.Plain)
        self.wynikY.setMidLineWidth(0)
        self.wynikY.setText("")
        self.wynikY.setWordWrap(False)
        self.wynikY.setObjectName("wynikY")
        self.stFI = QtWidgets.QLineEdit(program)
        self.stFI.setGeometry(QtCore.QRect(210, 30, 41, 41))
        self.stFI.setObjectName("stFI")
        self.mFI = QtWidgets.QLineEdit(program)
        self.mFI.setGeometry(QtCore.QRect(260, 30, 41, 41))
        self.mFI.setObjectName("mFI")
        self.secFI = QtWidgets.QLineEdit(program)
        self.secFI.setGeometry(QtCore.QRect(310, 30, 41, 41))
        self.secFI.setObjectName("secFI")
        self.stLAM = QtWidgets.QLineEdit(program)
        self.stLAM.setGeometry(QtCore.QRect(210, 90, 41, 41))
        self.stLAM.setObjectName("stLAM")
        self.secLAM = QtWidgets.QLineEdit(program)
        self.secLAM.setGeometry(QtCore.QRect(310, 90, 41, 41))
        self.secLAM.setObjectName("secLAM")
        self.mLAM = QtWidgets.QLineEdit(program)
        self.mLAM.setGeometry(QtCore.QRect(260, 90, 41, 41))
        self.mLAM.setObjectName("mLAM")
        self.Oblicz92 = QtWidgets.QPushButton(program)
        self.Oblicz92.setGeometry(QtCore.QRect(200, 450, 75, 23))
        self.Oblicz92.setObjectName("Oblicz92")
        self.h = QtWidgets.QLineEdit(program)
        self.h.setGeometry(QtCore.QRect(210, 150, 41, 41))
        self.h.setObjectName("h")
        self.podajH = QtWidgets.QLabel(program)
        self.podajH.setGeometry(QtCore.QRect(80, 150, 121, 41))
        self.podajH.setTextFormat(QtCore.Qt.RichText)
        self.podajH.setObjectName("podajH")
        self.wynikZ = QtWidgets.QLabel(program)
        self.wynikZ.setGeometry(QtCore.QRect(320, 530, 111, 31))
        self.wynikZ.setFrameShape(QtWidgets.QFrame.Box)
        self.wynikZ.setFrameShadow(QtWidgets.QFrame.Plain)
        self.wynikZ.setMidLineWidth(0)
        self.wynikZ.setText("")
        self.wynikZ.setWordWrap(False)
        self.wynikZ.setObjectName("wynikZ")
        self.ObliczXYZ = QtWidgets.QPushButton(program)
        self.ObliczXYZ.setGeometry(QtCore.QRect(290, 450, 75, 23))
        self.ObliczXYZ.setObjectName("ObliczXYZ")
        self.sA = QtWidgets.QRadioButton(program)
        self.sA.setGeometry(QtCore.QRect(320, 380, 101, 31))
        self.sA.setAcceptDrops(False)
        self.sA.setToolTipDuration(-1)
        self.sA.setAutoFillBackground(False)
        self.sA.setChecked(False)
        self.sA.setAutoRepeat(False)
        self.sA.setAutoExclusive(True)
        self.sA.setAutoRepeatDelay(300)
        self.sA.setObjectName("sA")
        self.buttonGroup.addButton(self.sA)

        self.retranslateUi(program)
        QtCore.QMetaObject.connectSlotsByName(program)

    def retranslateUi(self, program):
        _translate = QtCore.QCoreApplication.translate
        program.setWindowTitle(_translate("program", "Aplikacja"))
        self.podajFI.setText(_translate("program", "<html><head/><body><p align=\"center\">Podaj wartość φ:<br/>(stopnie minuty sekundy)</p></body></html>"))
        self.podajLAM.setText(_translate("program", "<html><head/><body><p align=\"center\">Podaj wartość λ:<br/>(stopnie minuty sekundy)</p></body></html>"))
        self.label.setText(_translate("program", "Elipsoida odniesienia:"))
        self.s5.setText(_translate("program", "5"))
        self.s6.setText(_translate("program", "6"))
        self.label_2.setText(_translate("program", "Strefa dla układu 2000:"))
        self.s7.setText(_translate("program", "7"))
        self.s8.setText(_translate("program", "8"))
        self.label_4.setText(_translate("program", "Wynik:"))
        self.Oblicz00.setText(_translate("program", "2000"))
        self.GRS80.setText(_translate("program", "GRS80"))
        self.WGS84.setText(_translate("program", "WGS84"))
        self.Oblicz92.setText(_translate("program", "1992"))
        self.podajH.setText(_translate("program", "<html><head/><body><p>Podaj wysokość h[m]:</p></body></html>"))
        self.ObliczXYZ.setText(_translate("program", "XYZ"))
        self.sA.setText(_translate("program", "Automatyczne \n"
"dopasowanie"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    program = QtWidgets.QDialog()
    ui = Ui_program()
    ui.setupUi(program)
    program.show()
    sys.exit(app.exec_())
