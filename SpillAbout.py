#!/usr/bin/python
# -*- coding: utf-8 -*-

from PyQt5 import QtCore, QtWidgets
import SpillRes

class UiDialog(QtWidgets.QWidget):

    def __init__(self):
        super(UiDialog, self).__init__()    # Python2
        # super().__init__()                # Python3
        self.setupUi()

    def setupUi(self):
        self.setFixedSize(600, 250)
        self.setWindowTitle("About")

        self.textBrowser = QtWidgets.QTextBrowser(self)
        self.textBrowser.setOpenExternalLinks(True)
        self.textBrowser.setGeometry(QtCore.QRect(10, 10, 580, 230))

        self.textBrowser.setSource(QtCore.QUrl.fromLocalFile(":/about/about.html"))
        frameGm = self.frameGeometry()
        screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
        frameGm.moveCenter(centerPoint)
        self.move(frameGm.topLeft())