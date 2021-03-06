# This file is part of QuTiP.
#
#    QuTiP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    QuTiP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with QuTiP.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2013, Paul D. Nation & Robert J. Johansson
#
###########################################################################
import sys
import os

try:
    # python 2
    from urllib2 import urlopen
except:
    # python 3
    from urllib.request import urlopen

if os.environ['QUTIP_GUI'] == "PYSIDE":
    from PySide import QtGui, QtCore

elif os.environ['QUTIP_GUI'] == "PYQT4":
    from PyQt4 import QtGui, QtCore

import numpy
import scipy
import matplotlib
from qutip import version2int, __version__ as qutip_version

CD_BASE = os.path.dirname(__file__)


style_quit = """
QPushButton {
    font-family: Arial;
    border-width: 2px;
    border-color: #666666;
    border-style: solid;
    border-radius: 7;
    font-size: 16px;
    background-color: QLinearGradient(x1: 0, y1: 0, x2: 0, y2: 1,
                            stop: 0 #FFAAAA, stop: 0.1 #FF9999,
                            stop: 0.49 #FF8888, stop: 0.5 #FF7777,
                            stop: 1 #FF6666)
}
"""


class AboutBox(QtGui.QWidget):
    def __init__(self):
        from qutip import _version
        if _version.release:
            version = _version.short_version
        else:
            version = 'HEAD'
        QtGui.QWidget.__init__(self)

        self.setGeometry(0, 0, 360, 480)
        self.setWindowTitle("About QuTiP")
        self.setWindowIcon(QtGui.QIcon(CD_BASE + "/logo.png"))
        self.resize(360, 480)
        self.setMinimumSize(360, 480)
        self.setMaximumSize(360, 480)
        self.center()
        self.setFocus()

        logo = QtGui.QLabel(self)
        logo.setGeometry((self.width() - 250) / 2, 0, 250, 163)
        logo.setPixmap(QtGui.QPixmap(CD_BASE + "/logo.png"))

        tlabel = QtGui.QLabel(self)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setBold(True)
        if sys.platform == 'darwin':
            font.setPointSize(17)
        else:
            font.setPointSize(13)
        fm = QtGui.QFontMetrics(font)
        tstring = "QuTiP: The Quantum Toolbox in Python"
        pixelswide = fm.width(tstring)
        tlabel.setFont(font)
        tlabel.setText(tstring)
        tlabel.move((self.width() - pixelswide) / 2.0, 165)

        # set tab widget and tabs
        tab_widget = QtGui.QTabWidget(self)
        tab_widget.move(10, 200)
        tab_widget.resize(340, 220)

        tab1 = QtGui.QWidget(self)
        tab_widget.addTab(tab1, "Version Info")
        tab1_vert = QtGui.QVBoxLayout(tab1)
        tab2 = QtGui.QWidget()
        tab_widget.addTab(tab2, "Developers")
        t2_vert = QtGui.QVBoxLayout(tab2)

        # first tab text
        # call _set_strings function to get label and label2 widgets
        label, label2 = _set_strings()
        tab1_vert.addWidget(label)
        tab1_vert.addWidget(label2)

        # tab2 text
        t2_font = QtGui.QFont()
        t2_font.setFamily("Arial")
        t2_font.setBold(True)
        t2_font.setUnderline(True)
        if sys.platform == 'darwin':
            t2_font.setPointSize(15)
        else:
            t2_font.setPointSize(12)

        t2_font2 = QtGui.QFont()
        t2_font2.setFamily("Arial")
        t2_font2.setBold(False)
        t2_font2.setUnderline(False)
        if sys.platform == 'darwin':
            t2_font2.setPointSize(15)
        else:
            t2_font2.setPointSize(12)

        tab2_text_1 = QtGui.QLabel()
        dev_string = "Lead Developers:"
        tab2_text_1.setFont(t2_font)
        tab2_text_1.setText(dev_string)
        t2_vert.addWidget(tab2_text_1)

        tab2_text_2 = QtGui.QLabel()
        dev_string2 = (
            "<a href=http://dml.riken.jp/~rob>Robert Johansson</a>" +
            " & <a href=http://dml.riken.jp/~paul>Paul Nation</a>")
        tab2_text_2.setOpenExternalLinks(True)
        tab2_text_2.setFont(t2_font2)
        tab2_text_2.setText(dev_string2)
        t2_vert.addWidget(tab2_text_2)

        tab2_text_3 = QtGui.QLabel()
        contrib_string = "Contributors:"
        tab2_text_3.setFont(t2_font)
        tab2_text_3.setText(contrib_string)
        t2_vert.addWidget(tab2_text_3)

        tab2_text_4 = QtGui.QLabel()
        contrib_string2 = "Arne Grimsmo, Markus Baden\n"
        tab2_text_4.setFont(t2_font2)
        tab2_text_4.setText(contrib_string2)
        t2_vert.addWidget(tab2_text_4)

        tab2_text_5 = QtGui.QLabel()
        bug_string = "For a list of bug hunters and other"
        tab2_text_5.setFont(t2_font2)
        tab2_text_5.setText(bug_string)
        t2_vert.addWidget(tab2_text_5)

        tab2_text_6 = QtGui.QLabel()
        bug_string2 = ("supporters, see the " +
                       "<a href=http://qutip.googlecode.com/svn/doc/"
                       + version +
                       "/html/contributors.html>QuTiP documentation</a>.")
        tab2_text_6.setOpenExternalLinks(True)
        tab2_text_6.setFont(t2_font2)
        tab2_text_6.setText(bug_string2)
        t2_vert.addWidget(tab2_text_6)

        # copyright text
        p1_font = QtGui.QFont()
        p1_font.setFamily("Arial")
        p1_font.setBold(False)
        if sys.platform == 'darwin':
            p1_font.setPointSize(14)
        else:
            p1_font.setPointSize(11)
        alabel = QtGui.QLabel(self)
        astring = "Copyright (c) 2011-2013,\nP. D. Nation & J. R. Johansson"
        pixelswide = fm.width(astring)
        alabel.setFont(p1_font)
        alabel.setText(astring)
        alabel.move(10, 430)

        # QUIT BUTTON-----------------
        quit = HoverExit('Close', self)
        quit.setStyleSheet(style_quit)
        quit.setGeometry((self.width() - 80), 430, 70, 35)
        quit.clicked.connect(self.close)

    def center(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size = self.frameSize()
        self.move((screen.width() - size.width()) / 2,
                  (screen.height() - size.height()) / 2)


def _set_strings():
    t1_font = QtGui.QFont()
    t1_font.setFamily("Arial")
    t1_font.setBold(False)
    if sys.platform == 'darwin':
        t1_font.setPointSize(15)
    else:
        t1_font.setPointSize(12)
    # qutip text
    lstring = "QuTiP Version:           " + qutip_version
    label = QtGui.QLabel()
    label.setText(lstring)

    try:
        url = "http://qutip.googlecode.com/svn/doc/current_version.txt"
        current = urlopen(url).read()
        current = str(current, "utf-8")  # for python 3
    except:
        try:
            current = str(current)  # for python 2
        except:
            current = None

    if current and version2int(current) > version2int(qutip_version):
        label.setOpenExternalLinks(True)
        lstring += " (<a href=http://code.google.com/p/qutip/wiki/Download>Update</a>)"

    t1_font.setBold(True)
    label.setFont(t1_font)
    label.setText(lstring)
    # dependencies text
    label2 = QtGui.QLabel()
    t1_font.setBold(False)
    label2.setFont(t1_font)
    try:
        import PySide
        pyside_ver = PySide.__version__
    except:
        pyside_ver = 'None'
    try:
        import PyQt4.QtCore as qt4Core
        pyqt4_ver = qt4Core.PYQT_VERSION_STR
    except:
        pyqt4_ver = 'None'
    if sys.platform == 'darwin':
        try:
            import Foundation
            pyobjc = 'Yes'
        except:
            pyobjc = 'No'
    lstring2 = "NumPy Version:           " + str(numpy.__version__) + "\n"
    lstring2 += "SciPy Version:              " + str(scipy.__version__) + "\n"
    lstring2 += "Matplotlib Version:       " + str(
        matplotlib.__version__) + "\n\n"
    lstring2 += "PySide Version:           " + str(pyside_ver) + "\n"
    lstring2 += "PyQt4 Version:             " + str(pyqt4_ver) + "\n"
    if sys.platform == 'darwin':
        lstring2 += "PyObjc Installed:          " + str(pyobjc)
    label2.setText(lstring2)
    return label, label2


style_hover_enter = """
QPushButton {
    font-family: Arial;
    border-width: 3px;
    border-color:#111111;
    border-style: solid;
    border-radius: 7;
    font-size: 16px;
    background: QLinearGradient(x1: 0, y1: 0, x2: 0, y2: 1,
                    stop: 0 #FFAAAA, stop: 0.1 #FF9999,
                    stop: 0.49 #FF8888, stop: 0.5 #FF7777,
                    stop: 1 #FF6666)
}
"""

style_hover_exit = """
QPushButton {
    font-family: Arial;
    border-width: 2px;
    border-color: #666666;
    border-style: solid;
    border-radius: 7;
    font-size: 16px;
    background: QLinearGradient(x1: 0, y1: 0, x2: 0, y2: 1,
                        stop: 0 #FFAAAA, stop: 0.1 #FF9999,
                        stop: 0.49 #FF8888, stop: 0.5 #FF7777,
                        stop: 1 #FF6666)
}
"""


class HoverExit(QtGui.QPushButton):
    def enterEvent(self, event):
        self.setStyleSheet(style_hover_enter)

    def leaveEvent(self, event):
        self.setStyleSheet(style_hover_exit)
