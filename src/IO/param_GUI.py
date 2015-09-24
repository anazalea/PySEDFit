# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 19:13:43 2015

@author: diego
"""

import sys, glob, shutil,os,subprocess
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import param5 as p

class Form(QMainWindow):
    """
    GUI class. Takes the parameter intialized by init_params from param5.py
    It follows the different options inside the parameter specification class
    in order to create the appropriate input format.

    """
    def __init__(self, parent=None):
        super(Form, self).__init__(parent)
        self.setWindowTitle('PySEDfit - Parameter selection')
        self.createUI()

    def createUI(self):
        self.fitsedParams,self.makesedParams = p.init_params()

        self.main_frame = QWidget()
        
        pbox = QVBoxLayout()
        self.labels = [] # Will hold parameter names for reference
        self.fields = [] # Will hold the actual values (default and modified) for each parameter
        fbsize = 350
        for i in range(len(self.fitsedParams)):
            param = self.fitsedParams[i]
            tempHbox = QHBoxLayout()
            templabel = QLabel(param.key)
            self.labels.append(templabel)
 
            if isinstance(param.allowedValues,list) and param.multipleKey==False and param.isList==False:
                temptbox = QComboBox()
                temptbox.setFixedWidth(fbsize-2)
                temptbox.addItems(list(param.allowedValues))
            
            elif param.multipleKey==True and param.isList==True:
                multibox = QHBoxLayout()
                
                nf = param.maxSize
                for j in range(nf):
                    if len(param.allowedValues)==nf:
                        if isinstance(param.allowedValues[j],list):
                            temptbox = QComboBox()
                            temptbox.setFixedWidth(fbsize/nf-5)
                            temptbox.addItems(list(param.allowedValues[j]))
                            #print "array",param.allowedValues[j]
                        else:
                            if not isinstance(param.dataType,list) and param.dataType==bool:
                                temptbox = QComboBox()
                                temptbox.setFixedWidth(fbsize-2)
                                temptbox.addItems(list([True,False]))

                            else:
                                temptbox = QLineEdit()
                                temptbox.setText(str(param.defaultValue[j]))
                                temptbox.setFixedWidth(fbsize/nf-5)
                                #print "single",param.allowedValues[j]
                                
                    else:
                        temptbox = QLineEdit()
                        temptbox.setText(str(param.defaultValue[j]))
                        temptbox.setFixedWidth(fbsize/4-5)
                        
                    multibox.addWidget(temptbox)
            else:
                if param.dataType != bool: 
                    temptbox = QLineEdit()
                    temptbox.setText(str(param.defaultValue))
                    temptbox.setFixedWidth(fbsize-2)
                elif param.dataType == bool:
                    temptbox = QComboBox()
                    temptbox.setFixedWidth(fbsize-2)
                    temptbox.addItems(["False","True"])
            
            self.fields.append(temptbox)
            
            if param.multipleKey==False and param.isList==False:
                tempHbox.addWidget(templabel)
                tempHbox.addStretch(1)
                tempHbox.addWidget(temptbox)
                pbox.addLayout(tempHbox)
            elif param.multipleKey==False and param.isList==True:
                tempHbox.addWidget(templabel)
                tempHbox.addStretch(1)
                tempHbox.addWidget(temptbox)
                pbox.addLayout(tempHbox)
            else:
                tempHbox.addWidget(templabel)
                tempHbox.addStretch(1)
                tempHbox.addLayout(multibox)
                pbox.addLayout(tempHbox)
        
        self.main_frame.setLayout(pbox)
        self.setCentralWidget(self.main_frame)        
        self.setMinimumSize(400, 50)
        print "Oh yeah!"
        


app = QApplication(sys.argv)
form = Form()
form.show()
app.exec_()