#!/usr/bin/env python
#
# -*- coding: utf-8 -*-
#
#  SpillMain.py(w)
#  
#  Copyright 2022 Farrokh A. Ghavanini <ghavanini@gmail.com>
#  
#  This file is part of Persian Gulf Pollutant Tracker program.
#  The program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  Network Common Data Form (NetCDF) software from UCAR/Unidata was
#  used in the manipulation process of the data.
#  https://www.unidata.ucar.edu/software/netcdf/
#
#  Qt software was used in creating GUI through using PyQt bindings
#  http://www.riverbankcomputing.com/software/pyqt/

from PyQt5 import QtGui, QtCore, QtWidgets
from SpillUI import UI_MainWindow
import sys
import math
# import pickle
import os
import subprocess
import shutil
import SpillRes
import SpillAbout
import GnomeModule, TamocModule
import datetime as dt
import netCDF4

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()  # Python2
        # super().__init__()                # Python3
        self.ui = UI_MainWindow()
        self.ui.setupUi(self)

        # Display MainWindow at the center of the screen
        self.frame_geometry = self.frameGeometry()
        self.center_position = QtWidgets.QDesktopWidget().availableGeometry().center()
        self.frame_geometry.moveCenter(self.center_position)
        self.move(self.frame_geometry.topLeft())

        # Setting app icon from .rcc resource file
        app_icon = QtGui.QIcon(":/icons/logo.png")
        self.setWindowIcon(app_icon)
        self.setWindowTitle("Pollution Tracking")

        """
        Initiating variables
        """
        self.reset_and_init_vars()

        """
        UI initial condition
        """        
        self.ui.wind_mag_lnEdt.setEnabled(False)
        self.ui.wind_dir_lnEdt.setEnabled(False)
        self.ui.current_mag_lnEdt.setEnabled(False)
        self.ui.current_dir_lnEdt.setEnabled(False)
        self.ui.current_sal_lnEdt.setEnabled(False)
        self.ui.current_temp_lnEdt.setEnabled(False)
        self.ui.release_pollutant_lnEdt.setEnabled(False)
        self.ui.release_pollutant_file_Btn.setEnabled(False)
        self.ui.release_pollutant_comboBox.addItems(["Condensate", "MEG", "Natural Gas"])
        self.ui.release_amount_comboBox.addItems(["BBL", "kg"])

        """
        Buttons and menu actions
        """  
        self.ui.wind_file_radioButton.toggled.connect(self.wind_rdBtn_toggled)
        self.ui.current_file_radioButton.toggled.connect(self.current_rdBtn_toggled)
        # For the following four lines refer to https://stackoverflow.com/questions/35819538/using-lambda-expression-to-connect-slots-in-pyqt
        self.ui.wind_file_select_Btn.clicked.connect(lambda state : self.input_file_btn_clicked(self.ui.wind_file_select_lnEd))
        self.ui.current_file_select_Btn.clicked.connect(lambda state : self.input_file_btn_clicked(self.ui.current_file_select_lnEd))
        self.ui.release_output_Btn.clicked.connect(lambda state : self.output_folder_btn_clicked(self.ui.release_output_lnEdt))
        self.ui.release_pollutant_file_Btn.clicked.connect(lambda state : self.output_file_btn_clicked(self.ui.release_pollutant_lnEdt))
        self.ui.release_pollutant_comboBox.currentTextChanged.connect(self.release_pollutant_cbx_changed)
        self.ui.start_button.clicked.connect(self.start_btn_clicked)
        self.ui.close_button.clicked.connect(self.close_btn_clicked)
        self.ui.action_new.triggered.connect(self.new_session)
        self.ui.action_about.triggered.connect(self.about_dialog)
        self.ui.statusbar.setStyleSheet("QStatusBar{background-color: #f5f5f5;color:black}")
        self.ui.statusbar.showMessage("Ready")

        """
        Input validators
        """ 
        float_validator = QtGui.QRegExpValidator(QtCore.QRegExp('[-+]?([0-9]*\.[0-9]+|[0-9]+)'), self)
        positive_float_validator = QtGui.QRegExpValidator(QtCore.QRegExp('[+]?([0-9]*\.[0-9]+|[0-9]+)'), self)
        self.ui.wind_mag_lnEdt.setValidator(positive_float_validator)
        self.ui.wind_dir_lnEdt.setValidator(positive_float_validator)
        self.ui.current_mag_lnEdt.setValidator(positive_float_validator)
        self.ui.current_dir_lnEdt.setValidator(positive_float_validator)
        self.ui.current_sal_lnEdt.setValidator(positive_float_validator)
        self.ui.current_temp_lnEdt.setValidator(float_validator)
        self.ui.release_amount_lnEdt.setValidator(positive_float_validator)
        self.ui.release_latitude_lnEdt.setValidator(float_validator)
        self.ui.release_longitude_lnEdt.setValidator(float_validator)
        self.ui.release_depth_lnEdt.setValidator(positive_float_validator)

        """
        Initial values (For testing purposes)
        """
        self.ui.start_strtDateTime.setDateTime(dt.datetime(1992, 5, 10, 12, 0, 0, 0))
        self.ui.start_endDateTime.setDateTime(dt.datetime(1992, 5, 11, 12, 0, 0, 0))
        self.ui.wind_file_select_lnEd.setText(u'/media/iwf/Archive/Farrokh/Personal/Programming/Python_Projects/GUI/New_Project/withQFileDialog/wind_1992.nc')
        self.ui.wind_mag_lnEdt.setText(str(10))
        self.ui.wind_dir_lnEdt.setText(str(90))
        self.ui.current_mag_lnEdt.setText(str(1.2))
        self.ui.current_dir_lnEdt.setText(str(135))
        self.ui.current_sal_lnEdt.setText(str(30))
        self.ui.current_temp_lnEdt.setText(str(20))
        self.ui.current_file_select_lnEd.setText(u'/media/iwf/Archive/Farrokh/Personal/Programming/Python_Projects/GUI/New_Project/withQFileDialog/his_00005.nc')
        self.ui.release_output_lnEdt.setText(u'/media/iwf/Archive/Farrokh/Personal/Programming/Python_Projects/GUI/New_Project/withQFileDialog')
        self.ui.release_pollutant_lnEdt.setText(u'/media/iwf/Archive/Farrokh/Personal/Programming/Python_Projects/GUI/New_Project/withQFileDialog/aercoare.orig/input.csv')
        self.ui.release_start_DateTime.setDateTime(dt.datetime(1992, 5, 10, 12, 30, 0, 0))
        self.ui.release_end_DateTime.setDateTime(dt.datetime(1992, 5, 10, 16, 0, 0, 0))
        self.ui.release_longitude_lnEdt.setText(str(52.0))
        self.ui.release_latitude_lnEdt.setText(str(27.0))
        self.ui.release_depth_lnEdt.setText(str(5.0))
        self.ui.release_amount_lnEdt.setText(str(1000))

    def start_btn_clicked(self):
        """ 
        Starts the main process which creates a binary file containing all required variables to be exported
        (spill_vars_file) 
        """
        # while True:

        try:
            script_location = os.path.dirname(os.path.abspath(__file__))
            self.spill_vars_file = open(os.path.join(script_location, 'spill_vars.bin'), 'wb')
        except PermissionError:
            self.error_messageBox("Error!", "I/O Error!\n Permission Denied!")
            return
        
        try:
            self.model_strt = self.ui.start_strtDateTime.dateTime().toPyDateTime()
            self.model_end  = self.ui.start_endDateTime.dateTime().toPyDateTime()
            self.model_strt = self.model_strt.replace(second=0, microsecond=0)
            self.model_end  = self.model_end.replace(second=0, microsecond=0)
            if self.model_strt >= self.model_end:
                raise InvalidTimeRange
            else:
                self.model_timeRange = self.model_end - self.model_strt
                self.spill_vars["self.model_strt"] = self.model_strt
                self.spill_vars["self.model_end"] = self.model_end
                self.spill_vars["self.model_timeRange"] = self.model_timeRange
        except InvalidTimeRange:
                self.error_messageBox("Error!", "Model's end date is earlier than the start date!")
                return

        try:
            if self.ui.wind_file_radioButton.isChecked():
                if self.ui.wind_file_select_lnEd.text() == "":
                    raise EmptyPath
                else:
                    self.wind_file_path = self.ui.wind_file_select_lnEd.text()
                    self.spill_vars["self.wind_file_path"] = self.wind_file_path
            else:
                if (self.ui.wind_mag_lnEdt.text() == "" or
                    self.ui.wind_dir_lnEdt.text() == "" or
                    float(self.ui.wind_dir_lnEdt.text()) >= 360):
                    raise EmptyValue
                else:
                    self.wind_cnst_values = (float(self.ui.wind_mag_lnEdt.text()), float(self.ui.wind_dir_lnEdt.text()))
                    self.spill_vars["self.wind_cnst_values"] = self.wind_cnst_values                        
        except EmptyPath:                    
                self.error_messageBox("Error!", "Path to the input wind file is empty!")
                return
        except EmptyValue:                    
                self.error_messageBox("Error!", "Wind constant values are not set properly!")
                return

        try:
            if self.ui.current_file_radioButton.isChecked():
                if self.ui.current_file_select_lnEd.text() == "":
                    raise EmptyPath
                else:
                    self.current_file_path = self.ui.current_file_select_lnEd.text()
                    self.spill_vars["self.current_file_path"] = self.current_file_path
            else:
                if (self.ui.current_mag_lnEdt.text() == "" or
                    self.ui.current_dir_lnEdt.text() == "" or
                    self.ui.current_sal_lnEdt.text() == "" or
                    self.ui.current_temp_lnEdt.text() == "" or
                    float(self.ui.current_dir_lnEdt.text()) >= 360):
                    raise EmptyValue
                else:
                    self.current_u_comp = float(self.ui.current_mag_lnEdt.text()) * math.sin(math.radians(float(self.ui.current_dir_lnEdt.text())))
                    self.current_v_comp = float(self.ui.current_mag_lnEdt.text()) * math.cos(math.radians(float(self.ui.current_dir_lnEdt.text())))
                    self.current_components = (self.current_u_comp, self.current_v_comp)
                    self.current_temp = float(self.ui.current_temp_lnEdt.text())
                    self.current_sal = float(self.ui.current_sal_lnEdt.text())
                    self.spill_vars["self.current_cnst_values"] = self.current_components                    
        except EmptyPath:                    
                self.error_messageBox("Error!", "Path to the input current file is empty!")
                return
        except EmptyValue:                    
                self.error_messageBox("Error!", "Current constant values are not set properly!")
                return                                         
        
        try:
            if self.ui.release_output_lnEdt.text() == "":
                raise EmptyPath
            else:
                self.output_directory_path = self.ui.release_output_lnEdt.text()
                self.spill_vars["self.output_file_path"] = self.output_directory_path
        except EmptyPath:
                self.error_messageBox("Error!", "Output directory is not set!")
                return

        try:
            if self.ui.release_pollutant_comboBox.currentText() == "Natural Gas":
                if self.ui.release_pollutant_lnEdt.text() == "":
                    raise EmptyPath
                else:
                    self.aercoare_input_file = self.ui.release_pollutant_lnEdt.text()
                    self.spill_vars["aercoare_input_file"] = self.aercoare_input_file            
        except EmptyPath:
                self.error_messageBox("Error!", "AERCOARE input file is not set!")
                return
        
        try:
            if self.ui.release_amount_lnEdt.text() == "":
                raise EmptyValue
            else:
                self.pollutant_amount = (float(self.ui.release_amount_lnEdt.text()))
                self.pollutant_material = self.ui.release_pollutant_comboBox.currentText()
                if self.ui.release_amount_comboBox.currentText()  == "BBL":
                    self.pollutant_unit = "bbl"
                else:
                    self.pollutant_unit = "kg"
                self.spill_vars["self.pollutant_amount"] = self.pollutant_amount
                self.spill_vars["self.pollutant_material"] = self.pollutant_material
                self.spill_vars["self.pollutant_unit"] = self.pollutant_unit
        except EmptyValue:                    
                self.error_messageBox("Error!", "Pollutant amount is not set!")
                return

        try:
            self.release_strt = self.ui.release_start_DateTime.dateTime().toPyDateTime()
            self.release_end  = self.ui.release_end_DateTime.dateTime().toPyDateTime()
            self.release_strt = self.release_strt.replace(second=0, microsecond=0)
            self.release_end  = self.release_end.replace(second=0, microsecond=0)
            if self.release_strt > self.release_end:
                raise InvalidTimeRange
            else:
                if self.release_strt < self.model_strt or self.release_end > self.model_end:
                    raise InvalidRlsTimeRange
                else:
                    self.release_timeRange = self.release_end - self.release_strt
                    self.spill_vars["self.release_strt"] = self.release_strt
                    self.spill_vars["self.release_end"] = self.release_end
                    self.spill_vars["self.release_timeRange"] = self.release_timeRange
        except InvalidTimeRange:
                self.error_messageBox("Error!", "Release end date is earlier than the start date!")
                return
        except InvalidRlsTimeRange:
                self.error_messageBox("Error!", "Release dates are not in the range of model dates!")
                return

        try:
            if (self.ui.release_latitude_lnEdt.text() == "" or
                float(self.ui.release_latitude_lnEdt.text()) > 90 or
                float(self.ui.release_latitude_lnEdt.text()) < -90):
                raise EmptyValue
            else:
                self.release_latitude = float(self.ui.release_latitude_lnEdt.text())
        except EmptyValue:                    
                self.error_messageBox("Error!", "Release latitude value is not set properly!")
                return

        try:
            if (self.ui.release_longitude_lnEdt.text() == "" or
                float(self.ui.release_longitude_lnEdt.text()) > 180 or
                float(self.ui.release_longitude_lnEdt.text()) <= -180):
                raise EmptyValue
            else:
                self.release_longitude = float(self.ui.release_longitude_lnEdt.text())
        except EmptyValue:                    
                self.error_messageBox("Error!", "Release longitude value is not set properly!")
                return                          

        try:
            if self.ui.release_depth_lnEdt.text() == "":
                raise EmptyValue
            else:
                self.release_depth = float(self.ui.release_depth_lnEdt.text())
                self.release_location = (self.release_longitude, self.release_latitude, self.release_depth)
                self.spill_vars["self.release_location"] = self.release_location
        except EmptyValue:                    
                self.error_messageBox("Error!", "Release depth value is not set properly!")
                return

        # try:
        #     script_location = os.path.dirname(os.path.abspath(__file__))
        #     self.spill_vars_file = open(os.path.join(script_location, 'spill_vars.bin'), 'wb')         
        #     pickle.dump(self.spill_vars, self.spill_vars_file)
        #     self.spill_vars_file.close()

        #     # for testing pickle:
        #     # self.spill_vars_file = open(os.path.join(script_location, 'spill_vars.bin'), 'rb')
        #     # self.spill_vars = None
        #     # self.spill_vars = pickle.load(self.spill_vars_file)
        #     # print(len(self.spill_vars))
        #     # print(self.spill_vars)

        # except PermissionError:
        #     self.messageBox("Error!", "I/O Error!\n Permission Denied!")
        #     break                
        # break
        # self.ui.statusbar.setStyleSheet("QStatusBar{background-color:#f5f5f5;color:green}")
        # self.ui.statusbar.showMessage("Processing, please wait...")

        # try:
            # print(self.model_strt,
            #       self.model_timeRange,
            #       self.wind_file_path,
            #       self.wind_cnst_values,
            #       self.current_file_path,
            #       self.current_components,
            #       self.output_file_path,
            #       self.release_strt,
            #       self.release_end,
            #       self.release_location,
            #       self.pollutant_amount,
            #       self.pollutant_unit,
            #       self.pollutant_material)
        self.ui.statusbar.setStyleSheet("QStatusBar{background-color:#f5f5f5;color:green}")
        self.ui.statusbar.showMessage("Processing, please wait...")
        QtCore.QCoreApplication.processEvents()

        if self.ui.current_file_radioButton.isChecked() == True:
            current_dataset = netCDF4.Dataset(self.current_file_path, 'r')
            units = current_dataset.variables['ocean_time'].units
            calendar = current_dataset.variables['ocean_time'].calendar            
            self.timestep = (netCDF4.num2date(current_dataset.variables['ocean_time'][1], units=units, calendar=calendar) -
                                netCDF4.num2date(current_dataset.variables['ocean_time'][0], units=units, calendar=calendar))

            self.start_timestep_index = netCDF4.date2index(self.release_strt, current_dataset.variables['ocean_time'], calendar = calendar, select=u'nearest')
            self.end_timestep_index = netCDF4.date2index(self.release_end, current_dataset.variables['ocean_time'], calendar = calendar, select=u'nearest')

            time = self.release_strt
            for time_index in range(self.start_timestep_index, self.end_timestep_index + 1):
                (self.tamoc_x_final, self.tamoc_y_final) = TamocModule.process(method = 'file',
                                                                               current_file_path = self.current_file_path,
                                                                               release_location = self.release_location,
                                                                               release_timestep_index = time_index,
                                                                               release_timestep_datetime = time,
                                                                               release_time = (self.release_end - self.release_strt).total_seconds(),                                                                               
                                                                               release_pollutant = self.ui.release_pollutant_comboBox.currentText(),
                                                                               pollutant_amount = self.pollutant_amount,
                                                                               pollutant_unit = self.pollutant_unit)
                time = time + self.timestep

        else:
            (self.tamoc_x_final, self.tamoc_y_final) = TamocModule.process(method = 'constant',
                                                                           release_location = self.release_location,
                                                                           release_timestep_datetime = self.release_strt,
                                                                           release_time = (self.release_end - self.release_strt).total_seconds(),
                                                                           release_pollutant = self.ui.release_pollutant_comboBox.currentText(),
                                                                           pollutant_amount = self.pollutant_amount,
                                                                           pollutant_unit = self.pollutant_unit,
                                                                           current_components = self.current_components,
                                                                           current_sal = self.current_sal,
                                                                           current_temp = self.current_temp)

        if self.ui.release_pollutant_comboBox.currentText() == 'Condensate':
            GnomeModule.process(self.model_strt,
                                self.model_timeRange,
                                self.wind_file_path,
                                self.wind_cnst_values,
                                self.current_file_path,
                                self.current_components,
                                self.output_directory_path,
                                self.release_strt,
                                self.release_end,
                                self.release_location,
                                self.pollutant_amount,
                                self.pollutant_unit)

        elif self.ui.release_pollutant_comboBox.currentText() == 'Natural Gas': ### Needs some error handling ###
            shutil.rmtree('aercoare.temp', ignore_errors=True)
            shutil.copytree('./aercoare/','./aercoare.temp')
            shutil.copy2(self.aercoare_input_file, './aercoare.temp/')
            os.chdir('aercoare.temp')

            with open('aercoare.inp', 'w') as aercoare_inp:
                aercoare_inp.write("'" + os.path.split(self.aercoare_input_file)[1] + "'\n")
                aercoare_inp.write(str(self.release_latitude) + "\n")
                aercoare_inp.write(str(self.release_longitude))
            subprocess.call("./aercoare.exe")

            with open('utm.txt', 'r') as file:
                data = file.readline()
                lon_utm, lat_utm = data.split()[0], data.split()[1]

            os.chdir('../.')
            shutil.rmtree('aermod.temp', ignore_errors=True)
            shutil.copytree('./aermod/','./aermod.temp')
            shutil.copy2('./aercoare.temp/filepfl.pfl', './aermod.temp/')
            shutil.copy2('./aercoare.temp/filesfc.sfc', './aermod.temp/')
            os.chdir('aermod.temp')            

            with open(self.aercoare_input_file, 'r') as file:
                line_count = len(file.readlines())
                file.seek(0)
                data = file.readlines()
            
            start_line = data[1].split(',')
            end_line = data[line_count - 1].split(',')

            with open('aermod.inp', 'r') as file:
                data = file.readlines()
            if self.pollutant_unit == "kg":
                pol_discharge_aermod = (self.pollutant_amount * 1000) / (self.release_end - self.release_strt).total_seconds()
            else:
                pol_discharge_aermod = (self.pollutant_amount * 1000 * 0.09) / (self.release_end - self.release_strt).total_seconds()

            start_time_aermod = " ".join(start_line[:4])
            end_time_aermod = " ".join(end_line[:4])
            data[23] = "   LOCATION AREA1        AREA       " + str(lon_utm) + "    " + str(lat_utm) + "        2.000\n"
            data[25] = "   SRCPARAM AREA1            " + str(pol_discharge_aermod) + "     2.000    " + str(self.tamoc_x_final) + "    " +  str(self.tamoc_y_final) + "     0.000\n"
            data[38] = "                    XYINC " + str(float(lon_utm) - (25*500)) + " 200 50 " + str(float(lat_utm) - (25*500)) + " 200 50\n"
            data[13651] = "   SURFDATA " + str(self.release_strt.year) + " " + str(self.release_strt.year) + " Ali\n"
            data[13652] = "   UAIRDATA " + str(self.release_strt.year) + " " + str(self.release_strt.year) + " Ali\n"
            data[13654] = "   STARTEND " + start_time_aermod + " " + end_time_aermod + "\n"


            with open('aermod.inp', 'w') as file:
                file.writelines(data)

            subprocess.call("./aermod.exe")
            os.chdir('../.')            

        else:
            # Only TAMOC runs for 'MEG'
            pass

        # except Exception as e:
        #     self.error_messageBox("Error!", "There was a fatal error with spill module!")
        #     self.ui.statusbar.setStyleSheet("QStatusBar{background-color:#f5f5f5;color:black}")
        #     self.ui.statusbar.showMessage("Ready")
        #     QtCore.QCoreApplication.processEvents()
        #     return
        self.info_messageBox("Info!", "Process is done!")
        self.ui.statusbar.setStyleSheet("QStatusBar{background-color:#f5f5f5;color:black}")
        self.ui.statusbar.showMessage("Ready")
        QtCore.QCoreApplication.processEvents()
        return

    def error_messageBox(self, title, message):
        """ 
        Creates and shows a critical messagebox with the given title and message
        """        
        QtWidgets.QMessageBox.critical(self, title, message)

    def info_messageBox(self, title, message):
        """ 
        Creates and shows an information messagebox with the given title and message
        """        
        QtWidgets.QMessageBox.information(self, title, message)

    def wind_rdBtn_toggled(self):
        """ 
        Enable/Disable related widgets according to the state of the wind frame radiobutton
        """ 

        if self.ui.wind_file_radioButton.isChecked():
            self.ui.wind_file_select_lnEd.setEnabled(True)
            self.ui.wind_file_select_Btn.setEnabled(True)
            self.ui.wind_mag_lnEdt.setEnabled(False)
            self.ui.wind_dir_lnEdt.setEnabled(False)
        else:
            self.ui.wind_file_select_lnEd.setEnabled(False)
            self.ui.wind_file_select_Btn.setEnabled(False)
            self.ui.wind_mag_lnEdt.setEnabled(True)
            self.ui.wind_dir_lnEdt.setEnabled(True)            

    def current_rdBtn_toggled(self):
        """ 
        Enable/Disable related widgets according to the state of the current frame radiobutton
        """         
        if self.ui.current_file_radioButton.isChecked():
            self.ui.current_file_select_lnEd.setEnabled(True)
            self.ui.current_file_select_Btn.setEnabled(True)
            self.ui.current_mag_lnEdt.setEnabled(False)
            self.ui.current_dir_lnEdt.setEnabled(False)
            self.ui.current_sal_lnEdt.setEnabled(False)
            self.ui.current_temp_lnEdt.setEnabled(False)
        else:
            self.ui.current_file_select_lnEd.setEnabled(False)
            self.ui.current_file_select_Btn.setEnabled(False)
            self.ui.current_mag_lnEdt.setEnabled(True)
            self.ui.current_dir_lnEdt.setEnabled(True)
            self.ui.current_sal_lnEdt.setEnabled(True)
            self.ui.current_temp_lnEdt.setEnabled(True)            

    def input_file_btn_clicked(self, QLineEdit):
        """ 
        Opens up a file dialog to specify an input file for reading
        """
        QLineEdit.setText(QtWidgets.QFileDialog.getOpenFileName()[0])  

    def output_file_btn_clicked(self, QLineEdit):
        """ 
        Opens up a file dialog to specify an output file for writing
        """
        QLineEdit.setText(QtWidgets.QFileDialog.getOpenFileName()[0])

    def output_folder_btn_clicked(self, QLineEdit):
        """ 
        Opens up a file dialog to specify an output folder
        """
        QLineEdit.setText(QtWidgets.QFileDialog.getExistingDirectory(self, "Choose Directory", "/home"))

    def release_pollutant_cbx_changed(self):
        """ 
        Sets the state of corresponding lineEdit and button according to combobox
        """
        if self.ui.release_pollutant_comboBox.currentText() == 'Natural Gas':
            self.ui.release_pollutant_lnEdt.setEnabled(True)
            self.ui.release_pollutant_file_Btn.setEnabled(True)
        else:
            self.ui.release_pollutant_lnEdt.setEnabled(False)
            self.ui.release_pollutant_file_Btn.setEnabled(False)

    def close_btn_clicked(self):
        """ 
        Close the application window and exit
        """         
        self.close()

    def new_session(self):
        """ 
        Clear fields values and starts over
        """
        self.ui.wind_file_radioButton.setChecked(True)
        self.ui.wind_file_select_lnEd.setText("")
        self.ui.wind_mag_lnEdt.setText("")
        self.ui.wind_dir_lnEdt.setText("")
        self.ui.current_file_radioButton.setChecked(True)
        self.ui.current_file_select_lnEd.setText("")
        self.ui.current_mag_lnEdt.setText("")
        self.ui.current_dir_lnEdt.setText("")
        self.ui.current_sal_lnEdt.setText("")
        self.ui.current_temp_lnEdt.setText("")
        self.ui.release_output_lnEdt.setText("")
        self.ui.release_pollutant_lnEdt.setText("")
        self.ui.release_amount_lnEdt.setText("")
        self.ui.release_latitude_lnEdt.setText("")
        self.ui.release_longitude_lnEdt.setText("")
        self.ui.release_depth_lnEdt.setText("")        
        self.reset_and_init_vars()
        self.ui.statusbar.setStyleSheet("QStatusBar{background-color:#f5f5f5;color:black}")
        self.ui.statusbar.showMessage("Ready")
        QtCore.QCoreApplication.processEvents()      

    def about_dialog(self):
        """ 
        Shows the about dialog
        """
        self.dlg = SpillAbout.UiDialog()
        self.dlg.show()
    
    def reset_and_init_vars(self):
        """ 
        Resets and initializes variables
        """
        self.model_strt = None
        self.model_end = None
        self.model_timeRange = None
        self.wind_file_path = None
        self.wind_cnst_values = None
        self.current_file_path = None
        self.current_u_comp = None
        self.current_v_comp = None
        self.current_components = None
        self.current_temp = None
        self.current_sal = None
        self.output_directory_path = None
        self.pollutant_material = None
        self.aercoare_input_file = None
        self.pollutant_amount = None
        self.pollutant_unit = None
        self.release_strt = None
        self.release_end = None
        self.release_timeRange = None
        self.release_latitude = None
        self.release_longitude = None
        self.release_depth = None
        self.release_location = None
        self.spill_vars = {}
        self.current_dataset = None
        self.start_timestep_index = None
        self.end_timestep_index = None
        self.tamoc_x_final = None
        self.tamoc_y_final = None

# define Python user-defined exceptions
class Error(Exception):
    """Base class for other exceptions"""
    pass

class InvalidTimeRange(Error):
    """Raised when model date range is invalid"""
    pass

class InvalidRlsTimeRange(Error):
    """Raised when release date range is invalid"""
    pass

class EmptyPath(Error):
    """Raised when input/output path is invalid"""
    pass

class EmptyValue(Error):
    """Raised when a field has invalid value"""
    pass

if __name__ == '__main__':

    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.setFixedSize(760, 710)
    window.show()
    sys.exit(app.exec_())