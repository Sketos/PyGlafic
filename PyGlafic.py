import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt
import pickle as pickle
from astropy.io import fits


class PyGlafic(object):

    def __init__(self):
        super(PyGlafic, self).__init__()

        # FILENAME
        self.FILENAME = "PyGlafic.INPUT"

        # INITIALIZE
        self.INPUTFILE_LINES = []

        self.INPUTFILE_LINES.append("\n")
        self.INPUTFILE_LINES.append("# PRI_PARAMETERS \n")
        self.INPUTFILE_LINES.append("\n")

        self.INPUTFILE_LINES.append("\n")
        self.INPUTFILE_LINES.append("# SEC_PARAMETERS \n")
        self.INPUTFILE_LINES.append("\n")

        self.INPUTFILE_LINES.append("\n")
        self.INPUTFILE_LINES.append("startup \n")
        self.INPUTFILE_LINES.append("\n")

        self.INPUTFILE_LINES.append("\n")
        self.INPUTFILE_LINES.append("# LEN_MODEL \n")
        self.INPUTFILE_LINES.append("\n")

        self.INPUTFILE_LINES.append("\n")
        self.INPUTFILE_LINES.append("# SRC_MODEL \n")
        self.INPUTFILE_LINES.append("\n")

        self.INPUTFILE_LINES.append("\n")
        self.INPUTFILE_LINES.append("end_startup \n")
        self.INPUTFILE_LINES.append("\n")

        # self.INPUTFILE_LINES.append("\n")
        # self.INPUTFILE_LINES.append("start_setopt \n")
        # self.INPUTFILE_LINES.append("\n")
        #
        # self.INPUTFILE_LINES.append("\n")
        # self.INPUTFILE_LINES.append("end_setopt \n")
        # self.INPUTFILE_LINES.append("\n")

        self.INPUTFILE_LINES.append("\n")
        self.INPUTFILE_LINES.append("start_command \n")
        self.INPUTFILE_LINES.append("\n")

        self.INPUTFILE_LINES.append("\n")
        self.INPUTFILE_LINES.append("quit \n")
        self.INPUTFILE_LINES.append("\n")

        # INITIALIZE MODELS
        self.LEN_MODELS = {}
        self.SRC_MODELS = {}


    def NAME(self):
        pass


    def COSMOLOGY(self, COSMOLOGICAL_PARAMETERS = None, \
        UPDATE = False, *args, **kwards):

        for i, line_i in enumerate(self.INPUTFILE_LINES):
            if line_i.split(" ")[0] == "#" and line_i.split(" ")[1] == "PRI_PARAMETERS":
                for j in range(i+1, len(self.INPUTFILE_LINES)):
                    if self.INPUTFILE_LINES[j].split(" ")[0] == "\n":
                        self.INPUTFILE_LINES.insert(j, "")
                        break

        if COSMOLOGICAL_PARAMETERS is None:
            COSMOLOGICAL_PARAMETERS = {"Omega_m" : 0.3, "Omega_l" : 0.7, \
                "w" : -1.00, "h" : 0.67 }
        else:
            pass # TODO :

        for key, value in COSMOLOGICAL_PARAMETERS.iteritems():
            if key == "Omega_m":
                self.INPUTFILE_LINES[j] += "omega" + "\t" + str(value) + "\n"
            if key == "Omega_l":
                self.INPUTFILE_LINES[j] += "lambda" + "\t" + str(value) + "\n"
            if key == "w":
                self.INPUTFILE_LINES[j] += "weos" + "\t" + str(value) + "\n"
            if key == "h":
                self.INPUTFILE_LINES[j] += "hubble" + "\t" + str(value) + "\n"


    def LENSING_CONFIGURATION(self, LEN_REDSHIFT = None, SRC_REDSHIFT = None, \
        UPDATE = False, *args, **kwards):

        self.LEN_REDSHIFT = LEN_REDSHIFT
        self.SRC_REDSHIFT = SRC_REDSHIFT

        if UPDATE:
            for i, line_i in enumerate(self.INPUTFILE_LINES):
                if line_i.split(" ")[0] == "#" and line_i.split(" ")[1] == "PRI_PARAMETERS":
                    for j in range(i+1, len(self.INPUTFILE_LINES)):
                        if self.INPUTFILE_LINES[j].split(" ")[0] == "zl":
                            self.INPUTFILE_LINES[j] = ""
                            self.INPUTFILE_LINES[j] += "zl" + " " + "\t" + \
                                str(self.LEN_REDSHIFT) + " " + "\n"
                            break
        else:
            for i, line_i in enumerate(self.INPUTFILE_LINES):
                if line_i.split(" ")[0] == "#" and line_i.split(" ")[1] == "PRI_PARAMETERS":
                    for j in range(i+1, len(self.INPUTFILE_LINES)):
                        if self.INPUTFILE_LINES[j].split(" ")[0] == "\n":
                            self.INPUTFILE_LINES.insert(j, "")
                            break

            self.INPUTFILE_LINES[j] += "zl" + " " + "\t" + \
                str(self.LEN_REDSHIFT) + " " + "\n"


    def GRID(self, ANGULAR_DISTANCE = None, LEN_PLANE_RESOLUTION = None, SRC_PLANE_RESOLUTION = None, \
        UPDATE = False, *args, **kwards):

        self.ANGULAR_DISTANCE = ANGULAR_DISTANCE

        for i, line_i in enumerate(self.INPUTFILE_LINES):
            if line_i.split(" ")[0] == "#" and line_i.split(" ")[1] == "PRI_PARAMETERS":
                for j in range(i+1, len(self.INPUTFILE_LINES)):
                    if self.INPUTFILE_LINES[j].split(" ")[0] == "\n":
                        self.INPUTFILE_LINES.insert(j, "")
                        break

        if ANGULAR_DISTANCE is None:
            raise ValueError
        else:
            self.INPUTFILE_LINES[j] += "xmin" + " " + "\t" + str(-ANGULAR_DISTANCE) + " " + "\n" + "ymin" + " " + "\t" + str(-ANGULAR_DISTANCE) + " " + "\n" + \
                "xmax" + " " + "\t" + str(+ANGULAR_DISTANCE) + " " + "\n" + "ymax" + " " + "\t" + str(+ANGULAR_DISTANCE) + " " + "\n"

        if SRC_PLANE_RESOLUTION is None:
            self.INPUTFILE_LINES[j] += "pix_poi" + " " + "\t" + str(1.0) + " " + "\n" + \
                "pix_ext" + " " + "\t" + str(1.0) + " " + "\n"
        else:
            self.INPUTFILE_LINES[j] += "pix_poi" + " " + "\t" + str(SRC_PLANE_RESOLUTION) + " " + "\n" + \
                "pix_ext" + " " + "\t" + str(SRC_PLANE_RESOLUTION) + " " + "\n"


    def LEN_MODEL(self, MODEL_TYPE = None, MODEL_PARAMETERS = None, \
        UPDATE = False, *args, **kwards):

        if UPDATE:
            for i, line_i in enumerate(self.INPUTFILE_LINES):
                if line_i.split(" ")[0] == "#" and line_i.split(" ")[1] == "LEN_MODEL":
                    for j in range(i+1, len(self.INPUTFILE_LINES)):
                        if self.INPUTFILE_LINES[j].split(" ")[0] == "lens" and self.INPUTFILE_LINES[j].split(" ")[2] == MODEL_TYPE:
                            self.INPUTFILE_LINES[j] = ""

                            self.INPUTFILE_LINES[j] += "lens" + " " + "\t" + \
                                " " + MODEL_TYPE + " " + "\t" + " "

                            for PARAMETER in MODEL_PARAMETERS:
                                self.INPUTFILE_LINES[j] += str(PARAMETER) + " " + "\t" + " "
                            self.INPUTFILE_LINES[j] += "\n"
                            break
        else:
            self.LEN_MODELS[str(len(self.LEN_MODELS) + 1)] = MODEL_TYPE
            for i, line_i in enumerate(self.INPUTFILE_LINES):
                if line_i.split(" ")[0] == "#" and line_i.split(" ")[1] == "LEN_MODEL":
                    for j in range(i+1, len(self.INPUTFILE_LINES)):
                        if self.INPUTFILE_LINES[j].split(" ")[0] == "\n":
                            self.INPUTFILE_LINES.insert(j, "")
                            break
            self.INPUTFILE_LINES[j] += "lens" + " " + "\t" + \
                " " + MODEL_TYPE + " " + "\t" + " "
            for PARAMETER in MODEL_PARAMETERS:
                self.INPUTFILE_LINES[j] += str(PARAMETER) + " " + "\t" + " "
            self.INPUTFILE_LINES[j] += "\n"


    def SRC_MODEL(self, MODEL_TYPE = None, MODEL_PARAMETERS = None, \
        UPDATE = False, *args, **kwards):

        self.SRC_MODELS[str(len(self.SRC_MODELS) + 1)] = MODEL_TYPE
        for i, line_i in enumerate(self.INPUTFILE_LINES):
            if line_i.split(" ")[0] == "#" and line_i.split(" ")[1] == "SRC_MODEL":
                for j in range(i+1, len(self.INPUTFILE_LINES)):
                    if self.INPUTFILE_LINES[j].split(" ")[0] == "\n":
                        self.INPUTFILE_LINES.insert(j, "")
                        break

        if MODEL_TYPE == "point":
            self.INPUTFILE_LINES[j] += MODEL_TYPE + " " + "\t" + " "
        else:
            self.INPUTFILE_LINES[j] += "extend" + " " + "\t" + " " + \
                MODEL_TYPE + " " + "\t" + " "

        if self.SRC_REDSHIFT:
            self.INPUTFILE_LINES[j] += str(self.SRC_REDSHIFT) + " " + "\t" + " "
        else:
            pass # TODO :

        for PARAMETER in MODEL_PARAMETERS:
            self.INPUTFILE_LINES[j] += str(PARAMETER) + " " + "\t" + " "
        self.INPUTFILE_LINES[j] += "\n"


    def MODELS(self):

        self.N_SRC_POI_MODELS = 0
        self.N_SRC_EXT_MODELS = 0

        for SRC_MODEL_KEY in self.SRC_MODELS.keys():
            if self.SRC_MODELS[SRC_MODEL_KEY] == "point":
                self.N_SRC_POI_MODELS += 1
            else:
                self.N_SRC_EXT_MODELS += 1

        for i, line_i in enumerate(self.INPUTFILE_LINES):
            if line_i.split(" ")[0] == "startup":
                break

        self.INPUTFILE_LINES[i] = "startup" + " " + "\t" + str(len(self.LEN_MODELS)) + " "  + \
            "\t" + str(self.N_SRC_EXT_MODELS) + " " + "\t" + str(self.N_SRC_POI_MODELS) + " "


    # def COMMANDS(self, COMMAND = None, UPDATE = None, COMMAND_ARGS = {}, *args, **kwards):
    #
    #     # NOTE : FIX
    #     if UPDATE:
    #         for i, line_i in enumerate(self.INPUTFILE_LINES):
    #             if line_i.split(" ")[0] == "start_command":
    #                 for j in range(i+1, len(self.INPUTFILE_LINES)):
    #                     if self.INPUTFILE_LINES[j].split(" ")[0] == COMMAND and COMMAND is not None:
    #                         self.INPUTFILE_LINES[j] = ""
    #                         if COMMAND == "writeimage":
    #                             self.INPUTFILE_LINES[j] += COMMAND + " " + "\t" + str(0.0) + "\t" + str(0.0) + " " + "\n"
    #                         break
    #     else:
    #         for i, line_i in enumerate(self.INPUTFILE_LINES):
    #             if line_i.split(" ")[0] == "start_command":
    #                 for j in range(i+1, len(self.INPUTFILE_LINES)):
    #                     if self.INPUTFILE_LINES[j].split(" ")[0] == "\n":
    #                         self.INPUTFILE_LINES.insert(j, "")
    #                         break
    #
    #         if COMMAND is None:
    #             pass # TODO :
    #         elif COMMAND == "writelens":
    #             self.INPUTFILE_LINES[j] += COMMAND + " " + "\t" + \
    #                 str(self.SRC_REDSHIFT) + " " + "\n"
    #         elif COMMAND == "writeimage":
    #             if "SKY" and "NOISE" in COMMAND_ARGS.keys():
    #                 self.INPUTFILE_LINES[j] += COMMAND + " " + "\t" + \
    #                     str(COMMAND_ARGS["SKY"]) + "\t" + str(COMMAND_ARGS["NOISE"]) + " " + "\n"
    #             else:
    #                 pass # TODO :


    def CREATE_PyGlafic_INPUTFILE(self):

        self.MODELS()

        with open("PyGlafic.INPUT", "w") as INPUTFILE:
            INPUTFILE.writelines(self.INPUTFILE_LINES)
        INPUTFILE.close()



if __name__ == "__main__":
    pass


    
