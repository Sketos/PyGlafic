import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from PyGlafic import PyGlafic


class DEFLECTION_ANGLE(PyGlafic):

    def __init__(self, COSMOLOGICAL_PARAMETERS = None, LEN_REDSHIFT = None, SRC_REDSHIFT = None, ANGULAR_DISTANCE = None, \
        LEN_PLANE_RESOLUTION = None, SRC_PLANE_RESOLUTION = None, MODEL_TYPE = None, MODEL_PARAMETERS = None):
        super(DEFLECTION_ANGLE, self).__init__()

        # COSMOLOGY
        self.COSMOLOGY(COSMOLOGICAL_PARAMETERS = COSMOLOGICAL_PARAMETERS, \
            UPDATE = False)

        if LEN_REDSHIFT >= SRC_REDSHIFT:
            raise ValueError
        else:
            self.LENSING_CONFIGURATION(LEN_REDSHIFT = LEN_REDSHIFT, \
                SRC_REDSHIFT = SRC_REDSHIFT)

        self.SRC_PLANE_RESOLUTION = SRC_PLANE_RESOLUTION

        if self.SRC_PLANE_RESOLUTION is None:
            raise ValueError
        else:
            self.GRID(ANGULAR_DISTANCE = ANGULAR_DISTANCE, \
                SRC_PLANE_RESOLUTION = SRC_PLANE_RESOLUTION)

        self.LEN_MODEL(MODEL_TYPE = MODEL_TYPE, MODEL_PARAMETERS = MODEL_PARAMETERS, \
            UPDATE = False)

        # COMMAND
        self.COMMAND()


    def COMMAND(self, COMMAND = "writelens", COMMAND_ARGS = None, UPDATE = False):

        if UPDATE:
            pass # TODO :
        else:
            for i, line_i in enumerate(self.INPUTFILE_LINES):
                if line_i.split(" ")[0] == "start_command":
                    for j in range(i+1, len(self.INPUTFILE_LINES)):
                        if self.INPUTFILE_LINES[j].split(" ")[0] == "\n":
                            self.INPUTFILE_LINES.insert(j, "")
                            break

            if COMMAND is None:
                raise ValueError
            elif COMMAND == "writelens":
                self.INPUTFILE_LINES[j] += COMMAND + \
                    " " + "\t" + str(self.SRC_REDSHIFT) + " " + "\n"
            else:
                raise ValueError

    def COMPUTE(self, BINSIZE = 5.0, DISPLAY = True):

        self.CREATE_PyGlafic_INPUTFILE()

        if os.path.isfile("./PyGlafic.INPUT"):
            os.system("./glafic PyGlafic.INPUT")
        else:
            IOError

        N_PIXELS = int(2.0 * self.ANGULAR_DISTANCE / \
            self.SRC_PLANE_RESOLUTION)

        LENSING_PROPERTIES_HDU = fits.open("out_lens.fits")

        # if CLEAN:
        #     os.system("rm out_lens.fits")

        LENSING_POTENTIAL_GRAD_X, LENSING_POTENTIAL_GRAD_Y = \
            np.gradient(LENSING_PROPERTIES_HDU["PRIMARY"].data[2])

        DEFLECTION_ANGLE_2D = np.sqrt(LENSING_POTENTIAL_GRAD_X**2.0 + \
            LENSING_POTENTIAL_GRAD_Y**2.0) / self.SRC_PLANE_RESOLUTION

        if DISPLAY:
            figure, axes = plt.subplots()

            IMAGE = plt.imshow(DEFLECTION_ANGLE_2D, cmap = "inferno", extent = (-self.ANGULAR_DISTANCE, self.ANGULAR_DISTANCE, \
                -self.ANGULAR_DISTANCE, self.ANGULAR_DISTANCE), vmin = 0.0, vmax = 10.0)

            plt.text(-0.85 * self.ANGULAR_DISTANCE, -0.85 * self.ANGULAR_DISTANCE, "$z_l = $" + str(self.LEN_REDSHIFT) + \
                "$,$" + "$z_s = $" + str(self.SRC_REDSHIFT), fontsize = 15)

            plt.xticks([-0.75 * self.ANGULAR_DISTANCE, -0.5 * self.ANGULAR_DISTANCE, -0.25 * self.ANGULAR_DISTANCE, 0.0, \
                0.25 * self.ANGULAR_DISTANCE, 0.5 * self.ANGULAR_DISTANCE, 0.75 * self.ANGULAR_DISTANCE])
            plt.yticks([-0.75 * self.ANGULAR_DISTANCE, -0.5 * self.ANGULAR_DISTANCE, -0.25 * self.ANGULAR_DISTANCE, 0.0, \
                0.25 * self.ANGULAR_DISTANCE, 0.5 * self.ANGULAR_DISTANCE, 0.75 * self.ANGULAR_DISTANCE])
            plt.xlabel(r"$x \,\, (arcsec)$", fontsize = 15)
            plt.ylabel(r"$y \,\, (arcsec)$", fontsize = 15)

            COLORBAR = figure.colorbar(IMAGE, ticks=[1.0, 3.0, 5.0, 7.0, 9.0], \
                cax = figure.add_axes([0.225, 0.90, 0.575, 0.05]), orientation = "horizontal")

            COLORBAR.set_label(r"$Deflection \,\, Angle$", \
                color = "black", fontsize = 15)


        N_PIXELS_Y, N_PIXELS_X = \
            DEFLECTION_ANGLE_2D.shape

        if (N_PIXELS_X != N_PIXELS) or (N_PIXELS_Y != N_PIXELS):
            raise ValueError

        Y_idx, X_idx = np.indices(DEFLECTION_ANGLE_2D.shape)
        CENTER = np.array([(X_idx.max() - X_idx.min()) / 2.0, \
            (Y_idx.max() - Y_idx.min()) / 2.0])
        R_2D = np.hypot(X_idx - CENTER[0], Y_idx - CENTER[1]) * \
            self.SRC_PLANE_RESOLUTION

        R_1D = np.linspace(R_2D.min(), R_2D.max(), \
            int(np.sqrt(N_PIXELS_X**2.0 + N_PIXELS_Y**2.0) / BINSIZE))

        DEFLECTION_ANGLE_1D = np.zeros((len(R_1D) - 1))
        for i in range(0, len(R_1D) - 1):
            idx = np.where((R_2D > R_1D[i]) & \
                (R_2D <= R_1D[i+1]))
            DEFLECTION_ANGLE_1D[i] = np.mean(DEFLECTION_ANGLE_2D[idx])

        R_1D = 0.5 * (R_1D[1:] + R_1D[:-1])

        # plt.figure()
        # plt.loglog(R_1D, DEFLECTION_ANGLE_1D, \
        #     marker = "o", linewidth = 4, alpha = 0.25)
        # print N_PIXELS * self.SRC_PLANE_RESOLUTION
        # plt.axvline(0.5 * N_PIXELS * self.SRC_PLANE_RESOLUTION)

        if DISPLAY:
            plt.show()


if __name__ == "__main__":

    # DEFLECTION_ANGLE_obj = DEFLECTION_ANGLE(LEN_REDSHIFT = 0.5, SRC_REDSHIFT = 2.0, ANGULAR_DISTANCE = 5.0, SRC_PLANE_RESOLUTION = 0.01, \
    #     MODEL_TYPE = "sie", MODEL_PARAMETERS = [300.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    DEFLECTION_ANGLE_obj = DEFLECTION_ANGLE(LEN_REDSHIFT = 0.5, SRC_REDSHIFT = 2.0, ANGULAR_DISTANCE = 100.0, SRC_PLANE_RESOLUTION = 1.0, \
        MODEL_TYPE = "nfw", MODEL_PARAMETERS = [10**14.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0])

    DEFLECTION_ANGLE_obj.COMPUTE()



