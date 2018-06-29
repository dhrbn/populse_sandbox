# Trait import
from nipype.interfaces.base import traits

# Other import
import os
from skimage.transform import resize
import nibabel as nib
import numpy as np
import subprocess

from PipelineManager.Process_mia import Process_mia


class Write_results(Process_mia):

    def __init__(self):
        super(Write_results, self).__init__()

        # Inputs
        self.add_trait("parametric_maps", traits.List(traits.File(exists=True), output=False))
        self.add_trait("roi_list", traits.List(output=False))

        # Outputs
        self.add_trait("mean_out_files", traits.List(traits.File(), output=True))
        self.add_trait("std_out_files", traits.List(traits.File(), output=True))

    def list_outputs(self):
        if not self.roi_list:
            return {}
        if not self.parametric_maps:
            return {}

        roi_dir = os.path.join(os.path.dirname(self.parametric_maps[0]), 'roi')
        conv_dir = os.path.join(roi_dir, 'convROI_BOLD')
        analysis_dir = os.path.join(roi_dir, 'ROI_analysis')

        if not os.path.isdir(conv_dir):
            print("No 'convROI_BOLD' folder in the working directory {0}.".
                  format(os.path.dirname(self.parametric_maps[0])))
            return {}

        if not os.path.isdir(analysis_dir):
            os.mkdir(analysis_dir)

        mean_out_files = []
        std_out_files = []

        for parametric_map in self.parametric_maps:
            for roi in self.roi_list:
                map_name = os.path.basename(parametric_map)[0:4] + '_BOLD'  # spmT_BOLD or beta_BOLD
                mean_out_files.append(os.path.join(analysis_dir, roi[0] + roi[1] + '_mean' + map_name + '.txt'))
                std_out_files.append(os.path.join(analysis_dir, roi[0] + roi[1] + '_std' + map_name + '.txt'))

        return {"mean_out_files": mean_out_files, "std_out_files": std_out_files}, {}

    def _run_process(self):

        roi_dir = os.path.join(os.path.dirname(self.parametric_maps[0]), 'roi')
        conv_dir = os.path.join(roi_dir, 'convROI_BOLD')
        analysis_dir = os.path.join(roi_dir, 'ROI_analysis')

        if not os.path.isdir(conv_dir):
            print("No 'convROI_BOLD' folder in the working directory {0}.".
                  format(os.path.dirname(self.parametric_maps[0])))
            return

        for parametric_map in self.parametric_maps:

            roi_1 = self.roi_list[0]
            roi_file = os.path.join(roi_dir, roi_1[0] + roi_1[1] + '.nii')
            roi_img = nib.load(roi_file)
            roi_data = roi_img.get_data()
            roi_size = roi_data.shape[:3]

            # Reading parametric map
            map_img = nib.load(parametric_map)
            map_data = map_img.get_data()

            # Setting the NaN to 0
            map_data = np.nan_to_num(map_data)
            map_name = os.path.basename(parametric_map)[0:4] + '_BOLD'  # spmT_BOLD or beta_BOLD

            # Making sure that the image are at the same size
            if roi_size != map_data.shape[:3]:
                map_data_max = max(map_data.max(), -map_data.min())
                map_data = resize(map_data / map_data_max, roi_size) * map_data_max

            for roi in self.roi_list:

                # Reading ROI file
                roi_file = os.path.join(conv_dir, 'conv' + roi[0] + roi[1] + '.nii')
                roi_img = nib.load(roi_file)
                roi_data = roi_img.get_data()

                # Making sure that the image are at the same size
                roi_data = np.nan_to_num(roi_data)

                # Convolution of the parametric map with the ROI
                roi_thresh = (roi_data > 0).astype(float)
                result = map_data * roi_thresh

                # Calculating mean and standard deviation
                if np.size(result[result.nonzero()]) == 0:
                    print("Warning: No data found in ROI {0}".format(roi[0] + roi[1]))
                    mean_result = 0
                    std_result = 0
                else:
                    mean_result = result[result.nonzero()].mean()
                    std_result = result[result.nonzero()].std()

                # Writing the value in the corresponding file

                mean_out_file = os.path.join(analysis_dir, roi[0] + roi[1] + '_mean' + map_name + '.txt')
                with open(mean_out_file, 'w') as f:
                    f.write("%.3f" % mean_result)

                std_out_file = os.path.join(analysis_dir, roi[0] + roi[1] + '_std' + map_name + '.txt')
                with open(std_out_file, 'w') as f:
                    f.write("%.3f" % std_result)

            print('{0} ROI files saved'.format(map_name))


class Grattefile(Process_mia):

    def __init__(self):
        super(Grattefile, self).__init__()

        patient_info_dict = {"name": "alej170316_testMIA24052018", "patho": "ACMD", "age": 64, "sex": "M",
                             "mr": "3T", "gas": "BACTAL", "admin": "MASK"}

        # Inputs
        self.add_trait("parametric_maps", traits.List(traits.File(exists=True), output=False))
        self.add_trait("data", traits.String("BOLD", output=False, optional=True))
        self.add_trait("calculs", traits.List(["mean", "std", "IL_mean", "IL_std"], output=False, optional=True))
        self.add_trait("mean_in_files", traits.List(traits.File(), output=False))
        self.add_trait("std_in_files", traits.List(traits.File(), output=False))
        self.add_trait("roi_list", traits.List(output=False))
        self.add_trait("patient_info", traits.Dict(patient_info_dict, output=False, optional=True))

        # Outputs
        self.add_trait("out_files", traits.List(traits.File(), output=True))

    def list_outputs(self):
        if not self.calculs:
            return {}
        if not self.parametric_maps:
            return {}

        roi_dir = os.path.join(os.path.dirname(self.parametric_maps[0]), 'roi')
        analysis_dir = os.path.join(roi_dir, 'ROI_analysis')

        if not os.path.isdir(analysis_dir):
            print("No 'ROI_analysis' folder in the working directory {0}.".
                  format(os.path.dirname(self.parametric_maps[0])))
            return {}

        out_files = []
        for parametric_map in self.parametric_maps:
            for calcul in self.calculs:
                out_files.append(os.path.join
                                 (analysis_dir,
                                  "{0}_{1}_{2}.xls".format(self.data, calcul,
                                                           os.path.basename(parametric_map)[0:9])))

        return {"out_files": out_files}, {}

    def _run_process(self):

        # Getting the list of all positions (roi_list without hemisphere)
        pos_list = []
        for roi in self.roi_list:
            pos = roi[0]
            if pos not in pos_list:
                pos_list.append(pos)

        roi_dir = os.path.join(os.path.dirname(self.parametric_maps[0]), 'roi')
        analysis_dir = os.path.join(roi_dir, 'ROI_analysis')

        if not os.path.isdir(analysis_dir):
            print("No 'ROI_analysis' folder in the working directory {0}.".
                  format(os.path.dirname(self.parametric_maps[0])))
            return {}

        for parametric_map in self.parametric_maps:
            map_name_file = os.path.basename(parametric_map)[0:9]
            map_name = map_name_file[0:4]
            for calcul in self.calculs:
                out_file = os.path.join(analysis_dir,
                                        "{0}_{1}_{2}.xls".format(self.data, calcul, map_name_file))

                with open(out_file, 'w') as f:
                    f.write("{0}\t".format('subjects'))
                    f.write("{0}\t".format('patho'))
                    f.write("{0}\t".format('age'))
                    f.write("{0}\t".format('sex'))
                    f.write("{0}\t".format('MR'))
                    f.write("{0}\t".format('Gaz'))
                    f.write("{0}\t".format('Admin'))

                    if calcul not in ['IL_mean', 'IL_std']:
                        for roi in self.roi_list:
                            f.write("{0}_{1}\t".format(map_name, roi[0] + roi[1]))
                    else:
                        for pos in pos_list:
                            f.write("{0}_{1}\t".format(map_name, pos))

                    # We should iterate on each patient here
                    f.write("\n{0}\t".format(self.patient_info["name"]))
                    if "patho" in self.patient_info.keys():
                        f.write("{0}\t".format(self.patient_info["patho"]))
                    else:
                        f.write("\t")
                    if "age" in self.patient_info.keys():
                        f.write("%3.1f\t" % self.patient_info["age"])
                    else:
                        f.write("\t")
                    if "sex" in self.patient_info.keys():
                        f.write("{0}\t".format(self.patient_info["sex"]))
                    else:
                        f.write("\t")
                    if "mr" in self.patient_info.keys():
                        f.write("{0}\t".format(self.patient_info["mr"]))
                    else:
                        f.write("\t")
                    if "gas" in self.patient_info.keys():
                        f.write("{0}\t".format(self.patient_info["gas"]))
                    else:
                        f.write("\t")
                    if "admin" in self.patient_info.keys():
                        f.write("{0}\t".format(self.patient_info["admin"]))
                    else:
                        f.write("\t")

                    if calcul == 'mean':
                        for roi in self.roi_list:
                            roi_file = os.path.join(analysis_dir,
                                                    "{0}_mean{1}_{2}.txt".format(roi[0] + roi[1], map_name, self.data))
                            with open(roi_file, 'r') as f_read:
                                final_res = float(f_read.read())
                            f.write("{0}\t".format(final_res))

                    elif calcul == 'IL_mean':
                        roi_checked = []
                        for roi in self.roi_list:
                            if roi[0] in roi_checked:
                                continue
                            roi_file = os.path.join(analysis_dir,
                                                    "{0}_mean{1}_{2}.txt".format(roi[0] + roi[1], map_name, self.data))
                            with open(roi_file, 'r') as f_read:
                                roi_value = float(f_read.read())

                            # Searching the roi that has the same first element
                            roi_2 = [s for s in self.roi_list if roi[0] in s[0] and roi[1] != s[1]][0]
                            roi_file_2 = os.path.join(analysis_dir,
                                                      "{0}_mean{1}_{2}.txt".format(roi_2[0] + roi_2[1],
                                                                                   map_name, self.data))
                            with open(roi_file_2, 'r') as f_read:
                                roi_value_2 = float(f_read.read())

                            if roi[1] == '_L':
                                sub_1 = roi_value
                                sub_2 = roi_value_2
                            else:
                                sub_1 = roi_value_2
                                sub_2 = roi_value

                            final_res = (sub_1 - sub_2) / (sub_1 + sub_2)
                            f.write("{0}\t".format(final_res))

                            roi_checked.append(roi[0])

                    elif calcul == "std":
                        for roi in self.roi_list:
                            roi_file = os.path.join(analysis_dir,
                                                    "{0}_std{1}_{2}.txt".format(roi[0] + roi[1], map_name, self.data))
                            with open(roi_file, 'r') as f_read:
                                final_res = float(f_read.read())
                            f.write("{0}\t".format(final_res))

                    elif calcul == 'IL_std':
                        roi_checked = []
                        for roi in self.roi_list:
                            if roi[0] in roi_checked:
                                continue
                            roi_file = os.path.join(analysis_dir,
                                                    "{0}_std{1}_{2}.txt".format(roi[0] + roi[1], map_name, self.data))
                            with open(roi_file, 'r') as f_read:
                                roi_value = float(f_read.read())

                            # Searching the roi that has the same first element
                            roi_2 = [s for s in self.roi_list if roi[0] in s[0] and roi[1] != s[1]][0]
                            roi_file_2 = os.path.join(analysis_dir,
                                                      "{0}_std{1}_{2}.txt".format(roi_2[0] + roi_2[1],
                                                                                  map_name, self.data))
                            with open(roi_file_2, 'r') as f_read:
                                roi_value_2 = float(f_read.read())

                            if roi[1] == '_L':
                                sub_1 = roi_value
                                sub_2 = roi_value_2
                            else:
                                sub_1 = roi_value_2
                                sub_2 = roi_value

                            final_res = (sub_1 - sub_2) / (sub_1 + sub_2)
                            f.write("{0}\t".format(final_res))

                            roi_checked.append(roi[0])


class BOLD_disp(Process_mia):
    """
    BOLD_disp(Patient, plane, tranche, tresh, native, tranche_native)
    Patient: Patient cell array.
    plane: axial, coronal or sagittal (mandatory parameter).
    tranche: First_plane:step:last_plane (mandatory parameter)
    tresh: Y (to change threshold in the parametric maps) or N (mandatory parameter)
    native: Y or N. If Y: generation of the 1st dynamic native image (mandatory parameter)
    tranche_native: If native = Y, tranche_native is used for the imagegeneration (only mandatory if native = Y)
    Example:
    BOLD_disp(Pat, 'axial', -50:5:65, Pat{1}.Thresholding, 'Y',-50:5:70)
    BOLD_disp(Pat,'coronal',-80:5:30, 'N', 'N')
    BOLD_disp(Pat,'sagittal',-60:2:-32, Pat{1}.Thresholding, 'N')
    """
    def __init__(self):
        super(BOLD_disp, self).__init__()

        # Inputs have to be .mat files
        self.add_trait("matlab_function", traits.File(output=False))
        self.add_trait("Patient", traits.File(output=False))
        self.add_trait("plane", traits.File(output=False))
        self.add_trait("tranche", traits.File(output=False))
        self.add_trait("tresh", traits.File(output=False))
        self.add_trait("native", traits.File(output=False))
        self.add_trait("tranche_native", traits.File(output=False))
        self.add_trait("output_directory", traits.Directory(output=False, optional=True))
        self.add_trait("dir_data", traits.Directory(output=False, optional=True))
        self.add_trait("dir_result", traits.Directory(output=False, optional=True))
        self.add_trait("dir_jpg", traits.Directory(output=False, optional=True))
        self.add_trait("todo", traits.File(output=False, optional=True))

    def list_outputs(self):
        pass

    def _run_process(self):
        verbose = False
        function_inputs = ["Patient", "plane", "tranche", "tresh", "native", "tranche_native",
                           "dir_data", "dir_result", "dir_jpg", "todo"]
        # Loading mat files
        matlab_script = ""
        for attribute in function_inputs:
            file_name = getattr(self, attribute)
            matlab_script += 'load("{0}","{1}");'.format(file_name, attribute)
            if attribute in ["dir_data", "dir_result", "dir_jpg", "todo"]:
                matlab_script += 'global {0};'.format(attribute)
            if verbose:
                matlab_script += 'disp("{0} loaded");disp({1});'.format(attribute, attribute)

        # Checking if there is an output directory
        if hasattr(self, "output_directory"):
            if self.output_directory:
                matlab_script += 'cd("{0}");'.format(self.output_directory)

        # Adding the function path
        head, tail = os.path.split(self.matlab_function)
        if verbose:
            matlab_script += 'disp(pwd);'
        matlab_script += 'addpath("{0}");'.format(head)

        # Adding the real path for display functions
        eric_path = '/home/david/Resultats_Pipeline_Eric/FICHIERS_FINAUX/IRMAGE_matlab_scripts/working_batchs/display/display_slices'
        matlab_script += 'addpath("{0}");'.format(eric_path)

        # Adding spm to Matlab path
        spm_path = "/home/david/code_matlab/spm12"
        matlab_script += 'addpath("{0}");'.format(spm_path)

        function_name = os.path.splitext(tail)[0]

        # Calling the function
        matlab_script += '{0}({1},{2},{3},{4},{5},{6});'.format(function_name,
                                                                "Patient",
                                                                "plane",
                                                                "tranche",
                                                                "tresh",
                                                                "native",
                                                                "tranche_native")

        # Exiting Matlab
        matlab_script += 'exit'

        # Running the function
        test = subprocess.run(['matlab', '-nodisplay', '-r', matlab_script])


class ANAT_disp(Process_mia):

    def __init__(self):
        super(ANAT_disp, self).__init__()

        # Inputs have to be .mat files
        self.add_trait("matlab_function", traits.File(output=False))
        self.add_trait("Patient", traits.File(output=False))
        self.add_trait("plane", traits.File(output=False))
        self.add_trait("tranche_native", traits.File(output=False))
        self.add_trait("output_directory", traits.Directory(output=False, optional=True))
        self.add_trait("dir_data", traits.Directory(output=False, optional=True))
        self.add_trait("dir_result", traits.Directory(output=False, optional=True))
        self.add_trait("dir_jpg", traits.Directory(output=False, optional=True))
        self.add_trait("todo", traits.File(output=False, optional=True))

    def list_outputs(self):
        pass

    def _run_process(self):
        verbose = False
        function_inputs = ["Patient", "plane", "tranche_native", "dir_data", "dir_result", "dir_jpg", "todo"]
        # Loading mat files
        matlab_script = ""
        for attribute in function_inputs:
            file_name = getattr(self, attribute)
            matlab_script += 'load("{0}","{1}");'.format(file_name, attribute)
            if attribute in ["dir_data", "dir_result", "dir_jpg", "todo"]:
                matlab_script += 'global {0};'.format(attribute)
            if verbose:
                matlab_script += 'disp("{0} loaded");disp({1});'.format(attribute, attribute)

        # Checking if there is an output directory
        if hasattr(self, "output_directory"):
            if self.output_directory:
                matlab_script += 'cd("{0}");'.format(self.output_directory)

        # Adding the function path
        head, tail = os.path.split(self.matlab_function)
        if verbose:
            matlab_script += 'disp(pwd);'
        matlab_script += 'addpath("{0}");'.format(head)

        # Adding the real path for display functions
        eric_path = '/home/david/Resultats_Pipeline_Eric/FICHIERS_FINAUX/IRMAGE_matlab_scripts/working_batchs/display/display_slices'
        matlab_script += 'addpath("{0}");'.format(eric_path)

        # Adding spm to Matlab path
        spm_path = "/home/david/code_matlab/spm12"
        matlab_script += 'addpath("{0}");'.format(spm_path)

        function_name = os.path.splitext(tail)[0]

        # Calling the function
        matlab_script += '{0}({1},{2},{3});'.format(function_name,
                                                    "Patient",
                                                    "plane",
                                                    "tranche_native")

        # Exiting Matlab
        matlab_script += 'exit'

        # Running the function
        test = subprocess.run(['matlab', '-nodisplay', '-r', matlab_script])


class Timecourse_fullTask(Process_mia):

    def __init__(self):
        super(Timecourse_fullTask, self).__init__()

        # Inputs have to be .mat files
        self.add_trait("matlab_function", traits.File(output=False))
        self.add_trait("Patient", traits.File(output=False))
        self.add_trait("output_directory", traits.Directory(output=False, optional=True))
        self.add_trait("dir_data", traits.Directory(output=False, optional=True))
        self.add_trait("dir_result", traits.Directory(output=False, optional=True))
        self.add_trait("dir_jpg", traits.Directory(output=False, optional=True))
        self.add_trait("todo", traits.File(output=False, optional=True))

    def list_outputs(self):
        pass

    def _run_process(self):
        verbose = True
        function_inputs = ["Patient", "dir_data", "dir_result", "dir_jpg", "todo"]
        # Loading mat files
        matlab_script = ""
        for attribute in function_inputs:
            file_name = getattr(self, attribute)
            matlab_script += 'load("{0}","{1}");'.format(file_name, attribute)
            if attribute in ["dir_data", "dir_result", "dir_jpg", "todo", "Patient"]:
                matlab_script += 'global {0};'.format(attribute)
            if verbose:
                matlab_script += 'disp("{0} loaded");disp({1});'.format(attribute, attribute)

        # Checking if there is an output directory
        if hasattr(self, "output_directory"):
            if self.output_directory:
                matlab_script += 'cd("{0}");'.format(self.output_directory)

        # Adding the function path
        head, tail = os.path.split(self.matlab_function)
        if verbose:
            matlab_script += 'disp(pwd);'
        matlab_script += 'addpath("{0}");'.format(head)

        # Adding the real path for display functions
        eric_path = '/home/david/Resultats_Pipeline_Eric/FICHIERS_FINAUX/IRMAGE_matlab_scripts/working_batchs/display/display_slices'
        matlab_script += 'addpath("{0}");'.format(eric_path)

        # Adding spm to Matlab path
        spm_path = "/home/david/code_matlab/spm12"
        matlab_script += 'addpath("{0}");'.format(spm_path)

        function_name = os.path.splitext(tail)[0]

        # Calling the function
        matlab_script += '{0};'.format(function_name)

        # Exiting Matlab
        matlab_script += 'exit'

        # Running the function
        test = subprocess.run(['matlab', '-nodisplay', '-r', matlab_script])
        # test = subprocess.run(['matlab', '-r', matlab_script])
