# Capsul import
from capsul.api import Process
from capsul.api import StudyConfig, get_process_instance

# Trait import
from traits.api import Float, TraitListObject, Undefined
from nipype.interfaces.base import OutputMultiPath, TraitedSpec, isdefined, InputMultiPath, File, Str, traits
from nipype.interfaces.spm.base import ImageFileSPM

# MIA2 import
from Project.Filter import Filter

# Other import
import os
from nipype.interfaces import spm
from skimage.transform import resize
import nibabel as nib
import numpy as np
import subprocess







def check_inputs(input_value):
    if input_value is not None:
        return input_value
    else:
        return {}


class Addition(Process):

    def __init__(self):
        super(Addition, self).__init__()

        self.add_trait("in_1", Float(output=False))
        self.add_trait("in_2", Float(output=False))
        self.add_trait("out", Float(output=True))

    def _run_process(self):
        self.out = self.in_1 + self.in_2
        print('Addition\n...\nInputs: {', self.in_1, ', ',
              self.in_2, '}\nOutput: ', self.out, '\n...\n')


class Substraction(Process):

    def __init__(self):
        super(Substraction, self).__init__()

        self.add_trait("in_1", Float(output=False))
        self.add_trait("in_2", Float(output=False))
        self.add_trait("out", Float(output=True))

    def _run_process(self):
        self.out = self.in_1 - self.in_2
        print('Substraction\n...\nInputs: {', self.in_1, ', ',
              self.in_2, '}\nOutput: ', self.out, '\n...\n')


class FSL_Smooth(Process):

    def __init__(self):
        super(FSL_Smooth, self).__init__()


        self.add_trait("in_file", File(output=False))
        self.add_trait("fwhm", Float(output=False))
        #self.add_trait(node_name + "_sigma", Float(output=False, optional=True))
        self.add_trait("out_file", File(output=True))

    def _run_process(self):

        import subprocess
        study_config = StudyConfig(modules=StudyConfig.default_modules + ['NipypeConfig'])

        # Process
        if study_config.use_nipype:
            #try:
            if study_config.use_nipype:
                smooth_process = get_process_instance("nipype.interfaces.fsl.Smooth")
                smooth_process.output_type = 'NIFTI'
                smooth_process.in_file = self.in_file

                # To resolve the sform/qform bug
                #subprocess.check_output(['fslorient', '-deleteorient', '1', self.in_file])
                #subprocess.check_output(['fslorient', '-setqformcode', '1', self.in_file])
                print("FWHM:", self.fwhm)

                if self.fwhm > 0:
                    smooth_process.sigma = self.fwhm / 2.355
                    #smooth_process.fwhm = self.fwhm
                elif self.sigma > 0:
                    smooth_process.sigma = self.sigma
                else:
                    smooth_process.fwhm = 2.0

                smooth_process.output_directory = os.path.split(self.out_file)[0]
                #smooth_process.output_directory = '/home/david/Nifti_data/'

            """except:
                smooth_process = None
                print('Smooth module of FSL is not present.')"""

        else:
            smooth_process = None
            print('NiPype is not present.')

        if smooth_process:
            study_config.reset_process_counter()
            study_config.run(smooth_process, verbose=1)

            # Display
            print('Smoothing with FSL\n...\nInputs: {', self.in_file, ', ',
                  self.fwhm, '}\nOutput: ', self.out_file, '\n...\n')


            #subprocess.check_output(['fslview', '/home/david/Nifti_data/1103/3/NIFTI/1103_3.nii'])
            #subprocess.check_output(['fslview', '/home/david/Nifti_data/1103_3_smooth.nii'])
            out_file = os.path.join(smooth_process.output_directory, os.path.basename(self.in_file)[:-4] + '_smooth.nii')

            subprocess.check_output(['fslview', self.in_file])
            subprocess.check_output(['fslview', out_file])


class Smooth(Process):

    in_file = File(optional=False, output=False)
    fwhm = Float(optional=False, output=False)
    out_file = File(optional=False, output=True)

    def _run_process(self):
        from capsul.subprocess import fsl
        #fsl.check_call(self.study_config, ['smooth', self.in_file, self.fwhm, self.out_file])
        #fsl.check_call(self.study_config, [self.in_file, "-s", self.fwhm, self.out_file])

        import subprocess
        subprocess.call(["fslmaths", self.in_file, "-s", str(self.fwhm), self.out_file])






class Populse_Filter(Process):

    def __init__(self, project, scans_list):
        super(Populse_Filter, self).__init__()

        self.project = project
        self.database = self.project.database
        if scans_list:
            self.scans_list = scans_list
        else:
            self.scans_list = self.project.database.get_documents_names()
        self.filter = Filter(None, [], [], [], [], [], "")

        self.add_trait("input", traits.List(traits.File, output=False))
        #self.add_trait("input", traits.Either(traits.List(File(exists=True)), File(exists=True), output=False))
        self.add_trait("output", traits.List(traits.File, output=True))

    def list_outputs(self):
        if self.input:
            self.scans_list = self.input
        else:
            self.scans_list = self.project.database.get_documents_names()

        filt = self.filter
        output = self.database.get_documents_matching_advanced_search(filt.links, filt.fields, filt.conditions,
                                                                      filt.values, filt.nots,
                                                                      self.scans_list)
        for idx, element in enumerate(output):
            full_path = os.path.join(self.project.folder, element)
            output[idx] = full_path

        return {'output': output}

    def _run_process(self):
        if self.input:
            self.scans_list = self.input
        else:
            self.scans_list = self.project.database.get_documents_names()
        self.output = []
        filt = self.filter
        # TODO: WHAT FUNCTION TO CALL
        output = self.database.get_documents_matching_advanced_search(filt.links, filt.fields, filt.conditions,
                                                                      filt.values, filt.nots,
                                                                      self.scans_list)

        for idx, element in enumerate(output):
            full_path = os.path.join(self.project.folder, element)
            output[idx] = full_path

        self.output = output



class Normalize_Spatial_Mask(Process):

    def __init__(self):
        super(Normalize_Spatial_Mask, self).__init__()

        # Inputs
        """self.add_trait("apply_to_files", InputMultiPath(traits.Either(
            ImageFileSPM(), traits.List(ImageFileSPM()), output=False)))"""
        self.add_trait("apply_to_files", traits.List(ImageFileSPM(), output=False))
        self.add_trait("deformation_file", ImageFileSPM(output=False))

        self.add_trait("jobtype", traits.String('write',
                                                usedefault=True, output=False, optional=True))
        # self.add_trait("write_bounding_box", traits.List(traits.List(traits.Float()), output=False, optional=True))
        self.add_trait("write_bounding_box", traits.List([[-78, -112, -50], [78, 76, 85]],
                                                         output=False, optional=True))
        # self.add_trait("write_voxel_sizes", traits.List(traits.Float(), output=False, optional=True))
        self.add_trait("write_voxel_sizes", traits.List([1, 1, 1], output=False, optional=True))
        # self.add_trait("write_interp", traits.Range(low=0, high=7, output=False, optional=True))
        self.add_trait("write_interp", traits.Int(1, output=False, optional=True))

        # Output
        self.add_trait("normalized_files", OutputMultiPath(File(), output=True))

    def _check_file_names(self):
        # Need to take only the three first ROIs
        files = []
        for file_name in self.apply_to_files:
            if type(file_name) in [list, TraitListObject]:
                file_name = file_name[0]
            if "c1" not in file_name:
                if "c2" not in file_name:
                    if "c3" not in file_name:
                        continue
                    else:
                        files.append(file_name)
                else:
                    files.append(file_name)
            else:
                files.append(file_name)

        return files

    def list_outputs(self):
        process = spm.Normalize12()
        if not self.apply_to_files:
            return {}
        else:
            process.inputs.apply_to_files = self._check_file_names()

        if not self.deformation_file:
            return {}
        else:
            process.inputs.deformation_file = self.deformation_file
        if not self.jobtype:
            return {}
        else:
            process.inputs.jobtype = self.jobtype

        outputs = process._list_outputs()

        return outputs

    def _run_process(self):

        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

        process = spm.Normalize12()
        process.inputs.apply_to_files = self._check_file_names()
        process.inputs.deformation_file = self.deformation_file
        process.inputs.jobtype = self.jobtype
        process.inputs.write_bounding_box = self.write_bounding_box
        process.inputs.write_voxel_sizes = self.write_voxel_sizes
        process.inputs.write_interp = self.write_interp

        process.run()


class Threshold(Process):

    def __init__(self):
        super(Threshold, self).__init__()

        # Inputs
        self.add_trait("in_files", InputMultiPath(traits.Either(
            ImageFileSPM(), traits.List(ImageFileSPM()), output=False)))
        self.add_trait("threshold", traits.Float(0.3, output=False, optional=True))
        self.add_trait("suffix", traits.String("_002", output=False, optional=True))
        self.add_trait("prefix", traits.String("", output=False, optional=True))

        # Output
        self.add_trait("out_files", OutputMultiPath(File(), output=True))

    def _check_file_names(self):
        # Need to take only the first ROI
        for file_name in self.in_files:
            if type(file_name) in [list, TraitListObject]:
                file_name = file_name[0]
            if "c1" not in file_name:
                continue
            else:
                return file_name

    def list_outputs(self):

        if not self.in_files:
            return {}

        if not self.suffix:
            self.suffix = ""

        if not self.prefix:
            self.prefix = ""

        file_name = self._check_file_names()
        path, file_name = os.path.split(file_name)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        out_file = os.path.join(path, self.prefix + file_name_no_ext + self.suffix + file_extension)

        d = {'out_files': out_file}
        return d

    def _run_process(self):

        if not self.suffix:
            self.suffix = ""

        if not self.prefix:
            self.prefix = ""

        file_name = self._check_file_names()
        if type(file_name) in [list, TraitListObject]:
            file_name = file_name[0]

        # Image processing
        img_final = threshold(file_name, self.threshold)

        # Image save
        path, file_name = os.path.split(file_name)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        out_file = os.path.join(path, self.prefix + file_name_no_ext + self.suffix + file_extension)
        nib.save(img_final, out_file)


class Resize(Process):

    def __init__(self):
        super(Resize, self).__init__()

        # Inputs
        self.add_trait("reference_image", InputMultiPath(traits.Either(
            ImageFileSPM(), traits.List(ImageFileSPM()), output=False)))
        self.add_trait("mask_to_resize", InputMultiPath(traits.Either(
            ImageFileSPM(), traits.List(ImageFileSPM()), output=False)))
        self.add_trait("suffix", traits.String("_003", output=False, optional=True))
        self.add_trait("prefix", traits.String(" ", output=False, optional=True))
        self.add_trait("interp", traits.Int(1, output=False, optional=True))

        # Output
        self.add_trait("out_file", ImageFileSPM(output=True))

    def list_outputs(self):

        if not self.mask_to_resize:
            return {}

        if not self.suffix or self.suffix in [Undefined, "<undefined>"]:
            self.suffix = " "

        if not self.prefix or self.prefix in [Undefined, "<undefined>"]:
            self.prefix = " "

        mask_name = self.mask_to_resize
        if type(mask_name) in [list, TraitListObject]:
            mask_name = mask_name[0]

        path, file_name = os.path.split(mask_name)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        if file_name_no_ext[-4:] == "_002":
            file_name_no_ext = file_name_no_ext[:-4]
        out_file = os.path.join(path, self.prefix.strip() + file_name_no_ext + self.suffix.strip() + file_extension)

        d = {'out_file': out_file}
        return d

    def _check_interp(self):
        if self.interp == 1:
            # Trilinear TODO: IS IT THE SAME AS BILINEAR?
            return "bilinear"
        elif self.interp == 0:
            # Nearest neighbour
            return "nearest"
        else:
            # TODO: WHAT ABOUT THE N-SINC INTERPOLATION?
            return None

    def _run_process(self):

        mask_name = self.mask_to_resize
        if type(mask_name) in [list, TraitListObject]:
            mask_name = mask_name[0]

        ref_name = self.reference_image
        if type(ref_name) in [list, TraitListObject]:
            ref_name = ref_name[0]

        if not self.suffix or self.suffix in [Undefined, "<undefined>"]:
            self.suffix = " "

        if not self.prefix or self.prefix in [Undefined, "<undefined>"]:
            self.prefix = " "

        # Image processing
        import nibabel as nib

        mask = nib.load(mask_name)
        mask_data = mask.get_data()

        ref = nib.load(ref_name)
        ref_data = ref.get_data()
        # Taking the first volume
        ref_size = ref_data.shape[:3]

        interp = self._check_interp()
        if not interp:
            raise ValueError("interp value of a Resize process has to be 0 (Nearest neighbour) or 1 (Trilinear).")

        # TODO: no info about the interp
        resized_mask = resize(mask_data, ref_size)
        # TODO: Taking info of mask's or ref's header?
        mask_final = nib.Nifti1Image(resized_mask, ref.affine, ref.header)

        # Image save
        path, file_name = os.path.split(mask_name)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        if file_name_no_ext[-4:] == "_002":
            file_name_no_ext = file_name_no_ext[:-4]
        out_file = os.path.join(path, self.prefix.strip() + file_name_no_ext + self.suffix.strip() + file_extension)
        nib.save(mask_final, out_file)




class ROI_List_Generator(Process):

    def __init__(self):
        super(ROI_List_Generator, self).__init__()

        list_pos = ['ROI_OCC', 'ROI_PAR', 'ROI_TEMP', 'ROI_INSULA', 'ROI_FRON', 'ROI_CING', 'ROI_THA',
                    'ROI_STR', 'ACP', 'ACA', 'ACM', 'PICA', 'SCA']

        # Inputs
        self.add_trait("pos", traits.List(list_pos, output=False))
        self.add_trait("hemi", traits.List(['_L', '_R'], output=False))

        # Output
        self.add_trait("roi_list", traits.List(output=True))

    def list_outputs(self):
        out_list = []
        for elm_pos in self.pos:
            for elm_hemi in self.hemi:
                out_list.append([elm_pos, elm_hemi])
        return {"roi_list": out_list}, {}

    def _run_process(self):
        out_list = []
        for elm_pos in self.pos:
            for elm_hemi in self.hemi:
                out_list.append([elm_pos, elm_hemi])

        self.roi_list = out_list


class Conv_ROI(Process):

    def __init__(self):
        super(Conv_ROI, self).__init__()

        # Inputs
        self.add_trait("roi_list", traits.List(output=False))
        self.add_trait("mask", File(output=False))

        # Outputs
        self.add_trait("out_masks", OutputMultiPath(output=True))

    def list_outputs(self):
        if not self.roi_list:
            return {}
        if not self.mask:
            return {}

        roi_dir = os.path.join(os.path.dirname(self.mask), 'roi')
        if not os.path.isdir(os.path.join(roi_dir, 'convROI_BOLD')):
            os.mkdir(os.path.join(roi_dir, 'convROI_BOLD'))

        conv_dir = os.path.join(roi_dir, 'convROI_BOLD')
        list_out = []

        for roi in self.roi_list:
            list_out.append(os.path.join(conv_dir, 'conv' + roi[0] + roi[1] + '.nii'))

        return {"out_masks": list_out}, {}

    def _run_process(self):

        roi_dir = os.path.join(os.path.dirname(self.mask), 'roi')
        conv_dir = os.path.join(roi_dir, 'convROI_BOLD')

        # Resizing the mask to the size of the ROIs
        roi_1 = self.roi_list[0]
        roi_file = os.path.join(roi_dir, roi_1[0] + roi_1[1] + '.nii')
        roi_img = nib.load(roi_file)
        roi_data = roi_img.get_data()
        roi_size = roi_data.shape[:3]

        mask_thresh = threshold(self.mask, 0.5).get_data()
        resized_mask = resize(mask_thresh, roi_size)

        for roi in self.roi_list:
            roi_file = os.path.join(roi_dir, roi[0] + roi[1] + '.nii')
            roi_img = nib.load(roi_file)
            roi_data = roi_img.get_data()
            mult = (roi_data * resized_mask).astype(float)
            # TODO: Should we take ino from ROI or from the mask ?
            mult_img = nib.Nifti1Image(mult, roi_img.affine, roi_img.header)

            # Image save
            out_file = os.path.join(conv_dir, 'conv' + roi[0] + roi[1] + '.nii')
            nib.save(mult_img, out_file)

            print('{0} saved'.format(os.path.basename(out_file)))


class Conv_ROI2(Process):

    def __init__(self):
        super(Conv_ROI2, self).__init__()

        # Inputs
        self.add_trait("roi_list", traits.List(output=False))
        self.add_trait("mask", File(output=False))

        # Outputs
        self.add_trait("out_masks2", OutputMultiPath(output=True))

    def list_outputs(self):
        if not self.roi_list:
            return {}
        if not self.mask:
            return {}

        roi_dir = os.path.join(os.path.dirname(self.mask), 'roi')
        if not os.path.isdir(os.path.join(roi_dir, 'convROI_BOLD2')):
            os.mkdir(os.path.join(roi_dir, 'convROI_BOLD2'))

        conv_dir = os.path.join(roi_dir, 'convROI_BOLD2')
        list_out = []

        for roi in self.roi_list:
            list_out.append(os.path.join(conv_dir, 'conv' + roi[0] + roi[1] + '2.nii'))

        return {"out_masks2": list_out}, {}

    def _run_process(self):

        roi_dir = os.path.join(os.path.dirname(self.mask), 'roi')
        conv_dir = os.path.join(roi_dir, 'convROI_BOLD')
        conv_dir2 = os.path.join(roi_dir, 'convROI_BOLD2')

        # Setting ROIs to the resolution of the functional
        mask = nib.load(self.mask).get_data()
        mask_size = mask.shape[:3]

        for roi in self.roi_list:
            roi_file = os.path.join(conv_dir, 'conv' + roi[0] + roi[1] + '.nii')
            roi_img = nib.load(roi_file)
            roi_data = roi_img.get_data()
            resized_roi = resize(roi_data, mask_size)
            resized_img = nib.Nifti1Image(resized_roi, roi_img.affine, roi_img.header)

            # Image save
            out_file = os.path.join(conv_dir2, 'conv' + roi[0] + roi[1] + '2.nii')
            nib.save(resized_img, out_file)

            print('{0} saved'.format(os.path.basename(out_file)))


class Write_results(Process):

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
                map_data = resize(map_data/map_data_max, roi_size) * map_data_max

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


class Grattefile(Process):

    def __init__(self):
        super(Grattefile, self).__init__()

        # Inputs
        self.add_trait("parametric_maps", traits.List(traits.File(exists=True), output=False))
        self.add_trait("data", traits.String("BOLD", output=False, optional=True))
        self.add_trait("calculs", traits.List(["mean", "std", "IL_mean", "IL_std"], output=False, optional=True))
        self.add_trait("mean_in_files", traits.List(traits.File(), output=False))
        self.add_trait("std_in_files", traits.List(traits.File(), output=False))
        self.add_trait("roi_list", traits.List(output=False))
        self.add_trait("patient_info", traits.Dict(output=False))

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

class BOLD_disp(Process):
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
        print(matlab_script)
        test = subprocess.run(['matlab', '-nodisplay', '-r', matlab_script])


class ANAT_disp(Process):

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
        print(matlab_script)
        test = subprocess.run(['matlab', '-nodisplay', '-r', matlab_script])


class Timecourse_fullTask(Process):

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


def threshold(file_name, thresh):
    img = nib.load(file_name)
    img_data = img.get_data()
    img_thresh = (img_data > thresh).astype(float)
    img_final = nib.Nifti1Image(img_thresh, img.affine, img.header)
    return img_final
