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


# To change to 'run_spm12.sh_location MCR_folder script"
from SoftwareProperties.Config import Config

config = Config()
spm_path = config.get_spm_path()
matlab_path = config.get_matlab_path()
matlab_cmd = '{0}/run_spm12.sh {1}/ script'.format(spm_path, matlab_path)

def refresh_matlab_command():
    """
    Refresh SPM and MatLab paths
    """

    global  matlab_cmd
    config = Config()
    spm_path = config.get_spm_path()
    matlab_path = config.get_matlab_path()
    matlab_cmd = '{0}/run_spm12.sh {1}/ script'.format(spm_path, matlab_path)

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

            #To resolve the sform/qform bug
            #subprocess.check_output(['fslorient', '-deleteorient', '1', out_file])
            #subprocess.check_output(['fslorient', '-setqformcode', '1', out_file])

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


class FSL_Smooth_Real(Process):


    def __init__(self):
        super(FSL_Smooth_Real, self).__init__()

        self.add_trait("in_file", File(output=False))
        self.add_trait("fwhm", Float(output=False))
        #self.add_trait(node_name + "_sigma", Float(output=False, optional=True))
        self.add_trait("out_file", File(output=True))

    def _run_process(self):

        import subprocess
        study_config = StudyConfig(
            modules=["FSLConfig"],
            use_fsl=True,
            fsl_config="/usr/share/fsl/5.0/etc/fslconf/fsl.sh",
            use_smart_caching=True,
            output_directory="/tmp/capsul_demo")
        #study_config = StudyConfig(modules=StudyConfig.default_modules + ['NipypeConfig'])

        # Process
        if study_config.use_fsl:
            smooth_process = get_process_instance(Smooth)
            #smooth_process.output_type = 'NIFTI'
            smooth_process.in_file = self.in_file
            smooth_process.out_file = self.out_file
            #smooth_process.output_directory = "/tmp/capsul_demo"
            smooth_process.output_type = 'NIFTI'

            # To resolve the sform/qform bug
            #subprocess.check_output(['fslorient', '-deleteorient', '1', self.in_file])
            #subprocess.check_output(['fslorient', '-setqformcode', '1', self.in_file])

            if self.fwhm > 0:
                smooth_process.fwhm = self.fwhm
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

            #subprocess.check_output(["fslchfiletype", "NIFTI", self.out_file + ".gz", self.out_file])

            import nibabel as nib
            img = nib.load(self.out_file + ".gz")
            nib.save(img, self.out_file)

            # Display
            print('Smoothing with FSL\n...\nInputs: {', self.in_file, ', ',
                  self.fwhm, '}\nOutput: ', self.out_file, '\n...\n')


            #subprocess.check_output(['fslview', '/home/david/Nifti_data/1103/3/NIFTI/1103_3.nii'])
            #subprocess.check_output(['fslview', '/home/david/Nifti_data/1103_3_smooth.nii'])
            #out_file = os.path.join(smooth_process.output_directory, os.path.basename(self.in_file)[:-4] + '_smooth.nii')
            #out_file = self.out_file

            #To resolve the sform/qform bug
            #subprocess.check_output(['fslorient', '-deleteorient', '1', out_file])
            #subprocess.check_output(['fslorient', '-setqformcode', '1', out_file])

            #subprocess.check_output(['fslview', self.in_file])
            #subprocess.check_output(['fslview', out_file])


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


class SPM_Smooth(Process):

    def __init__(self):
        super(SPM_Smooth, self).__init__()

        # Inputs
        self.add_trait("in_files", InputMultiPath(ImageFileSPM(), copyfile=False, output=False))
        self.add_trait("fwhm", traits.List([6, 6, 6], output=False, optional=True))
        #self.add_trait("fwhm", traits.List(traits.Float, output=False, optional=True))
        self.add_trait("data_type", traits.Int(output=False, optional=True))
        self.add_trait("implicit_masking", traits.Bool(output=False, optional=True))
        self.add_trait("out_prefix", traits.String('s', usedefault=True, output=False, optional=True))

        # Output
        self.add_trait("smoothed_files", OutputMultiPath(File(), output=True))

    def list_outputs(self):
        process = spm.Smooth()

        if not self.in_files:
            return {}
        else:
            process.inputs.in_files = self.in_files

        if self.out_prefix:
            process.inputs.out_prefix = self.out_prefix

        outputs = process._list_outputs()

        """raw_data_folder = os.path.join("data", "raw_data")
        derived_data_folder = os.path.join("data", "derived_data")
        for out_name, out_value in outputs.items():
            if type(out_value) is list:
                for idx, element in enumerate(out_value):
                    # To change the raw_data_folder to the derived_data_folder
                    element = element.replace(raw_data_folder, derived_data_folder)
                    outputs[out_name][idx] = element

        self.smoothed_files = outputs["smoothed_files"]"""

        inheritance_dict = {}
        for key, values in outputs.items():
            if key == "smoothed_files":
                for fullname in values:
                    path, filename = os.path.split(fullname)
                    if self.out_prefix:
                        filename_without_prefix = filename[len(self.out_prefix):]
                    else:
                        filename_without_prefix = filename[len('s'):]

                    if os.path.join(path, filename_without_prefix) in self.in_files:
                        inheritance_dict[fullname] = os.path.join(path, filename_without_prefix)

        return outputs, inheritance_dict

    def _run_process(self):

        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

        process = spm.Smooth()
        for idx, element in enumerate(self.in_files):
            full_path = os.path.relpath(element)
            self.in_files[idx] = full_path
        process.inputs.in_files = self.in_files
        process.inputs.fwhm = self.fwhm
        process.inputs.data_type = self.data_type
        process.inputs.implicit_masking = self.implicit_masking
        process.inputs.out_prefix = self.out_prefix

        process.run()


class SPM_NewSegment(Process):

    def __init__(self):
        super(SPM_NewSegment, self).__init__()

        # Inputs
        self.add_trait("affine_regularization", traits.Enum('mni', 'eastern', 'subj', 'none',
                                                            output=False, optional=True))
        self.add_trait("channel_files", InputMultiPath(output=False))
        """self.add_trait("channel_info", traits.Tuple(traits.Float(), traits.Float(),
                                                    traits.Tuple(traits.Bool, traits.Bool(True)),
                                                    output=False, optional=True))"""
        self.add_trait("channel_info", traits.Tuple((0.0001, 60, (False, True)), output=False, optional=True))
        self.add_trait("sampling_distance", Float(3.0, output=False, optional=True))
        self.add_trait("tissues", traits.List(traits.Tuple(
            traits.Tuple(ImageFileSPM(exists=True), traits.Int()),
            traits.Int(), traits.Tuple(traits.Bool, traits.Bool),
            traits.Tuple(traits.Bool, traits.Bool)), output=False, optional=True))
        """self.add_trait("warping_regularization", traits.Either(
            traits.List(traits.Float(), minlen=5, maxlen=5),
            traits.Float(), output=False))"""

        """self.add_trait("warping_regularization",
            traits.List(traits.Float(), output=False, optional=True))"""
        self.add_trait("warping_regularization",
            traits.List([0, 0.001, 0.5, 0.05, 0.2], output=False, optional=True))
        """self.add_trait("write_deformation_fields", traits.List(
            traits.Bool(), output=False, optional=True))"""
        self.add_trait("write_deformation_fields", traits.List(
            [False, True], output=False, optional=True))

        # Output
        self.add_trait("forward_deformation_field", File(output=True))
        self.add_trait("bias_field_images", File(output=True))
        self.add_trait("native_class_images", traits.List(traits.List(File()), output=True))

    def list_outputs(self):
        process = spm.NewSegment()

        if not self.channel_files:
            return {}
        else:
            process.inputs.channel_files = self.channel_files

        if not self.write_deformation_fields:
            return {}
        else:
            process.inputs.write_deformation_fields = self.write_deformation_fields

        if not self.channel_info:
            return {}
        else:
            process.inputs.channel_info = self.channel_info

        outputs = process._list_outputs()

        """raw_data_folder = os.path.join("data", "raw_data")
        derived_data_folder = os.path.join("data", "derived_data")
        for out_name in list(outputs):
            out_value = outputs[out_name]
            if out_name not in ["forward_deformation_field"]:
                del outputs[out_name]
            else:
                if type(out_value) is list:
                    for idx, element in enumerate(out_value):
                        if not element:
                            continue
                        if type(element) is list:
                            for idx_2, element_2 in enumerate(element):
                                element_2 = element_2.replace(raw_data_folder, derived_data_folder)
                                outputs[out_name][idx][idx_2] = element
                        else:

                            print("element: ", element)
                            element = element.replace(raw_data_folder, derived_data_folder)
                            outputs[out_name][idx] = element"""

        return outputs

    def _run_process(self):

        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

        process = spm.NewSegment()
        process.inputs.affine_regularization = self.affine_regularization
        process.inputs.channel_files = self.channel_files
        process.inputs.channel_info = self.channel_info
        process.inputs.sampling_distance = self.sampling_distance
        process.inputs.tissues = self.tissues
        process.inputs.warping_regularization = self.warping_regularization
        process.inputs.write_deformation_fields = self.write_deformation_fields

        process.run()


class SPM_Normalize(Process):

    def __init__(self):
        super(SPM_Normalize, self).__init__()

        # Inputs
        # self.add_trait("apply_to_files", InputMultiPath(traits.Either(
        #    ImageFileSPM(exists=True), traits.List(ImageFileSPM(exists=True)), output=False)))
        self.add_trait("apply_to_files", InputMultiPath(traits.Either(
            ImageFileSPM(), traits.List(ImageFileSPM()), output=False)))
        self.add_trait("deformation_file", ImageFileSPM(output=False))

        """self.add_trait("jobtype", traits.Enum('write', 'est', 'estwrite',
                                              usedefault=True, output=False, optional=True))"""
        self.add_trait("jobtype", traits.String('write',
                                                usedefault=True, output=False, optional=True))
        # self.add_trait("write_bounding_box", traits.List(traits.List(traits.Float()), output=False, optional=True))
        self.add_trait("write_bounding_box", traits.List(traits.List([[-78, -112, -50], [78, 76, 85]]),
                                                         output=False, optional=True))
        # self.add_trait("write_voxel_sizes", traits.List(traits.Float(), output=False, optional=True))
        self.add_trait("write_voxel_sizes", traits.List([1, 1, 1], output=False, optional=True))
        # self.add_trait("write_interp", traits.Range(low=0, high=7, output=False, optional=True))
        self.add_trait("write_interp", traits.Int(1, output=False, optional=True))

        # Output
        self.add_trait("normalized_files", OutputMultiPath(File(), output=True))

    def list_outputs(self):
        process = spm.Normalize12()
        if not self.apply_to_files:
            return {}
        else:
            process.inputs.apply_to_files = self.apply_to_files
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
        process.inputs.apply_to_files = self.apply_to_files
        process.inputs.deformation_file = self.deformation_file
        process.inputs.jobtype = self.jobtype
        process.inputs.write_bounding_box = self.write_bounding_box
        process.inputs.write_voxel_sizes = self.write_voxel_sizes
        process.inputs.write_interp = self.write_interp

        process.run()


class SPM_Realign(Process):

    def __init__(self):
        super(SPM_Realign, self).__init__()

        # Inputs
        #self.add_trait("in_files", InputMultiPath(traits.Either(
        #    ImageFileSPM(exists=True), traits.List(ImageFileSPM(exists=True)), output=False, copyfile=True)))
        self.add_trait("in_files", InputMultiPath(traits.Either(
            ImageFileSPM(), traits.List(ImageFileSPM()), output=False, copyfile=True)))

        """self.add_trait("fwhm",
                       traits.Range(low=0.0, output=False, optional=True))"""
        self.add_trait("fwhm", traits.Float(5.0, output=False, optional=True))

        # self.add_trait("interp", traits.Range(low=0, high=7, output=False, optional=True))
        self.add_trait("interp", traits.Int(2, output=False, optional=True))

        self.add_trait("jobtype", traits.Enum('estwrite', 'estimate', 'write', usedefault=True,
                                              output=False, optional=True))
        self.add_trait("out_prefix", traits.String('r', output=False, optional=True))
        # self.add_trait("quality", traits.Range(low=0.0, high=1.0, output=False, optional=True))

        # self.add_trait("quality", traits.Float(output=False, optional=True))
        self.add_trait("quality", traits.Float(0.9, output=False, optional=True))

        # self.add_trait("register_to_mean", traits.Bool(output=False, optional=True))
        self.add_trait("register_to_mean", traits.Bool(True, output=False, optional=True))

        # self.add_trait("separation", traits.Range(low=0.0, output=False, optional=True))
        self.add_trait("separation", traits.Float(4.0, output=False, optional=True))

        # self.add_trait("separation", traits.Float(output=False, optional=True))
        # self.add_trait("weight_img", File(output=False, optional=True))

        # self.add_trait("wrap", traits.List(traits.Int(), output=False, optional=True))
        self.add_trait("wrap", traits.List([0, 0, 0], output=False, optional=True))

        # self.add_trait("write_interp", traits.Range(low=0, high=7, output=False, optional=True))
        self.add_trait("write_interp", traits.Int(4, output=False, optional=True))

        # self.add_trait("write_interp", traits.Int(output=False, optional=True))
        # self.add_trait("write_mask", traits.Bool(output=False, optional=True))
        self.add_trait("write_mask", traits.Bool(True, output=False, optional=True))

        self.add_trait("write_which", traits.ListInt([2, 1], usedefault=True, output=False, optional=True))
        #self.add_trait("write_wrap", traits.Range(traits.Int(), output=False, optional=True))
        self.add_trait("write_wrap", traits.ListInt([0, 0, 0], output=False, optional=True))

        # Outputs
        self.add_trait("realigned_files", OutputMultiPath(
            traits.Either(traits.List(File()), File()), output=True))
        self.add_trait("mean_image", File(output=True))
        self.add_trait("realignment_parameters", File(output=True)) #rp_

    def list_outputs(self):
        process = spm.Realign()
        if not self.in_files:
            return {}
        else:
            process.inputs.in_files = self.in_files
        outputs = process._list_outputs()
        return outputs

    def _run_process(self):

        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

        process = spm.Realign()
        process.inputs.in_files = self.in_files
        process.inputs.fwhm = self.fwhm
        process.inputs.interp = self.interp
        process.inputs.jobtype = self.jobtype
        process.inputs.out_prefix = self.out_prefix
        process.inputs.quality = self.quality
        process.inputs.register_to_mean = self.register_to_mean
        process.inputs.separation = self.separation
        #seg.inputs.weight_img = self.weight_img
        process.inputs.wrap = self.wrap
        process.inputs.write_interp = self.write_interp
        process.inputs.write_mask = self.write_mask
        process.inputs.write_which = self.write_which
        process.inputs.write_wrap = self.write_wrap

        process.run()


class SPM_Coregister(Process):

    def __init__(self):
        super(SPM_Coregister, self).__init__()

        # Inputs
        self.add_trait("target", ImageFileSPM(output=False, copyfile=False))
        self.add_trait("source", InputMultiPath(ImageFileSPM(), output=False, copyfile=True))
        self.add_trait("jobtype", traits.Enum('estimate', 'estwrite', 'write', usedefault=True,
                                              output=False, optional=True))
        self.add_trait("apply_to_files", InputMultiPath(File(), output=False, copyfile=True)) # Optional for Nipype
        self.add_trait("cost_function", traits.Enum('nmi', 'mi', 'ecc', 'ncc', output=False, optional=True))

        # self.add_trait("fwhm", traits.List(traits.Float(), output=False, optional=True))
        self.add_trait("fwhm", traits.List([7, 7], output=False, optional=True))

        # self.add_trait("separation", traits.List(traits.Float(), output=False, optional=True))
        self.add_trait("separation", traits.List([4, 2], output=False, optional=True))

        # self.add_trait("tolerance", traits.List(traits.Float(), output=False, optional=True))
        self.add_trait("tolerance",
                       traits.List([.02, .02, .02, 0.001, 0.001, 0.001, .01, .01, .01, 0.001, 0.001, 0.001],
                                   output=False, optional=True))

        # Outputs
        self.add_trait("coregistered_files", OutputMultiPath(File(), output=True))

    def list_outputs(self):
        process = spm.Coregister()
        if not self.target:
            return {}
        else:
            process.inputs.target = self.target
        if not self.source:
            return {}
        else:
            process.inputs.source = self.source
        if not self.apply_to_files:
            return {}
        else:
            process.inputs.apply_to_files = self.apply_to_files
        if not self.jobtype:
            return {}
        else:
            process.inputs.jobtype = self.jobtype
        outputs = process._list_outputs()

        return outputs

    def _run_process(self):

        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

        process = spm.Coregister()
        process.inputs.target = self.target
        process.inputs.source = self.source
        process.inputs.jobtype = self.jobtype
        process.inputs.apply_to_files = self.apply_to_files
        process.inputs.cost_function = self.cost_function
        process.inputs.fwhm = self.fwhm
        process.inputs.separation = self.separation
        process.inputs.tolerance = self.tolerance

        process.run()


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


class SPM_Level1Design(Process):

    def __init__(self):
        super(SPM_Level1Design, self).__init__()

        # Inputs
        self.add_trait("timing_units", traits.Enum('scans', 'secs', use_default=True,
                                                   output=False, copyfile=False, optional=True))
        self.add_trait("interscan_interval", traits.Float(3.0, output=False, usedefault=True))

        self.add_trait("microtime_resolution", traits.Int(16, usedefault=True, output=False, optional=True))
        self.add_trait("microtime_onset", traits.Float(1.0, usedefault=True, output=False, optional=True))
        self.add_trait("session_info", traits.Any(output=False, optional=True)) # TODO: Find the value to add
        self.add_trait("factor_info", traits.List(traits.Dict(traits.Enum('name', 'levels')),
                                                   output=False, optional=True))
        #self.add_trait("factor_info", traits.Dict(traits.Enum('name', 'levels'), output=False, optional=True))
        self.add_trait("bases", traits.Dict(traits.Enum('hrf', 'fourier', 'fourier_han', 'gamma', 'fir'),
                                            output=False))
        self.add_trait("volterra_expansion_order", traits.Int(1, output=False, optional=True))
        self.add_trait("global_intensity_normalization", traits.Enum('none', 'scaling', output=False, optional=True))
        self.add_trait("mask_image", traits.File(output=False, optional=True))
        self.add_trait("mask_threshold", traits.Float(0.8, output=False, optional=True))
        self.add_trait("model_serial_correlations", traits.Enum('AR(1)', 'FAST', 'none', output=False, optional=True))

        # Output
        self.add_trait("spm_mat_file", traits.File(output=True))

    def _run_process(self):

        # Removing the spm_mat_file to avoid a bug
        cur_dir = os.getcwd()
        out_file = os.path.join(cur_dir, 'SPM.mat')
        if os.path.isfile(out_file):
            os.remove(out_file)

        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

        process = spm.Level1Design()
        process.inputs.timing_units = self.timing_units
        process.inputs.interscan_interval = self.interscan_interval
        process.inputs.microtime_resolution = self.microtime_resolution
        process.inputs.microtime_onset = self.microtime_onset
        process.inputs.factor_info = self.factor_info
        process.inputs.bases = self.bases
        process.inputs.volterra_expansion_order = self.volterra_expansion_order
        process.inputs.global_intensity_normalization = self.global_intensity_normalization
        process.inputs.mask_image = self.mask_image
        process.inputs.mask_threshold = self.mask_threshold
        process.inputs.model_serial_correlations = self.model_serial_correlations

        #process.inputs.session_info = self.get_session_info()
        process.inputs.session_info = self.session_info

        process.run()

        # Copying the generated SPM.mat file in the data directory
        mask_image_folder, mask_image_name = os.path.split(self.mask_image)
        from shutil import copy2
        copy2(out_file, mask_image_folder)

        if os.path.isfile(out_file):
            os.remove(out_file)

    def get_session_info(self):
        from nipype.algorithms import modelgen
        from nipype.interfaces.base import Bunch
        s = modelgen.SpecifyModel()

        for key, value in self.session_info.items():
            if key == 'subject_info' and value == {}:
                value = None
            setattr(s.inputs, key, value)

        return s

    def list_outputs(self):
        # Copying the generated SPM.mat file in the data directory
        mask_image_folder, mask_image_name = os.path.split(self.mask_image)
        out_file = os.path.join(mask_image_folder, 'SPM.mat')

        d = {'spm_mat_file': out_file}
        return d


class SPM_EstimateModel(Process):

    def __init__(self):
        super(SPM_EstimateModel, self).__init__()

        # Inputs
        self.add_trait("spm_mat_file", File(output=False, copyfile=True))
        self.add_trait("estimation_method", traits.Dict(traits.Enum('Classical', 'Bayesian2', 'Bayesian'),
                                                        output=False))
        self.add_trait("write_residuals", traits.Bool(output=False, optional=True))
        self.add_trait("flags", traits.Dict(output=False, optional=True))
        self.add_trait("version", traits.String("spm12", output=False, optional=True))

        # Outputs
        self.add_trait("mask_image", ImageFileSPM(output=True, optional=True))
        self.add_trait("beta_images", OutputMultiPath(output=True, optional=True))
        self.add_trait("residual_image", ImageFileSPM(output=True, optional=True))
        self.add_trait("residual_images", OutputMultiPath(output=True, optional=True))
        self.add_trait("RPVimage", ImageFileSPM(output=True, optional=True))
        self.add_trait("out_spm_mat_file", File(output=True))
        self.add_trait("labels", ImageFileSPM(output=True, optional=True))
        self.add_trait("SDerror", OutputMultiPath(ImageFileSPM(), output=True, optional=True))
        self.add_trait("ARcoef", OutputMultiPath(ImageFileSPM(), output=True, optional=True))
        self.add_trait("Cbetas", OutputMultiPath(ImageFileSPM(), output=True, optional=True))
        self.add_trait("SDbetas", OutputMultiPath(ImageFileSPM(), output=True, optional=True))

    def list_outputs(self):
        process = spm.EstimateModel()
        if not self.spm_mat_file:
            return {}
        else:
            process.inputs.spm_mat_file = self.spm_mat_file
        if not self.estimation_method:
            return {}
        else:
            process.inputs.estimation_method = self.estimation_method
        process.inputs.write_residuals = self.write_residuals
        process.inputs.flags = self.flags

        outputs = process._list_outputs()
        print("OUTPUT ESTIMATE MODEL", outputs)

        outputs["out_spm_mat_file"] = outputs.pop("spm_mat_file")

        return outputs

    def _run_process(self):

        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

        process = spm.EstimateModel()
        process.inputs.spm_mat_file = self.spm_mat_file

        # Removing the image files to avoid a bug
        outputs = self.list_outputs()
        for key, value in outputs.items():
            if key not in ["out_spm_mat_file"]:
                if value not in ["<undefined>", Undefined]:
                    if os.path.isfile(value):
                        os.remove(value)

        process.inputs.estimation_method = self.estimation_method
        process.inputs.write_residuals = self.write_residuals
        process.inputs.flags = self.flags

        process.run()


class SPM_EstimateContrast(Process):

    def __init__(self):
        super(SPM_EstimateContrast, self).__init__()

        # Inputs
        self.add_trait("spm_mat_file", File(output=False, copyfile=True))
        self.add_trait("contrasts", traits.List(output=False))
        self.add_trait("beta_images", InputMultiPath(File(), output=False, copyfile=False))
        self.add_trait("residual_image", File(output=False, copyfile=False))
        self.add_trait("use_derivs", traits.Bool(output=False, optional=True, xor=['group_contrast']))
        self.add_trait("group_contrast", traits.Bool(output=False, optional=True, xor=['use_derivs']))

        # Outputs
        self.add_trait("con_images", OutputMultiPath(File(), optional=True, output=True))
        self.add_trait("spmT_images ", OutputMultiPath(File(), optional=True, output=True))
        self.add_trait("spmF_images", OutputMultiPath(File(), optional=True, output=True))
        self.add_trait("out_spm_mat_file", File(output=True, copyfile=False))

    def list_outputs(self):
        process = spm.EstimateContrast()
        if not self.spm_mat_file:
            return {}
        else:
            process.inputs.spm_mat_file = self.spm_mat_file
        if not self.contrasts:
            return {}
        else:
            process.inputs.contrasts = self.contrasts
        if not self.beta_images:
            return {}
        else:
            process.inputs.beta_images = self.beta_images
        if not self.residual_image:
            return {}
        else:
            process.inputs.residual_image = self.residual_image

        outputs = process._list_outputs()
        outputs["out_spm_mat_file"] = outputs.pop("spm_mat_file")
        return outputs

    def _run_process(self):
        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
        process = spm.EstimateContrast()
        process.inputs.spm_mat_file = self.spm_mat_file
        process.inputs.contrasts = self.contrasts
        process.inputs.beta_images = self.beta_images
        process.inputs.residual_image = self.residual_image
        if self.use_derivs is not None:
            process.inputs.use_derivs = self.use_derivs
        else:
            if self.group_contrast is not None:
                process.inputs.group_contrast = self.group_contrast
            else:
                process.inputs.use_derivs = False
        process.run()


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
                                                           os.path.basename(parametric_map)[0:4])))
                
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
            map_name = os.path.basename(parametric_map)[0:4]
            for calcul in self.calculs:
                out_file = os.path.join(analysis_dir,
                                        "{0}_{1}_{2}.xls".format(self.data, calcul, map_name))

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
                            roi_file = os.path.join(analysis_dir, "{0}_mean{1}.txt".format(roi[0] + roi[1], map_name))
                            with open(roi_file, 'r') as f_read:
                                final_res = float(f_read.read())
                            f.write("{0}\t".format(final_res))

                    elif calcul == 'IL_mean':
                        roi_checked = []
                        for roi in self.roi_list:
                            if roi[0] in roi_checked:
                                continue
                            roi_file = os.path.join(analysis_dir, "{0}_mean{1}.txt".format(roi[0] + roi[1], map_name))
                            with open(roi_file, 'r') as f_read:
                                roi_value = float(f_read.read())

                            # Searching the roi that has the same first element
                            roi_2 = [s for s in self.roi_list if roi[0] in s[0] and roi[1] != s[1]][0]
                            roi_file_2 = os.path.join(analysis_dir,
                                                      "{0}_mean{1}.txt".format(roi_2[0] + roi_2[1], map_name))
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
                            roi_file = os.path.join(analysis_dir, "{0}_std{1}.txt".format(roi[0] + roi[1], map_name))
                            with open(roi_file, 'r') as f_read:
                                final_res = float(f_read.read())
                            f.write("{0}\t".format(final_res))

                    elif calcul == 'IL_std':
                        roi_checked = []
                        for roi in self.roi_list:
                            if roi[0] in roi_checked:
                                continue
                            roi_file = os.path.join(analysis_dir, "{0}_std{1}.txt".format(roi[0] + roi[1], map_name))
                            with open(roi_file, 'r') as f_read:
                                roi_value = float(f_read.read())

                            # Searching the roi that has the same first element
                            roi_2 = [s for s in self.roi_list if roi[0] in s[0] and roi[1] != s[1]][0]
                            roi_file_2 = os.path.join(analysis_dir,
                                                      "{0}_std{1}.txt".format(roi_2[0] + roi_2[1], map_name))
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


def threshold(file_name, thresh):
    img = nib.load(file_name)
    img_data = img.get_data()
    img_thresh = (img_data > thresh).astype(float)
    img_final = nib.Nifti1Image(img_thresh, img.affine, img.header)
    return img_final
