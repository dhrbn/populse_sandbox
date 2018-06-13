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

# To change to 'run_spm12.sh_location MCR_folder script"
matlab_cmd = '/home/david/spm12/run_spm12.sh /usr/local/MATLAB/MATLAB_Runtime/v93/ script'


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

        return outputs

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
        import nibabel as nib
        img = nib.load(file_name)
        img_data = img.get_data()
        img_thresh = (img_data > self.threshold).astype(float)
        img_final = nib.Nifti1Image(img_thresh, img.affine, img.header)

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
        self.add_trait("prefix", traits.String("", output=False, optional=True))
        self.add_trait("interp", traits.Int(1, output=False, optional=True))

        # Output
        self.add_trait("out_file", ImageFileSPM(output=True))

    def list_outputs(self):

        if not self.mask_to_resize:
            return {}

        if not self.suffix:
            self.suffix = ""

        if not self.prefix:
            self.prefix = ""

        mask_name = self.mask_to_resize
        if type(mask_name) in [list, TraitListObject]:
            mask_name = mask_name[0]

        path, file_name = os.path.split(mask_name)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        if file_name_no_ext[-4:] == "_002":
            file_name_no_ext = file_name_no_ext[:-4]
        out_file = os.path.join(path, self.prefix + file_name_no_ext + self.suffix + file_extension)

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

        if not self.suffix:
            self.suffix = ""

        if not self.prefix:
            self.prefix = ""

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
        out_file = os.path.join(path, self.prefix + file_name_no_ext + self.suffix + file_extension)
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

        outputs["out_spm_mat_file"] = outputs.pop("spm_mat_file")

        return outputs

    def _run_process(self):

        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

        process = spm.EstimateModel()
        process.inputs.spm_mat_file = self.spm_mat_file
        # TODO: TEST

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
        print("ALLO")
        if not self.spm_mat_file:
            return {}
        else:
            process.inputs.spm_mat_file = self.spm_mat_file
        print("contrasts", self.contrasts)
        if not self.contrasts:
            print('YOLOOO')
            return {}
        else:
            process.inputs.contrasts = self.contrasts
        print("beta_images", self.beta_images)
        if not self.beta_images:
            return {}
        else:
            process.inputs.beta_images = self.beta_images
        print("ALLO")
        if not self.residual_image:
            return {}
        else:
            process.inputs.residual_image = self.residual_image

        outputs = process._list_outputs()
        print("OUTPUTS", outputs)
        outputs["out_spm_mat_file"] = outputs.pop("spm_mat_file")
        return outputs

    def _run_process(self):
        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
        process = spm.EstimateContrast()
        process.inputs.spm_mat_file = self.spm_mat_file
        process.inputs.contrasts = self.contrasts
        process.inputs.beta_images = self.beta_images
        process.inputs.residual_image = self.residual_image
        process.inputs.use_derivs = self.use_derivs
        process.inputs.group_contrast = self.group_contrast
        process.run()


