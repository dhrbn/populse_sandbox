# Capsul import
from capsul.api import Process
from capsul.api import StudyConfig, get_process_instance

# Trait import
from traits.api import Float, File
from nipype.interfaces.base import OutputMultiPath, TraitedSpec, isdefined, InputMultiPath, File, Str, traits
from nipype.interfaces.spm.base import ImageFileSPM

# Other imports
import os

from nipype.interfaces import spm

matlab_cmd = '/home/david/spm12/run_spm12.sh /usr/local/MATLAB/MATLAB_Runtime/v93/ script'
spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)


class AvailableProcesses(list):
    # List of all the Processes in this file
    def __init__(self):
        super(AvailableProcesses, self).__init__()
        self.append(Addition)
        self.append(Substraction)
        self.append(FSL_Smooth)
        self.append(SPM_Smooth)


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


class SPM_Smooth(Process):

    def __init__(self):
        super(SPM_Smooth, self).__init__()

        self.add_trait("in_files", InputMultiPath(ImageFileSPM(), copyfile=False, output=False))
        self.add_trait("fwhm", traits.List(traits.Float, output=False, optional=True))
        self.add_trait("data_type", traits.Int(output=False, optional=True))
        self.add_trait("implicit_masking", traits.Bool(output=False, optional=True))
        self.add_trait("out_prefix", traits.String('s', usedefault=True, output=False, optional=True))

        self.add_trait("smoothed_files", OutputMultiPath(File(), output=True))

    def _run_process(self):

        seg = spm.Smooth()
        seg.inputs.in_files = self.in_files
        seg.inputs.fwhm = self.fwhm
        seg.inputs.data_type = self.data_type
        seg.inputs.implicit_masking = self.implicit_masking
        seg.inputs.out_prefix = self.out_prefix

        seg.run()


class SPM_NewSegment(Process):

    def __init__(self):
        super(SPM_NewSegment, self).__init__()

        self.add_trait("affine_regularization", traits.Enum('mni', 'eastern', 'subj', 'none',
                                                            output=False, optional=True))
        self.add_trait("channel_files", InputMultiPath(output=False))
        self.add_trait("channel_info", traits.Tuple(traits.Float(), traits.Float(), traits.Tuple(traits.Bool, traits.Bool),
                                                    output=False, optional=True))
        self.add_trait("sampling_distance", Float(output=False, optional=True))
        self.add_trait("tissues", traits.List(traits.Tuple(
            traits.Tuple(ImageFileSPM(exists=True), traits.Int()),
            traits.Int(), traits.Tuple(traits.Bool, traits.Bool),
            traits.Tuple(traits.Bool, traits.Bool)), output=False, optional=True))
        """self.add_trait("warping_regularization", traits.Either(
            traits.List(traits.Float(), minlen=5, maxlen=5),
            traits.Float(), output=False))"""
        self.add_trait("warping_regularization",
            traits.List(traits.Float(), output=False, optional=True))
        self.add_trait("write_deformation_fields", traits.List(
            traits.Bool(), output=False, optional=True))

        #self.add_trait(node_name + "_sigma", Float(output=False, optional=True))
        self.add_trait("forward_deformation_field", File(output=True))

    def _run_process(self):

        seg = spm.NewSegment()
        seg.inputs.affine_regularization = self.affine_regularization
        seg.inputs.channel_files = self.channel_files
        seg.inputs.channel_info = self.channel_info
        seg.inputs.sampling_distance = self.sampling_distance
        seg.inputs.tissues = self.tissues
        seg.inputs.warping_regularization = self.warping_regularization
        seg.inputs.write_deformation_fields = self.write_deformation_fields

        seg.run()


class SPM_Normalize(Process):

    def __init__(self):
        super(SPM_Normalize, self).__init__()

        self.add_trait("apply_to_files", InputMultiPath(traits.Either(
            ImageFileSPM(exists=True), traits.List(ImageFileSPM(exists=True)), output=False)))
        self.add_trait("deformation_file",
                       ImageFileSPM(output=False))
        self.add_trait("jobtype", traits.Enum('estwrite', 'est', 'write',
                                              usedefault=True, output=False, optional=True))
        self.add_trait("write_bounding_box", traits.List(traits.List(traits.Float()), output=False, optional=True))
        self.add_trait("write_voxel_sizes", traits.List(traits.Float(), output=False, optional=True))
        self.add_trait("write_interp", traits.Range(low=0, high=7, output=False, optional=True))

        self.add_trait("normalized_files", OutputMultiPath(File(), output=True))

    def _run_process(self):

        seg = spm.Normalize12()
        seg.inputs.apply_to_files = self.apply_to_files
        seg.inputs.deformation_file = self.deformation_file
        seg.inputs.jobtype = self.jobtype
        seg.inputs.write_bounding_box = self.write_bounding_box
        seg.inputs.write_voxel_sizes = self.write_voxel_sizes
        seg.inputs.write_interp = self.write_interp

        seg.run()


class SPM_Realign(Process):

    def __init__(self):
        super(SPM_Realign, self).__init__()

        self.add_trait("in_files", InputMultiPath(traits.Either(
            ImageFileSPM(exists=True), traits.List(ImageFileSPM(exists=True)), output=False, copyfile=True)))
        self.add_trait("fwhm",
                       traits.Range(low=0.0, output=False, optional=True))
        self.add_trait("interp", traits.Range(low=0, high=7, output=False, optional=True))
        self.add_trait("jobtype", traits.Enum('estwrite', 'estimate', 'write', usedefault=True,
                                              output=False, optional=True))
        self.add_trait("out_prefix", traits.String('r', output=False, optional=True))
        self.add_trait("quality", traits.Range(low=0.0, high=1.0, output=False, optional=True))
        #self.add_trait("quality", traits.Float(output=False, optional=True))
        self.add_trait("register_to_mean", traits.Bool(output=False, optional=True))
        self.add_trait("separation", traits.Range(low=0.0, output=False, optional=True))
        #self.add_trait("separation", traits.Float(output=False, optional=True))
        #self.add_trait("weight_img", File(output=False, optional=True))
        self.add_trait("wrap", traits.List(traits.Int(), output=False, optional=True))
        self.add_trait("write_interp", traits.Range(low=0, high=7, output=False, optional=True))
        #self.add_trait("write_interp", traits.Int(output=False, optional=True))
        self.add_trait("write_mask", traits.Bool(output=False, optional=True))
        self.add_trait("write_which", traits.ListInt([2, 1], usedefault=True, output=False, optional=True))
        #self.add_trait("write_wrap", traits.Range(traits.Int(), output=False, optional=True))
        self.add_trait("write_wrap", traits.ListInt([0, 0, 0], output=False, optional=True))

        self.add_trait("realigned_files", OutputMultiPath(
            traits.Either(traits.List(File()), File()), output=True))
        self.add_trait("mean_image", File(output=True))
        self.add_trait("realignment_parameters", File(output=True)) #rp_


    def _run_process(self):

        seg = spm.Realign()
        seg.inputs.in_files = self.in_files
        seg.inputs.fwhm = self.fwhm
        seg.inputs.interp = self.interp
        seg.inputs.jobtype = self.jobtype
        seg.inputs.out_prefix = self.out_prefix
        seg.inputs.quality = self.quality
        seg.inputs.register_to_mean = self.register_to_mean
        seg.inputs.separation = self.separation
        #seg.inputs.weight_img = self.weight_img
        seg.inputs.wrap = self.wrap
        seg.inputs.write_interp = self.write_interp
        seg.inputs.write_mask = self.write_mask
        seg.inputs.write_which = self.write_which
        seg.inputs.write_wrap = self.write_wrap

        seg.run()


class SPM_Coregister(Process):

    def __init__(self):
        super(SPM_Coregister, self).__init__()

        self.add_trait("target", ImageFileSPM(output=False, copyfile=False))
        self.add_trait("source", InputMultiPath(ImageFileSPM(), output=False, copyfile=True))
        self.add_trait("jobtype", traits.Enum('estwrite', 'estimate', 'write', usedefault=True,
                                              output=False, optional=True))
        self.add_trait("apply_to_files", InputMultiPath(File(), output=False, copyfile=True)) # Optional for Nipype
        self.add_trait("cost_function", traits.Enum('mi', 'nmi', 'ecc', 'ncc', output=False, optional=True))
        self.add_trait("fwhm", traits.List(traits.Float(), output=False, optional=True))
        self.add_trait("separation", traits.List(traits.Float(), output=False, optional=True))
        self.add_trait("tolerance", traits.List(traits.Float(), output=False, optional=True))

    def _run_process(self):

        seg = spm.Coregister()
        seg.inputs.target = self.target
        seg.inputs.source = self.source
        seg.inputs.jobtype = self.jobtype
        seg.inputs.apply_to_files = self.apply_to_files
        seg.inputs.cost_function = self.cost_function
        seg.inputs.fwhm = self.fwhm
        seg.inputs.separation = self.separation
        seg.inputs.tolerance = self.tolerance

        seg.run()


'''
class FSL_Smooth(Process):

    def __init__(self, node_name):
        super(FSL_Smooth, self).__init__()

        self.node_name = node_name

        self.add_trait(node_name + "_in_file", File(output=False))
        self.add_trait(node_name + "_fwhm", Float(output=False))
        #self.add_trait(node_name + "_sigma", Float(output=False, optional=True))
        self.add_trait(node_name + "_out_file", File(output=True))

        self.in_file = getattr(self, self.node_name + "_in_file")
        self.fwhm = getattr(self, self.node_name + "_fwhm")
        #self.sigma = getattr(self, self.node_name + "_sigma")
        self.out_file = getattr(self, self.node_name + "_out_file")

    def _run_process(self):
        self.get_plugs_values()

        study_config = StudyConfig(
            modules=["FSLConfig"],
            fsl_config="/etc/fsl/5.0/fsl.sh",
            use_smart_caching=True,
            output_directory="/tmp/capsul_demo",
            use_fsl=True)

        # Process
        if study_config.use_fsl:
            try:

                """smooth_process = get_process_instance("nipype.interfaces.fsl.Smooth")
                smooth_process.output_type = 'NIFTI'
                #smooth_process.in_file = self.in_file
                smooth_process.in_file = "/home/david/Nifti_data/1103/3/NIFTI/1103_3.nii"
                if self.fwhm > 0:
                    smooth_process.fwhm = self.fwhm
                elif self.sigma > 0:
                    smooth_process.sigma = self.sigma
                else:
                    smooth_process.fwhm = 2.0
                #smooth_process.output_directory = os.path.split(self.out_file)[0]
                smooth_process.output_directory = '/home/david/Nifti_data'"""

                fsl.check_call(study_config, ['fslmaths', "/home/david/Nifti_data/1103/3/NIFTI/1103_3.nii",
                                              "-kernel", str(self.fwhm), "-fmean",
                                              '/home/david/Nifti_data/test_fsl.nii'])

            except:
                smooth_process = None
                print('Smooth module of FSL is not present.')

        import subprocess
        subprocess.check_output(['fslview', '/home/david/Nifti_data/1103/3/NIFTI/1103_3.nii'])
        subprocess.check_output(['fslview', '/home/david/Nifti_data/test_fsl.nii'])

        self.set_plugs_values()

    def get_plugs_values(self):
        self.in_file = getattr(self, self.node_name + "_in_file")
        self.fwhm = getattr(self, self.node_name + "_fwhm")
        #self.sigma = getattr(self, self.node_name + "_sigma")
        self.out_file = getattr(self, self.node_name + "_out_file")

    def set_plugs_values(self):
        setattr(self, self.node_name + "_in_file", self.in_file)
        setattr(self, self.node_name + "_fwhm", self.fwhm)
        #setattr(self, self.node_name + "_sigma", self.sigma)
        setattr(self, self.node_name + "_out_file", self.out_file)


'''
