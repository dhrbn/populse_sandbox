# Capsul import
from capsul.api import Process
from capsul.api import StudyConfig, get_process_instance
from capsul.subprocess import fsl
from soma.path import find_in_path
from nilearn.plotting import plot_anat, show

# Trait import
from traits.api import Float, File

# Other imports
import os


class AvailableProcesses(list):
    # List of all the processes in this file
    def __init__(self):
        super(AvailableProcesses, self).__init__()
        self.append(Addition)
        self.append(Substraction)
        self.append(FSL_Smooth)


class Addition(Process):

    def __init__(self, node_name):
        super(Addition, self).__init__()

        self.node_name = node_name

        self.add_trait(node_name + "_in_1", Float(output=False))
        self.add_trait(node_name + "_in_2", Float(output=False))
        self.add_trait(node_name + "_out", Float(output=True))

        self.in_1 = getattr(self, node_name + "_in_1")
        self.in_2 = getattr(self, node_name + "_in_2")
        self.out = getattr(self, node_name + "_out")

    def _run_process(self):
        self.get_plugs_values()

        self.out = self.in_1 + self.in_2
        print('Addition - ', self.node_name, '\n...\nInputs: {', self.in_1, ', ',
              self.in_2, '}\nOutput: ', self.out, '\n...\n')

        self.set_plugs_values()

    def get_plugs_values(self):
        self.in_1 = getattr(self, self.node_name + "_in_1")
        self.in_2 = getattr(self, self.node_name + "_in_2")
        self.out = getattr(self, self.node_name + "_out")

    def set_plugs_values(self):
        setattr(self, self.node_name + "_in_1", self.in_1)
        setattr(self, self.node_name + "_in_2", self.in_2)
        setattr(self, self.node_name + "_out", self.out)


class Substraction(Process):

    def __init__(self, node_name):
        super(Substraction, self).__init__()

        self.node_name = node_name

        self.add_trait(node_name + "_in_1", Float(output=False))
        self.add_trait(node_name + "_in_2", Float(output=False))
        self.add_trait(node_name + "_out", Float(output=True))

        self.in_1 = getattr(self, node_name + "_in_1")
        self.in_2 = getattr(self, node_name + "_in_2")
        self.out = getattr(self, node_name + "_out")

    def _run_process(self):
        self.get_plugs_values()

        self.out = self.in_1 - self.in_2
        print('Substraction - ', self.node_name, '\n...\nInputs: {', self.in_1, ', ',
              self.in_2, '}\nOutput: ', self.out, '\n...\n')

        self.set_plugs_values()

    def get_plugs_values(self):
        self.in_1 = getattr(self, self.node_name + "_in_1")
        self.in_2 = getattr(self, self.node_name + "_in_2")
        self.out = getattr(self, self.node_name + "_out")

    def set_plugs_values(self):
        setattr(self, self.node_name + "_in_1", self.in_1)
        setattr(self, self.node_name + "_in_2", self.in_2)
        setattr(self, self.node_name + "_out", self.out)


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

        study_config = StudyConfig(modules=StudyConfig.default_modules + ['NipypeConfig'])

        # Process
        if study_config.use_nipype:
            try:
                smooth_process = get_process_instance("nipype.interfaces.fsl.Smooth")
                smooth_process.output_type = 'NIFTI'
                print(self.in_file)
                smooth_process.in_file = self.in_file
                print(smooth_process.in_file)
                #smooth_process.in_file = "/home/david/Nifti_data/1103/3/NIFTI/1103_3.nii"
                if self.fwhm > 0:
                    smooth_process.fwhm = self.fwhm
                elif self.sigma > 0:
                    smooth_process.sigma = self.sigma
                else:
                    smooth_process.fwhm = 2.0
                print(os.path.split(self.out_file)[0])
                smooth_process.output_directory = os.path.split(self.out_file)[0]
                print(smooth_process.output_directory)
                #smooth_process.output_directory = '/home/david/Nifti_data'

            except:
                smooth_process = None
                print('Smooth module of FSL is not present.')

        else:
            smooth_process = None
            print('NiPype is not present.')

        if smooth_process:
            study_config.reset_process_counter()
            study_config.run(smooth_process, verbose=1)

            # Display
            print('Smoothing with FSL - ', self.node_name, '\n...\nInputs: {', self.in_file, ', ',
                  self.fwhm, '}\nOutput: ', self.out_file, '\n...\n')

            import subprocess
            #subprocess.check_output(['fslview', '/home/david/Nifti_data/1103/3/NIFTI/1103_3.nii'])
            #subprocess.check_output(['fslview', '/home/david/Nifti_data/1103_3_smooth.nii'])
            out_file = os.path.join(smooth_process.output_directory, os.path.basename(self.in_file)[:-4] + '_smooth.nii')
            print(out_file)
            subprocess.check_output(['fslview', self.in_file])
            subprocess.check_output(['fslview', out_file])

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
