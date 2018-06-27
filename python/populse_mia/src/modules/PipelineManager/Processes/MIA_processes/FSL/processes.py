# Capsul import
from capsul.api import StudyConfig, get_process_instance

# Trait import
from traits.api import Float
from nipype.interfaces.base import File

# Other import
import os
from nipype.interfaces import fsl

from PipelineManager.Process_mia import Process_mia


class Smooth(Process_mia):

    def __init__(self):
        super(Smooth, self).__init__()

        self.add_trait("in_file", File(output=False))
        self.add_trait("fwhm", Float(output=False))
        self.add_trait("out_file", File(output=True))

    def _run_process(self):

        smooth_process = get_process_instance(fsl.Smooth)
        smooth_process.in_file = self.in_file
        smooth_process.out_file = self.out_file
        smooth_process.output_type = 'NIFTI'


        if self.fwhm > 0:
            smooth_process.fwhm = self.fwhm
        else:
            smooth_process.fwhm = 2.0

        smooth_process.output_directory = os.path.split(self.out_file)[0]
        smooth_process.run()

        import nibabel as nib
        img = nib.load(self.out_file + ".gz")
        nib.save(img, self.out_file)

        # Display
        print('Smoothing with FSL\n...\nInputs: {', self.in_file, ', ',
              self.fwhm, '}\nOutput: ', self.out_file, '\n...\n')