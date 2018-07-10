from nipype.interfaces import spm
from nipype.interfaces.spm.base import (SPMCommand, SPMCommandInputSpec, ImageFileSPM)
from nipype.interfaces.base import (OutputMultiPath, TraitedSpec, traits, InputMultiPath, File)


class NewSegmentInputSpec(SPMCommandInputSpec):
    channel_files = InputMultiPath(
        ImageFileSPM(exists=True),
        mandatory=True,
        desc="A list of files to be segmented",
        field='channel',
        copyfile=False)
    channel_info = traits.Tuple(
        traits.Float(),
        traits.Float(),
        traits.Tuple(traits.Bool, traits.Bool),
        desc="""A tuple with the following fields:
            - bias reguralisation (0-10)
            - FWHM of Gaussian smoothness of bias
            - which maps to save (Corrected, Field) - a tuple of two boolean values""",
        field='channel')
    tissues = traits.List(
        traits.Tuple(
            traits.Tuple(ImageFileSPM(exists=True), traits.Int()),
            traits.Int(), traits.Tuple(traits.Bool, traits.Bool),
            traits.Tuple(traits.Bool, traits.Bool)),
        desc="""A list of tuples (one per tissue) with the following fields:
            - tissue probability map (4D), 1-based index to frame
            - number of gaussians
            - which maps to save [Native, DARTEL] - a tuple of two boolean values
            - which maps to save [Unmodulated, Modulated] - a tuple of two boolean values""",
        field='tissue')
    affine_regularization = traits.Enum(
        'mni',
        'eastern',
        'subj',
        'none',
        field='warp.affreg',
        desc='mni, eastern, subj, none ')
    warping_regularization = traits.Either(
        traits.List(traits.Float(), minlen=5, maxlen=5),
        traits.Float(),
        field='warp.reg',
        desc=('Warping regularization '
              'parameter(s). Accepts float '
              'or list of floats (the '
              'latter is required by '
              'SPM12)'))
    sampling_distance = traits.Float(
        field='warp.samp',
        desc=('Sampling distance on data for '
              'parameter estimation'))
    write_deformation_fields = traits.List(
        traits.Bool(),
        minlen=2,
        maxlen=2,
        field='warp.write',
        desc=("Which deformation fields to "
              "write:[Inverse, Forward]"))


class NewSegmentOutputSpec(TraitedSpec):
    native_class_images = traits.List(
        traits.List(File(exists=True)), desc='native space probability maps')
    dartel_input_images = traits.List(
        traits.List(File(exists=True)), desc='dartel imported class images')
    normalized_class_images = traits.List(
        traits.List(File(exists=True)), desc='normalized class images')
    modulated_class_images = traits.List(
        traits.List(File(exists=True)),
        desc=('modulated+normalized class '
              'images'))
    transformation_mat = OutputMultiPath(
        File(exists=True), desc='Normalization transformation')
    bias_corrected_images = OutputMultiPath(
        File(exists=True), desc='bias corrected images')
    bias_field_images = OutputMultiPath(
        File(exists=True), desc='bias field images')
    forward_deformation_field = OutputMultiPath(File(exists=True))
    inverse_deformation_field = OutputMultiPath(File(exists=True))


class NewSegmentIrmage(spm.NewSegment):
    """ Derived class to avoid to instance SPMCommand that takes around 6 seconds. """

    input_spec = NewSegmentInputSpec
    output_spec = NewSegmentOutputSpec

    def __init__(self, **inputs):
        self._jobtype = 'spatial'
        self._jobname = 'preproc'
        SPMCommand.__init__(self, **inputs)

