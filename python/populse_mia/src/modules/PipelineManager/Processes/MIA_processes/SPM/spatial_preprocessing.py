# Trait import
from traits.api import Float
from nipype.interfaces.base import OutputMultiPath, InputMultiPath, File, traits
from nipype.interfaces.spm.base import ImageFileSPM

# Other import
import os
from nipype.interfaces import spm

# MIA import
from PipelineManager.Process_mia import Process_mia
from SoftwareProperties.Config import Config
from .nipype_extension import NewSegmentIrmage

config = Config()


class Smooth(Process_mia):

    def __init__(self):
        super(Smooth, self).__init__()

        # Inputs
        self.add_trait("in_files", InputMultiPath(ImageFileSPM(), copyfile=False, output=False))
        self.add_trait("fwhm", traits.List([6, 6, 6], output=False, optional=True))
        self.add_trait("data_type", traits.Int(output=False, optional=True))
        self.add_trait("implicit_masking", traits.Bool(output=False, optional=True))
        self.add_trait("out_prefix", traits.String('s', usedefault=True, output=False, optional=True))

        # Output
        self.add_trait("smoothed_files", OutputMultiPath(File(), output=True))

    def list_outputs(self):
        Process_mia.list_outputs(self)
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

        # Can be removed if the SPMMCRCMD and FORCE_SPMMCR environment variable are set correctly.
        # spm.SPMCommand.set_mlab_paths(matlab_cmd=config.get_matlab_command(), use_mcr=True)

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


class NewSegment(Process_mia):

    def __init__(self):
        super(NewSegment, self).__init__()

        tissues_list = [(('/home/david/spm12/spm12_mcr/spm/spm12/tpm/TPM.nii', 1), 2, (True, False), (False, False)),
                        (('/home/david/spm12/spm12_mcr/spm/spm12/tpm/TPM.nii', 2), 2, (True, False), (False, False)),
                        (('/home/david/spm12/spm12_mcr/spm/spm12/tpm/TPM.nii', 3), 2, (True, False), (False, False)),
                        (('/home/david/spm12/spm12_mcr/spm/spm12/tpm/TPM.nii', 4), 3, (True, False), (False, False)),
                        (('/home/david/spm12/spm12_mcr/spm/spm12/tpm/TPM.nii', 5), 4, (True, False), (False, False)),
                        (('/home/david/spm12/spm12_mcr/spm/spm12/tpm/TPM.nii', 6), 2, (True, False), (False, False))]

        # Inputs
        self.add_trait("affine_regularization", traits.Enum('mni', 'eastern', 'subj', 'none',
                                                            output=False, optional=True))
        self.add_trait("channel_files", InputMultiPath(output=False))
        """self.add_trait("channel_info", traits.Tuple(traits.Float(), traits.Float(),
                                                    traits.Tuple(traits.Bool, traits.Bool(True)),
                                                    output=False, optional=True))"""
        self.add_trait("channel_info", traits.Tuple((0.0001, 60, (False, True)), output=False, optional=True))
        self.add_trait("sampling_distance", Float(3.0, output=False, optional=True))
        """self.add_trait("tissues", traits.List(traits.Tuple(
            traits.Tuple(ImageFileSPM(exists=True), traits.Int()),
            traits.Int(), traits.Tuple(traits.Bool, traits.Bool),
            traits.Tuple(traits.Bool, traits.Bool)), output=False, optional=True))"""

        self.add_trait("tissues", traits.List(tissues_list, output=False, optional=True))

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
        Process_mia.list_outputs(self)
        process = NewSegmentIrmage()

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

        # Can be removed if the SPMMCRCMD and FORCE_SPMMCR environment variable are set correctly.
        # spm.SPMCommand.set_mlab_paths(matlab_cmd=config.get_matlab_command(), use_mcr=True)

        process = spm.NewSegment()
        process.inputs.affine_regularization = self.affine_regularization
        process.inputs.channel_files = self.channel_files
        process.inputs.channel_info = self.channel_info
        process.inputs.sampling_distance = self.sampling_distance
        process.inputs.tissues = self.tissues
        process.inputs.warping_regularization = self.warping_regularization
        process.inputs.write_deformation_fields = self.write_deformation_fields

        process.run()


class Normalize(Process_mia):

    def __init__(self):
        super(Normalize, self).__init__()

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
        self.add_trait("write_bounding_box", traits.List([[-78, -112, -50], [78, 76, 85]],
                                                         output=False, optional=True))
        # self.add_trait("write_voxel_sizes", traits.List(traits.Float(), output=False, optional=True))
        self.add_trait("write_voxel_sizes", traits.List([1, 1, 1], output=False, optional=True))
        # self.add_trait("write_interp", traits.Range(low=0, high=7, output=False, optional=True))
        self.add_trait("write_interp", traits.Int(1, output=False, optional=True))

        # Output
        self.add_trait("normalized_files", OutputMultiPath(File(), output=True))

    def list_outputs(self):
        Process_mia.list_outputs(self)
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

        # Can be removed if the SPMMCRCMD and FORCE_SPMMCR environment variable are set correctly.
        # spm.SPMCommand.set_mlab_paths(matlab_cmd=config.get_matlab_command(), use_mcr=True)

        process = spm.Normalize12()
        process.inputs.apply_to_files = self.apply_to_files
        process.inputs.deformation_file = self.deformation_file
        process.inputs.jobtype = self.jobtype
        process.inputs.write_bounding_box = self.write_bounding_box
        process.inputs.write_voxel_sizes = self.write_voxel_sizes
        process.inputs.write_interp = self.write_interp

        process.run()


class Realign(Process_mia):

    def __init__(self):
        super(Realign, self).__init__()

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
        Process_mia.list_outputs(self)
        process = spm.Realign()
        if not self.in_files:
            return {}
        else:
            process.inputs.in_files = self.in_files
        outputs = process._list_outputs()
        return outputs

    def _run_process(self):

        # Can be removed if the SPMMCRCMD and FORCE_SPMMCR environment variable are set correctly.
        # spm.SPMCommand.set_mlab_paths(matlab_cmd=config.get_matlab_command(), use_mcr=True)

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


class Coregister(Process_mia):

    def __init__(self):
        super(Coregister, self).__init__()

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
        Process_mia.list_outputs(self)
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

        if self.jobtype == "estimate":
            outputs["coregistered_files"] = self.apply_to_files

        return outputs, {}

    def _run_process(self):

        # Can be removed if the SPMMCRCMD and FORCE_SPMMCR environment variable are set correctly.
        # spm.SPMCommand.set_mlab_paths(matlab_cmd=config.get_matlab_command(), use_mcr=True)

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
