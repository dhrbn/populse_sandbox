# Trait import
from traits.api import Undefined
from nipype.interfaces.base import OutputMultiPath, InputMultiPath, File, traits
from nipype.interfaces.spm.base import ImageFileSPM

# Other import
import os
from nipype.interfaces import spm

# MIA import
from PipelineManager.Process_mia import Process_mia
from SoftwareProperties.Config import Config

config = Config()


class Level1Design(Process_mia):

    def __init__(self):
        super(Level1Design, self).__init__()

        session_info_dict = {'cond': {}, 'hpf': 427.2, 'multi': {}, 'regress': {}}

        # Inputs
        self.add_trait("timing_units", traits.Enum('scans', 'secs', use_default=True,
                                                   output=False, copyfile=False, optional=True))
        self.add_trait("interscan_interval", traits.Float(3.0, output=False, usedefault=True, optional=True))

        self.add_trait("microtime_resolution", traits.Int(16, usedefault=True, output=False, optional=True))
        self.add_trait("microtime_onset", traits.Float(1.0, usedefault=True, output=False, optional=True))
        self.add_trait("session_info", traits.Dict(session_info_dict, output=False, optional=True))
        self.add_trait("factor_info", traits.List([{}], output=False, optional=True))
        self.add_trait("bases", traits.Dict({'hrf': {'derivs': [0, 0]}}, output=False, optional=True))
        self.add_trait("scans", traits.List(traits.File(), output=False))
        self.add_trait("multi_reg", traits.List(traits.File(), output=False))
        self.add_trait("volterra_expansion_order", traits.Int(1, output=False, optional=True))
        self.add_trait("global_intensity_normalization", traits.Enum('none', 'scaling', output=False, optional=True))
        self.add_trait("mask_image", InputMultiPath(output=False))
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

        spm.SPMCommand.set_mlab_paths(matlab_cmd=config.get_matlab_command(), use_mcr=True)

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

        process.inputs.session_info = self.get_session_info()
        # process.inputs.session_info = self.session_info

        process.run()

        # Copying the generated SPM.mat file in the data directory
        mask_image = self.mask_image[0]
        mask_image_folder, mask_image_name = os.path.split(mask_image)
        from shutil import copy2
        copy2(out_file, mask_image_folder)

        if os.path.isfile(out_file):
            os.remove(out_file)

    def get_session_info(self):
        session_info = self.session_info
        session_info['scans'] = self.scans
        session_info['multi_reg'] = []
        for reg in self.multi_reg:
            reg_to_add = [{reg}]
            session_info['multi_reg'].append(reg_to_add)

        return session_info

    def list_outputs(self):

        if not self.mask_image or self.mask_image in ['<undefined>', Undefined]:
            return {}, {}

        if not self.scans or self.scans in ['<undefined>', Undefined]:
            return {}, {}

        if not self.multi_reg or self.multi_reg in ['<undefined>', Undefined]:
            return {}, {}

        # Copying the generated SPM.mat file in the data directory
        if self.scans and self.scans not in ['<undefined>', Undefined] \
                and self.scans[0] not in ['<undefined>', Undefined]:
            mask_image_folder, mask_image_name = os.path.split(self.scans[0])
            out_file = os.path.join(mask_image_folder, 'SPM.mat')

            d = {'spm_mat_file': out_file}
        else:
            d = {}
        return d


class EstimateModel(Process_mia):

    def __init__(self):
        super(EstimateModel, self).__init__()

        # Inputs
        self.add_trait("spm_mat_file", File(output=False, copyfile=True))
        self.add_trait("estimation_method", traits.Dict({'Classical': 1}, output=False, optional=True))
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
            return {}, {}
        else:
            process.inputs.spm_mat_file = self.spm_mat_file
        if not self.estimation_method:
            return {}, {}
        else:
            process.inputs.estimation_method = self.estimation_method

        process.inputs.write_residuals = self.write_residuals
        process.inputs.flags = self.flags

        outputs = process._list_outputs()
        outputs["out_spm_mat_file"] = outputs.pop("spm_mat_file")

        return outputs

    def _run_process(self):

        spm.SPMCommand.set_mlab_paths(matlab_cmd=config.get_matlab_command(), use_mcr=True)

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


class EstimateContrast(Process_mia):

    def __init__(self):
        super(EstimateContrast, self).__init__()

        # Inputs
        self.add_trait("spm_mat_file", File(output=False, copyfile=True))
        self.add_trait("contrasts", traits.List([('+', 'T', ['R1_1'], [1])], output=False, optional=True))
        self.add_trait("beta_images", InputMultiPath(File(), output=False, copyfile=False))
        self.add_trait("residual_image", File(output=False, copyfile=False))
        self.add_trait("use_derivs", traits.Bool(output=False, optional=True, xor=['group_contrast']))
        self.add_trait("group_contrast", traits.Bool(output=False, optional=True, xor=['use_derivs']))

        # Outputs
        self.add_trait("con_images", OutputMultiPath(File(), optional=True, output=True))
        self.add_trait("spmT_images", OutputMultiPath(File(), optional=True, output=True))
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
        return outputs, {}

    def _run_process(self):
        spm.SPMCommand.set_mlab_paths(matlab_cmd=config.get_matlab_command(), use_mcr=True)
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