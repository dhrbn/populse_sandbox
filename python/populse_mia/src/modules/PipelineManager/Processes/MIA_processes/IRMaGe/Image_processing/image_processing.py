# Trait import
from nipype.interfaces.base import OutputMultiPath, InputMultiPath, File, traits, TraitListObject, Undefined
from nipype.interfaces.spm.base import ImageFileSPM

# Other import
import os
from nipype.interfaces import spm
from skimage.transform import resize
import nibabel as nib
from distutils.dir_util import copy_tree

# MIA import
from PipelineManager.Process_mia import Process_mia
from SoftwareProperties.Config import Config

config = Config()


class Normalize_Spatial_Mask(Process_mia):

    def __init__(self):
        super(Normalize_Spatial_Mask, self).__init__()

        # Inputs
        """self.add_trait("apply_to_files", InputMultiPath(traits.Either(
            ImageFileSPM(), traits.List(ImageFileSPM()), output=False)))"""
        self.add_trait("apply_to_files", traits.List(output=False))
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

        return outputs, {}

    def _run_process(self):

        # Can be removed if the SPMMCRCMD and FORCE_SPMMCR environment variable are set correctly.
        # spm.SPMCommand.set_mlab_paths(matlab_cmd=config.get_matlab_command(), use_mcr=True)

        process = spm.Normalize12()
        process.inputs.apply_to_files = self._check_file_names()
        process.inputs.deformation_file = self.deformation_file
        process.inputs.jobtype = self.jobtype
        process.inputs.write_bounding_box = self.write_bounding_box
        process.inputs.write_voxel_sizes = self.write_voxel_sizes
        process.inputs.write_interp = self.write_interp

        process.run()


class Threshold(Process_mia):

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
            return {}, {}

        if not self.suffix:
            self.suffix = ""

        if not self.prefix:
            self.prefix = ""

        file_name = self._check_file_names()
        path, file_name = os.path.split(file_name)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        out_file = os.path.join(path, self.prefix.strip() + file_name_no_ext + self.suffix.strip() + file_extension)

        d = {'out_files': out_file}
        return d, {}

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
        out_file = os.path.join(path, self.prefix.strip() + file_name_no_ext + self.suffix.strip() + file_extension)
        nib.save(img_final, out_file)


class Resize(Process_mia):

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

        if file_name_no_ext.strip()[-4:] == "_002":
            file_name_no_ext = file_name_no_ext.strip()[:-4]
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


class Conv_ROI(Process_mia):

    def __init__(self):
        super(Conv_ROI, self).__init__()

        # Inputs
        self.add_trait("roi_list", traits.List(output=False))
        self.add_trait("mask", InputMultiPath(output=False))

        # Outputs
        self.add_trait("out_masks", OutputMultiPath(output=True))

    def list_outputs(self):
        if not self.roi_list:
            return {}
        if not self.mask:
            return {}
        mask = self.mask[0]

        roi_dir = os.path.join(os.path.dirname(mask), 'roi')
        if not os.path.isdir(roi_dir):
            os.mkdir(roi_dir)
            # Copying the ROIs from the ressources folder
            ref_dir = os.path.join('..', '..', 'ressources', 'reference_data', 'roi')
            copy_tree(ref_dir, roi_dir)

        if not os.path.isdir(os.path.join(roi_dir, 'convROI_BOLD')):
            os.mkdir(os.path.join(roi_dir, 'convROI_BOLD'))

        conv_dir = os.path.join(roi_dir, 'convROI_BOLD')
        list_out = []

        for roi in self.roi_list:
            list_out.append(os.path.join(conv_dir, 'conv' + roi[0] + roi[1] + '.nii'))

        return {"out_masks": list_out}, {}

    def _run_process(self):

        mask = self.mask[0]

        roi_dir = os.path.join(os.path.dirname(mask), 'roi')
        if not os.path.isdir(roi_dir):
            os.mkdir(roi_dir)
        conv_dir = os.path.join(roi_dir, 'convROI_BOLD')

        if not os.path.isdir(os.path.join(roi_dir, 'convROI_BOLD')):
            os.mkdir(os.path.join(roi_dir, 'convROI_BOLD'))

        # Resizing the mask to the size of the ROIs
        roi_1 = self.roi_list[0]
        roi_file = os.path.join(roi_dir, roi_1[0] + roi_1[1] + '.nii')
        roi_img = nib.load(roi_file)
        roi_data = roi_img.get_data()
        roi_size = roi_data.shape[:3]

        mask_thresh = threshold(mask, 0.5).get_data()
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


class Conv_ROI2(Process_mia):

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
        if not os.path.isdir(roi_dir):
            os.mkdir(roi_dir)
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


def threshold(file_name, thresh):
    img = nib.load(file_name)
    img_data = img.get_data()
    img_thresh = (img_data > thresh).astype(float)
    img_final = nib.Nifti1Image(img_thresh, img.affine, img.header)
    return img_final
