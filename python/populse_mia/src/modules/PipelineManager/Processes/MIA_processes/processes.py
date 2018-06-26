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

