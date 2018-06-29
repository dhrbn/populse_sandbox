# Trait import
from nipype.interfaces.base import traits

# MIA import
from PipelineManager.Process_mia import Process_mia


class Duplicate_File(Process_mia):
    """
    From a file name, generating a list containing this file name and the file name itself.
    """

    def __init__(self):
        super(Duplicate_File, self).__init__()

        # Inputs
        self.add_trait("file1", traits.File(output=False))

        # Outputs
        self.add_trait("out_file1", traits.File(output=True))
        self.add_trait("out_file2", traits.File(output=True))

    def list_outputs(self):
        print('#' * 100)
        print('#' * 100)
        print('FILE_NAME DUPLICATE FILE', self.file1)
        print('#' * 100)
        print('#' * 100)
        if not self.file1 or self.file1 in ["<undefined>", traits.Undefined]:
            return {}, {}
        return {"out_file1": self.file1, "out_file2": self.file1}, {}

    def _run_process(self):
        self.out_file1 = self.file1
        self.out_file2 = self.file1


class Find_In_List(Process_mia):
    """
    From a list, returning the element that contains the pattern.
    """

    def __init__(self):
        super(Find_In_List, self).__init__()

        # Inputs
        self.add_trait("in_list", traits.List(output=False))
        self.add_trait("pattern", traits.String("0001", output=False, optional=True))

        # Outputs
        self.add_trait("out_file", traits.File(output=True))

    def list_outputs(self):
        if not self.in_list or self.in_list in ["<undefined>", traits.Undefined]:
            return {}, {}
        for file in self.in_list:
            if self.pattern in file:
                return {'out_file': file}, {}

    def _run_process(self):
        out_file = None
        for file in self.in_list:
            if self.pattern in file:
                out_file = file
                break

        self.out_file = out_file


class Files_To_List(Process_mia):
    """
    From several file names, generating a list containing all theses file names.
    """

    def __init__(self):
        super(Files_To_List, self).__init__()

        # Inputs
        self.add_trait("file1", traits.File(output=False))
        self.add_trait("file2", traits.File(output=False))

        # Outputs
        self.add_trait("file_list", traits.List(output=True))

    def list_outputs(self):
        out_list = [self.file1, self.file2]
        return {'file_list': out_list}, {}

    def _run_process(self):
        out_list = [self.file1, self.file2]
        self.file_list = out_list


class List_Duplicate(Process_mia):
    """
    From a file name, generating a list containing this file name and the file name itself.
    """

    def __init__(self):
        super(List_Duplicate, self).__init__()

        # Inputs
        self.add_trait("file_name", traits.File(output=False))

        # Outputs
        self.add_trait("out_file", traits.File(output=True))
        self.add_trait("out_list", traits.List(output=True))

    def list_outputs(self):
        return {"out_list": [self.file_name], "out_file": self.file_name}, {}

    def _run_process(self):
        self.out_list = [self.file_name]
        self.out_file = self.file_name


class ROI_List_Generator(Process_mia):

    def __init__(self):
        super(ROI_List_Generator, self).__init__()

        list_pos = ['ROI_OCC', 'ROI_PAR', 'ROI_TEMP', 'ROI_INSULA', 'ROI_FRON', 'ROI_CING', 'ROI_THA',
                    'ROI_STR', 'ACP', 'ACA', 'ACM', 'PICA', 'SCA']

        # Inputs
        self.add_trait("pos", traits.List(list_pos, output=False, optional=True))
        self.add_trait("hemi", traits.List(['_L', '_R'], output=False, optional=True))

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


"""
class Populse_Filter(Process_mia):

    def __init__(self, scans_list):
        super(Populse_Filter, self).__init__()

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
"""


