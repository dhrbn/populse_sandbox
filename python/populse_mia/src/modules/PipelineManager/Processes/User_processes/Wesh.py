from capsul.api import Pipeline
import traits.api as traits


class Wesh(Pipeline):

    def pipeline_definition(self):
        # nodes
        self.add_process("realign1", "MIA_processes.SPM.spatial_preprocessing.Realign")
        self.nodes["realign1"].process.in_files = ['../../projects/test64/data/raw_data/Guerbet-C6-2014-Rat-K53-Tube20-2014-02-13_10-28-03-02-G1_Guerbet_Anat_modified_-RARE__pvm_-00-02-20.000.nii', '../../projects/test64/data/raw_data/Guerbet-C6-2014-Rat-K54-Tube21-2014-02-13_11-13-20-02-G1_Guerbet_Anat-RARE__pvm_-00-02-20.000.nii', '../../projects/test64/data/raw_data/Guerbet-C6-2014-Rat-K55-Tube22-2014-02-13_12-27-59-02-G1_Guerbet_Anat-RARE__pvm_-00-02-20.000.nii', '../../projects/test64/data/raw_data/sGuerbet-C6-2014-Rat-K53-Tube20-2014-02-13_10-28-03-02-G1_Guerbet_Anat_modified_-RARE__pvm_-00-02-20.000.nii', '../../projects/test64/data/raw_data/sGuerbet-C6-2014-Rat-K54-Tube21-2014-02-13_11-13-20-02-G1_Guerbet_Anat-RARE__pvm_-00-02-20.000.nii', '../../projects/test64/data/raw_data/sGuerbet-C6-2014-Rat-K55-Tube22-2014-02-13_12-27-59-02-G1_Guerbet_Anat-RARE__pvm_-00-02-20.000.nii']
        self.add_process("smooth1", "MIA_processes.SPM.spatial_preprocessing.Smooth")
        self.nodes["smooth1"].process.in_files = traits.Undefined

        # links
        self.add_link("realign1.realigned_files->smooth1.in_files")

        # nodes positions
        self.node_position = {
            "smooth1": (136.0, 100.0),
            "realign1": (-227.0, -159.0),
        }

        self.do_autoexport_nodes_parameters = False
