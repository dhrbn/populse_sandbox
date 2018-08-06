from PyQt5 import QtWidgets, QtCore
from SoftwareProperties.Config import Config
from PopUps.Ui_Filter_Selection import Ui_Filter_Selection

class Ui_Select_Filter(Ui_Filter_Selection):
    """
    Is called when the user wants to open a filter saved before
    """

    def __init__(self, project, databrowser):
        super(Ui_Select_Filter, self).__init__(project)
        self.project = project
        self.databrowser = databrowser
        self.config = Config()
        self.setWindowTitle("Open a filter")

        # Filling the filter list
        for filter in self.project.filters:
            item = QtWidgets.QListWidgetItem()
            self.list_widget_filters.addItem(item)
            item.setText(filter.name)

    def ok_clicked(self):
        for item in self.list_widget_filters.selectedItems():

            # Current filter updated
            filter_name = item.text()
            filter_object = self.project.getFilter(filter_name)
            self.project.setCurrentFilter(filter_object)
            break

        self.databrowser.open_filter_infos()

        self.accept()
        self.close()
