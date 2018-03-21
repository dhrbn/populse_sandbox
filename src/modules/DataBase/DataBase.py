from sqlalchemy import create_engine, cast, VARCHAR
from sqlalchemy.orm import sessionmaker
from sqlalchemy.interfaces import PoolListener
from DataBase.DataBaseModel import Tag, Scan, Value, Base, createDatabase
import os
import tempfile
import yaml
from SoftwareProperties.Config import Config
from DataBase.DataBaseModel import TAG_TYPE_STRING, TAG_ORIGIN_USER, TAG_ORIGIN_RAW
import pickle

class ForeignKeysListener(PoolListener):
    """
    Class to activate the pragma case_sensitive_like, that makes the like and contains functions case sensitive
    """
    def connect(self, dbapi_con, con_record):
        db_cursor = dbapi_con.execute('pragma case_sensitive_like=ON')

class DataBase:

    def __init__(self, project_root_folder, new_project):
        """ Database constructor
            New instance when switching project
            :param project_root_folder: projet root folder
                   If None, a temporary folder will be created

            :param new_project: Boolean to know if we have to create the database
                                True when New Project pop up or at the beginning for unnamed project
                                False when Open Project or Save as pop up

            self.isTempProject: To know if it's still unnamed project or not
            self.folder: Project root folder, relative
            self.properties: Properties of the project: Name, date of creation, sorted tag, and sort order (Not for Unnamed project)
            self.unsavedModifications: To know if there are unsaved modifications
            self.undos: Stack of undo actions we can do
            self.redos: Stack of redo actions we can do

            Memory approximation: the database file takes approximately 26 000 octets (1 bytes, 8 bits) per scan

        """

        # We don't have a project root folder at the opening of the software (Unnamed project), we generate a temporary folder
        if(project_root_folder == None):
            self.isTempProject = True
            self.folder = os.path.relpath(tempfile.mkdtemp())
        # We have a project root folder(New, Open, Save As)
        else:
            self.isTempProject = False
            self.folder = project_root_folder
            self.properties = self.loadProperties()
        # We create the database if it does not exists yet (for Unnamed project and New project)
        if(new_project):
            createDatabase(self.folder)
        # We open the database
        engine = create_engine('sqlite:///' + os.path.join(self.folder, 'database', 'mia2.db'), listeners=[ForeignKeysListener()])
        Base.metadata.bind = engine
        # We create a session
        DBSession = sessionmaker(bind=engine)
        self.session = DBSession()
        # If it's a new project, we refresh the list of tags, in case of changes in MIA2 preferences
        if new_project:
            self.refreshTags()
        # Initialisation
        self.unsavedModifications = False
        self.undos = []
        self.redos = []

    """ FROM properties/properties.yml """

    def refreshTags(self):
        """ Refreshes the tag if the project is new
            Makes sense if the user just changed the list of default tags in MIA2 preferences
        """

        # Tags cleared
        tags = self.session.query(Tag).filter().all()
        for tag in tags:
            self.session.delete(tag)

        # New tags added
        config = Config()
        if config.getDefaultTags() != None:
            for default_tag in config.getDefaultTags():
                if not self.hasTag(default_tag):
                    # Tags by default set as visible
                    self.addTag(default_tag, True, TAG_ORIGIN_USER, TAG_TYPE_STRING, None, None, None)
                    self.saveModifications()

        # FileName as raw tag
        self.setTagOrigin("FileName", TAG_ORIGIN_RAW)

    def loadProperties(self):
        """ Loads the properties file (Unnamed project does not have this file) """
        with open(os.path.join(self.folder, 'properties', 'properties.yml'), 'r') as stream:
            try:
                return yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)

    def getName(self):
        """ Returns the name of the project if it's not Unnamed project, otherwise empty string """
        if(self.isTempProject):
            return ""
        else:
            return self.properties["name"]

    def setName(self, name):
        """ Sets the name of the project if it's not Unnamed project, otherwise does nothing
            :param name: new name of the project
        """
        if not self.isTempProject:
            self.properties["name"] = name
            self.saveConfig()

    def saveConfig(self):
        """ Save the changes in the properties file """
        with open(os.path.join(self.folder, 'properties', 'properties.yml'), 'w', encoding='utf8') as configfile:
            yaml.dump(self.properties, configfile, default_flow_style=False, allow_unicode=True)

    def getDate(self):
        """ Returns the date of creation of the project if it's not Unnamed project, otherwise empty string """
        if(self.isTempProject):
            return ""
        else:
            return self.properties["date"]

    def getSortedTag(self):
        """ Returns the sorted tag of the project if it's not Unnamed project, otherwise empty string """
        if (self.isTempProject):
            return ""
        else:
            return self.properties["sorted_tag"]

    def setSortedTag(self, tag):
        """ Sets the sorted tag of the project if it's not Unnamed project, otherwise does nothing
            :param tag: new sorted tag of the project
        """
        if not self.isTempProject:
            self.properties["sorted_tag"] = tag
            self.unsavedModifications = True

    def getSortOrder(self):
        """ Returns the sort order of the project if it's not Unnamed project, otherwise empty string """
        if (self.isTempProject):
            return ""
        else:
            return self.properties["sort_order"]

    def setSortOrder(self, order):
        """ Sets the sort order of the project if it's not Unnamed project, otherwise does nothing
            :param order: new sort order of the project (ascending or descending)
        """
        if not self.isTempProject:
            self.properties["sort_order"] = order
            self.unsavedModifications = True


    """ FROM SQlite Database """

    def addScan(self, name, checksum):
        """
        To add a new scan to the project
        :param name: FileName of the scan
        :param checksum: MD5 checksum of the scan file
        """
        scan = Scan(scan=name, checksum=checksum)
        self.session.add(scan)
        self.unsavedModifications = True

    def addTag(self, tag, visible, origin, type, unit, default, description):
        """
        To add a new tag to the project
        :param tag: Tag name
        :param visible: True or False
        :param origin: Raw or user
        :param type: String, Integer, Float, or List
        :param unit:
        :param default:
        :param description:
        """
        tag = Tag(tag=tag, visible=visible, origin=origin, type=type, unit=unit, default=default, description=description)
        self.session.add(tag)
        self.unsavedModifications = True

    def addValue(self, scan, tag, current_value, raw_value):
        """
        To add a new value (a new cell in the databrowser)
        :param scan: FileName of the scan
        :param tag: Tag name
        :param current_value: Current value displayed in the DataBrowser (Raw_value for raw tags and default_value for user tags)
        :param raw_value: Raw value of the cell (None of user tags)
        """
        value = Value(scan=scan, tag=tag, current_value=current_value, raw_value=raw_value)
        self.session.add(value)
        self.unsavedModifications = True

    def getScans(self):
        """
        To get all Scan objects
        :return: All Scan objects
        """
        scans = self.session.query(Scan).filter().all()
        return scans

    def getValues(self):
        """
        To get all values
        :return: All Value objects
        """
        values = self.session.query(Value).filter().all()
        return values

    def getScan(self, scan):
        """
        To get the Scan object of the scan
        :param scan: FileName of the scan
        :return: The Scan object of the scan
        """
        scans = self.session.query(Scan).filter(Scan.scan == scan).all()
        if len(scans) == 1:
            return scans[0]

    def getValuesGivenTag(self, tag):
        """
        To get all the values of the tag
        :param tag: Tag name
        :return: A list of Value objects for the tag
        """
        values = self.session.query(Value).filter(Value.tag == tag).all()
        return values

    def getValuesGivenScan(self, scan):
        """
        To get all the values of the scan
        :param scan: FileName of the scan
        :return: A list of Value objects for the scan
        """
        values = self.session.query(Value).filter(Value.scan == scan).all()
        return values

    def getValue(self, scan, tag):
        """
        To get the Value object of the cell <Scan, Tag>
        :param scan: FileName of the scan
        :param tag: Tag name
        :return: The Value object of the cell
        """
        values = self.session.query(Value).filter(Value.tag == tag).filter(Value.scan == scan).all()
        if len(values) == 1:
            return values[0]

    def scanHasTag(self, scan, tag):
        """
        To know if the scan has a value for the tag
        :param scan: FileName of the scan
        :param tag: Tag name
        :return: True if the scan has a value for the tag, False otherwise
        """
        values = self.session.query(Value).filter(Value.tag == tag).filter(Value.scan == scan).all()
        return len(values) == 1

    def hasScan(self, scan):
        """
        To know if the scan exists
        :param scan: FileName of the scan
        :return: True if the scan exists, False otherwise
        """
        scans = self.session.query(Scan).filter(Scan.scan == scan).all()
        return len(scans) == 1

    def getTags(self):
        """
        Gives a list of the tags of the project (Tag objects)
        :return: A list of Tag objects of the project
        """
        tags = self.session.query(Tag).filter().all()
        return tags

    def getVisualizedTags(self):
        """
        Gives a list of the visualized tags of the project
        :return: A list of the visualized tags of the project
        """
        tags = self.session.query(Tag).filter(Tag.visible == True).all()
        return tags

    def hasTag(self, tag):
        """
        To know if the tag is in the project
        :param tag: Tag name
        :return: True if the tag is in the project, False otherwise
        """
        tags = self.session.query(Tag).filter(Tag.tag == tag).all()
        return len(tags) == 1

    def getTagOrigin(self, tag):
        """
        Gives the origin of the tag
        :param tag: Tag name
        :return: The origin of the tag (raw or user)
        """
        return self.getTag(tag).origin

    def getTagVisibility(self, tag):
        """
        Gives the visibility of the tag
        :param tag: Tag name
        :return: The visibility of the tag (True or False)
        """
        return self.getTag(tag).visible

    def getTagDefault(self, tag):
        """
        Gives the default value of the tag
        :param tag: Tag name
        :return: Default value of the tag
        """
        return self.getTag(tag).default

    def getTagUnit(self, tag):
        """
        Gives the unit of the tag
        :param tag: Tag name
        :return: The unit of the tag
        """
        return self.getTag(tag).unit

    def getTagDescription(self, tag):
        """
        Gives the description of the tag
        :param tag: Tag name
        :return: The description of the tag (tooltip displayed when putting the mouse on the headers in the DataBrowser)
        """
        return self.getTag(tag).description

    def getTagType(self, tag):
        """
        Gives the type of the tag
        :param tag: Tag name
        :return: The type of the tag (String, Integer, Float, or List)
        """
        return self.getTag(tag).type

    def getTag(self, tag):
        """
        Gives the Tag object of the tag
        :param tag: Tag name
        :return: The Tag object associated to the tag
        """
        tags = self.session.query(Tag).filter(Tag.tag == tag).all()
        if len(tags) == 1:
            return tags[0]

    def getUserTags(self):
        """
        Gives all the user tags of the project
        :return: The list of user tags of the project
        """
        tags = self.session.query(Tag).filter(Tag.origin == TAG_ORIGIN_USER).all()
        return tags

    def setTagVisibility(self, name, visibility):
        """
        Sets the visibility of the tag
        :param name: Tag name
        :param visibility: New visibility (True or False)
        """
        tag = self.getTag(name)
        tag.visible = visibility
        self.unsavedModifications = True

    def setTagOrigin(self, name, origin):
        """
        Sets the origin of the tag
        :param name: Tag name
        :param origin: New origin of the tag (raw or user)
        """
        tag = self.getTag(name)
        tag.origin = origin
        self.unsavedModifications = True

    def setTagDescription(self, name, description):
        """
        Sets the description of the tag
        :param name: Tag name
        :param description: New description of the tag
        """
        tag = self.getTag(name)
        tag.description = description
        self.unsavedModifications = True

    def setTagUnit(self, name, unit):
        """
        Sets the unit of the tag
        :param name: Tag name
        :param unit: New unit of the tag
        """
        tag = self.getTag(name)
        tag.unit = unit
        self.unsavedModifications = True

    def setTagType(self, name, type):
        """
        Sets the type of the tag
        :param name: Tag name
        :param unit: New type of the tag
        """
        tag = self.getTag(name)
        tag.type = type
        self.unsavedModifications = True

    def resetAllVisibilities(self):
        """
        Puts the visibility of all tags of the project to False
        """
        tags = self.session.query(Tag).filter().all()
        for tag in tags:
            tag.visible = False
        self.unsavedModifications = True

    def setTagValue(self, scan, tag, new_value):
        """ Sets the value of the cell asked
            :param scan: The FileName of the scan to reset
            :param tag: The tag name to reset
            :param new_value: New value of the cell
        """
        # We only change the cell if the tag is not FileName
        if not tag == "FileName":
            tags = self.session.query(Value).filter(Value.scan==scan).filter(Value.tag==tag).all()
            # There is already a value
            if len(tags) == 1:
                tag = tags[0]
                tag.current_value = new_value
            self.unsavedModifications = True

    def resetTag(self, scan, tag):
        """
        Resets the value of the cell asked, only done on raw tags, does not make sense on user tags
        :param scan: The FileName of the scan to reset
        :param tag: The tag name to reset
        """
        tag = self.getValue(scan, tag)
        if(tag.raw_value != None):
            tag.current_value = tag.raw_value
            self.unsavedModifications = True
            return True
        else:
            return False

    def removeScan(self, scan):
        """
        Removes a scan from the project (corresponds to a row in the databrowser)
        Removes the scan from the list of scans (Scan table), and all the values corresponding to this scan (Value table)
        :param scan: FileName of the scan to remove
        """
        # All the values of the scan are removed
        tags = self.session.query(Value).filter(Value.scan == scan).all()
        for tag in tags:
            self.session.delete(tag)
        # The scan is removed from the list of scans
        scans = self.session.query(Scan).filter(Scan.scan == scan).all()
        if len(scans) == 1:
            self.session.delete(scans[0])
            self.unsavedModifications = True

    def removeTag(self, tag):
        """
        Removes a tag from the project (corresponds to a column in the databrowser)
        We can only remove user tags from the software
        :param tag: Name of the tag to remove
        """
        # Values associated to the tags removed
        values = self.session.query(Value).filter(Value.tag == tag).all()
        for value in values:
            self.session.delete(value)
        # Tag removed
        tags = self.session.query(Tag).filter(Tag.tag == tag).all()
        # TODO return error if len(tags) != 1
        self.session.delete(tags[0])
        self.unsavedModifications = True

    def removeValue(self, scan, tag):
        """
        Removes the value of the tuple <scan, tag> (corresponds to a cell in the databrowser)
        :param scan: FileName of the scan
        :param tag: Name of the tag
        """
        value = self.getValue(scan, tag)
        self.session.delete(value)
        self.unsavedModifications = True

    def saveModifications(self):
        """
        Saves the pending operations of the project (actions still not saved)
        """
        self.session.commit()
        if not self.isTempProject:
            self.saveConfig()
        self.unsavedModifications = False

    def unsaveModifications(self):
        """
        Unsaves the pending operations of the project (actions still not saved)
        """
        self.session.rollback()
        self.unsavedModifications = False

    def getScansMissingTags(self):
        """
        Checks for scans that have missing values in the visualized tags
        :return: The list of scans with missing values in the visualized tags
        """
        return_list = []
        for scan in self.getScans():
            for tag in self.getVisualizedTags():
                if not self.scanHasTag(scan.scan, tag.tag) and not scan.scan in return_list:
                    return_list.append(scan.scan)
        return return_list

    def getScansSimpleSearch(self, search):
        """
        Executes the rapid search (search bar) and returns the list of scans that contain the search in their tag values

        :param search: string corresponding to the search
        :return: The list of scans containing the search
        """

        # We take all the scans that contain the search in at least one of their visible tag values
        values = self.session.query(Value.scan).filter(Value.current_value.like("%" + search + "%"), Value.tag == Tag.tag, Tag.visible == True).distinct().all()

        scans = [] # List of scans to return
        for value in values:
            # We add the scan to the list
            scans.append(value.scan)
        return scans

    def getScansAdvancedSearch(self, links, fields, conditions, values, nots):
        """
        Executes the advanced search and returns the list of scans corresponding to the criterias

        :param links: list of links (AND/OR)
        :param fields: list of fields (tag name or all visualized tags)
        :param conditions: list of conditions (=, !=, <=, >=, <, >, IN, BETWEEN, CONTAINS)
        :param values: list of values (tag values)
        :param nots: list of nots (NOT?)
        :return: list of scans corresponding to those criterias
        """
        result = [] # list of scans to return
        queries = [] # list of scans of each query (row)
        i = 0
        while i < len(conditions):
            rowResult = self.session.query(Value.scan) # Every scan to start with
            if(fields[i] != "All visualized tags"):
                # Filter on the tag if necessary
                rowResult = rowResult.filter(Value.tag == fields[i])
            # Filter on the condition on the tag value (condition type + value)
            if(conditions[i] == "="):
                rowResult = rowResult.filter(Value.current_value == values[i])
            elif(conditions[i] == "!="):
                rowResult = rowResult.filter(Value.current_value != values[i])
            elif(conditions[i] == ">="):
                rowResult = rowResult.filter(Value.current_value >= values[i])
            elif (conditions[i] == "<="):
                rowResult = rowResult.filter(Value.current_value <= values[i])
            elif (conditions[i] == ">"):
                rowResult = rowResult.filter(Value.current_value > values[i])
            elif (conditions[i] == "<"):
                rowResult = rowResult.filter(Value.current_value < values[i])
            elif (conditions[i] == "CONTAINS"):
                rowResult = rowResult.filter(Value.current_value.contains(values[i]))
            elif (conditions[i] == "IN"):
                rowResult = rowResult.filter(Value.current_value.in_(values[i].split(', ')))
            elif (conditions[i] == "BETWEEN"):
                borders = values[i].split(', ')
                rowResult = rowResult.filter(Value.current_value.between(borders[0], borders[1]))
            if(nots[i] == "NOT"):
                # If NOT, we take the opposite: All scans MINUS the result of the query
                rowResult = self.session.query(Value.scan).except_(rowResult)
            # We add the list of scans to the subresults
            queries.append(rowResult)
            i = i + 1
        # We start with the first row to put the link between the conditions
        # Links are made row by row, there is no priority like in SQL where AND is stronger than OR
        finalQuery = queries[0]
        i = 0
        while i < len(links):
            if(links[i] == "AND"):
                # If the link is AND, we do an intersection between the current result and the next row
                finalQuery = finalQuery.intersect(queries[i + 1])
            else:
                # If the link is OR, we do an union between the current result and the next row
                finalQuery = finalQuery.union(queries[i + 1])
            i = i + 1
        # We take the list of all the scans that respect every conditions
        finalQuery = finalQuery.distinct().all()
        # We create the return list with the name of the scans (FileName)
        for value in finalQuery:
            result.append(value.scan)
        return result

    def check_count_table(self, values):
        """
        Checks if a scan contains all the couples <tag, value> given in parameter
        :param values: List of couple <tag, value> to check
        :return: True if a scan has all the couples <tag, value> given in parameter, False otherwise
        """
        queries = []
        for value in values:
            tag = value[0]
            value = value[1]
            cellResult = self.session.query(Value.scan)  # Every scan to start with
            cellResult = cellResult.filter(Value.tag == tag).filter(Value.current_value == value) # We apply the tag and value filter
            queries.append(cellResult)
        # We put all the cells together, with intersections
        # We only want a scan that has all the values
        finalResult = queries[0]
        i = 0
        while i < len(queries) - 1:
            finalResult = finalResult.intersect(queries[i + 1])
            i = i + 1
        finalResult = finalResult.distinct().all()
        # If there is a scan, we return True, and False otherwise
        return finalResult

    def undo(self):
        """
        Undo the last action made by the user on the project
        """

        # We can undo if we have an action to revert
        if len(self.undos) > 0:
            toUndo = self.undos.pop()
            self.redos.append(toUndo) # We pop the undo action in the redo stack
            # The first element of the list is the type of action made by the user (add_tag, remove_tags, add_scans, remove_scans, or modified_values)
            action = toUndo[0]
            if(action == "add_tag"):
                # For removing the tag added, we just have to memorize the tag name, and remove it
                tagToRemove = toUndo[1]
                self.removeTag(tagToRemove)
            if (action == "remove_tags"):
                # To reput the removed tags, we need to reput the tag in the tag list, and all the tags values associated to this tag
                tagsRemoved = toUndo[1] # The second element is a list of the removed tags (Tag class)
                i = 0
                while i < len(tagsRemoved):
                    # We reput each tag in the tag list, keeping all the tags params
                    tagToReput = tagsRemoved[i]
                    self.addTag(tagToReput.tag, tagToReput.visible, tagToReput.origin, tagToReput.type, tagToReput.unit, tagToReput.default, tagToReput.description)
                    i += 1
                valuesRemoved = toUndo[2] # The third element is a list of tags values (Value class)
                self.reput_values(valuesRemoved)
            if (action == "add_scans"):
                # To remove added scans, we just need their file name
                scansAdded = toUndo[1] # The second element is a list of added scans to remove
                i = 0
                while i < len(scansAdded):
                    # We remove each scan added
                    scanToRemove = scansAdded[i][0]
                    self.removeScan(scanToRemove)
                    i += 1
            if(action == "remove_scans"):
                # To reput a removed scan, we need the scans names, and all the values associated
                scansRemoved = toUndo[1] # The second element is the list of removed scans (Scan class)
                i = 0
                while i < len(scansRemoved):
                    # We reput each scan, keeping the same values
                    scanToReput = scansRemoved[i]
                    self.addScan(scanToReput.scan, scanToReput.checksum)
                    i += 1
                valuesRemoved = toUndo[2] # The third element is the list of removed values (Value class)
                self.reput_values(valuesRemoved)

            if (action == "modified_values"):
                # To revert a value changed in the databrowser, we need two things: the cell (scan and tag, and the old value)
                modifiedValues = toUndo[1] # The second element is a list of modified values (reset, or value changed)
                i = 0
                while i < len(modifiedValues):
                    # Each modified value is a list of 3 elements: scan, tag, and old_value
                    valueToRestore = modifiedValues[i]
                    scan = valueToRestore[0]
                    tag = valueToRestore[1]
                    old_value = valueToRestore[2]
                    if(old_value == None):
                        # If the cell was not defined before, we reput it
                        self.removeValue(scan, tag)
                    else:
                        # If the cell was there before, we just set it to the old value
                        self.setTagValue(scan, tag, str(old_value))
                    i += 1
            if (action == "modified_visibilities"):
                # To revert the modifications of the visualized tags
                visibles = toUndo[1] # List of the tags visibles before the modification (Tag objects)
                self.resetAllVisibilities() # Reset of the visibilities
                for visible in visibles:
                    # We reput each old tag visible
                    self.setTagVisibility(visible.tag, True)
            if (action == "modified_sort"):
                # To revert a sort change
                old_sorted_tag = toUndo[1]
                old_sort_order = toUndo[2]
                self.update_sort(old_sorted_tag, old_sort_order)

        #print(len(pickle.dumps(self.history, -1))) # Memory approximation in number of bits

        """
        Approximate results
        
        - Changing 1 cell value: 180 bits
        - Changing 10 cell values: 1350 bits
        - Removing 1 scan: 20000 bits
        - Removing 10 scans: 200000 bits
        - Removing 1 tag: From 300 bits to 4000 bits => Depends on the number of values defined
        - Adding 1 tag: 40 bits
        => By storing the minimum of data to be able to revert the actions, we take way less memory than by storing the whole database after each action
        """

    def update_sort(self, tag, order):
        self.setSortedTag(tag)
        self.setSortOrder(order)

    def reput_values(self, values):
        """
        To reput the Value objects in the database
        :param values: List of Value objects
        """
        i = 0
        while i < len(values):
            # We reput each value, exactly the same as it was before
            valueToReput = values[i]
            self.addValue(valueToReput.scan, valueToReput.tag, valueToReput.current_value, valueToReput.raw_value)
            i += 1

    def redo(self):
        """
        Redo the last action made by the user on the project
        """
        # We can redo if we have an action to make again
        if len(self.redos) > 0:
            toRedo = self.redos.pop()
            self.undos.append(toRedo)  # We pop the redo action in the undo stack
            # The first element of the list is the type of action made by the user (add_tag, remove_tags, add_scans, remove_scans, or modified_values)
            action = toRedo[0]
            if(action == "add_tag"):
                # For adding the tag, we need the tag name, and all its attributes
                tagToAdd = toRedo[1]
                tagType = toRedo[2]
                tagUnit = toRedo[3]
                tagDefaultValue = toRedo[4]
                tagDescription = toRedo[5]
                values = toRedo[6] # List of values stored
                # Adding the tag
                self.addTag(tagToAdd, True, TAG_ORIGIN_USER, tagType, tagUnit, tagDefaultValue, tagDescription)
                # Adding all the values associated
                for value in values:
                    self.addValue(value[0], value[1], value[2], value[3])
            if (action == "remove_tags"):
                # To remove the tags, we need the names
                tagsRemoved = toRedo[1]  # The second element is a list of the removed tags (Tag class)
                i = 0
                while i < len(tagsRemoved):
                    # We reput each tag in the tag list, keeping all the tags params
                    tagToRemove = tagsRemoved[i].tag
                    self.removeTag(tagToRemove)
                    i += 1
            if (action == "add_scans"):
                # To add the scans, we need the FileNames and the values associated to the scans
                scansAdded = toRedo[1]  # The second element is a list of the scans to add
                # We add all the scans
                i = 0
                while i < len(scansAdded):
                    # We remove each scan added
                    scanToAdd = scansAdded[i]
                    self.addScan(scanToAdd[0], scanToAdd[1])
                    i += 1
                # We add all the values
                i = 0
                valuesAdded = toRedo[2] # The third element is a list of the values to add
                while i < len(valuesAdded):
                    valueToAdd = valuesAdded[i]
                    self.addValue(valueToAdd[0], valueToAdd[1], valueToAdd[2], valueToAdd[2])
                    i += 1
            if(action == "remove_scans"):
                # To remove a scan, we only need the FileName of the scan
                scansRemoved = toRedo[1]  # The second element is the list of removed scans (Scan class)
                i = 0
                while i < len(scansRemoved):
                    # We reput each scan, keeping the same values
                    scanToRemove = scansRemoved[i].scan
                    self.removeScan(scanToRemove)
                    i += 1
            if (action == "modified_values"):
                # To modily the values, we need the cells, and the updated values
                modifiedValues = toRedo[1]  # The second element is a list of modified values (reset, or value changed)
                i = 0
                while i < len(modifiedValues):
                    # Each modified value is a list of 3 elements: scan, tag, and old_value
                    valueToRestore = modifiedValues[i]
                    scan = valueToRestore[0]
                    tag = valueToRestore[1]
                    # valueToRestore[2] is the old value of the cell
                    new_value = valueToRestore[3]
                    self.setTagValue(scan, tag, new_value)
                    i += 1
            if (action == "modified_visibilities"):
                # To revert the modifications of the visualized tags
                visibles = toRedo[2]  # List of the tags visibles after the modification (Tag objects)
                self.resetAllVisibilities()  # Reset of the visibilities
                for visible in visibles:
                    # We reput each new tag visible
                    self.setTagVisibility(visible, True)
            if (action == "modified_sort"):
                # To revert a sort change
                # toUndo[1] and toUndo[2] are old values
                new_sorted_tag = toRedo[3]
                new_sort_order = toRedo[4]
                self.update_sort(new_sorted_tag, new_sort_order)

        #print(len(pickle.dumps(self.history, -1))) # Memory approximation in number of bits