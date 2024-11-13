#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, flush=True, **kwargs)

import os
from datetime import datetime

from sqlalchemy import Column, ForeignKey, Integer, Float, String, DateTime, Text, PickleType, LargeBinary
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import desc
from sqlalchemy import inspect

Base = declarative_base()

debug = False

#### Define the database tables as classes
class LibraryRecord(Base):
    __tablename__ = 'library_record'
    library_record_id = Column(Integer, primary_key=True)
    status = Column(String(25), nullable=False)
    id_name = Column(String(25), nullable=False)
    source = Column(String(100), nullable=True)
    species = Column(String(255), nullable=True)
    keywords = Column(String(100), nullable=True)
    version_tag = Column(String(255), nullable=False)
    original_filename = Column(String(255), nullable=False)
    original_md5_checksum = Column(String(50), nullable=True)
    local_filename = Column(String(255), nullable=True)
    local_md5_checksum = Column(String(50), nullable=True)
    converted_to_mzSpecLib = Column(String(5), nullable=True)
    metadata_quality_score = Column(String(25), nullable=True)
    QC_score = Column(String(25), nullable=True)
    QC_report_url = Column(String(255), nullable=True)
    title = Column(String(255), nullable=True)
    description = Column(Text, nullable=True)
    submitter_full_name = Column(String(255), nullable=True)
    submitter_email = Column(String(255), nullable=True)
    submitter_affiliation = Column(String(255), nullable=True)
    lab_head_full_name = Column(String(255), nullable=True)
    lab_head_email = Column(String(255), nullable=True)
    lab_head_affiliation = Column(String(255), nullable=True)
    library_type = Column(String(255), nullable=True)
    library_building_software = Column(String(255), nullable=True)
    library_building_protocol = Column(Text, nullable=True)
    instruments = Column(String(255), nullable=True)
    fragmentation_type = Column(String(255), nullable=True)
    mass_modifications = Column(String(255), nullable=True)
    intended_workflow = Column(String(255), nullable=True)
    intended_sample_type = Column(String(255), nullable=True)
    publication = Column(String(255), nullable=True)
    documentation_url = Column(String(255), nullable=True)
    source_url = Column(String(255), nullable=True)
    provenance_information = Column(String(255), nullable=True)
    n_entries = Column(Integer, nullable=True)
    record_created_datetime = Column(DateTime, nullable=False)
    record_human_updated_datetime = Column(DateTime, nullable=True)
    record_automation_updated_datetime = Column(DateTime, nullable=True)
    changelog_comments = Column(Text, nullable=True)


class SpectrumLibraryCollection:
    """
    SpectrumLibraryCollection - Class for a collection of spectrum libraries

    Attributes
    ----------
    filename : string
        Filename of the SQLite database file that contains information about the collection of libraries available.

    Methods
    -------
    create - Create a new spectral library collection
    show - Return a string that summarizes the state of the collection
    get_libraries - Return a list of available libraries
    get_library - Return attributes of a specific library
    add_library - Add a new library
    create_index - Create a master index from all the constituent library indexes to be able to find spectra in any library
    find_spectra - Return a list of spectra given query constraints

    """


    def __init__(self, filename=None):
        """
        __init__ - SpectrumLibraryCollection constructor

        Parameters
        ----------
        filename : string
            Filename of the SQLite database file that contains information about the collection of libraries available.

        """

        self.filename = filename
        if os.path.exists(self.filename):
            if debug: eprint(f"DEBUG: Library collection {self.filename} exists")
            self.connect()
        else:
            if debug: eprint(f"DEBUG: Library collection {self.filename} not found. Will create new.")
            self.createDatabase()


    #### Destructor
    def __del__(self):
        self.disconnect()


    #### Define getter/setter for attribute filename
    @property
    def filename(self):
        return(self._filename)
    @filename.setter
    def filename(self, filename):
        self._filename = filename


    #### Define attribute session
    @property
    def session(self) -> str:
        return self._session

    @session.setter
    def session(self, session: str):
        self._session = session


    #### Define attribute engine
    @property
    def engine(self) -> str:
        return self._engine

    @engine.setter
    def engine(self, engine: str):
        self._engine = engine



    #### Delete and create the database. Careful!
    def createDatabase(self):
        if os.path.exists(self.filename):
            eprint("INFO: Deleting previous database file " + self.filename)
            os.remove(self.filename)
        eprint("INFO: Creating database " + self.filename)
        engine = create_engine("sqlite:///"+self.filename)
        Base.metadata.create_all(engine)
        self.connect()



    #### Create and store a database connection
    def connect(self):
        if debug:
            eprint("DEBUG: Connecting to database " + self.filename)
        engine = create_engine("sqlite:///"+self.filename)
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        self.session = session
        self.engine = engine



    #### Destroy the database connection
    def disconnect(self):
        if debug: eprint("DEBUG: Disconnecting from database " + self.filename)
        session = self.session
        engine = self.engine
        session.close()
        engine.dispose()



    def create(self, overwrite_existing=False):
        """
        create - Create a new spectral library collection

        Extended description of function.

        Parameters
        ----------
        overwrite_existing : boolean
            Set to true in order to write over the previous file if it exists

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        if debug:
            eprint("DEBUG: Creating database " + self.filename)
        if os.path.exists(self.filename):
            os.remove(self.filename)
        engine = create_engine("sqlite:///"+self.filename)
        Base.metadata.create_all(engine)
        self.connect()
        return()



    def get_libraries(self):
        """
        get_libraries - Return a list of available libraries

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        if debug:
            eprint("DEBUG: Fetching all libraries")
        session = self.session
        libraries = session.query(LibraryRecord).all()
        return(libraries)



    def get_settable_library_attribute_names(self):
        """
        get_settable_library_attribute_names - Return the names of library attributes

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        list
            List of the available attributes
        """

        attribute_names = [ 'source',
                            'species',
                            'keywords',
                            'version_tag',
                            'original_filename',
                            'original_md5_checksum',
                            'local_filename',
                            'local_md5_checksum',
                            'converted_to_mzSpecLib',
                            'metadata_quality_score',
                            'QC_score',
                            'QC_report_url',
                            'title',
                            'description',
                            'submitter_full_name',
                            'submitter_email',
                            'submitter_affiliation',
                            'lab_head_full_name',
                            'lab_head_email',
                            'lab_head_affiliation',
                            'library_type',
                            'library_building_software',
                            'library_building_protocol',
                            'instruments',
                            'fragmentation_type',
                            'mass_modifications',
                            'intended_workflow',
                            'intended_sample_type',
                            'publication',
                            'documentation_url',
                            'source_url',
                            'provenance_information',
                            'n_entries',
                            'changelog_comments' ]

        return attribute_names


    def get_all_library_attribute_names(self):
        """
        get_all_library_attribute_names - Return the names of library attributes

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        list
            List of the available attributes
        """

        attribute_names = [ 'library_record_id', 'status', 'id_name' ]
        attribute_names.extend(self.get_settable_library_attribute_names())
        attribute_names.extend([ 'record_created_datetime', 'record_human_updated_datetime', 'record_automation_updated_datetime' ])

        return attribute_names



    def get_library(self,identifier=None,version_tag=None,filename=None):
        """
        get_library - Return attributes of a specific library

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        session = self.session
        if identifier is not None and identifier > "":
            libraries = session.query(LibraryRecord).filter(LibraryRecord.id_name==identifier).order_by(desc(LibraryRecord.version_tag)).all()
            if len(libraries) == 0:
                return
            else:
                libraryListStr = ""
                for library in libraries:
                    libraryListStr += library.version_tag+","
                    if version_tag is not None and version_tag > "":
                        if version_tag == library.version_tag:
                            return(library)
                if version_tag is not None and version_tag > "":
                    raise Exception(f"Unable to find version {version_tag} of library {identifier}")
            if len(libraries) == 1:
                return(library)
            raise Exception(f"There are several version of this library ({libraryListStr}). Please specify a version_tag")
            return()

        elif filename is not None and filename > "":
            raise Exception("Search by filename not implemented")
        else:
            raise Exception("Not enough information to find library")


    def add_library(self, attributes):
        """
        add_library - Add a new library

        Adds a new library to the library collection database.

        Parameters
        ----------
        attributes - Dict of library attributes that correspond to the table definition
        
        Returns
        -------
        None
        """

        if debug:
            eprint("DEBUG: Adding a library entry")

        #### Create a sanitized version of the attributes
        checked_attributes = attributes.copy()
        attribute_names = self.get_settable_library_attribute_names()
        for attribute_name in attribute_names:
            if attribute_name not in checked_attributes:
                checked_attributes[attribute_name] = None

        session = self.session
        library_record = LibraryRecord(
            status='Pending',
            id_name="PXL000000",
            source=checked_attributes['source'],
            species=checked_attributes['species'],
            keywords=checked_attributes['keywords'],
            version_tag=checked_attributes['version_tag'],
            original_filename=checked_attributes['original_filename'],
            original_md5_checksum=checked_attributes['original_md5_checksum'],
            local_filename=checked_attributes['local_filename'],
            local_md5_checksum=checked_attributes['local_md5_checksum'],
            converted_to_mzSpecLib=checked_attributes['converted_to_mzSpecLib'],
            metadata_quality_score=checked_attributes['metadata_quality_score'],
            QC_score=checked_attributes['QC_score'],
            QC_report_url=checked_attributes['QC_report_url'],
            title=checked_attributes['title'],
            description=checked_attributes['description'],
            submitter_full_name=checked_attributes['submitter_full_name'],
            submitter_email=checked_attributes['submitter_email'],
            submitter_affiliation=checked_attributes['submitter_affiliation'],
            lab_head_full_name=checked_attributes['lab_head_full_name'],
            lab_head_email=checked_attributes['lab_head_email'],
            lab_head_affiliation=checked_attributes['lab_head_affiliation'],
            library_type=checked_attributes['library_type'],
            library_building_software=checked_attributes['library_building_software'],
            library_building_protocol=checked_attributes['library_building_protocol'],
            instruments=checked_attributes['instruments'],
            fragmentation_type=checked_attributes['fragmentation_type'],
            mass_modifications=checked_attributes['mass_modifications'],
            intended_workflow=checked_attributes['intended_workflow'],
            intended_sample_type=checked_attributes['intended_sample_type'],
            publication=checked_attributes['publication'],
            documentation_url=checked_attributes['documentation_url'],
            source_url=checked_attributes['source_url'],
            provenance_information=checked_attributes['provenance_information'],
            n_entries=checked_attributes['n_entries'],
            record_created_datetime=datetime.now(),
            changelog_comments=checked_attributes['changelog_comments']
            )
        session.add(library_record)
        session.flush()

        #### Get the primary key identifier
        assert(library_record.library_record_id)
        idstr = str(library_record.library_record_id)
        if debug:
            eprint(f"DEBUG: Returned id={idstr}")
        idstr_length = len(idstr)
        assert(idstr_length)

        #### Create and store the PXL identifier
        padding = "000000"
        new_idstr = "PXL" + padding[0:len(padding)-idstr_length] + idstr
        library_record.id_name = new_idstr
        library_record.status = 'OK'
        session.flush()
        session.commit()
        if debug:
            eprint(f"DEBUG: Record for {new_idstr} created")
        return new_idstr


    def update_library_metadata(self, id, attributes, update_type, changelog_comments):
        """
        update_library_metadata - Update the metadata associated with a library

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        if debug: eprint(f"INFO: Updating the library entry with id {id}")
        session = self.session
        library = session.query(LibraryRecord).filter(LibraryRecord.library_record_id==id).first()
        if library is not None:
            made_change = False
            result = ''
            if changelog_comments is None:
                changelog_comments = ''
            else:
                changelog_comments += ', '
            attribute_names = self.get_settable_library_attribute_names()
            for attribute_name in attribute_names:
                if attribute_name in attributes:
                    if attributes[attribute_name] is not None:
                        previous_value = getattr(library, attribute_name)
                        if attributes[attribute_name] != previous_value:
                            print(f"    Updating {attribute_name} to {attributes[attribute_name]}")
                            setattr(library, attribute_name, attributes[attribute_name])
                            made_change = True
                            result += f"{attribute_name}: '{previous_value}' --> '{attributes[attribute_name]}', "
                            changelog_comments += f"{attribute_name}: '{previous_value}' --> '{attributes[attribute_name]}', "

            if made_change:
                if update_type == 'human':
                    library.record_human_updated_datetime = datetime.now()
                elif update_type == 'automation':
                    library.record_automation_updated_datetime = datetime.now()
                else:
                    print(f"ERROR: Illegal update_type 'update_type', assuming 'automation'")
                    library.record_automation_updated_datetime = datetime.now()
                library.changelog_comments = changelog_comments
                session.flush()
                session.commit()
            else:
                result = 'Nothing to change'
            return result

        else:
            print(f"ERROR: Library entry with id {id} not found")



    def create_index(self):
        """
        create_index - Create a master index from all the constituent library indexes to be able to find spectra in any library

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here

        return()



    def find_spectra(self):
        """
        find_spectra - Return a list of spectra given query constraints

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here

        return()



#### Example using this class
def example():

    #### Create a new or attach to existing library collection
    spec_lib_collection = SpectrumLibraryCollection("zzTest.sqlite")

    #### Add a library
    spec_lib_collection.add_library()

    #### Show all libraries
    libraries = spec_lib_collection.get_libraries()
    if ( libraries is None ):
        print("The library collection is empty")
    else:
        for library in libraries:
            print(library.library_record_id,library.id_name,library.original_name)

    return()



#### Example using this class
def example2():

    #### Create a new or attach to existing library collection
    spec_lib_collection = SpectrumLibraryCollection("../spectralLibraries/SpectrumLibraryCollection.sqlite")

    try:
        library = spec_lib_collection.get_library(identifier="PXL000003", version="05-24-2011")
    except Exception as error:
        print("ERROR:",error)
        return()
    print("\t".join([str(library.library_record_id),library.id_name,library.version,library.original_name]))
    return()



#### If this class is run from the command line, perform a short little test to see if it is working correctly
def main():

    #### Run an example
    example2()
    return()

if __name__ == "__main__": main()
