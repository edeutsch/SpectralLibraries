#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, flush=True, **kwargs)

import os
import os.path
import argparse
import re
import timeit

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../lib")
from SpectrumLibraryCollection import SpectrumLibraryCollection
from SpectrumLibrary import SpectrumLibrary

import io
import tarfile
import urllib.request
import shutil
from zipfile import ZipFile


#### Open a URL for a gzipped doument as a file
def open_url(url: str, buffer_size: int=4 ** 20):
    buffer = io.BytesIO(
        urllib.request.urlopen(url).read(buffer_size)
    )

    if url.endswith('.zip'):
        zipfile = ZipFile(buffer)
        filename = list(zipfile.namelist())[0]
        fh = zipfile.open(filename)
        return fh

    elif url.endswith('.tar.gz'):
        arc = tarfile.open(
            fileobj=buffer,
            mode='r:gz')

        ti = arc.firstmember
        fh = arc.extractfile(ti)
        return fh

    else:
        print(f"ERROR: Don't know how to fetch and decompress 'url'")


#### Main
def main():

    argparser = argparse.ArgumentParser(description='Front-end for maintenance of a spectral library collection')

    argparser.add_argument('--list', action='store_true', help="List the available libraries in the collection")
    argparser.add_argument('--id', action='store', help="Integer identifier (e.g. 1) to access")
    argparser.add_argument('--show', action='store_true', help="Show the library idenrtified with --id")
    argparser.add_argument('--update', action='store_true', help="Rescan the spectralLibraries directory and update the collection")
    argparser.add_argument('--DROP', action='store_true', help="DROP the current collection database and re-create it. Careful!")

    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    params = argparser.parse_args()

    collection_dir = os.path.dirname(os.path.abspath(__file__))+"/../spectralLibraries"
    database_file = 'SpectrumLibraryCollection.sqlite'
    database_filepath = f"{collection_dir}/{database_file}"

    #### If --list, then list existing libraries
    if params.DROP is True:
        if os.path.exists(database_filepath):
            print("INFO: DELETING the collections database and re-creating")
            os.remove(database_filepath)
        else:
            eprint(f"ERROR: Database {database_filepath} does not exist")
        return


    #### Get a listing of what entries are already in the collection
    spectrum_library_collection = SpectrumLibraryCollection(database_filepath)
    libraries = spectrum_library_collection.get_libraries()


    #### If --list, then list existing libraries
    if params.list is True:
        if ( len(libraries) == 0 ):
            print("The library collection is empty")
        else:
            for library in libraries:
                print(f"{library.library_record_id:4d}  {library.id_name}  {library.version_tag:15s}  {library.title:40s}  {library.original_filename}")
        return


    #### If --show --id PXDnnnnnn, then show that record
    if params.show:
        if not params.id:
            eprint(f"ERROR: Please specify an --id PXLnnnnnn to show")
            return
        library = spectrum_library_collection.get_library(identifier=params.id)
        if library is None:
            return
        attribute_names = spectrum_library_collection.get_all_library_attribute_names()
        for attribute_name in attribute_names:
            print(f"  {attribute_name:27s} = {getattr(library,attribute_name)}")
        return


    #### If --update, then sync the files in spectralLibraries with the collection
    if params.update is True:

        #### Also create a dict for the filenames
        libraries_dict = {}
        counter = 0
        for library in libraries:
            libraries_dict[library.original_filename] = counter
            counter += 1

        #### Read the metadata file for the collection
        metadatafile = collection_dir + "/SpectrumLibraryCollection.tsv"
        metadata_dict = {}
        metadata_files = []
        metadata_entries = []
        line_counter = 0
        expected_column_names = spectrum_library_collection.get_all_library_attribute_names()
        with open(metadatafile, 'r') as infile:
            print(f"INFO: Reading {metadatafile}")
            for line in infile:
                line_counter += 1
                line = line.rstrip()
                columns = line.split("\t")
                if line_counter == 1:
                    icolumn = 0
                    for expected_column_name in expected_column_names:
                        if columns[icolumn] != expected_column_name:
                            eprint(f"ERROR: Expected fourth column to be '{expected_column_name}' but it was '{columns[icolumn]}'")
                            exit()
                        icolumn += 1
                    print("INFO: Column headers in input file are correct")
                    continue

                icolumn = 0
                attributes = {}
                for expected_column_name in expected_column_names:
                    try:
                        value = columns[icolumn]
                        if value == '':
                            value = None
                        attributes[expected_column_name] = value
                    except:
                        attributes[expected_column_name] = None
                    icolumn += 1

                if attributes['keywords'] is not None:
                    attributes['keywords']= attributes['keywords'].replace('"','')
                if attributes['version_tag'] is not None:
                    attributes['version_tag']= attributes['version_tag'].replace('"','')
                original_filename = columns[7]
                metadata_entries.append(attributes)
                metadata_dict[original_filename] = attributes
                metadata_files.append(original_filename)

            print(f"INFO: Read {len(metadata_files)} entries")

        #### Loop over all files in the directory and put them in a hash
        local_file_dict = {}
        for filename in os.listdir(collection_dir):
            match = re.fullmatch("(.+)\.(msp|sptxt)",filename)
            if match is not None:
                local_file_dict[filename] = 1
        print(f"INFO: Found {len(local_file_dict)} local library files")

        #### Loop over all entries in the metadata files and check against the database
        i_library = 0
        for library_entry in metadata_entries:

            print('========================================')
            i_library += 1
            library_record_id = library_entry['library_record_id']
            id_name = library_entry['id_name']
            if library_record_id is None or library_record_id == '':
                print(f"INFO: Row {i_library} has no library_record_id. Ending")
                break

            print(f"Processing library {id_name} ")
            existing_libraries = spectrum_library_collection.get_library(identifier=id_name)
            if existing_libraries is None:
                print(f"  - This library is not yet in the database. Adding it.")
                result_id_name = spectrum_library_collection.add_library(library_entry)
                print(f"  - Library {result_id_name} has been added")
            else:
                print(f"  - This library is already in the database. Updating metadata")
                result = spectrum_library_collection.update_library_metadata(id=library_record_id, attributes=library_entry, update_type='automation', changelog_comments='Updated metadata from master reference tsv')
                print(f"  - Update result was: {result}")

            filename = library_entry['original_filename']

            #### If the file doesn't exist here, then maybe we can fetch it
            if filename not in local_file_dict:
                print(f"  - Did not find file {filename} locally. Try to fetch it")
                try:
                    url = metadata_dict[filename]['source_url']
                except:
                    print(f"ERROR: No URL for {filename}")
                    return

                print(f"  - URL for {filename} is {url}")
                print(f"  - Fetching library..")
                enable_library_downloading = True
                if enable_library_downloading:
                    if True:
                        with open_url(url) as infile:
                            with open(filename, 'wb') as outfile:
                                shutil.copyfileobj(infile, outfile)
                        local_file_dict[filename] = 1
                    else:
                        print(f"ERROR: Unable to fetch {url}")
                        return
                else:
                    print(f"  - (Remote file fetch turned off for the moment)")

            if filename in local_file_dict:
                print(f"  - Processing {filename}")

                #### Check to see if there is an index file yet
                if os.path.exists(f"{collection_dir}/{filename}.splindex"):
                    print("  - Index already exists")
                else:
                    print("  - Need to create an index")
                    spectrum_library = SpectrumLibrary()
                    spectrum_library.filename = f"{collection_dir}/{filename}"
                    t0 = timeit.default_timer()
                    spectrum_library.create_index()
                    t1 = timeit.default_timer()
                    print(f"\n  - Elapsed time: {t1-t0}")

            else:
                print(f"ERROR: filename '{filename}' in the metadata file not found locally")
                return

        return

    #### If we got here, there was no recognized command
    print("ERROR: No valid set of parameters. See --help for more information")

    return

if __name__ == "__main__": main()


'''
                #### Check to see if the file is already registered
                if filename in libraries_dict:
                    print(f" - This library is already in the collection as {libraries[libraries_dict[filename]].id_name}")
                    print(f"    (stored version={libraries[libraries_dict[filename]].version_tag}")
                    version_tag = metadata_dict[filename]['version_tag']
                    print(f"    (metadata sheet version={version_tag}")
                    spectrum_library_collection.update_library_metadata(id=libraries[libraries_dict[filename]].library_record_id, attributes=metadata_dict[filename])

                else:
                    print(f"  Need to create a record for filename")
                    if filename in metadata_dict:
                        attributes = metadata_dict[filename]
                        attributes['local_filename'] = attributes['original_filename']
                        spectrum_library_collection.add_library(attributes)
                    else:
                        print(f"ERROR: There is no metadata entry for {filename} yet. Please create one and update again")
                        return

'''
