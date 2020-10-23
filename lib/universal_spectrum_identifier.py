#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

import re
import json

from proforma_peptidoform import ProformaPeptidoform
from response import Response

# Define a subset of useful atomic masses and the proton
atomic_masses = {
    'proton': 1.00727646688,
    'H': 1.007825035,
    'C': 12.0000000,
    'N': 14.0030740,
    'O': 15.99491463,
    'P': 30.973762,
    'S': 31.9720707,
}


class UniversalSpectrumIdentifier(object):

    # usi object takes usi_string an automatically parses it and stores attributes
    # usi objects can still exist even if the usi str is incorrect.
    # it will simply show where the error in the string is

    def __init__(self, usi=None):

        self.usi = usi

        self.is_valid = False
        self.identifier_type = None

        self.collection_identifier = None
        self.collection_type = None
        self.dataset_subfolder = None
        self.ms_run_name = None
        self.index_type = None
        self.index = None
        self.interpretation = None
        self.peptidoform_string = None
        self.peptidoform = None
        self.charge = None
        self.provenance_identifier = None

        self.error = 0
        self.error_code = None
        self.error_message = None
        self.warning_message = None

        if usi:
            self.parse(usi,verbose=None)
        

    # Attributes:
    #   usi
    #   collection_identifier
    #   dataset_subfolder
    #   ms_run_name
    #   index_type
    #   index
    #   interpretation
    #   peptidoform
    #   charge


    #### Set the error state with supplied information
    def set_error(self, error_code, error_message):
        self.error_code = error_code
        self.error_message = error_message
        self.is_valid = False


    #### Parse the USI string
    def parse(self, usi, verbose=False):

        # Reset all destintion values for a fresh parse
        self.is_valid = False
        self.identifier_type = None

        self.collection_identifier = None
        self.collection_type = None
        self.dataset_subfolder = None
        self.ms_run_name = None
        self.index_type = None
        self.index = None
        self.interpretation = None
        self.peptidoform = None
        self.peptidoform_string = None
        self.charge = None
        self.provenance_identifier = None

        self.error = 0
        self.error_code = None
        self.error_message = None
        self.warning_message = None

        # Get or set the usi string and ensure it is not not None
        if usi is None:
            usi = self.usi
        else:
            self.usi = usi
        if usi is None:
            self.set_error("NullUSI","USI is NULL")
            return self

        # Ensure that the usi is a string
        usi = str(usi)

        # Handle verbose mode
        verboseprint = print if verbose else lambda *a, **k: None
        verboseprint(f"INFO: Parsing USI string '{usi}'")

        # Ensure that the string starts with 'mzspec:' else we can stop right here
        if usi.startswith("mzspec:"):
            usi_body = usi[len("mzspec:"):]
        else:
            self.set_error("MissingPrefix","USI string does not begin with prefix 'mszpec:'")
            return self

        # creates list of potential usi components
        elements = usi_body.split(":")
        n_elements = len(elements)
        offset = 0

        # checks if usi has at least 2 colon-separated fields after mzspec:
        if n_elements < 2:
            self.set_error("InsufficientComponents","USI string does not have the minimum required 2 colon-separated components after mzspec:")
            return self

        # Extract the collection identifier field
        self.collection_identifier = elements[offset]
        if self.collection_identifier == '':
            self.set_error("EmptyDatasetIdentifier","USI component collection identifier is empty. Not permitted.")
            return self
        verboseprint(f"INFO: collection identifier is {self.collection_identifier}")

        # Extract the MS run name field
        offset += 1
        self.ms_run_name = elements[offset]
        verboseprint(f"INFO: MS run name so far is {self.ms_run_name}")

        # Scan the rest of the elements for one of the permitted index_flags
        # Be forgiving and allow improper case
        permitted_index_flags = { 'SCAN': 'scan', 'INDEX': 'index', 'NATIVEID': 'nativeId', 'TRACE': 'trace' }
        offset += 1
        while offset < n_elements:
            if elements[offset].upper() in permitted_index_flags:
                self.index_type = permitted_index_flags[elements[offset].upper()]
                verboseprint(f"INFO: Found index type '{self.index_type}'")
                break
            else:
                self.ms_run_name += ':' + elements[offset]
                verboseprint(f"INFO: MS run name is now {self.ms_run_name}")
            offset += 1

        if self.index_type is None:
            verboseprint(f"INFO: Did not detect an index flag of 'scan', 'index', nativeId', or 'trace'. Assuming this is just a Universal MS Run Identifier. Not a true USI, but maybe okay.")
            return self

        # Parse the index number
        offset += 1
        if offset < n_elements:
            self.index = elements[offset]
            if self.index:
                verboseprint("INFO Index is " + self.index)
            else:
                self.set_error("MissingIndex","Index number empty! Not permitted.")
                self.error += 1
        else:
            self.set_error("MissingIndex",f"There is no component after '{self.index_type}'")
            return self

        # If we got to here, it is at least a USI
        self.identifier_type = 'USI'

        # Extract the interpretation string
        offset += 1
        while offset < n_elements:
            if self.interpretation is None:
                self.interpretation = elements[offset]
                self.identifier_type = 'UPSMI'
            else:
                self.interpretation += ':' + elements[offset]
            offset += 1

        # Try to decompose the interpretation string
        if self.interpretation is not None:

            # First check to see if this is PSM provenance identifier
            match = re.match(r'(.+)/(\d+):(.+)$', self.interpretation)
            if match:
                self.interpretation = match.group(1) + '/' + match.group(2)
                self.provenance_identifier = match.group(3)
                self.peptidoform_string = match.group(1)
                self.charge = int(match.group(2))
                self.identifier_type = 'UPSMPI'

            # Otherwise try to decompose the interpretation
            else:
                match = re.match(r'(.+)/(\d+)$', self.interpretation)
                if match:
                    self.peptidoform_string = match.group(1)
                    self.charge = int(match.group(2))
                else:
                    self.set_error("MissingCharge",f"There is no charge number (e.g. '/2') provided in the interpretation '{self.interpretation}'")
                    return self

        # If the MS run name begins with a [ then try to extract a subfolder
        if self.ms_run_name.startswith('['):
            split_ms_run_name = self.ms_run_name.split(']')
            if len(split_ms_run_name) == 1:
                pass
            elif len(split_ms_run_name) == 2 :
                self.dataset_subfolder = split_ms_run_name[0][1:]
                self.ms_run_name = split_ms_run_name[1]
            else:
                self.dataset_subfolder = split_ms_run_name[0][1:]
                subfolder_complete = False
                open_count = 0
                close_count = 0
                self.ms_run_name = ''
                i = 1
                while i < len(split_ms_run_name):
                    verboseprint(f"Looping: {i} {subfolder_complete}: '{split_ms_run_name[i]}'")
                    if not subfolder_complete:
                        open_count = self.dataset_subfolder.count('[')
                        close_count = self.dataset_subfolder.count(']')
                        if open_count == close_count:
                            self.ms_run_name += split_ms_run_name[i]
                            subfolder_complete = True
                        else:
                            self.dataset_subfolder += split_ms_run_name[i] + ']'
                    else:
                        self.ms_run_name += split_ms_run_name[i] + ']'
                    i += 1

        # Validate the collection identifier against the currently allowed set
        if self.collection_identifier is not None:
            possible_templates = { 'PXD': r'PXD\d{6}$', 'PXL': r'PXL\d{6}$', 'MSV': r'MSV\d{9}$', 'placeholder': r'USI000000' }
            for template_type,template in possible_templates.items():
                match = re.match(template,self.collection_identifier)
                if match:
                    self.collection_type = template_type
                    break
            if self.collection_type is None:
                self.set_error("UnsupportedCollection",f"The collection identifier does not match a supported template")
                return self

        #### Try to parse the peptidoform and extract useful information from it
        if self.peptidoform_string is not None:
            peptidoform = ProformaPeptidoform(self.peptidoform_string, verbose=verbose)
            self.peptidoform = peptidoform.to_dict()

            if peptidoform.response.n_errors == 0:
                if self.peptidoform['neutral_mass'] and self.charge > 0:
                    self.peptidoform['ion_mz'] = ( self.peptidoform['neutral_mass'] + atomic_masses['proton'] * self.charge ) / self.charge

            self.error += peptidoform.response.n_errors

        # If there are no recorded errors, then we're in good shape
        if self.error == 0:
           self.is_valid = True

        # But if there were errors found, then USI is not valid
        else:
            verboseprint("Number of errors: " + str(self.error))
            self.is_valid = False
            verboseprint("ERROR: Invalid USI " + self.usi)

        return self


    # prints out USI attributes
    def show(self):
        if self.error_code is not None:
            print(f"error_code: {self.error_code}")
            print(f"error_message: {self.error_message}")
        print("USI: " + str(self.usi))
        print("is_valid: " + str(self.is_valid))
        print("identifier type: " + str(self.identifier_type))
        print("Collection Identifier: " + str(self.collection_identifier))
        print("Collection Type: " + str(self.collection_type))
        print("Dataset Subfolder: " + str(self.dataset_subfolder))
        print("MS run name: " + str(self.ms_run_name))
        print("Index flag: " + str(self.index_type))
        print("Index number: " + str(self.index))
        print("Interpretation: " + str(self.interpretation))
        print("Peptidoform string: " + str(self.peptidoform_string))
        print("Charge: " + str(self.charge))
        print("Provenance identifier: " + str(self.provenance_identifier))
        #print("Peptidoform: " + str(self.peptidoform))



############################################################################################
#### Define example peptidoforms to parse
def define_examples():
    return [
        [   "valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951" ],
        [ "invalid", "PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951" ],
        [   "valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951:PEPT[Phospho]IDELVISK/2" ],
        [ "invalid", "mzspec:PASS002437:00261_A06_P001564_B00E_A00_R1:scan:10951:PEPT[Phospho]IDELVISK/2" ],
        [ "invalid", None ],
        [ "invalid", 3 ],
        [ "invalid", "mzspec" ],
        [ "invalid", "mzspec:" ],
        [ "invalid", "mzspec:PXD001234" ],
        [ "invalid", "mzspec:PXD001234:00261_A06_P001564_B00E_A00_R1:scan" ],
        [   "valid", "mzspec:PXD001234:00261_A06_P001564_B00E_A00_R1:index:10951" ],
        [   "valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951:PEPT[Phospho]IDELVISK/2" ],
        [   "valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951:PEPT[+79]IDELVISK/2" ],
        [   "valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951:PEPT[UNIMOD:34]IDELVISK/2" ],
        [   "valid", "mzspec:PXD001234:Dilution1:4:scan:10951"],
        [   "valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:test1:scan:10951:PEPT[Phospho]IDELVISK/2" ],
        [   "valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1\\:test1:scan:10951:PEPT[Phospho]IDELVISK/2" ],
        [   "valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951:PEPT[Phospho]IDELVISK/2:PA-28732" ],
        [   "valid", "mzspec:PXD001234:[Control]fr10:scan:10951" ],
        [   "valid", "mzspec:PXD001234:[Control[2]]fr10:scan:10951" ],
        [   "valid", "mzspec:PXD001234:[Control]fr10[7]:scan:10951" ],
        [   "valid", "mzspec:PXD001234:[Control[2]]fr10[7]:scan:10951" ],
        [   "valid", "mzspec:MSV000086127:[Control[2]]fr10[7]:scan:10951" ],
        [   "valid", "mzspec:PXD002437:fr5:scan:10951:[UNIMOD:214]PEPT[Phospho]IDEL[+12.0123]VIS[UNIMOD:1]K[iTRAQ4plex]/2" ],
        [   "valid", "mzspec:USI000000:a:scan:1:EM[Oxidation]EVEES[UNIMOD:21]PEK/2"], # 24
        [   "valid", "mzspec:USI000000:a:scan:1:ELV[+11.9784|info:suspected frobinylation]IS/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:ELV[INFO:suspected frobinylation|+11.9784]IS/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:ELV[+11.9784|info:suspected frobinylation | INFO:confidence:high]IS/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:ELV[Oxidation | INFO:confidence:high]IS/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:ELV[info:AnyString]IS/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:EM[L-methionine sulfoxide]EVEES[MOD:00046]PEK/2"], # 30
        [   "valid", "mzspec:USI000000:a:scan:1:[iTRAQ4plex]-EMEVNESPEK[UNIMOD:214]-[Methyl]/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:[iTRAQ4plex]-EM[Oxidation]EVNESPEK[UNIMOD:214]-[Methyl]/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:SEQUEN[Glycan:HexNAc1Hex 2]CE/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:SEQUEN[Formula:C12H20O2]CE/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:YPVLN[GNO:G62765YT]VTMPN[GNO:G02815KT]NSNGKFDK/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:SN{Hex|INFO:completely labile}ACK/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:{Hex|INFO:completely labile}DINNER/2"],
        [   "valid", "mzspec:PXD000000:a:scan:1:{Hex|INFO:completely labile}[iTRAQ4plex]-EM[Oxidation]EVNESPEK[UNIMOD:214]-[Methyl]/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:EM[Oxidation]EVEES[U:Phospho]PEK/2"],
        [   "valid", "mzspec:USI000000:a:scan:1:EM[Carboxyamidomethylation]EVEES[U:homoarginine]PEK/2"], # 40
        [ "invalid", "mzspec:USI000000:a:scan:1:EM[U:L-methionine sulfoxide]E[U:Phospho]V[P:Phospho]EES[P:L-methionine sulfoxide]P[UNIMOD:99999]E[MOD:99999]K/2"],
    ]
 

 # Run all examples through the parser and see if the results are as expected
def run_tests():
    examples = define_examples()
    validity_map = { 'valid': True, 'invalid': False }

    # Loop over each test USI, parse it, and determine if it is valid or not, and print the index number
    print("Testing example USIs:")
    i_counter = 0
    for example in examples:
        expected_status = example[0]
        usi_string = example[1]

        usi = UniversalSpectrumIdentifier()
        usi.parse(usi_string)

        status = 'OK'
        if usi.is_valid != validity_map[expected_status]:
            status = 'FAIL'

        print(f"{i_counter}\t{status}\texpected {expected_status}\t{usi_string}")
        i_counter += 1


#### A very simple example of using this class
def run_one_example(example_number):
    examples = define_examples()
 
    usi_string = examples[example_number][1]

    usi = UniversalSpectrumIdentifier()
    usi.parse(usi_string, verbose=1)

    print('==== USI:')
    print(json.dumps(usi.__dict__, sort_keys=True, indent=2))


#### If class is run directly
def main():

    #### Parse command line options
    import argparse
    argparser = argparse.ArgumentParser(description='Command line interface to the UniversalSpectrumIdentifier class')
    argparser.add_argument('--verbose', action='count', help='If set, print out messages to STDERR as they are generated' )
    argparser.add_argument('--example', type=int, help='Specify an example to run instead of unit tests (use --example=1)')
    argparser.add_argument('--test', action='count', help='If set, run all tests')
    params = argparser.parse_args()

    #### Set verbosity of the Response class
    if params.verbose is not None:
        Response.output = 'STDERR'

    #### If --test is specified, run the full test suite adn return
    if params.test is not None:
        run_tests()
        return

    #### If --example is specified, run that example number
    example_number = 23
    if params.example is not None:
        example_number = params.example
    run_one_example(example_number)


if __name__ == "__main__": main()
