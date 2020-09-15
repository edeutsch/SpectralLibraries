#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

import re

class UniversalSpectrumIdentifier(object):

    # usi object takes usi_string an automatically parses it and stores attributes
    # usi objects can still exist even if the usi str is incorrect.
    # it will simply show where the error in the string is

    def __init__(self, usi=None):

        self.usi = usi

        self.is_valid = False
        self.identifier_type = None

        self.collection_identifier = None
        self.dataset_subfolder = None
        self.ms_run_name = None
        self.index_type = None
        self.index = None
        self.interpretation = None
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
        self.collection_identifier = None
        self.dataset_subfolder = None
        self.ms_run_name = None
        self.index_type = None
        self.index = None
        self.interpretation = None
        self.peptidoform = None
        self.charge = None
        self.provenance_identifier = None

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
                self.peptidoform = match.group(1)
                self.charge = int(match.group(2))
                self.identifier_type = 'UPSMPI'

            # Otherwise try to decompose the interpretation
            else:
                match = re.match(r'(.+)/(\d+)$', self.interpretation)
                if match:
                    self.peptidoform = match.group(1)
                    self.charge = int(match.group(2))

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
        print("Dataset Identifier: " + str(self.collection_identifier))
        print("Dataset Subfolder: " + str(self.dataset_subfolder))
        print("MS run name: " + str(self.ms_run_name))
        print("Index flag: " + str(self.index_type))
        print("Index number: " + str(self.index))
        print("Interpretation: " + str(self.interpretation))
        print("Peptidoform: " + str(self.peptidoform))
        print("Charge: " + str(self.charge))
        print("Provenance identifier: " + str(self.provenance_identifier))


# If this class is run from the command line, perform a short little test to see if it is working correctly
def run_tests():
    test_usis = [
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
    ]
    test_usisValid = []

    # Loop over each test USI, parse it, and determine if it is valid or not, and print the index number
    print("Testing example USIs:")
    for usiSet in test_usis:
        expectedStatus = usiSet[0]
        usi_string = usiSet[1]

        # Create a new UniversalSpectrumIdentifier object
        # made the USI object itself take a string so that parse does not need to be called explicitly
        usi = UniversalSpectrumIdentifier(usi_string)
        expected_validity = True
        if expectedStatus == 'invalid':
            expected_validity = False
        else:
            expectedStatus = 'valid  '

        status = 'PASS'
        if usi.is_valid is not expected_validity:
            status = 'FAIL'

        response = usi.is_valid
        test_usisValid.append(response)
        print(f"{status}\texpected {expectedStatus}\t{usi_string}")
    # check to see if parsing is correct
    #print(test_usisValid)


#### A very simple example of using this class
def run_one_test():
    test_usis = [
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
    ]
 
    usi_string = test_usis[21][1]

    usi = UniversalSpectrumIdentifier()
    usi.parse(usi_string, verbose=1)
    print('==== Result:')
    usi.show()


#### A very simple example of using this class
def example():
    usi_string = "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951"
    usi_string = "mzspec:PXL000001:01-09-2015:index:500"
    usi_string = "mzspec:PXL000001::index:500"
    usi = UniversalSpectrumIdentifier(usi_string)
    #usi.parse(verbose=1)
    usi.show()


#### If class is run directly
def main():
    #example()
    run_tests()
    #run_one_test()


if __name__ == "__main__": main()
