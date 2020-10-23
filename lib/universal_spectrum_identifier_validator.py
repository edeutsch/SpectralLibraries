#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

import re
import ast
import json

from universal_spectrum_identifier import UniversalSpectrumIdentifier


class UniversalSpectrumIdentifierValidator(object):

    def __init__(self, usi_list=None, options=None):

        self.usi_list = usi_list
        self.options = options

        self.error_code = None
        self.error_message = None
        self.response = { 'error_code': 'OK', 'error_message': '', 'validation_results': {} }

        if usi_list is not None:
            self.validate_usi_list()


    #### Set the error state with supplied information
    def set_error(self, error_code, error_message):
        self.error_code = error_code
        self.error_message = error_message
        self.response['error_code'] = error_code
        self.response['error_message'] = error_message


    #### Parse the USI string
    def validate_usi_list(self, usi_list=None, options=None, verbose=False):

        # Get or set the usi_list and ensure that it is valid
        if usi_list is None:
            usi_list = self.usi_list
        else:
            self.usi_list = usi_list
        if usi_list is None:
            self.set_error("NullUSIList","USI List is NULL")
            return self
        if not isinstance(usi_list, list):
            self.set_error("NotListOfUSIs","The input list of USIs is not a list")
            return self

        # Handle verbose mode
        verboseprint = print if verbose else lambda *a, **k: None
        verboseprint(f"INFO: Validating list of {len(usi_list)} USIs")

        self.response = { 'error_code': 'OK', 'error_message': '', 'validation_results': {} }
        validation_results = self.response['validation_results']

        #### Loop over the list of USIs and validate them
        for usi_str in usi_list:

            #### Check to see if we've already done this one
            if usi_str in validation_results:
                verboseprint(f"INFO: Skipping USI already validated '{usi_str}'")
                continue

            verboseprint(f"INFO: Validating USI '{usi_str}'")

            usi = UniversalSpectrumIdentifier()
            usi.parse(usi_str, verbose=verbose)

            result = usi.__dict__
            validation_results[usi_str] = result

        return self


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


#### Validate one of several testing USIs for testing purposes
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
 
    usi_string = test_usis[1][1]

    usi_validator = UniversalSpectrumIdentifierValidator()
    options = {}

    usi_validator.validate_usi_list( [usi_string], options, verbose=1)
    print('==== Result:')
    print(json.dumps(usi_validator.response,sort_keys=True,indent=2))


#### Validate the whole list of testing USIs in one batch
def run_list_test():
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
 
    usi_list = []
    for row in test_usis:
        if row[1] is not None and isinstance(row[1],str):
            usi_list.append(row[1])

    usi_validator = UniversalSpectrumIdentifierValidator()
    options = {}

    usi_validator.validate_usi_list( usi_list, options)
    print('==== Result:')
    print(json.dumps(usi_validator.response,sort_keys=True,indent=2))


#### A very simple example of using this class
def example():
    usi_string = "mzspec:PXD000000:a:scan:1:{Hex|INFO:completely labile}[iTRAQ4plex]-EM[Oxidation]EVNESPEK[UNIMOD:214]-[Methyl]/2"
    usi_validator = UniversalSpectrumIdentifierValidator([usi_string])
    print(json.dumps(usi_validator.response,sort_keys=True,indent=2))


#### If class is invoked directly
def main():
    example()
    #run_one_test()
    #run_list_test()


if __name__ == "__main__": main()
