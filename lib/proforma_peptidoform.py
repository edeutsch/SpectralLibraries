#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, flush=True, **kwargs)

import re
import json
import os

from ontology import Ontology
from response import Response

#### Class-level structure to hold all the ontology data. Only load it once
ontologies = {}

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

# Define the basic set of amino acids
# From https://proteomicsresource.washington.edu/protocols06/masses.php
amino_acid_masses = {
    'G': 57.021463735,
    'A': 71.037113805,
    'S': 87.032028435,
    'P': 97.052763875,
    'V': 99.068413945,
    'T': 101.047678505,
    'C': 103.009184505,
    'L': 113.084064015,
    'I': 113.084064015,
    'N': 114.042927470,
    'D': 115.026943065,
    'Q': 128.058577540,
    'K': 128.094963050,
    'E': 129.042593135,
    'M': 131.040484645,
    'O': 132.089877680,  # probably not useful to have this here. But consider UNIMOD:372 Arg->Orn
    'H': 137.058911875,
    'F': 147.068413945,
    'U': 150.953633405,  # selenocysteine
    'R': 156.101111050,
    'Y': 163.063328575,
    'W': 186.079312980,
    'X': 0.0,            # Used to denote a gap of one or more residues like RTAAX[+367.0537]WT
    #'B': 114.534935,
    #'Z': 128.550585,
}


class ProformaPeptidoform(object):

    ########################################################################################
    def __init__(self, peptidoform_string=None, verbose=False):

        self.peptidoform_string = None
        self.peptide_sequence = None
        self.residue_modifications = {}
        self.terminal_modifications = {}
        self.unlocalized_mass_modifications = {}
        self.neutral_mass = None
        self.residues = []

        self.response = Response()
        self.is_valid = False

        if len(ontologies) == 0:
            if verbose:
                eprint("Loading ontologies...")
            possible_locations = [ '.', 'C:/local/Repositories/SVN/proteomics/var/CV', '/net/dblocal/wwwspecial/proteomecentral/extern/CVs' ]
            ontology_filenames = {
                'UNIMOD': 'unimod.obo',
                'PSI-MOD': 'PSI-MOD.obo'
            }
            for ontology_key, ontology_filename in ontology_filenames.items():
                for possible_location in possible_locations:
                    ontology_path = possible_location + '/' + ontology_filename
                    if os.path.exists( ontology_path ):
                        if verbose:
                            eprint(f" - Loading {ontology_key} from {ontology_path}")
                        ontologies[ontology_key] = Ontology(filename=ontology_path)
                        break
                if ontology_key not in ontologies:
                    eprint(f"ERROR: Unable to locate {ontology_key} with filename {ontology_filename}")
            if verbose:
                eprint(" - Done")

        if peptidoform_string:
            self.parse(peptidoform_string, verbose=None)


    ########################################################################################
    def to_dict(self):
        result = {}
        for key in [ 'peptidoform_string', 'peptide_sequence', 'residue_modifications', 'unlocalized_mass_modifications',
            'terminal_modifications', 'neutral_mass', 'is_valid', 'residues' ]:
            result[key] = self.__getattribute__(key)
        return result


    ########################################################################################
    #### Set the error state with supplied information
    def set_error(self, error_code, error_message):
        self.error_code = error_code
        self.error_message = error_message
        self.is_valid = False


    ########################################################################################
    def parse(self, peptidoform_string=None, verbose=False):

        # Verify that we were either passed a peptidoform_string or there was already one available
        if peptidoform_string is None:
            if self.peptidoform_string is None:
                self.response.error(f"No peptidoform string is available to parse", error_code="NoPeptidoform", http_status=400)
                return self.response
            peptidoform_string = self.peptidoform_string
        else:
            self.peptidoform_string = peptidoform_string

        # Set up a stack counters for special paired characters
        character_stack = { 'square_brackets': 0, 'curly_brackets': 0 }

        # Create a blank data structure to contain the parsed information
        residue_modifications = {}
        terminal_modifications = {}
        residues = []
        unlocalized_mass_modifications = []
        self.residue_modifications = residue_modifications
        self.terminal_modifications = terminal_modifications
        self.residues = residues
        self.unlocalized_mass_modifications = unlocalized_mass_modifications

        # Split the input string into a list of characters to process
        characters = list(peptidoform_string)

        # Prepare for looping
        i_residue = 0
        current_residue = ''
        i_char = 0

        # Loop over all characters, processing what we find
        for char in characters:
            if char == '[':
                character_stack['square_brackets'] += 1
                current_residue += char
            elif char == ']':
                if character_stack['square_brackets'] == 0:
                    self.response.error(f"Unmatched square bracket at position {i_char}", error_code="ErrorInPeptidoform", http_status=400)
                elif character_stack['square_brackets'] == 1:
                    current_residue += char
                    character_stack['square_brackets'] -= 1
                else:
                    current_residue += char
                    character_stack['square_brackets'] -= 1
            elif char == '{':
                character_stack['curly_brackets'] += 1
                current_residue += char
            elif char == '}':
                if character_stack['curly_brackets'] == 0:
                    self.response.error(f"Unmatched curly bracket at position {i_char}", error_code="ErrorInPeptidoform", http_status=400)
                elif character_stack['curly_brackets'] == 1:
                    current_residue += char
                    character_stack['curly_brackets'] -= 1
                    #### If we're still at the starting position, then move this to a labile mod pile
                    if i_residue == 0:
                        self.unlocalized_mass_modifications.append(( { 'residue_string': current_residue, 'index': -1 } ))
                        current_residue = ''
                else:
                    current_residue += char
                    character_stack['curly_brackets'] -= 1
            else:
                if character_stack['square_brackets'] > 0 or character_stack['curly_brackets'] > 0:
                    current_residue += char
                else:
                    if char == '-' and i_residue == 0:
                        current_residue += char
                    else:
                        residues.append( { 'residue_string': current_residue, 'index': i_residue } )
                        i_residue += 1
                        current_residue = char

            #print(f"{i_char}=={char}. {i_residue}=={current_residue}")
            i_char += 1

        if current_residue != '':
            residues.append( { 'residue_string': current_residue, 'index': i_residue } )

        if residues[0]['residue_string'] == '':
            residues = residues[1:]
        #print(residues)

        #### Set the peptidoform mass to H2O and build the mass
        self.neutral_mass = atomic_masses['H'] * 2 + atomic_masses['O']
        self.peptide_sequence = ''

        for residue in residues:
            #print(residue)
            amino_acid = ''
            if len(residue['residue_string']) > 1:
                if residue['index'] == 0:
                    residue['base_residue'] = 'nterm'
                    offset = 1
                else:
                    amino_acid = residue['residue_string'][0]
                    self.peptide_sequence += amino_acid
                    residue['base_residue'] = amino_acid
                    offset = 2

                trimmed_string = residue['residue_string'].rstrip('-')
                residue['modification_string'] = trimmed_string[offset:-1]
                residue['modification_type'] = ''

                if residue['base_residue'] == '-':
                    residue['base_residue'] = 'cterm'

                #### Further parse and interpret the modification string
                self.parse_modification_string(residue)
                if 'errors' in residue:
                    self.response.error(f"Error parsing residue {residue['index']}: {residue['residue_string']}", error_code="ErrorInPeptidoform", http_status=400)

                #### Add the delta mass to the peptidoform mass
                if 'delta_mass' in residue and residue['delta_mass'] is not None:
                    self.neutral_mass += residue['delta_mass']

                #### If this is a curly bracket modification, then try to adjust the mass classification
                if residue['residue_string'][1] == '{':
                    if 'delta_mass' in residue:
                        residue['labile_delta_mass'] = residue['delta_mass']
                        residue['delta_mass'] = 0.0

                if residue['base_residue'].endswith('term'):
                    terminal_modifications[residue['base_residue']] = residue
                else:
                    residue_modifications[residue['index']] = residue

            else:
                amino_acid = residue['residue_string'][0]
                self.peptide_sequence += amino_acid

            # Add to the peptidoform mass
            if amino_acid > '' and amino_acid != '-':
                if amino_acid in amino_acid_masses:
                    self.neutral_mass += amino_acid_masses[amino_acid]
                else:
                    self.response.error(f"Unable to determine mass of amino acid {amino_acid}", error_code="ErrorInPeptidoform", http_status=400)

        #### Handle any unlocalized mass modifications
        for unlocalized_mass_modification in self.unlocalized_mass_modifications:
            unlocalized_mass_modification['base_residue'] = ''
            unlocalized_mass_modification['modification_string'] = unlocalized_mass_modification['residue_string'][1:-1]
            self.parse_modification_string(unlocalized_mass_modification)
            if 'errors' in unlocalized_mass_modification:
                self.response.error(f"Error parsing unlocalized mass modification {unlocalized_mass_modification['residue_string']}", error_code="ErrorInPeptidoform", http_status=400)
            #### Add the delta mass to the peptidoform mass
            if 'delta_mass' in unlocalized_mass_modification and unlocalized_mass_modification['delta_mass'] is not None:
                self.neutral_mass += unlocalized_mass_modification['delta_mass']
            del(unlocalized_mass_modification['base_residue'])
            del(unlocalized_mass_modification['index'])

        #### If the peptide sequence ends with a -, then just clip it off:
        if self.peptide_sequence[-1] == '-':
            self.peptide_sequence = self.peptide_sequence.rstrip('-')

        #### Set final status
        if self.response.n_errors == 0:
            self.is_valid = True

        #### Set the modification structures to None if empty
        if len(self.unlocalized_mass_modifications) == 0:
            self.unlocalized_mass_modifications = None
        if len(self.terminal_modifications) == 0:
            self.terminal_modifications = None
        if len(self.residue_modifications) == 0:
            self.residue_modifications = None


    ########################################################################################
    #### Parse the modification string of a residue and extract the necessary information
    def parse_modification_string(self, residue, verbose=False):

        modification_string = residue['modification_string']
        if verbose:
            print(residue)

        #### try to handle the case of two consecutive mods
        multiple_mods = modification_string.split('][')
        is_first_mod = True
        if len(multiple_mods) > 1:
            for one_mod in multiple_mods:
                one_residue = residue.copy()
                one_residue['modification_string'] = one_mod
                self.parse_modification_string(one_residue, verbose=verbose)
                if is_first_mod is True:
                    residue['modification_string'] = one_mod
                    residue['modification_type'] = one_residue['modification_type']
                    residue['modification_name'] = one_residue['modification_name']
                    residue['modification_curie'] = one_residue['modification_curie']
                    residue['delta_mass'] = one_residue['delta_mass']
                else:
                    residue['modification_string'] += ' & ' + one_mod
                    residue['modification_type'] += ' & ' + one_residue['modification_type']
                    residue['modification_name'] += ' & ' + one_residue['modification_name']
                    residue['modification_curie'] += ' & ' + one_residue['modification_curie']
                    residue['delta_mass'] += one_residue['delta_mass']

                is_first_mod = False
            return


        #### Process the single mod
        modification_components = modification_string.split('|')

        for component in modification_components:

            component = component.strip()
            found_match = False

            if not found_match:
                match = re.match(r'Glycan:(.+)$',component,re.IGNORECASE)
                if match:
                    if 'warnings' not in residue:
                        residue['warnings'] = []
                    residue['warnings'].append("The 'Glycan:' prefix is recognized and legal but is not yet supported by this system")
                    residue['modification_type'] = 'unsupported'
                    found_match = True

            if not found_match:
                match = re.match(r'Formula:(.+)$',component,re.IGNORECASE)
                if match:
                    if 'warnings' not in residue:
                        residue['warnings'] = []
                    residue['warnings'].append("The 'Formula:' prefix is recognized and legal but is not yet supported by this system")
                    residue['modification_type'] = 'unsupported'
                    found_match = True

            if not found_match:
                match = re.match(r'GNO:(.+)$',component,re.IGNORECASE)
                if match:
                    if 'warnings' not in residue:
                        residue['warnings'] = []
                    residue['warnings'].append("The 'GNO:' prefix is recognized and legal but is not yet supported by this system")
                    residue['modification_type'] = 'unsupported'
                    found_match = True

            if not found_match:
                match = re.match(r'RESID:(.+)$',component,re.IGNORECASE)
                if match:
                    if 'warnings' not in residue:
                        residue['warnings'] = []
                    residue['warnings'].append("The 'RESID:' prefix is recognized and legal but is not yet supported by this system")
                    residue['modification_type'] = 'unsupported'
                    found_match = True

            if not found_match:
                match = re.match(r'info:(.+)$',component,re.IGNORECASE)
                if match:
                    if 'custom_info' not in residue:
                        residue['custom_info'] = []
                    residue['custom_info'].append(match.group(1))
                    if 'modification_type' not in residue or residue['modification_type'] is None or residue['modification_type'] == '':
                        residue['modification_type'] = 'custom_info'
                    found_match = True

            if not found_match:
                match = re.match(r'[\+\-][\d\.]+$',component)
                if match:
                    residue['delta_mass'] = float(match.group(0))
                    residue['modification_type'] = 'delta_mass'
                    found_match = True

            if not found_match:
                match = re.match(r'UNIMOD:\d+$',component)
                if match:
                    identifier = match.group(0)
                    if identifier in ontologies['UNIMOD'].terms:
                        term = ontologies['UNIMOD'].terms[identifier]
                        residue['delta_mass'] = term.monoisotopic_mass
                        residue['modification_name'] = term.name
                        residue['modification_curie'] = term.curie
                    else:
                        if 'errors' not in residue:
                            residue['errors'] = []
                        residue['errors'].append(f"The curie '{component}' cannot be found in Unimod")
                        self.response.error(f"The curie '{component}' cannot be found in Unimod", error_code="ErrorInPeptidoform", http_status=400)
                    residue['modification_type'] = 'UNIMOD_identifier'
                    found_match = True

            if not found_match:
                if component.startswith('U:'):
                    if component[2:].upper() in ontologies['UNIMOD'].uc_names:
                        matching_terms = ontologies['UNIMOD'].uc_names[component[2:].upper()]
                        for identifier in matching_terms:
                            if identifier in ontologies['UNIMOD'].terms:
                                term = ontologies['UNIMOD'].terms[identifier]
                                residue['delta_mass'] = term.monoisotopic_mass
                                residue['modification_name'] = term.name
                                residue['modification_curie'] = term.curie
                                residue['modification_type'] = 'UNIMOD_name'
                        found_match = True
                    else:
                        if 'errors' not in residue:
                            residue['errors'] = []
                        residue['errors'].append(f"The name after the U: in '{component}' cannot be found in UNIMOD")
                        self.response.error(f"The name after the U: in '{component}' cannot be found in UNIMOD", error_code="ErrorInPeptidoform", http_status=400)
                        residue['modification_type'] = 'unknown'
                        found_match = True

            if not found_match:
                if component.upper() in ontologies['UNIMOD'].uc_names:
                    matching_terms = ontologies['UNIMOD'].uc_names[component.upper()]
                    for identifier in matching_terms:
                        if identifier in ontologies['UNIMOD'].terms:
                            term = ontologies['UNIMOD'].terms[identifier]
                            residue['delta_mass'] = term.monoisotopic_mass
                            residue['modification_name'] = term.name
                            residue['modification_curie'] = term.curie
                            residue['modification_type'] = 'UNIMOD_name'
                    found_match = True

            if not found_match:
                match = re.match(r'MOD:\d+$',component)
                if match:
                    identifier = match.group(0)
                    if identifier in ontologies['PSI-MOD'].terms:
                        term = ontologies['PSI-MOD'].terms[identifier]
                        residue['delta_mass'] = term.monoisotopic_mass
                        residue['modification_name'] = term.name
                        residue['modification_curie'] = term.curie
                    else:
                        if 'errors' not in residue:
                            residue['errors'] = []
                        residue['errors'].append(f"The curie '{component}' cannot be found in PSI-MOD")
                        self.response.error(f"The curie '{component}' cannot be found in PSI-MOD", error_code="ErrorInPeptidoform", http_status=400)
                    residue['modification_type'] = 'PSI-MOD_identifier'
                    found_match = True

            if not found_match:
                if component.startswith('P:'):
                    if component[2:].upper() in ontologies['PSI-MOD'].uc_names:
                        matching_terms = ontologies['PSI-MOD'].uc_names[component[2:].upper()]
                        for identifier in matching_terms:
                            if identifier in ontologies['PSI-MOD'].terms:
                                term = ontologies['PSI-MOD'].terms[identifier]
                                residue['delta_mass'] = term.monoisotopic_mass
                                residue['modification_name'] = term.name
                                residue['modification_curie'] = term.curie
                                residue['modification_type'] = 'PSI-MOD_name'
                        found_match = True
                    else:
                        if 'errors' not in residue:
                            residue['errors'] = []
                        residue['errors'].append(f"The name after the P: in '{component}' cannot be found in PSI-MOD")
                        self.response.error(f"The name after the P: in '{component}' cannot be found in PSI-MOD", error_code="ErrorInPeptidoform", http_status=400)
                        residue['modification_type'] = 'unknown'
                        found_match = True

            if not found_match:
                if component.upper() in ontologies['PSI-MOD'].uc_names:
                    matching_terms = ontologies['PSI-MOD'].uc_names[component.upper()]
                    for identifier in matching_terms:
                        if identifier in ontologies['PSI-MOD'].terms:
                            term = ontologies['PSI-MOD'].terms[identifier]
                            residue['delta_mass'] = term.monoisotopic_mass
                            residue['modification_name'] = term.name
                            residue['modification_curie'] = term.curie
                            residue['modification_type'] = 'PSI-MOD_name'
                    found_match = True

            if not found_match:
                if 'errors' not in residue:
                    residue['errors'] = []
                residue['errors'].append(f"The modification element '{component}' cannot be understood by the parser")
                self.response.error(f"The modification element '{component}' cannot be understood by the parser", error_code="ErrorInPeptidoform", http_status=400)
                residue['modification_type'] = 'unknown'

        return self.response


############################################################################################
#### Define example peptidoforms to parse
def define_examples():
    tests = [
        [ 0,    "valid", "PEPT[Phospho]IDELVISK" ],
        [ 1,    "valid", "PEPT[phospho]IDELVISK" ],
        [ 2,    "valid", "PEPT[phosphorylation]IDELVISK" ],
        [ 3,    "invalid", "PEPT[phos]IDELVISK" ],
        [ 4,    "valid", "PEPT[+79]IDELVISK" ],
        [ 5,    "valid", "PEPT[UNIMOD:34]IDELVISK" ],
        [ 6,    "valid", "[UNIMOD:214]PEPT[Phospho]IDEL[+12.0123]VIS[UNIMOD:1]K[iTRAQ4plex]" ],
        [ 7,    "valid", "EM[Oxidation]EVEES[UNIMOD:21]PEK"],
        [ 8,    "valid", "ELV[+11.9784|info:suspected frobinylation]IS"],
        [ 9,    "valid", "ELV[INFO:suspected frobinylation|+11.9784]IS"],
        [ 10,    "valid", "ELV[+11.9784|info:suspected frobinylation | INFO:confidence:high]IS"],
        [ 11,    "valid", "ELV[Oxidation | INFO:confidence:high]IS"],
        [ 12,    "valid", "ELV[info:AnyString]IS"],
        [ 13,   "valid", "EM[L-methionine sulfoxide]EVEES[MOD:00046]PEK"],
        [ 14,   "valid", "[iTRAQ4plex]-EMEVNESPEK[UNIMOD:214]-[Methyl]"],
        [ 15,   "valid", "[iTRAQ4plex]-EM[Oxidation]EVNESPEK[UNIMOD:214]-[Methyl]"],
        [ 16,   "valid", "SEQUEN[Glycan:HexNAc1Hex2]CE"],
        [ 17,   "valid", "SEQUEN[Formula:C12H20O2]CE"],
        [ 18,   "valid", "YPVLN[GNO:G62765YT]VTMPN[GNO:G02815KT]NSNGKFDK"],
        [ 19,   "valid", "SN{Hex|INFO:completely labile}ACK"],
        [ 20,   "valid", "{Hex|INFO:completely labile}DINNER"],
        [ 21,   "valid", "{Hex|INFO:completely labile}[iTRAQ4plex]-EM[Oxidation]EVNESPEK[UNIMOD:214]-[Methyl]"],
        [ 22,   "valid", "EM[Oxidation]EVEES[U:Phospho]PEK"],
        [ 23,   "valid", "EM[Carboxyamidomethylation]EVEES[U:homoarginine]PEK"],
        [ 24, "invalid", "EM[U:L-methionine sulfoxide]E[U:Phospho]V[P:Phospho]EES[P:L-methionine sulfoxide]PEK"],
        [ 25,   "valid", "HPDIY[Phospho][Oxidation]AVPIK"],    # two mods on the same residue
        [ 26,   "valid", "[Dimethyl_2H(4)13C(2)]-DRM[Oxidation]YQM[Oxidation]DIQQELQR"],    # two mods on the same residue
    ]
 
    index = 0
    for test in tests:
        test[0] = index
        index += 1

    return tests


# Run all examples through the parser and see if the results are as expected
def run_tests():
    examples = define_examples()
    validity_map = { 'valid': True, 'invalid': False }
    reverse_validity_map = { True: 'valid', False: 'invalid' }

    # Loop over each test USI, parse it, and determine if it is valid or not, and print the index number
    print("Testing example USIs:")
    for example in examples:
        expected_status = example[1]
        peptidoform_string = example[2]

        peptidoform = ProformaPeptidoform(peptidoform_string, verbose=1)

        result = 'OK'
        validity = reverse_validity_map[peptidoform.is_valid]
        if peptidoform.is_valid != validity_map[expected_status]:
            result = 'FAIL'

        print(f"{example[0]:4d}  {result:6s}{validity:9s} expected {expected_status:9s}{peptidoform_string}")


#### A very simple example of using this class
def run_one_example(example_number):
    examples = define_examples()
 
    peptidoform_string = examples[example_number][2]

    peptidoform = ProformaPeptidoform(peptidoform_string, verbose=1)

    print('==== Peptidoform:')
    print(json.dumps(peptidoform.to_dict(), sort_keys=True, indent=2))

    print('\n==== Response:')
    print(json.dumps(peptidoform.response.__dict__, sort_keys=True, indent=2))


##############################################################################################################################
#### If class is run directly
def main():

    #### Parse command line options
    import argparse
    argparser = argparse.ArgumentParser(description='Command line interface to the ProformaPeptidoform class')
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
    example_number = 1
    if params.example is not None:
        example_number = params.example
    run_one_example(example_number)


if __name__ == "__main__": main()
