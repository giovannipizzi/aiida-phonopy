# -*- coding: utf-8 -*-

from aiida.orm.data.folder import FolderData
from aiida.parsers.parser import Parser
from aiida.common.datastructures import calc_states
from aiida.parsers.exceptions import OutputParsingError
from aiida.common.exceptions import UniquenessError
import numpy
from aiida.orm.data.array import ArrayData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
import json
from aiida_phonopy.calculations.phonopy import PhonopyCalculation

__copyright__ = u"Copyright (c), 2014-2015, École Polytechnique Fédérale de Lausanne (EPFL), Switzerland, Laboratory of Theory and Simulation of Materials (THEOS). All rights reserved."
__license__ = "Non-Commercial, End-User Software License Agreement, see LICENSE.txt file"
__version__ = "0.4.1"


class PhonopyParser(Parser):
    """
    This class is the implementation of the Parser class for a Phonopy calculator.
    """
    
    _out_band_name = 'phonon_frequencies'
    _out_dos_name = 'phonon_dos'
    _out_thermal_name = 'thermal_properties'
        
    def __init__(self,calc):
        """
        Initialize the instance of PhonopyParser
        """
        # check for valid input
        if not isinstance(calc,PhonopyCalculation):
            raise OutputParsingError("Input calculation must be a PhonopyCalculation")
        self._calc = calc
        
    def parse_with_retrieved(self, retrieved):
        """
        Parses the datafolder, stores results.
        This parser for this simple code does simply store in the DB a node
        representing the file of forces in real space
        """
        from aiida.common.exceptions import InvalidOperation
        
        # suppose at the start that the job is successful
        successful = True
        
        # check that calculation is in the right state
#        state = self._calc.get_state()
#        if state != calc_states.PARSING:
#            raise InvalidOperation("Calculation not in {} state"
#                                   .format(calc_states.PARSING) )
        
        # select the folder object
        # Check that the retrieved folder is there 
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            self.logger.error("No retrieved folder found")
            return False, ()
        
        # check what is inside the folder
        list_of_files = out_folder.get_folder_list()
        
        # at least the stdout should exist
        if not self._calc._OUTPUT_FILE_NAME in list_of_files or not self._calc._RESULT_FILE_NAME:
            successful = False
            self.logger.error("Output/results not found",extra=logger_extra)
            return successful,()
        
        # load the results dictionary
        json_outfile = out_folder.get_abs_path( self._calc._RESULT_FILE_NAME )
        with open(json_outfile,'r') as f:
            json_params = json.load(f) 
        
        # look at warnings
        warnings = []
        with open(out_folder.get_abs_path( self._calc._SCHED_ERROR_FILE )) as f:
            errors = f.read()
        if errors:
            warnings = [errors]
        
        # I implicitly assume that all data inside the json are arrays
        # it should be very often the case for the phonon properties
        
        # ====================== prepare the output nodes ======================
        
        # save the outputs
        new_nodes_list= []
        
        # save dos
        try:
            frequencies_dos = json_params['frequencies_dos']
            total_dos = json_params['total_dos']
            
            array_dos = ArrayData()
            array_dos.set_array('frequency', frequencies_dos)
            array_dos.set_array('phonon_dos', total_dos)
                
            new_nodes_list.append( (self._out_dos_name, array_dos) )
        except KeyError: # keys not found in json
            pass
        
        # save thermodynamic quantities
        try:
            temperature = json_params['temperature']
            free_energy = json_params['free_energy']
            entropy = json_params['entropy']
            cv = json_params['cv']
            
            array_thermal = ArrayData()
            array_thermal.set_array('temperature', temperature)
            array_thermal.set_array('free_energy', free_energy)
            array_thermal.set_array('entropy', entropy)
            array_thermal.set_array('specific_heat', cv)
            # TODO: in which units am I storing stuff???
            
            new_nodes_list.append( (self._out_thermal_name, array_thermal) )
        except KeyError: # keys not found in json
            pass
        
        # save frequencies
        
        array_freq = BandsData()
        
        try:
            structure = self._calc.inp.structure
        except AttributeError:
            structure = self._calc.inp.force_constants.structure
            
        inp_kpoints = self._calc.inp.qpoints
        
        try:
            inp_kpoints.get_kpoints()
            array_freq.set_kpointsdata(inp_kpoints)
        except AttributeError: # it had a mesh of kpoints in input
            try:
                cell = inp_kpoints.cell
            except AttributeError:
                cell = structure.cell
            try:
                pbc = inp_kpoints.pbc
            except AttributeError:
                pbc = structure.pbc
            
            try:
                the_kpoints = json_params['q_points']
            except KeyError:
                the_kpoints = inp_kpoints.get_kpoints()    
            
            try:
                the_weights = json_params['weights']
            except KeyError:
                the_weights = None
                
            array_freq.cell = cell
            array_freq.pbc = pbc
            array_freq.set_kpoints(the_kpoints, weights=the_weights)
            array_freq.labels = inp_kpoints.labels
        
        try:
            frequencies = json_params['frequencies']
        except KeyError:
            warnings.append('Unable to read phonon frequencies')
            new_nodes_list.append((self.get_linkname_outparams(), ParameterData(dict={'warnings': warnings})))
            return False, new_nodes_list
        
        labels = 'frequencies'
        bands = frequencies
        
        try:
            group_velocities = json_params['group_velocities']
            vx = [ _[0] for _ in group_velocities ]
            vy = [ _[1] for _ in group_velocities ]
            vz = [ _[2] for _ in group_velocities ]
            bands = [frequencies]
            labels = ['frequencies']
            bands.append(vx)
            bands.append(vy)
            bands.append(vz)
            labels += ['vx','vy','vz']
        except KeyError:
            pass
        
        array_freq.set_bands(frequencies, units='THz', occupations=None, labels='frequencies')
        #TODO: verify the units
        
        try:
            eigenvectors = json_params['eigenvectors']
            array_freq.set_array('eigenvectors', eigenvectors)
        except KeyError:
            pass
        
        print 'here'
            
        new_nodes_list.append( (self._out_band_name, array_freq) )
        #except KeyError as e: # keys not found in json
        #    raise e
        
        # add the dictionary with warnings       
        new_nodes_list.append( (self.get_linkname_outparams(), ParameterData(dict={'warnings': warnings})))
        
        return successful, new_nodes_list
        
