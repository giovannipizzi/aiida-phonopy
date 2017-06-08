# -*- coding: utf-8 -*-
import json

from aiida.orm.calculation.job import JobCalculation
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.array.kpoints import KpointsData 
from aiida.orm.data.structure import StructureData
from aiida.common.exceptions import InputValidationError
from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida import get_file_header
from aiida.common.utils import classproperty
from aiida_phonopy.data.forceconstants2 import Forceconstants2Data

__copyright__ = u"Copyright (c), 2014, École Polytechnique Fédérale de Lausanne (EPFL), Switzerland, Laboratory of Theory and Simulation of Materials (THEOS). All rights reserved."
__license__ = "Non-Commercial, End-User Software License Agreement, see LICENSE.txt file"
__version__ = "0.3.0"


class PhonopyCalculation(JobCalculation):
    """
    A generic plugin for calculations based on Phonopy.
    
    Requirement: the node should be able to import phonopy
    """
    
    def _init_internal_params(self):
        super(PhonopyCalculation, self)._init_internal_params()
        
        self._INPUT_FILE_NAME = 'aiida_script.py'
        self._OUTPUT_FILE_NAME = 'aiida.out'
        self._RESULT_FILE_NAME = 'results.json'
        self._default_parser = "phonopy"
        self._INPUT_JSON = 'input.json'
        self._default_is_eigenvectors = False
        
        self._default_tstep = 25.
        self._default_tmin = 0.
        self._default_tmax = 1000.
        
        self._default_parser = 'phonopy'
        
        self._DEFAULT_INPUT_FILE = 'aiida_script.py'
        self._DEFAULT_OUTPUT_FILE = 'aiida.out'

    @classproperty
    def _use_methods(cls):
        """
        Additional use_* methods for the namelists class.
        """
        retdict = JobCalculation._use_methods
        retdict.update({
            "parameters": {
               'valid_types': ParameterData,
               'additional_parameter': None,
               'linkname': 'parameters',
               'docstring': ("Use a node that specifies the input parameters "
                             "for the namelists"),
               },
            "structure": {
               'valid_types': StructureData,
               'additional_parameter': None,
               'linkname': 'structure',
               'docstring': "Use a node for the structure",
               },
            "supercell": {
               'valid_types': StructureData,
               'additional_parameter': None,
               'linkname': 'supercell',
               'docstring': "Use a node for the supercell structure",
               },
            "settings": {
               'valid_types': ParameterData,
               'additional_parameter': None,
               'linkname': 'settings',
               'docstring': ("Use a node that specifies the extra information "
                             "to be used by the calculation"),
               },
            "qpoints": {
               'valid_types': KpointsData,
               'additional_parameter': None,
               'linkname': 'qpoints',
               'docstring': ("Use a node that specifies the qpoints"),
               },
            "force_constants": {
               'valid_types': Forceconstants2Data,
               'additional_parameter': None,
               'linkname': 'force_constants',
               'docstring': "Use a node for the force constants",
               },
            })
        return retdict
    
    def _prepare_for_submission(self,tempfolder, inputdict):        
        """
        This is the routine to be called when you want to create
        the input files and related stuff with a plugin.
        
        :param tempfolder: a aiida.common.folders.Folder subclass where
                           the plugin should put all its files.
        :param inputdict: a dictionary with the input nodes, as they would
                be returned by get_inputdata_dict (without the Code!)
        """
        
        try:
            parameters_data = inputdict.pop(self.get_linkname('parameters'))
        except KeyError:
            pass
            #raise InputValidationError("No parameters specified for this "
            #                           "calculation")
        if not isinstance(parameters_data, ParameterData):
            raise InputValidationError("parameters is not of type "
                                       "ParameterData")
        parameters = parameters_data.get_dict()
        
        try:
            force_constants = inputdict.pop(self.get_linkname('force_constants'))
        except KeyError:
            raise InputValidationError("No force constants specified for this "
                                       "calculation")
        if not isinstance(force_constants,Forceconstants2Data):
            raise InputValidationError("force_constants node is not of type"
                                       "Forceconstants2Data")
        
        try:
            structure = inputdict.pop(self.get_linkname('structure'))
        except KeyError:
            try:
                structure = force_constants.structure
            except AttributeError:
                raise AttributeError("No structure found in input to this calculation")
        if not isinstance(structure,StructureData):
            raise InputValidationError("structure node is not of type"
                                       "StructureData")
        
        try:
            supercell = inputdict.pop(self.get_linkname('supercell'))
        except KeyError:
            try:
                supercell = force_constants.supercell
            except AttributeError:
                raise AttributeError("No supercell found in input to this calculation")
        if not isinstance(supercell,StructureData):
            raise InputValidationError("Supercell node is not of type"
                                       "StructureData")
        
        try:
            qpoints = inputdict.pop(self.get_linkname('qpoints'),None)
        except KeyError:
            raise InputValidationError("Qpoints are missing")
        if qpoints is not None:
            if not isinstance(qpoints, KpointsData):
                raise InputValidationError("qpoints is not of type KpointsData")
        
        try:
            settings = inputdict.pop(self.get_linkname('settings'),None)
        except KeyError:
            pass
        if settings is not None:
            if not isinstance(settings, ParameterData):
                raise InputValidationError("Settings is not of type "
                                           "ParameterData")
            settings_dict = settings.get_dict()
        else:
            settings_dict = {}
            
        try:
            code = inputdict.pop(self.get_linkname('code'))
        except KeyError:
            raise InputValidationError("no code is specified for this calculation")

        ##############################
        # END OF INITIAL INPUT CHECK #
        ##############################
        
        # =================== prepare the python script ========================
        
        # Prepare a json where to dump a lot of information
        input_json = {}
        input_json['cell'] = structure.cell
        input_json['positions'] = [ _.position for _ in structure.sites]
        input_json['symbols'] = [ [k for k in structure.kinds if k.name==_.kind_name  ][0].symbol for _ in structure.sites]
        
        input_txt = "import json\n"
        input_txt += "from phonopy import Phonopy\n"
        input_txt += "from phonopy.structure.atoms import Atoms as PhonopyAtoms\n\n"
        
        # load json with input information
        input_txt += "with open('{}','r') as f:\n".format(self._INPUT_JSON)
        input_txt += "    input_json = json.load(f)\n"
        
        # read data for the structure
        input_txt += "symbols = input_json['symbols']\n"
        input_txt += "cell = input_json['cell']\n"
        input_txt += "positions = input_json['positions']\n"
        
        # Set the Phonopy Atoms
        input_txt += "bulk = PhonopyAtoms(symbols=symbols)\n"
        input_txt += "bulk.set_cell(cell)\n"
        input_txt += "bulk.set_positions(positions)\n"
        
        # Set the Phonopy object
        # read the supercell from input parameters
        phonopy_pre_input = parameters.pop('phonopy_pre_input',None)
        # or from the supercell
        if phonopy_pre_input is None:
            supercell_matrix = supercell.cell
            for i in range(3):
                for j in range(3):
                    ratio = supercell_cell[i][j] / self.structure.cell[i][j]
                    if abs(int(ratio) - ratio)>1e-2:
                        raise ValidationError("The supercell is not an integer "
                                              "multiple of the primitive cell")
                    supercell_matrix[i][j] = int(supercell_cell[i][j] / self.structure.cell[i][j])
            phonopy_pre_input = {'supercell_matrix':supercell_matrix}
        
        this_args = ", ".join( [ '{}={}'.format(k,v) for k,v in 
                                 phonopy_pre_input.iteritems() ])
        
        input_txt += "phonon = Phonopy(bulk, {})\n".format(this_args)
        
        # leave the possibility for this manual setting of displacements
        displacement_dataset = parameters.pop('displacement_dataset', None)
        if displacement_dataset is not None:
            input_json['displacement_dataset'] = displacement_dataset
            input_txt += "displacement_dataset = input_json['displacement_dataset']\n"
            input_txt += "phonon.set_displacement_dataset(displacement_dataset)\n"
        
        # Set the force_constants
        input_json['force_constants'] = force_constants.get_force_constants().tolist()
        input_txt += "force_constants = input_json['force_constants']\n"
        input_txt += "phonon.set_force_constants(force_constants)\n"
        
        
        # ================= set the qpoints and the eigenvectors ===============
        
        try:
            mesh = qpoints.get_kpoints_mesh()[0] # must be in scaled coordinates!
            qpoint_list = None
            input_txt += "phonon.set_mesh({}".format(mesh)
        except AttributeError:
            qpoint_list = qpoints.get_kpoints().tolist()
            mesh = None
            input_json['qpoints'] = qpoint_list
            input_txt += "qpoints = input_json['qpoints']\n"
#            input_txt += "phonon.set_band_structure(qpoints"
        # understand if the user wants the eigenvectors
        is_eigenvectors = parameters.pop('is_eigenvectors', 
                                         self._default_is_eigenvectors)
#        # close the lines left open right before
#        if is_eigenvectors:
#            input_txt += ", is_eigenvectors=True)\n"
#        else:
#            input_txt += ")\n"
        
        
        # ===================== non analytical correction ======================
        
        add_z = parameters.pop("add_nac", None)
        if add_z or add_z is None:
            try:
                input_json['born'] = force_constants.effective_charges
                input_json['dielectric'] = force_constants.dielectric_tensor
                #input_json['factors'] = 14.4 #parameter['factor'] # Here I should understand the conversion factor
                
                input_txt += "born = input_json['born']\n"
                input_txt += "dielectric = input_json['dielectric']\n"
                #input_txt += "factors = input_json['factors']\n"
                # the default value of the conversion factor should treat forces
                # as ev/angstrom and convert frequencies in THz.
                
                input_txt += "phonon.set_nac_params({'born': born, 'factor': factors, 'dielectric': dielectric})\n"
            except AttributeError:
                if add_z:
                    raise InputValidationError("Non analytical correction requested, but no effective "
                                               "charges/dielectric tensor have been found") 
                else:
                    pass
        
        # ======================================================================
        
        # initialize an output json
        input_txt += '\noutput_json = {}\n'
        
        methods = parameters.pop('methods', [])
        if isinstance(methods, basestring):
            methods = [methods]
        #if methods is None:
        #    raise InputValidationError("Must provide the flag 'methods' in the parameters")
        
        # Always get the frequencies
        if mesh is None:
            #input_txt += "q_points, distances, frequencies, eigvecs = phonon.get_band_structure()\n"
            input_txt += 'frequencies = []\n'
            if is_eigenvectors:
                input_txt += 'eigvecs = []\n'
            
            input_txt += "for q in qpoints:\n"
            if is_eigenvectors:
                input_txt += "    fr,ei = phonon.get_frequencies_with_eigenvectors(q)\n"
                input_txt += "    frequencies.append(fr.tolist())\n"
                input_txt += "    eigvecs.append(ei.tolist())\n"
            else:
                input_txt += "    frequencies.append(phonon.get_frequencies(q).tolist())\n"
        else:
            input_txt += "q_points, weights, frequencies, eigvecs = phonon.get_mesh()\n"
            input_txt += "output_json['q_points'] = [_.tolist() for _ in q_points]\n"
            input_txt += "output_json['weights'] = [_.tolist() for _ in weights]\n"
            input_txt += "frequencies = frequencies.tolist()\n"
        input_txt += "output_json['frequencies'] = frequencies\n"
        #input_txt += "output_json['qpoints'] = q_points\n" # I should get them from input kpoints
        if is_eigenvectors:
            input_txt += "output_json['eigenvectors'] = [_.tolist() for _ in eigvecs]\n"

        # get the other results
        for method in methods:
            
            if method == 'DOS':
                if mesh is None:
                    raise InputValidationError("Cannot compute phonon DOS without setting a mesh")
                input_txt += "frequencies_dos, total_dos = phonon.get_total_DOS()\n"
                input_txt += "output_json['frequencies_dos'] = frequencies_dos\n"
                input_txt += "output_json['total_dos'] = total_dos\n"
                
            elif method == 'thermal':
                if mesh is None:
                    raise InputValidationError("Cannot compute thermal properties without setting a mesh")
                
                t_step = parameters.pop('t_step', self._default_tstep)
                t_min = parameters.pop('t_min', self._default_tmin)
                t_max = parameters.pop('t_max', self._default_tmax)
                
                input_txt += "phonon.set_thermal_properties(t_step={}, t_max={}, t_min={})\n".format(t_step,
                                                                                                     t_min,
                                                                                                     t_max)
                input_txt += "temperature, free_energy, entropy, cv = phonon.get_thermal_properties()\n"
                input_txt += "output_json['temperature'] = temperature\n"
                input_txt += "output_json['free_energy'] = free_energy\n"
                input_txt += "output_json['entropy'] = entropy\n"
                input_txt += "output_json['cv'] = cv\n"
            
            elif method == 'velocity':    
                if mesh is None:
                    input_txt += "group_velocities = []\n"
                    input_txt += "for q in qpoints:\n"
                    input_txt += "    group_velocities.append( phonon.get_group_velocity_at_q(q).tolist() )\n"
                else:
                    input_txt += "phonon.set_group_velocity()\n"
                    input_txt += "group_velocities = phonon.get_group_velocity()\n"
                    input_txt += "group_velocities = [_.tolist() for _ in group_velocities]\n"
                input_txt += "output_json['group_velocities'] = group_velocities\n"
        
        input_txt += '\n'
        input_txt += "with open('{}', 'w') as fil:\n".format(self._RESULT_FILE_NAME)
        input_txt += "    json.dump(output_json, fil)\n"
        
        # =========================== dump to file =============================
        
        input_filename = tempfolder.get_abs_path(self._INPUT_FILE_NAME)
        with open(input_filename,'w') as infile:
            infile.write(input_txt)
        
        with open(tempfolder.get_abs_path(self._INPUT_JSON), 'w') as f:
            json.dump(input_json,f)
        
        # ============================ calcinfo ================================
        
        # TODO: look at the qmmm infoL: it might be necessary to put
        # some singlefiles in the directory.
        # right now it has to be taken care in the pre_lines
        
        local_copy_list = []
        remote_copy_list = [] 
        additional_retrieve_list = settings_dict.pop("ADDITIONAL_RETRIEVE_LIST",[])
        
        calcinfo = CalcInfo()
        
        calcinfo.uuid = self.uuid
        # Empty command line by default
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        
        # Retrieve files
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self._OUTPUT_FILE_NAME)
        calcinfo.retrieve_list.append(self._RESULT_FILE_NAME)
        calcinfo.retrieve_list += additional_retrieve_list
        
        codeinfo = CodeInfo()
        codeinfo.cmdline_params = [self._INPUT_FILE_NAME]
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
     #   codeinfo.stderr_name = self._DEFAULT_ERROR_FILE
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        return calcinfo

