# -*- coding: utf-8 -*-
from aiida.orm.calculation.inline import make_inline
from aiida.orm import DataFactory, CalculationFactory
from aiida.orm.workflow import Workflow

StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')

@make_inline
def create_supercells_inline(**kwargs):
    """
    1) export input structure to PhononpyAtoms 
    2) phonon = Phonopy(bulk,
                        supercell,
                        distance=0.01)
    3) phonon.get_supercells_with_displacements()
    4) phonon.get_supercell() (maybe)
    """
    from phonopy import Phonopy
    
    phonopy_pre_input = kwargs.pop('phonopy_pre_input')
    structure = kwargs.pop('structure')
    
    bulk = structure.convert(object_format='phonopyatoms')
    extra_args = phonopy_pre_input.get_dict()
    
    if kwargs:
        raise ValueError("Some kwargs have not been recognized")
    
    phonon = Phonopy(bulk, **extra_args )
    
    supercell_phonoatoms = phonon.get_supercells_with_displacements()
    
    results = {}
    for i,s in enumerate(supercell_phonoatoms):
        structuredata = StructureData(cell=s.get_cell())
        for symbol, position in zip( s.get_chemical_symbols(),
                                     s.get_positions() ):
            structuredata.append_atom(position=position, symbols=[symbol])
        
        results['structure_{}'.format(i)] = structuredata
    
    # save also the bare supercell, which is needed for the force_constants matrix
    phonopy_supercell = phonon.get_supercell()
    supercell = StructureData(cell=phonopy_supercell.get_cell())
    for symbol, position in zip( phonopy_supercell.get_chemical_symbols(),
                                 phonopy_supercell.get_positions() ):
        supercell.append_atom(position=position, symbols=[symbol])        
    results['supercell'] = supercell
    
    return results
    
    
@make_inline
def create_force_matrices_inline(**kwargs):
    """
    1) extract forces from the ArrayDatas
    2) phonon.set_forces(sets_of_forces)
       phonon.produce_force_constants()
    """
    from phonopy import Phonopy
    from aiida_phonopy.data.forceconstants2 import Forceconstants2Data
    arraydatas = []
    for k in sorted( [ _ for _ in kwargs.keys() 
                       if _.startswith('array_with_forces_') ]):
        arraydatas.append( kwargs.pop(k) )
        
    sets_of_forces = []
    for arr in arraydatas:
        try:
            sets_of_forces.append( arr.get_array('forces')[-1] )
        except KeyError:
            raise KeyError("Forces not found in the ArrayData")
    
    phonopy_pre_input = kwargs.pop('phonopy_pre_input').get_dict()
    
    structure = kwargs.pop('structure')
    supercell = kwargs.pop('supercell')
    bulk = structure.convert(object_format='phonopyatoms')
    
    if kwargs:
        raise ValueError("Some kwargs have not been recognized")
    
    phonon = Phonopy(bulk, **phonopy_pre_input)
    
    phonon.set_forces(sets_of_forces)
    phonon.produce_force_constants()
    
    force_constants = phonon.get_force_constants()
    force_data = Forceconstants2Data()
    force_data.set_force_constants(force_constants,
                                   structure=structure,
                                   supercell=supercell)
    
    return {'force_constants': force_data}

#(if I already had the matrices, simply phonon.set_force_constants())

#(I should create a force constants matrix object that stores a NxNx3x3 array,
#and that than exports this to either QuantumEspresso matrices or to phonopy matrices (simpler))


class FrozenphononWorkflow(Workflow):
    """
    General workflow to launch Phonon dispersion calculations with phonopy
    Requires:
    1) Phonopy installed on the AiiDA server (the supercell creations are run locally)
    2) Phonopy installed on a remote computer as an AiiDA code
    
    Docs missing
    """
    _clean_workdir = False
    _use_qgrid_parallelization = False
    _default_forces_subworkflow = 'quantumespresso.pw'
    
    def __init__(self,**kwargs):
        super(FrozenphononWorkflow, self).__init__(**kwargs)
    
    @Workflow.step
    def start(self):
        # check input
        wf_params = self.get_parameters()
        
        if 'force_constants' in wf_params.keys(): 
            self.next(self.postprocessing)
        else:
            self.next(self.supersize)
            
            
    @Workflow.step
    def supersize(self):
        from aiida.orm.data.array.kpoints import KpointsData
        from aiida.common.exceptions import MissingPluginError
        from aiida.orm.code import Code
        
        wf_params = self.get_parameters()
        
        # run the phonopy inline supercell creator
        inline_params = {"structure": wf_params['structure'],
                         "phonopy_pre_input": wf_params['phonopy_pre_input'],
                         } 
        
        _, ret_dict = create_supercells_inline(**inline_params)
        
        
        
        initial_kpoints_mesh,offset = wf_params['pw_kpoints'].get_kpoints_mesh()
        supercell_matrix = wf_params['phonopy_pre_input'].get_dict()['supercell_matrix']
        diag = [ supercell_matrix[0][0], supercell_matrix[1][1], supercell_matrix[2][2] ]
        mesh = [ max(1, _[0]/_[1] ) for _ in zip(initial_kpoints_mesh, diag) ]
        # check: the mesh of the supercell is commensurate with the primitive
        if any( [ abs(_[0]/_[1] - float(_[0])/float(_[1] ))>1.e-4 
                  for _ in zip(initial_kpoints_mesh, diag) ] ):
            self.append_to_report("Careful: the kpoint mesh of the primitive "
                                  "cell is not commensurate with the supercell")
        
        sup_kpoints = KpointsData()
        sup_kpoints.set_kpoints_mesh(mesh=mesh,offset=offset)
        sup_kpoints.store()
        
        # launch the pw calculation for the structures created in output
        
        for k in sorted( [ _ for _ in ret_dict.iterkeys() if _.startswith('structure_')] ):
            this_structure = ret_dict[k]
            
            pw_parameters = wf_params['pw_parameters']
            try:
                pw_parameters['CONTROL']
            except KeyError:
                pw_parameters['CONTROL'] = {}    
            pw_parameters['CONTROL']['calculation'] = 'scf'
            pw_parameters['CONTROL']['tprnfor'] = True
            
            # understand if I want to use a subworkflow to manage the force 
            # calculations
            try:
                subworkflow = wf_params['input'].get('forces_subworkflow',
                                                     self._default_forces_subworkflow)
                if not subworkflow:
                    subworkflow = None
            except KeyError:
                subworkflow = None
                
            if subworkflow is not None:
                try:
                    from aiida.orm import WorkflowFactory
                    PwWorkflow = WorkflowFactory(subworkflow)
                    
                    new_wf_params= {'structure': this_structure,
                                    'pseudo_family': wf_params['pseudo_family'],
                                    'codename': wf_params['pw_codename'],
                                    'parameters': pw_parameters,
                                    'calculation_set': wf_params['pw_calculation_set'],
                                    'kpoints': sup_kpoints, #wf_params['pw_kpoints'],
                                    }
                    try:
                        new_wf_params['settings'] = wf_params['pw_settings']
                    except KeyError:
                        pass
                    
                    wf_pw = PwWorkflow(params=new_wf_params)
                    
                    wf_pw.store()
                    self.attach_workflow(wf_pw)
                    wf_pw.start()
                except MissingPluginError: # couldn't load the workflow
                    self.append_to_report("Couldn't load the requested subworkflow "
                                          "({}) ,falling back to simple calcs".format(wf_params['input']['use_subworkflow']))
                    subworkflow = None
                
            if subworkflow is None: # don't put an elif!
                code = Code.get_from_string(wf_params['pw_codename'])
                calc = code.new_calc()
                
                PwCalculation = CalculationFactory('quantumespresso.pw')
                if not isinstance(calc,PwCalculation):
                    raise ValueError("The code found is not configured for the "
                                     "quantumespresso.pw plugin")
                
                calc.use_parameters(ParameterData(dict=pw_parameters))
                calc.use_structure(this_structure)
                calc.use_pseudos_from_family(wf_params['pseudo_family'])
                calc.use_kpoints(sup_kpoints)
                if not isinstance(wf_params['pw_calculation_set'],dict):
                    raise ValueError("set must be a dictionary, found instead a {}".format(wf_params['pw_calculation_set'].__class__))
                try:
                    calc.set(**wf_params['pw_calculation_set'])
                except ValueError as e:
                    raise ValueError("A key in params['pw_calculation_set'] does not correspond"
                                     " to any calc.set_[key] method\n{}".format(e.message))
                try:
                    calc.use_settings(ParameterData(dict=wf_params['pw_settings']))
                except KeyError:
                    pass
                
                calc.store_all()
                self.attach_calculation(calc)
                
        self.next(self.recollect)
        
        
    @Workflow.step
    def recollect(self):
        from aiida.orm.calculation.inline import InlineCalculation
        
        wf_params = self.get_parameters()
        
        all_supercells_wf = self.get_step(self.supersize).get_sub_workflows()
        if not all_supercells_wf:
            # case where we launched calculations
            calculations = self.get_step_calculations(self.supersize)
            structures = [ _.inp.structure for _ in calculations ]
        else:
            all_supercells_wf = self.get_step(self.supersize).get_sub_workflows()
            structures = [ _.get_parameter('structure') for _ in all_supercells_wf ]
            calculations = [ _.get_result('pw_calculation') for _ in all_supercells_wf ]
        
        arraydatas = {}
        for s,pw_calculation in zip(structures,calculations):
            inp_dict = s.get_inputs_dict()
            k =  [ _ for _ in inp_dict.keys() if _.startswith('structure') ]
            if len(k)>1:
                raise AttributeError('Do not know how to align this workflow with the phonopy')
            index_order = int(k[0].split('_')[-1])
            arraydatas['array_with_forces_{}'.format(index_order)] = pw_calculation.out.output_array
        
        # s is the last structure to come from the previous loop
        inline = s.get_inputs()[0]
        if not isinstance(inline, InlineCalculation):
            raise AttributeError("An expected inlinecalculation was not found")
        supercell = inline.out.supercell
        
        _, ret_dict = create_force_matrices_inline(structure=wf_params['structure'],
                                                   phonopy_pre_input=wf_params['phonopy_pre_input'],
                                                   supercell=supercell,
                                                   **arraydatas)
        self.add_result('force_constants', ret_dict['force_constants'])
        
        self.next(self.postprocessing)
    
    
    @Workflow.step
    def postprocessing(self):
        from aiida.orm.code import Code
        wf_params = self.get_parameters()
        
        try:
            force_constants = wf_params['force_constants']
        except KeyError:
            force_constants = self.get_result('force_constants')
        
        code = Code.get_from_string(wf_params['phonopy_codename']) 
        calc = code.new_calc()
        
        phonopy_parameters_from_input = wf_params.pop('phonopy_parameters',{})
        if isinstance(phonopy_parameters_from_input,ParameterData):
            phonopy_parameters = phonopy_parameters_from_input.get_dict()
        else:
            phonopy_parameters = phonopy_parameters_from_input
        
        phonopy_pre_input = wf_params['phonopy_pre_input']
        if isinstance(phonopy_pre_input,ParameterData):
            pre_input_dict = phonopy_pre_input.get_dict()
        else:
            pre_input_dict = phonopy_pre_input
        phonopy_parameters['phonopy_pre_input'] = pre_input_dict
        
        phonopy_parameters_data = ParameterData(dict=phonopy_parameters).store()
        
        calc.use_parameters(phonopy_parameters_data)
        calc.use_qpoints(wf_params['phonopy_qpoints'])
        calc.use_force_constants(force_constants)
        
        phonopy_settings = wf_params.pop('phonopy_settings', None)
        if phonopy_settings is not None:
            if not isinstance(phonopy_settings,ParameterData):
                phonopy_settings_data = ParameterData(dict=phonopy_settings).store()
            else:
                phonopy_settings_data = phonopy_settings
            
            calc.use_settings(phonopy_settings_data)
        
        try:
            calc.set(**wf_params['phonopy_calculation_set'])
        except ValueError as e:
            raise ValueError(
                "A key in 'phonopy_calculation_set'] does not correspond"
                " to any calc.set_[key] method\n{}".format(e.message))

        calc.set_withmpi(False)
        calc.store_all()
        
        self.attach_calculation(calc)
        
        self.next(self.final_step)
        
    @Workflow.step
    def final_step(self):
        post_calc = self.get_step_calculations(self.postprocessing)[0]
        
        something_found = False
        # here append the results in output
        try:
            self.add_result('phonon_dispersion', post_calc.out.phonon_frequencies)
            something_found = True
        except AttributeError:
            pass
        try:
            self.add_result('phonon_dos', post_calc.out.phonon_dos)
            something_found = True
        except AttributeError:
            pass
        try:
            self.add_result('thermal_properties', post_calc.out.thermal_properties)
            something_found = True
        except AttributeError:
            pass
        
        if something_found:
            self.append_to_report("Finished correctly")
        else:
            self.append_to_report("Nothing produced...")
        self.next(self.exit)
        
