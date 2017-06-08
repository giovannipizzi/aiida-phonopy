from aiida.orm.data.structure import StructureData
from aiida.orm.data.array import ArrayData
from aiida.orm.data import Data
import numpy
from aiida.common.exceptions import ValidationError

class Forceconstants2Data(Data):
    """
    Docs missing
    """
    _array_prefix = "array|"
    _array_name = 'force_constants' 
    
    def __init__(self,*args,**kwargs):
        super(Forceconstants2Data, self).__init__(*args, **kwargs)
        self._cached_arrays = {}
    
    def _validate_structure(self, struc, type='Structure'):
        if not isinstance(struc, StructureData):
            raise ValidationError("{} must be an AiiDA StructureData".format(type))
        if struc.pk is None:
            raise ValidationError("{} passed to Forceconstants must be stored".format(type))
        return struc
    
    def _set_structure(self, struc, type='structure'):
        if type not in ['structure','supercell']:
            raise ValidationError('Type not understood')
        validated_structure = self._validate_structure(struc, type=type.capitalize())
        self._set_attr('reference_{}_uuid'.format(type), validated_structure.uuid)
        
    def _get_reference_structure(self, type='structure'):
        if type not in ['structure','supercell']:
            raise ValidationError('Type not understood')
        struc_uuid = self.get_attr('reference_{}_uuid'.format(type))
        q = StructureData.query(uuid=struc_uuid)
        if len(q) != 1:
            raise AttributeError("Should have found one reference structure."
                                 "Found instead {} objects".format(len(q)))
        return q[0] 
    
    def _validate_force_constants(self, in_constants):
        out_constants = numpy.array(in_constants, dtype='float64')
        shape = out_constants.shape
        if any([shape[0]!=shape[1], shape[2]!=3]):
            raise ValidationError('force_constants must be an array NxNx3x3')
        return out_constants
        
    def set_force_constants(self, force_constants, 
                            structure=None,
                            supercell=None,
                            dielectric_tensor=None,
                            effective_charges=None):
        validated_constants = self._validate_force_constants(force_constants)
        
        self._set_array(validated_constants)
        
#        if supercell is not None:
#            self._set_supercell(supercell)
        
        if structure is not None:
            self._set_structure(structure)
            
        if supercell is not None:
            self._set_structure(supercell, type='supercell')
            
        if dielectric_tensor is not None:
            self.dielectric_tensor = dielectric_tensor
        
        # set the charges always after the structure
        # because it will check with the number of atoms
        if effective_charges is not None:
            self.dielectric_charges = dielectric_charges
    
    def _delete_array(self):
        """
        Delete an array from the node. Can only be called before storing.
        
        :param name: The name of the array to delete from the node.
        """
        import numpy
        name = self._array_name

        fname = '{}.npy'.format(name)
        if fname not in self.get_folder_list():
            raise KeyError("Array with name '{}' not found in node pk= {}".format(
                name, self.pk))
        
        # remove both file and attribute
        self.remove_path(fname)
        try:
            self._del_attr("{}{}".format(self._array_prefix, name))
        except (KeyError, AttributeError):
            # Should not happen, but do not crash if for some reason the 
            # property was not set.
            pass
        
    def get_force_constants(self):
        return self._get_array()
    
    @property
    def supercell(self):
        """
        The supercell
        :return: a StructureData
        """
        return self._get_reference_structure(type='supercell')
     
    @property
    def structure(self):
        """
        The reference primitive crystal structure
        :return: a StructureData
        """
        return self._get_reference_structure(type='structure')

#    @supercell.setter
#    def supercell(self, supercell):
#        validated_supercell = self._validate_structure(supercell, type='Supercell')
#        self._set_structure(self, validated_supercell, type='supercell')
        
    @property
    def dielectric_tensor(self):
        """
        The dielectric tensor matrix.
        :return: a 3x3 tuple
        """
        return self.get_attr('dielectric_tensor')
    
    @dielectric_tensor.setter
    def dielectric_tensor(self, dielectric_tensor):
        validated_tensor = self._validate_dielectric_tensor(dielectric_tensor)
        self._set_attr('dielectric_tensor', valildated_tensor)

    def _validate_dielectric_tensor(self, in_tensor):
        out_tensor = numpy.array(in_tensor, dtype='float64')
        if out_tensor.shape != (3,3):
            raise ValidationError("Dielectric tensor should be a 3x3 float tensor")
        return out_tensor
        
    @property
    def effective_charges(self):
        """
        The effective charges for each atom.
        :return: a list of number_of_atoms elements, each being a 3x3 tuple
        """
        return self.get_attr('effective_charges')

    @effective_charges.setter
    def effective_charges(self, effective_charges):
        validated_charges = self._validate_charges(effective_charges)
        self._set_attr('effective_charges', valildated_charges)

    def _validate_charges(self, in_charges):
        out_charges = numpy.array(in_charges, dtype='float64')
        shape = out_charges.shape
        if shape[1] != 3 or len(shape)!=2:
            raise ValidationError("Dielectric tensor should be a Nx3 float tensor")
        
        try:
            struc = self._get_reference_structure()
            if len(struc.sites) != shape[0]:
                raise ValidationError("The number of charges does not match the"
                                      "number of atoms of the structure")
        except AttributeError:
            # no structure was set, cannot know the number of atoms
            pass
        return out_charges
        
    def _get_array(self):
        """
        Return an array stored in the node
        
        :param name: The name of the array to return.
        """
        import numpy
        name = self._array_name
        
        # raw function used only internally
        def get_array_from_file(self, name):
            fname = '{}.npy'.format(name)
            if fname not in self.get_folder_list():
                raise KeyError(
                    "Array with name '{}' not found in node pk= {}".format(
                    name, self.pk))
        
            array = numpy.load(self.get_abs_path(fname))
            return array

        
        # Return with proper caching, but only after storing. Before, instead,
        # always re-read from disk
        if self._to_be_stored:
            return get_array_from_file(self,name)
        else:
            if name not in self._cached_arrays:
                self._cached_arrays[name] = get_array_from_file(self,name)
            return self._cached_arrays[name]
    
    def clear_internal_cache(self):
        """
        Clear the internal memory cache where the arrays are stored after being
        read from disk (used in order to reduce at minimum the readings from
        disk).
        This function is useful if you want to keep the node in memory, but you
        do not want to waste memory to cache the arrays in RAM.
        """
        self._cached_arrays = {}
    
    def _set_array(self, array):
        """
        Store a new numpy array inside the node. Possibly overwrite the array
        if it already existed.
        
        Internally, it stores a name.npy file in numpy format.
        
        :param name: The name of the array.
        :param array: The numpy array to store.
        """
        import re
        import tempfile
        import numpy
        
        name = self._array_name
    
        if not(isinstance(array, numpy.ndarray)):
            raise TypeError("ArrayData can only store numpy arrays. Convert "
                            "the object to an array first")
        
        # Check if the name is valid
        if not(name) or re.sub('[0-9a-zA-Z_]', '', name):
            raise ValueError("The name assigned to the array ({}) is not valid,"
                             "it can only contain digits, letters or underscores")
                
        fname = "{}.npy".format(name)
        
        with tempfile.NamedTemporaryFile() as f:
            # Store in a temporary file, and then add to the node
            numpy.save(f, array)
            f.flush() # Important to flush here, otherwise the next copy command
                      # will just copy an empty file
            self.add_path(f.name, fname)
        
        # Mainly for convenience, for querying purposes (both stores the fact
        # that there is an array with that name, and its shape)
        self._set_attr("{}{}".format(self._array_prefix, name), list(array.shape))
    
    def export(self, path='.', filename='force_constants.dat', format='qe'):
        """
        Requires ASE (just to get scaled coordinates)
        """
        if format=='qe':
            from aiida.common.constants import bohr_to_ang
            import os
            import numpy
             
            lines = []
            structure = self.structure
            num_atoms = len(structure.sites)
            num_species = len(structure.kinds)
            
# check the order of next line:
            lines.append( "{}  {}  0  {}  ".format(num_species, num_atoms, cell[0][0]/bohr_to_ang
                                                   ) + "  ".join(['0.0000000']*5) +'\n')
            for i in structure.cell:
                lines.append( "  {}  {}  {}\n".format(*i) )
             
            for i,k in enumerate(s.kinds):
                lines.append( "{}  '{3s}'  {}\n".format( i+1, k.name, k.mass*amu_to_ry ) )
            
            pos = s.get_ase().get_scaled_positions() # laziness!
            for i,site in s.sites:
                this_pos = pos[i]
                lines.append( "{}  {}  {}  {}  {}\n".format(i+1, 
                                                            [ _.name for _ in s.kinds ].index(site.kind_name)+1, 
                                                            *this_pos) )
            try:
                dielectric_tensor = self.get_dielectric_tensor()
                charges = self.get_effective_charges()
                lines.append('T\N')
                for i in dielectric_tensor:
                    lines.append( "  {}  {}  {}\n".format(*i) )
                     
                for n,z in enumerate(charges):
                    lines.append("{}\n".format(i+1))
                    for i in z:
                        lines.append( "  {}  {}  {}\n".format(*i) )
            
            except AttributeError:
                lines.append('F\N')
            
            supercell = self.supercell
            structure = self.structure
            
            size = [ nint(supercell.cell[i][i]/structure.cell[i][i]) for i in range(3) ]
            if any([ abs(size[i] - supercell.cell[i][i]/structure.cell[i][i])>0.1 for i in range(3) ]):
                raise ValidationError("The supercell cell is not LxMxN times larger than the structure cell")
            
            lines.append("{}  {}  {}\n".format(*size))
            
            # now, for every atom, recognize the basis index and to which cell index belongs
            
            # laziness!
            ase_super = supercell.get_ase()
            ase_structure = structure.get_ase()
            supercell_dict = {}
            for this_position in ase_super.get_scaled_positions():
                # the following relies on the fact that scaled coordinates are always between 0 and 1
                if any( [ _ >= 1. for _ in this_position] ):
                    raise ValueError("scaled coordinates should be between 0 and 1")
                multiplied_pos = [ i[0]*i[1] for i in zip(this_position, size) ]
                cell_index = [ numpy.int(_) for _ in multiplied_pos ] # index from 0 to N-1
                folded_pos = [i[0] - i[1] for i in zip(multiplied_pos, cell_index)]
                
                # get the basis index
                positions_primitive = ase_structure.get_scaled_positions()
                distances = [numpy.linalg.norm(folded_pos, _) for _ in positions_primitive]
                basis_index = distances.index(min(distances)) # index from 0 to N-1
                
                # for the ith atom in the list of supercell sites
                # tells the unit cell it belongs to 
                in_primitive = True if cell_index==(0,0,0) else False
                supercell_dict[i] = [ this_position, basis_index, in_primitive ]
            
            force_constants = get_force_constants()
            
            # Fortran style!
            
            id = [None]*4
            for cart1 in range(3):
                id[2] = cart1
                for cart2 in range(3):
                    id[3] = cart2
                    for b1 in range(num_atoms):
                        id[0] = [ _[0] for _ in supercell_dict.iteritems() if _[1][2]==True and _[1][1]==b1 ][0]
                        for b2 in range(num_atoms):
                            lines.append("{} {} {} {}\n".format(cart1,cart2,b1,b2))
                            for sup3 in range(supercell[2]):
                                for sup2 in range(supercell[1]):
                                    for sup1 in range(supercell[0]):
                                        
                                        id[1] = [ _[0] for _ in supercell_dict.iteritems() if _[1][1]==(sup1,sup2,sup3) ][0]
                                        lines.append("{} {} {} {}\n".format(sup1+1,sup2+1,sup3+1,
                                                                            force_constants[id[0],id[1],id[2],id[3]]))
            
            with open( os.path.join(path, filename), 'w') as f:
                f.write(lines)
            
        elif format=='phonopy':
            from phonopy import Phonopy
            import copy
            bulk = self.structure.export(format='phonopy')
            supercell_cell = self.supercell.cell
            supercell_matrix = copy.copy(supercell_cell)
            for i in range(3):
                for j in range(3):
                    supercell_matrix[i][j] = supercell_cell[i][j] / self.structure.cell[i][j]
            phonon = Phonopy(bulk, supercell_matrix) 
            
            phonon.set_force_constants(self.get_force_constants())
            
            return phonon
            
        else:
            raise ValueError('Format {} is not recognized'.format(format))
            
    def read_qe_matrix(self, filepath, structure=None):
        raise NotImplementedError   
        with open(filepath) as f:
            lines = f.readlines()
        
        for count, line in enumerate(lines):
        
            if line.strip().capitalize()=='F' or line.strip().capitalize()=='T':
                with_charge_line = count
                with_charges = bool(line.strip().capitalize())
                break
        
        starting_line = with_charge_line
        
        if with_charges:        
            for count, line in enumerate(lines[with_charge_line:]):
                if True:
                    pass
            starting_line += count
            
            
        
        
            