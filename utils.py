import Bio.PDB
import Bio.Align 
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
import py3Dmol
import matplotlib.pyplot as plt

def align_on_ref(model_file, ref_file, ligname = None):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    ref_structure = pdb_parser.get_structure("reference", ref_file)
    model_structure = pdb_parser.get_structure("model", model_file)

    ref_seq = ''.join([aa3to1.get(i.resname, None) for i in ref_structure.get_residues() if aa3to1.get(i.resname, None) is not None])
    model_seq = ''.join([aa3to1.get(i.resname, 'X') for i in model_structure.get_residues()])

    seq_aligner = Bio.Align.PairwiseAligner()
    seq_align = seq_aligner.align(ref_seq, model_seq)
    aligner =  Bio.PDB.StructureAlignment(seq_align[0], ref_structure, model_structure)
    map_0, map_1 = aligner.get_maps()

    ref_atoms = []
    model_atoms = []

    for k,v in map_0.items():
        if k is None or v is None:
            continue
        else:
            ref_atoms.append(k['CA'])
            model_atoms.append(v['CA'])

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, model_atoms)
    super_imposer.apply(model_structure.get_atoms())

    io = Bio.PDB.PDBIO()

    if ligname is not None:
        for res in ref_structure.get_residues():
            if res.resname == ligname:
                coords = res.center_of_mass()
                break

        add_ghost_atom(ref_structure, coords)
        add_ghost_atom(model_structure, coords)
        io.set_structure(ref_structure) 
        io.save(ref_file.split('.')[0]+'_aligned.pdb')

    io.set_structure(model_structure) 
    io.save(model_file.split('.')[0]+'_aligned.pdb')
    

def add_ghost_atom(structure, coords):
    models = [m for m in structure.get_models()]
    chains = [c for c in models[0].get_chains()]
    chains[0].add(Bio.PDB.Residue.Residue(('GST', 9999, ' '), 'GST', ''))
    for r in chains[0].get_residues():
        if r.resname == 'GST':
            r.add(Bio.PDB.Atom.Atom('G', coords, 0, 1, ' ', 'G', '999', element = 'C'))


def show_aligned(model_file,
                ref_file,
                ligname = None,
                model_cartoon_color = 'red',
                ref_cartoon_color = 'blue',
                ligand_color = 'magenta',
                model_sidecahin_color = 'OrangeCarbon',
                ref_sidecahin_color = 'OrangeCarbon'):
    view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js', width = 800, height = 600)
    view.addModel(open(model_file,'r').read(),'pdb')
    view.addModel(open(ref_file,'r').read(),'pdb')

    #find ligand center of mass
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    ref_structure = pdb_parser.get_structure("reference", ref_file)
    
            
    sele = {'resn': 'GST'}

    view.setStyle({'model':0},{'cartoon': {'color':model_cartoon_color}})
    view.setStyle({'model':1},{'cartoon': {'color':ref_cartoon_color}})

    view.zoomTo()

    if ligname:
        selection_0 = {'model': 0, 'byres': 'true', 'within':{'distance': 10, 'sel': sele}}
        selection_1 = {'model': 1, 'byres': 'true', 'within':{'distance': 10, 'sel': sele}}
        view.addStyle(selection_0,{'stick':{'colorscheme':model_sidecahin_color,'radius':0.3}})
        view.addStyle(selection_1,{'stick':{'colorscheme':ref_sidecahin_color,'radius':0.3}})
        view.setStyle({'resn': ligname}, {'stick':{'colorscheme':ligand_color,'radius':0.3}})
        view.setClickable({},
                        'true',
                        "function(atom, viewer, event, container){"\
                        "if(atom.label){viewer.removeLabel(atom.label);delete atom.label;}" \
                        "else{atom.label=viewer.addLabel(atom.resn+atom.resi+':'+atom.atom, {'position': atom});}\n}")
        view.zoomTo({'resn': ligname})

    return view
    

def plot_chain_legend(dpi=100,
                        name_model = 'AF',
                        name_ref = 'Ref.',
                        model_cartoon_color = 'red',
                        ref_cartoon_color = 'blue',
                        ligand_color = 'magenta',
                        model_sidecahin_color = 'OrangeCarbon',
                        ref_sidecahin_color = 'OrangeCarbon'):
  thresh = [name_model,name_ref,f'{name_model} sidechains',f'{name_ref} sidechains','Ligand']
  plt.figure(figsize=(1,0.1),dpi=dpi)
  ########################################
  for c in [model_cartoon_color,ref_cartoon_color,model_sidecahin_color,ref_sidecahin_color,ligand_color]:
    plt.bar(0, 0, color=c)
  plt.legend(thresh, frameon=False,
             loc='center', ncol=5,
             handletextpad=1,
             columnspacing=1,
             markerscale=0.5,)
  plt.axis(False)
  return plt