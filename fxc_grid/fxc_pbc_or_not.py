from pyscf.pbc import gto, scf
from pyscf.tools import cubegen

cell = gto.M(
    atom = '''H 10 10 10.37
              H 10 10 9.63''',
    a = '''20 0 0
           0 20 0
           0 0 20''',
    pseudo = 'gth-pade',
    basis = 'gth-szv'
    )
mf = scf.KUKS(cell, kpts = cell.make_kpts([1,1,1]))
mf.xc='lda'
mf.kernel()
fxc1=cubegen.fxc_ao(cell, "H2-fxc-PZ-PBC", mf.make_rdm1(), 180,180,180)
# max(fxc)= 0.0
# min(fxc)= -4179557448.3036504

from pyscf import gto, scf
mol = gto.Mole(
    atom = '''H 10 10 10.37
              H 10 10 9.63''',
    pseudo = 'gth-pade',
    basis = 'gth-szv'
    )
mf1=scf.UKS(mol)
mf1.xc='lda'
mf1.kernel()
fxc2=cubegen.fxc_ao(mol, "H2-fxc-PZ", mf.make_rdm1(), 80,80,80)
#max(fxc)= -0.22391758138188966
#min(fxc)= -4899.047740238133