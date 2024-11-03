# add this to pyscf/tools/cubegen.py

def fxc_gto(mol, outfile, dm, nx=80, ny=80, nz=80):
    from pyscf.pbc.gto import Cell
    cc = Cube(mol, nx, ny, nz)
    GTOval = 'GTOval'
    if isinstance(mol, Cell):
        GTOval = 'PBC' + GTOval   
    coords = cc.get_coords()
    ao_value = mol.eval_gto(GTOval, coords)
    rho_u = numint.eval_rho(mol, ao_value, dm[0], xctype='LDA')
    rho_d = numint.eval_rho(mol, ao_value, dm[1], xctype='LDA')
    cc.write(rho_u.reshape(cc.nx,cc.ny,cc.nz), outfile+"rhou_mol.cube", comment='rho_u')
    cc.write(rho_d.reshape(cc.nx,cc.ny,cc.nz), outfile+"rhod_mol.cube", comment='rho_d')
    fxc=dft.libxc.eval_xc('LDA,PZ', (rho_u, rho_d), spin=1, deriv=2)[2][0]
    sp_map=("aa","ab","bb")
    for ispin in range(3):
        cc.write(fxc[:,ispin].reshape(cc.nx,cc.ny,cc.nz),outfile+str(ispin)+".cube" ,comment="fxc_"+f'{sp_map[ispin]}')
    return fxc

def fxc_ao(mol, outfile, dm, nx=80, ny=80, nz=80):
    cc = Cube(mol, nx, ny, nz)
    # Compute density on the .cube grid
    coords = cc.get_coords()
    ao_value = numint.eval_ao(mol, coords)
    rho_u = numint.eval_rho(mol, ao_value, dm[0], xctype='LDA')
    rho_d = numint.eval_rho(mol, ao_value, dm[1], xctype='LDA')
    fxc=dft.libxc.eval_xc('LDA,PZ', (rho_u, rho_d), spin=1, deriv=2)[2][0]
    print("max(fxc)=",numpy.max(fxc))
    print("min(fxc)=",numpy.min(fxc))
    sp_map=("aa","ab","bb")
    for ispin in range(3):
        cc.write(fxc[:,ispin].reshape(cc.nx,cc.ny,cc.nz),outfile+str(ispin)+".cube" ,comment="fxc_"+f'{sp_map[ispin]}')
    return fxc