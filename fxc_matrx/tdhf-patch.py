
# If you only want to print only the Hartree or fxc term of A-matrix element, 
# you can comment out the other terms in pyscf/scf/_response_functions.py

# Add this to line 985 of pyscf/tdscf/rhf.py (before the return of `vind`)
        #hack: output the full matrix
        nvir_tmp=3 #H2, nocc=1
        vmat=numpy.zeros([nz, nocc*nvir_tmp, nocc*nvir_tmp])
        for io in range(nocc):
            for iv in range(nvir_tmp):
                x=numpy.zeros([nz,nocc,nvir_tmp])
                x[:, io, iv]=1.0
                dms_test=lib.einsum('xov,qv,po->xpq', x*2, orbv.conj()[:, :nvir_tmp], orbo) #transition density matrix
                v1ao_test = vresp(dms_test) #v=f*dm
                v1ov_test = lib.einsum('xpq,po,qv->xov', v1ao_test, orbo.conj(), orbv[:,:nvir_tmp]) 
                vmat[:,io*nvir_tmp+iv,:]=v1ov_test.reshape(nz,-1) #(ai|jb)矩阵的一行（所有jb）
        print("vmat=",vmat)
        
        
# Add this to line 809 of pyscf/tdscf/uhf.py (before the return of `vind`)
        #hack: output the full matrix
        nvir_tmp=3 #H2, nocc=1
        vmat=numpy.zeros([nz, (nocca+noccb)*nvir_tmp, (nocca+noccb)*nvir_tmp])
        for ispin in range(1):
            for io in range(nocca): #nacc(ispin)
                for iv in range(nvir_tmp):
                    xa=numpy.zeros([nz,nocca,nvir_tmp])
                    xb=numpy.zeros([nz,noccb,nvir_tmp])
                    if ispin == 0:
                        xa[:,io,iv]=1.0
                    else:
                        xb[:,io,iv]=1.0
                    dmsa_test=lib.einsum('xov,qv,po->xpq', xa, orbva.conj()[:, :nvir_tmp], orboa) #transition density matrix
                    dmsb_test=lib.einsum('xov,qv,po->xpq', xb, orbvb.conj()[:, :nvir_tmp], orbob) #transition density matrix
                    v1ao_test = vresp(numpy.asarray((dmsa_test,dmsb_test)))
                    
                    v1aov_test = lib.einsum('xpq, po,qv->xov', v1ao_test[0], orboa.conj(), orbva[:,:nvir_tmp])# a-part row of A
                    v1bov_test = lib.einsum('xpq, po,qv->xov', v1ao_test[1], orbob.conj(), orbvb[:,:nvir_tmp])# b-part row of A
                    vmat[:, ispin*nocca*nvir_tmp +io*nvir_tmp+iv,:nocca*nvir_tmp]=v1aov_test.reshape(nz,-1) #(ai|jb)矩阵的一行（所有jb）
                    vmat[:,  ispin*nocca*nvir_tmp+ io*nvir_tmp+iv, nocca*nvir_tmp:] = v1bov_test.reshape(nz, -1)  # (ai|jb)矩阵的一行（所有jb）
        print("vmat=",vmat)