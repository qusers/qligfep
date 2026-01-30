load PDB
cmd.mouse('forward')
cmd.edit('(COMPLEX`TER)',None,None,None,pkresi=0,pkbond=0)
PRE
editor.attach_amino_acid('pk1', 'nme')
POST
cmd.edit('(COMPLEX`1)',None,None,None,pkresi=0,pkbond=0)
editor.attach_amino_acid('pk1', 'ace')
save TRIPEPTIDE, COMPLEX
