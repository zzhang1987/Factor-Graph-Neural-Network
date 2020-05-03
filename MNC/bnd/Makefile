# MAKEFILE FOR `code' 
include ../ansi/Makefile
mncN other_programs: ALWAYS
	cd $(bin) ; make $@
tar:
	tar cvf ../pub/code/code.tar ./*.c ./*.h ./*_var ./*_var6 ./commander.p ./Makefile ./ansi/*.c ./ansi/*.h ./ansi/Makefile 
	compress -f  ../pub/code/code.tar


nn_var_usg.c:	nn_var commander.p
	commander.p nn_var
fe_var_usg.c:	fe_var commander.p
	commander.p fe_var
fe_var6_usg.c:	fe_var6 commander.p
	commander.p fe_var6
fed2_var_usg.c:	fed2_var commander.p
	commander.p fed2_var
dc_var_usg.c:	dc_var commander.p
	commander.p dc_var
dc_var6_usg.c:	dc_var6 commander.p
	commander.p dc_var6
code4_var_usg.c:	code4_var commander.p
	commander.p code4_var
code5_var_usg.c:	code5_var commander.p
	commander.p code5_var
code_ez_var_usg.c:	code_ez_var commander.p
	commander.p code_ez_var
bnd_var_usg.c:  bnd_var commander.p
	commander.p bnd_var
codemod_var_usg.c:	codemod_var commander.p
	commander.p codemod_var
bndi_var_usg.c: bndi_var commander.p
	commander.p bndi_var
