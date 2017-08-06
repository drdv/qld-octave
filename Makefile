all: matlab install

matlab: get_qld
	matlab -nosplash -nodesktop -nojvm -r 'make_matlab; quit' | tail +10

octave: get_qld
	octave --no-history --silent make_octave.m

install:
	mv *.mex* ~/local/bin/qp_mex
	cp qld.m ~/local/bin/qp_mex

get_qld:
	cp ~/local/opt/qld/fortran/qld-2.13.f ./qld.f

clean:
	rm -f *.o
	rm -f *.mex*
	rm -f qld.f
