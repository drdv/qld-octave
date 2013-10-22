octave: forceclean
	octave --no-history --silent make_octave.m

matlab: forceclean
	matlab -nosplash -nodesktop -nojvm -r 'make_matlab; quit' | tail +10

clean:
	rm -f *.o

forceclean: clean
	rm -f *.mex*
