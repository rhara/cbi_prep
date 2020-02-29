run:
	python prep.py data.csv data -p 5

clean:
	rm -r $$(find . -name '__pycache__') || exit 0
	rm $$(find data -name '*.pdb') || exit 0
	rm $$(find data -name '*.apo.pdb.gz') || exit 0
	rm $$(find data -name '*.sdf') || exit 0
	rm $$(find data -name '*.pair*.pkl') || exit 0
