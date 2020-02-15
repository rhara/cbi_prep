.SHELL = bash

default:
	@echo Please specify target
	@for f in $$(cat Makefile | grep -E '^\w.*:' | cut -d: -f1) ; do echo '-' $$f ; done

diff_all: diff_get_pdb diff_run diff_values diff_pairing

diff_get_pdb:
	python get.py --csv diff_2018_2019.csv --out pdb_diff

diff_run:
	python run1.py --pdb pdb_diff --csv diff_2018_2019.csv --out data_diff
	python run2.py --pdb pdb_diff --out data_diff

diff_values:
	python val.py --csv diff_2018_2019.csv --out data_diff

diff_pairing:
	python pairing.py --dir data_diff

diff_archive:
	tar zcvf pdb_diff.tar.gz pdb_diff/
	tar zcvf data_diff.tar.gz data_diff/

select_all: select_get_pdb select_run select_values select_pairing

select_get_pdb:
	python get.py --csv pdb2018_kd_select.csv --out pdb_select

select_run:
	python run1.py --pdb pdb_select --csv pdb2018_kd_select.csv --out data_select
	python run2.py --pdb pdb_select --out data_select

select_values:
	python val.py --csv pdb2018_kd_select.csv --out data_select

select_pairing:
	python pairing.py --dir data_select

select_archive:
	tar zcvf pdb_select.tar.gz pdb_select/
	tar zcvf data_select.tar.gz data_select/

clean:
	rm -f *.log
