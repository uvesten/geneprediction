# Gene Predictor

A simple prokaryotic gene predictor, based on shine-delgarno boxes, ORF:s, and GC:content. 

## The command line app

Since it makes heavy use of Python 3.5's type hinte (PEP 484), it needs Python >= 3.5

(However, older versions of Python 3 with `mypy` might also work.)

The new predictor also needs `numpy` installed. 


To run: 

    ./predictor.py <genome_file>.fasta

## The web app

The web app is built using cherrypy and jinja2. 

To install: 

    pip3 install -r requirements.txt


To run: 


    python3 webapp.py


## TODO

* Add Shine-Dalgarno sequences to GFF3 output
* Add more statistics 
