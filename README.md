# Gene Predictor

A simple prokaryotic gene predictor, based on shine-delgarno boxes, ORF:s, and GC:content. 

## The command line app

It should be executable on any modern Python 3 installation (though only tested with 3.5)

To run: 

    ./predictor.py <genome_file>.fasta

## The web app

The web app is built using cherrypy and jinja2. 

To install: 

    pip3 install -r requirements.txt


To run: 


    python3 webapp.py


## TODO

* Clean up the code, and document
* Add the possibility to download an output file from the web interface
* Add more statistics 
