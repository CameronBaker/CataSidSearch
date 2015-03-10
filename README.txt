Cameron Baker

processproteins is a bash script to automatically run css_script.py through a list
of proteins

details of css_script.py are contained here http://catsid.llnl.gov/catsid/css_script.html

the proteins must be in a unix file format (dos2unix on windows txt files) at 1 protein per line

css_script.py will query the website server for a given protein and store the motif 
results for that protein. The script then retrieves the file from the server