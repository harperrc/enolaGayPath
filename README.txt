here are my hacky scripts

mappm.line & fullCoasts.line are lat,lon maps of the globe

i use gnuplot for plotting (the .des files contain various commands)

path.raw      is the data i got from the IMAGES folder
processRay.py convert path.raw into a set of data files-> enolaGayPath.txt 6km.dat 10mi.dat shikoku.dat
evalManeu.py  takes the results of "The Math Of Saving The Enola Gay" by Dr. Jorge S. Diaz (youtuber)
              and generates a plot of the data
hardway.py    uses the equations Diaz came up with but solves the equations in a brute force
              manner (expands the polynomials and solves.  i like to keep my mind busy)
can.py        is a very crude integration of the bomb path (mainly used for plotting) 
              (i have a runge kutta that i checked it with (using Cd & mass from a document i found) 
               and its good enought for this)
eglib.py      are various routines and constants needed

if you have any questions send me an e-mail robertharper@comcast.net
