# Filterimage - a 2D Savitzky–Golay filtering tool

Filterimage implements 2D Savitzky–Golay filters for images. It allows 
reading and writing of various image formats and the construction and 
application of Savitzky–Golay convolution filters. Filterimage can 
treat the edges of the image in 2 ways. Either the image is 
periodically continued _or_ the edges are treated seperately by 
constructing a set of Savitzky–Golay filters for the edges of the 
image. In the latter case, the image edges are filtered with _similar_ 
filters as the rest of the image but with adapted kernel sizes such 
that the kernel does not cross image boundaries. 

Filterimage can parse a sequence of commands to load, show, filter, and 
save images. The commands may be passed directly on the commandline, in
a Filterimage script file, or on an interactive shell.
 
Usage:

`FilterImage [-i] [-f \<input file\>] [command] [command] [command] ...`


Arguments are processed sequentially in the order of appearance, i.e. 
I can write an input file and execute its commands and after that have 
FilterImage drop you to an interactive shell where you can continue 
working on the images loaded/generated. To get information on available 
commands type the command "help".

To get you started a small sample script can be found in test.fi. To 
try it out run:

`FilterImage -f test.fi` 
