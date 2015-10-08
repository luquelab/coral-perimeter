# coral-perimeter
A set of algorithms to calculate the fractal dimension of coral boundaries. 

Image requirements: 
 - The corals are stored as binary images; image names and the pixels per centimeter parameter are stored in "pc_values.csv"

Data flow: Driver --> Model --> Driver
 - Driver: Richardson_Masteroftheuniverse.m
	This driver: 
	 - Loads the "pc_values.csv" file to get the names of images to read along with the "pixels per centimeter" parameter corresponding to each image
	 - Loads the images from the current directory one at a time
	 - Creates the "Ruler Length" vector and "Boundary Points" vector
	 - Writes the vectors to text files ("stepsizes.dat" and "boundary.dat")
	 - Writes a "params.dat" file that gives the lengths of the vectors: boundary vector length followed by ruler vector length
	 - Calls "richDist" precompiled C code
	 - Reads the "richdist.dat" and "polypts.dat" files that are written by "richDist"
	 - Writes the calculated "ruler length," "perimeter length," and "area" values, adjusted to be in "cm," to disk under the name of the image that was analyzed
	 - Writes the linear fit model that corresponds to the calculated "1-D" (slope) and "log10(M)" (intercept) values to disk under the name of the image that was analyzed

 - Model: richDist (C Code: richardsonDistance.c)
	This model:
	 - Loads the "params.dat," "stepsizes.dat," and "boundary.dat" data into memory
	 - Iterates over the stepsizes from "stepsizes.dat" and calculates the perimeter of the boundary from "boundary.dat"
	 - Simultaneously calculates the area enclosed by the perimeter
	 - Writes the results to file as "ruler length," "perimeter length," and "area"

UNUSED: richardsonDistance.m
